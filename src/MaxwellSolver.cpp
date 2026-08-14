#include "MaxwellSolver.hpp"

#include "EigenSolve.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace waveguide
{
    namespace
    {
        constexpr double kC0 = 299792458.0;
        constexpr double kPi = 3.14159265358979323846;

        // Yee-grid vector operator with Helmholtz (divergence) cleaning.
        // A' E = curl-curl(E) + beta*grad(div E) is used only to make the CG
        // solves invertible; every iterate is then projected onto the
        // divergence-free subspace, so the subspace iteration converges to the
        // physical curl-curl modes without spurious-gradient contamination.
        struct YeeOp
        {
            int nx, ny, nz;
            double dx, dy, dz, beta;
            mutable std::vector<double> Ex, Ey, Ez, Hx, Hy, Hz, phi;
            std::vector<double> outEx, outEy, outEz;
            std::vector<int> dComp, dIdx;
            std::vector<std::uint8_t> interior; // interior nodes (div=0 imposed here only)

            int nEx() const { return nx * (ny + 1) * (nz + 1); }
            int nEy() const { return (nx + 1) * ny * (nz + 1); }
            int nEz() const { return (nx + 1) * (ny + 1) * nz; }
            int iEx(int i, int j, int k) const { return (k * (ny + 1) + j) * nx + i; }
            int iEy(int i, int j, int k) const { return (k * ny + j) * (nx + 1) + i; }
            int iEz(int i, int j, int k) const { return (k * (ny + 1) + j) * (nx + 1) + i; }
            int iHx(int i, int j, int k) const { return (k * ny + j) * (nx + 1) + i; }
            int iHy(int i, int j, int k) const { return (k * (ny + 1) + j) * nx + i; }
            int iHz(int i, int j, int k) const { return (k * ny + j) * nx + i; }
            int iNode(int i, int j, int k) const { return (k * (ny + 1) + j) * (nx + 1) + i; }

            double gEx(int i, int j, int k) const {
                if (i < 0 || i >= nx || j < 0 || j > ny || k < 0 || k > nz) return 0.0;
                return Ex[iEx(i, j, k)]; }
            double gEy(int i, int j, int k) const {
                if (i < 0 || i > nx || j < 0 || j >= ny || k < 0 || k > nz) return 0.0;
                return Ey[iEy(i, j, k)]; }
            double gEz(int i, int j, int k) const {
                if (i < 0 || i > nx || j < 0 || j > ny || k < 0 || k >= nz) return 0.0;
                return Ez[iEz(i, j, k)]; }
            double gHx(int i, int j, int k) const {
                if (i < 0 || i > nx || j < 0 || j >= ny || k < 0 || k >= nz) return 0.0;
                return Hx[iHx(i, j, k)]; }
            double gHy(int i, int j, int k) const {
                if (i < 0 || i >= nx || j < 0 || j > ny || k < 0 || k >= nz) return 0.0;
                return Hy[iHy(i, j, k)]; }
            double gHz(int i, int j, int k) const {
                if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k > nz) return 0.0;
                return Hz[iHz(i, j, k)]; }
            double gPhi(int i, int j, int k) const {
                if (i < 0 || i > nx || j < 0 || j > ny || k < 0 || k > nz) return 0.0;
                return phi[iNode(i, j, k)]; }

            void alloc()
            {
                Ex.assign(nEx(), 0.0); Ey.assign(nEy(), 0.0); Ez.assign(nEz(), 0.0);
                Hx.assign((nx + 1) * ny * nz, 0.0);
                Hy.assign(nx * (ny + 1) * nz, 0.0);
                Hz.assign(nx * ny * (nz + 1), 0.0);
                phi.assign((nx + 1) * (ny + 1) * (nz + 1), 0.0);
                outEx = Ex; outEy = Ey; outEz = Ez;
            }
            void zeroE() const
            {
                std::fill(Ex.begin(), Ex.end(), 0.0);
                std::fill(Ey.begin(), Ey.end(), 0.0);
                std::fill(Ez.begin(), Ez.end(), 0.0);
            }
            int Ndof() const { return int(dComp.size()); }

            void scatterDof(const std::vector<double>& x) const
            {
                zeroE();
                const int Nd = int(dComp.size());
                for (int d = 0; d < Nd; ++d)
                {
                    if (dComp[d] == 0) Ex[dIdx[d]] = x[d];
                    else if (dComp[d] == 1) Ey[dIdx[d]] = x[d];
                    else Ez[dIdx[d]] = x[d];
                }
            }
            void computeH() const
            {
                for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i)
                    Hx[iHx(i, j, k)] = (gEz(i, j + 1, k) - gEz(i, j, k)) / dy
                                     - (gEy(i, j, k + 1) - gEy(i, j, k)) / dz;
                for (int k = 0; k < nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i)
                    Hy[iHy(i, j, k)] = (gEx(i, j, k + 1) - gEx(i, j, k)) / dz
                                     - (gEz(i + 1, j, k) - gEz(i, j, k)) / dx;
                for (int k = 0; k <= nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
                    Hz[iHz(i, j, k)] = (gEy(i + 1, j, k) - gEy(i, j, k)) / dx
                                     - (gEx(i, j + 1, k) - gEx(i, j, k)) / dy;
            }
            void computeDiv() const
            {
                for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i)
                    phi[iNode(i, j, k)] =
                        (gEx(i, j, k) - gEx(i - 1, j, k)) / dx +
                        (gEy(i, j, k) - gEy(i, j - 1, k)) / dy +
                        (gEz(i, j, k) - gEz(i, j, k - 1)) / dz;
            }
            // out = curl H + beta*grad(div) ; needs H and phi already computed.
            void curlHpenalty()
            {
                for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i)
                    outEx[iEx(i, j, k)] = (gHz(i, j, k) - gHz(i, j - 1, k)) / dy
                                        - (gHy(i, j, k) - gHy(i, j, k - 1)) / dz
                                        + beta * (gPhi(i, j, k) - gPhi(i + 1, j, k)) / dx;
                for (int k = 0; k <= nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i)
                    outEy[iEy(i, j, k)] = (gHx(i, j, k) - gHx(i, j, k - 1)) / dz
                                        - (gHz(i, j, k) - gHz(i - 1, j, k)) / dx
                                        + beta * (gPhi(i, j, k) - gPhi(i, j + 1, k)) / dy;
                for (int k = 0; k < nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i)
                    outEz[iEz(i, j, k)] = (gHy(i, j, k) - gHy(i - 1, j, k)) / dx
                                        - (gHx(i, j, k) - gHx(i, j - 1, k)) / dy
                                        + beta * (gPhi(i, j, k) - gPhi(i, j, k + 1)) / dz;
            }
            void gatherDof(std::vector<double>& y) const
            {
                const int Nd = int(dComp.size());
                y.assign(Nd, 0.0);
                for (int d = 0; d < Nd; ++d)
                {
                    if (dComp[d] == 0) y[d] = outEx[dIdx[d]];
                    else if (dComp[d] == 1) y[d] = outEy[dIdx[d]];
                    else y[d] = outEz[dIdx[d]];
                }
            }
            void applyAprime(const std::vector<double>& x, std::vector<double>& y)
            {
                scatterDof(x); computeH(); computeDiv();
                // Penalize only the INTERIOR divergence: boundary nodes carry
                // legitimate surface charge (normal E at PEC walls), so
                // penalizing them would wrongly lift physical modes like TE101.
                if (!interior.empty())
                    for (std::size_t n = 0; n < phi.size(); ++n) if (!interior[n]) phi[n] = 0.0;
                curlHpenalty(); gatherDof(y);
            }
            // Set Ex/Ey/Ez to grad(psi) restricted to the free (DOF) edges.
            void gradDof(const std::vector<double>& psi) const
            {
                zeroE();
                const int Nd = int(dComp.size());
                for (int d = 0; d < Nd; ++d)
                {
                    const int idx = dIdx[d];
                    if (dComp[d] == 0) { // Ex edge (i,j,k) -> (psi(i+1)-psi(i))/dx
                        const int i = idx % nx; const int r = idx / nx;
                        const int j = r % (ny + 1); const int k = r / (ny + 1);
                        Ex[idx] = (phiVal(psi, i + 1, j, k) - phiVal(psi, i, j, k)) / dx;
                    } else if (dComp[d] == 1) {
                        const int i = idx % (nx + 1); const int r = idx / (nx + 1);
                        const int j = r % ny; const int k = r / ny;
                        Ey[idx] = (phiVal(psi, i, j + 1, k) - phiVal(psi, i, j, k)) / dy;
                    } else {
                        const int i = idx % (nx + 1); const int r = idx / (nx + 1);
                        const int j = r % (ny + 1); const int k = r / (ny + 1);
                        Ez[idx] = (phiVal(psi, i, j, k + 1) - phiVal(psi, i, j, k)) / dz;
                    }
                }
            }
            double phiVal(const std::vector<double>& psi, int i, int j, int k) const
            {
                if (i < 0 || i > nx || j < 0 || j > ny || k < 0 || k > nz) return 0.0;
                return psi[iNode(i, j, k)];
            }
            // Poisson operator M psi = -div(gradDof(psi)), restricted to interior
            // nodes (psi = 0 on the boundary; div=0 only where there is no charge).
            void applyM(const std::vector<double>& psiIn, std::vector<double>& out) const
            {
                std::vector<double> psi = psiIn;
                for (std::size_t n = 0; n < psi.size(); ++n) if (!interior[n]) psi[n] = 0.0;
                gradDof(psi); computeDiv();
                out.assign(phi.size(), 0.0);
                for (std::size_t n = 0; n < phi.size(); ++n) out[n] = interior[n] ? -phi[n] : 0.0;
            }
            // Remove the gradient part of the DOF field x so its divergence
            // vanishes at the interior nodes (physical, charge-free condition).
            void projectDivFree(std::vector<double>& x, int poissonIters)
            {
                scatterDof(x); computeDiv();
                std::vector<double> negRho(phi.size(), 0.0);
                for (std::size_t n = 0; n < phi.size(); ++n)
                    if (interior[n]) negRho[n] = -phi[n];
                eig::ApplyOp M = [this](const std::vector<double>& in, std::vector<double>& o)
                { applyM(in, o); };
                std::vector<double> psi;
                eig::cg(M, negRho, psi, poissonIters, 1e-7);
                for (std::size_t n = 0; n < psi.size(); ++n) if (!interior[n]) psi[n] = 0.0;
                gradDof(psi); // Ex/Ey/Ez now hold grad psi on DOF edges
                const int Nd = int(dComp.size());
                for (int d = 0; d < Nd; ++d)
                {
                    if (dComp[d] == 0) x[d] -= Ex[dIdx[d]];
                    else if (dComp[d] == 1) x[d] -= Ey[dIdx[d]];
                    else x[d] -= Ez[dIdx[d]];
                }
            }
        };
    } // namespace

    double MaxwellSolver::debugTE(const VoxelMask& mask, int m, int l, double* divRatio)
    {
        const int nx = mask.nx, ny = mask.ny, nz = mask.nz;
        YeeOp op;
        op.nx = nx; op.ny = ny; op.nz = nz;
        op.dx = mask.box.sizeX() / nx; op.dy = mask.box.sizeY() / ny;
        op.dz = mask.box.sizeZ() / nz; op.beta = 0.0;
        op.alloc();
        // Ey(i,j,k) = sin(m*pi*i/nx) * sin(l*pi*k/nz), uniform in j.
        for (int k = 0; k <= nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i)
            op.Ey[op.iEy(i, j, k)] = std::sin(m * kPi * i / nx) * std::sin(l * kPi * k / nz);
        std::vector<double> ccEx, ccEy, ccEz;
        op.computeH(); op.computeDiv(); op.curlHpenalty();
        double num = 0.0, den = 0.0;
        for (std::size_t n = 0; n < op.Ey.size(); ++n) { num += op.Ey[n] * op.outEy[n]; den += op.Ey[n] * op.Ey[n]; }
        double curlN = 0.0; for (double v : op.Hx) curlN += v*v; for (double v : op.Hy) curlN += v*v; for (double v : op.Hz) curlN += v*v;
        double divN = 0.0; for (double v : op.phi) divN += v*v;
        if (divRatio) *divRatio = std::sqrt(curlN) / (std::sqrt(curlN) + std::sqrt(divN) + 1e-300);
        return (den > 0) ? num / den : 0.0;
    }

    double MaxwellResult::frequency(int mode) const
    {
        if (mode < 0 || mode >= int(modes.size())) return 0.0;
        return std::sqrt(std::max(0.0, modes[mode].lambda)) * kC0 / (2.0 * kPi);
    }

    MaxwellResult MaxwellSolver::solveCavity(const VoxelMask& mask, int nModes,
                                             double penalty, int cgIters, int subIters)
    {
        MaxwellResult res;
        res.nx = mask.nx; res.ny = mask.ny; res.nz = mask.nz; res.box = mask.box;
        res.solid = mask.occ;
        const int nx = res.nx, ny = res.ny, nz = res.nz;
        res.dx = mask.box.sizeX() / std::max(1, nx);
        res.dy = mask.box.sizeY() / std::max(1, ny);
        res.dz = mask.box.sizeZ() / std::max(1, nz);
        if (res.dx <= 0 || res.dy <= 0 || res.dz <= 0) return res;

        auto solid = [&](int i, int j, int k) -> bool {
            if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return false;
            return mask.occ[res.idx(i, j, k)] != 0;
        };

        YeeOp op;
        op.nx = nx; op.ny = ny; op.nz = nz;
        op.dx = res.dx; op.dy = res.dy; op.dz = res.dz; op.beta = penalty;
        op.alloc();

        for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i < nx; ++i)
            if (solid(i, j - 1, k - 1) && solid(i, j, k - 1) && solid(i, j - 1, k) && solid(i, j, k))
            { op.dComp.push_back(0); op.dIdx.push_back(op.iEx(i, j, k)); }
        for (int k = 0; k <= nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i <= nx; ++i)
            if (solid(i - 1, j, k - 1) && solid(i, j, k - 1) && solid(i - 1, j, k) && solid(i, j, k))
            { op.dComp.push_back(1); op.dIdx.push_back(op.iEy(i, j, k)); }
        for (int k = 0; k < nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i)
            if (solid(i - 1, j - 1, k) && solid(i, j - 1, k) && solid(i - 1, j, k) && solid(i, j, k))
            { op.dComp.push_back(2); op.dIdx.push_back(op.iEz(i, j, k)); }

        // Interior nodes: all 8 surrounding cells solid (div=0 imposed only here;
        // boundary nodes carry surface charge / free normal E).
        op.interior.assign((nx + 1) * (ny + 1) * (nz + 1), 0);
        for (int k = 0; k <= nz; ++k) for (int j = 0; j <= ny; ++j) for (int i = 0; i <= nx; ++i)
        {
            bool inter = true;
            for (int dk = -1; dk <= 0 && inter; ++dk)
            for (int dj = -1; dj <= 0 && inter; ++dj)
            for (int di = -1; di <= 0 && inter; ++di)
                if (!solid(i + di, j + dj, k + dk)) inter = false;
            op.interior[op.iNode(i, j, k)] = inter ? 1 : 0;
        }

        const int Nd = op.Ndof();
        if (Nd == 0) return res;

        // ---- projected subspace inverse iteration (div-free) ----
        const int p = std::min(Nd, std::max(1, nModes) + 4);
        const int poissonIters = std::max(200, 4 * (nx + ny + nz));
        std::vector<std::vector<double>> X(p, std::vector<double>(Nd)), Y(p, std::vector<double>(Nd));
        std::mt19937 rng(20240714u);
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        for (int j = 0; j < p; ++j)
        {
            for (int i = 0; i < Nd; ++i) X[j][i] = uni(rng);
            op.projectDivFree(X[j], poissonIters);
        }
        eig::orthonormalize(X);

        eig::ApplyOp applyA = [&op](const std::vector<double>& x, std::vector<double>& y)
        { op.applyAprime(x, y); };

        std::vector<double> theta(p, 0.0), prev(p, 1e300);
        for (int outer = 0; outer < subIters; ++outer)
        {
            for (int j = 0; j < p; ++j)
            {
                eig::cg(applyA, X[j], Y[j], cgIters, 1e-6);
                op.projectDivFree(Y[j], poissonIters);
            }
            eig::orthonormalize(Y);
            std::vector<std::vector<double>> AY(p, std::vector<double>(Nd));
            for (int j = 0; j < p; ++j) op.applyAprime(Y[j], AY[j]);
            std::vector<double> H(std::size_t(p) * p, 0.0);
            for (int a = 0; a < p; ++a) for (int b = a; b < p; ++b)
            { const double h = eig::dot(Y[a], AY[b]); H[a * p + b] = h; H[b * p + a] = h; }
            std::vector<double> ev, V;
            eig::jacobiEig(H, p, ev, V);
            for (int j = 0; j < p; ++j)
            {
                std::fill(X[j].begin(), X[j].end(), 0.0);
                for (int k = 0; k < p; ++k)
                {
                    const double s = V[k * p + j];
                    if (s == 0.0) continue;
                    for (int i = 0; i < Nd; ++i) X[j][i] += s * Y[k][i];
                }
            }
            theta = ev;
            double maxrel = 0.0;
            const int want = std::min(p, std::max(1, nModes));
            for (int j = 0; j < want; ++j)
                maxrel = std::max(maxrel, std::abs(theta[j] - prev[j]) / (std::abs(theta[j]) + 1e-30));
            prev = theta;
            if (outer > 2 && maxrel < 1e-6) break;
        }

        const int emit = std::min(p, std::max(1, nModes));
        res.modes.resize(emit);
        for (int m = 0; m < emit; ++m)
        {
            op.scatterDof(X[m]);
            op.computeH(); op.computeDiv();
            double curlN = 0.0; for (double v : op.Hx) curlN += v*v; for (double v : op.Hy) curlN += v*v; for (double v : op.Hz) curlN += v*v;
            double divN = 0.0;  for (std::size_t n = 0; n < op.phi.size(); ++n) if (op.interior[n]) divN += op.phi[n]*op.phi[n];
            curlN = std::sqrt(curlN); divN = std::sqrt(divN);
            MaxwellMode& md = res.modes[m];
            md.lambda = theta[m];
            md.divRatio = curlN / (curlN + divN + 1e-300);
            md.ex.assign(std::size_t(nx) * ny * nz, 0.0);
            md.ey.assign(std::size_t(nx) * ny * nz, 0.0);
            md.ez.assign(std::size_t(nx) * ny * nz, 0.0);
            for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
            {
                const std::size_t g = res.idx(i, j, k);
                md.ex[g] = 0.25 * (op.gEx(i, j, k) + op.gEx(i, j + 1, k) + op.gEx(i, j, k + 1) + op.gEx(i, j + 1, k + 1));
                md.ey[g] = 0.25 * (op.gEy(i, j, k) + op.gEy(i + 1, j, k) + op.gEy(i, j, k + 1) + op.gEy(i + 1, j, k + 1));
                md.ez[g] = 0.25 * (op.gEz(i, j, k) + op.gEz(i + 1, j, k) + op.gEz(i, j + 1, k) + op.gEz(i + 1, j + 1, k));
            }
        }
        return res;
    }

} // namespace waveguide
