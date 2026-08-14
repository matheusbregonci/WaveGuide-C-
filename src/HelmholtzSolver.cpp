#include "HelmholtzSolver.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace waveguide
{
    namespace
    {
        constexpr double kC0 = 299792458.0;
        constexpr double kPi = 3.14159265358979323846;

        // Dirichlet 7-point Laplacian as a matrix-free operator over the solid
        // cells. Neighbors that are outside the grid or non-solid contribute 0
        // (psi = 0 walls). The diagonal is constant on a uniform grid.
        struct Laplacian
        {
            int N = 0;                         // number of solid dofs
            double diag = 0.0;                 // 2(1/dx^2+1/dy^2+1/dz^2)
            double ix2 = 0, iy2 = 0, iz2 = 0;  // 1/dx^2, 1/dy^2, 1/dz^2
            std::vector<int> nbr;              // 6*N: x-,x+,y-,y+,z-,z+ dof or -1

            void apply(const std::vector<double>& x, std::vector<double>& y) const
            {
                for (int d = 0; d < N; ++d)
                {
                    double v = diag * x[d];
                    const int* nb = &nbr[std::size_t(d) * 6];
                    if (nb[0] >= 0) v -= ix2 * x[nb[0]];
                    if (nb[1] >= 0) v -= ix2 * x[nb[1]];
                    if (nb[2] >= 0) v -= iy2 * x[nb[2]];
                    if (nb[3] >= 0) v -= iy2 * x[nb[3]];
                    if (nb[4] >= 0) v -= iz2 * x[nb[4]];
                    if (nb[5] >= 0) v -= iz2 * x[nb[5]];
                    y[d] = v;
                }
            }
        };

        double dot(const std::vector<double>& a, const std::vector<double>& b)
        {
            double s = 0.0;
            for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
            return s;
        }

        // Conjugate gradient solve A x = b (A SPD). x starts at 0.
        void cg(const Laplacian& A, const std::vector<double>& b,
                std::vector<double>& x, int maxIt, double tol)
        {
            const int N = A.N;
            x.assign(N, 0.0);
            std::vector<double> r = b, p = b, Ap(N);
            double rs = dot(r, r);
            const double rs0 = rs;
            if (rs0 == 0.0) return;
            for (int it = 0; it < maxIt; ++it)
            {
                A.apply(p, Ap);
                const double denom = dot(p, Ap);
                if (denom <= 0.0) break;
                const double alpha = rs / denom;
                for (int i = 0; i < N; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
                const double rsN = dot(r, r);
                if (rsN <= tol * tol * rs0) break;
                const double beta = rsN / rs;
                for (int i = 0; i < N; ++i) p[i] = r[i] + beta * p[i];
                rs = rsN;
            }
        }

        // Modified Gram-Schmidt orthonormalization of the columns of X (p cols,
        // each length N). Drops (zeros) columns that collapse.
        void orthonormalize(std::vector<std::vector<double>>& X)
        {
            for (std::size_t j = 0; j < X.size(); ++j)
            {
                for (std::size_t k = 0; k < j; ++k)
                {
                    const double d = dot(X[j], X[k]);
                    for (std::size_t i = 0; i < X[j].size(); ++i) X[j][i] -= d * X[k][i];
                }
                double nrm = std::sqrt(dot(X[j], X[j]));
                if (nrm < 1e-300) nrm = 1.0;
                const double inv = 1.0 / nrm;
                for (double& v : X[j]) v *= inv;
            }
        }

        // Cyclic Jacobi eigensolver for a small symmetric dense matrix H (p x p,
        // row-major). Returns eigenvalues (ascending) and eigenvectors as
        // columns of V (row-major p x p): eigenvector j = V[i*p+j].
        void jacobiEig(std::vector<double> H, int p,
                       std::vector<double>& eval, std::vector<double>& V)
        {
            V.assign(std::size_t(p) * p, 0.0);
            for (int i = 0; i < p; ++i) V[i * p + i] = 1.0;
            for (int sweep = 0; sweep < 100; ++sweep)
            {
                double off = 0.0;
                for (int i = 0; i < p; ++i)
                    for (int j = i + 1; j < p; ++j) off += H[i * p + j] * H[i * p + j];
                if (off < 1e-24) break;
                for (int q = 0; q < p; ++q)
                    for (int r = q + 1; r < p; ++r)
                    {
                        const double apq = H[q * p + r];
                        if (std::abs(apq) < 1e-300) continue;
                        const double app = H[q * p + q], arr = H[r * p + r];
                        const double phi = 0.5 * std::atan2(2.0 * apq, arr - app);
                        const double c = std::cos(phi), s = std::sin(phi);
                        for (int k = 0; k < p; ++k)
                        {
                            const double hkq = H[k * p + q], hkr = H[k * p + r];
                            H[k * p + q] = c * hkq - s * hkr;
                            H[k * p + r] = s * hkq + c * hkr;
                        }
                        for (int k = 0; k < p; ++k)
                        {
                            const double hqk = H[q * p + k], hrk = H[r * p + k];
                            H[q * p + k] = c * hqk - s * hrk;
                            H[r * p + k] = s * hqk + c * hrk;
                        }
                        for (int k = 0; k < p; ++k)
                        {
                            const double vkq = V[k * p + q], vkr = V[k * p + r];
                            V[k * p + q] = c * vkq - s * vkr;
                            V[k * p + r] = s * vkq + c * vkr;
                        }
                    }
            }
            eval.resize(p);
            for (int i = 0; i < p; ++i) eval[i] = H[i * p + i];
            // sort ascending, permuting eigenvector columns
            std::vector<int> ord(p);
            for (int i = 0; i < p; ++i) ord[i] = i;
            std::sort(ord.begin(), ord.end(), [&](int a, int b) { return eval[a] < eval[b]; });
            std::vector<double> ev(p);
            std::vector<double> Vp(std::size_t(p) * p);
            for (int j = 0; j < p; ++j)
            {
                ev[j] = eval[ord[j]];
                for (int i = 0; i < p; ++i) Vp[i * p + j] = V[i * p + ord[j]];
            }
            eval = ev;
            V = Vp;
        }
    } // namespace

    double HelmholtzResult::frequency(int mode) const
    {
        if (mode < 0 || mode >= int(modes.size())) return 0.0;
        const double k = std::sqrt(std::max(0.0, modes[mode].lambda));
        return k * kC0 / (2.0 * kPi);
    }

    HelmholtzResult HelmholtzSolver::solveDirichlet(const VoxelMask& mask,
                                                    int nModes, int cgIters,
                                                    int subIters)
    {
        HelmholtzResult res;
        res.nx = mask.nx; res.ny = mask.ny; res.nz = mask.nz;
        res.box = mask.box;
        res.solid = mask.occ;
        res.dx = mask.box.sizeX() / std::max(1, mask.nx);
        res.dy = mask.box.sizeY() / std::max(1, mask.ny);
        res.dz = mask.box.sizeZ() / std::max(1, mask.nz);
        if (res.dx <= 0 || res.dy <= 0 || res.dz <= 0) return res;

        const int nx = res.nx, ny = res.ny, nz = res.nz;

        // Map solid cells to dof indices and record grid coordinates.
        std::vector<int> dofOf(std::size_t(nx) * ny * nz, -1);
        std::vector<int> ci, cj, ck;
        for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                {
                    const std::size_t g = res.idx(i, j, k);
                    if (mask.occ[g]) { dofOf[g] = int(ci.size()); ci.push_back(i); cj.push_back(j); ck.push_back(k); }
                }
        const int N = int(ci.size());
        if (N == 0) return res;

        Laplacian A;
        A.N = N;
        A.ix2 = 1.0 / (res.dx * res.dx);
        A.iy2 = 1.0 / (res.dy * res.dy);
        A.iz2 = 1.0 / (res.dz * res.dz);
        A.diag = 2.0 * (A.ix2 + A.iy2 + A.iz2);
        A.nbr.assign(std::size_t(N) * 6, -1);
        auto solidDof = [&](int i, int j, int k) -> int {
            if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return -1;
            return dofOf[res.idx(i, j, k)];
        };
        for (int d = 0; d < N; ++d)
        {
            const int i = ci[d], j = cj[d], k = ck[d];
            int* nb = &A.nbr[std::size_t(d) * 6];
            nb[0] = solidDof(i - 1, j, k); nb[1] = solidDof(i + 1, j, k);
            nb[2] = solidDof(i, j - 1, k); nb[3] = solidDof(i, j + 1, k);
            nb[4] = solidDof(i, j, k - 1); nb[5] = solidDof(i, j, k + 1);
        }

        // Subspace inverse iteration for the lowest eigenpairs.
        const int p = std::min(N, std::max(1, nModes) + 2); // guard vectors
        std::vector<std::vector<double>> X(p, std::vector<double>(N));
        std::mt19937 rng(12345);
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        for (int j = 0; j < p; ++j)
            for (int i = 0; i < N; ++i) X[j][i] = uni(rng);
        orthonormalize(X);

        std::vector<std::vector<double>> Y(p, std::vector<double>(N));
        std::vector<double> theta(p, 0.0), prevTheta(p, 1e300);
        const double cgTol = 1e-6;
        for (int outer = 0; outer < subIters; ++outer)
        {
            // Y = A^{-1} X (block CG solves), then orthonormalize.
            for (int j = 0; j < p; ++j) cg(A, X[j], Y[j], cgIters, cgTol);
            orthonormalize(Y);

            // Rayleigh-Ritz: H = Y^T A Y (p x p), symmetric.
            std::vector<std::vector<double>> AY(p, std::vector<double>(N));
            for (int j = 0; j < p; ++j) A.apply(Y[j], AY[j]);
            std::vector<double> H(std::size_t(p) * p, 0.0);
            for (int a = 0; a < p; ++a)
                for (int b = a; b < p; ++b)
                {
                    const double h = dot(Y[a], AY[b]);
                    H[a * p + b] = h; H[b * p + a] = h;
                }
            std::vector<double> eval, V;
            jacobiEig(H, p, eval, V);

            // New Ritz vectors X_j = sum_k Y_k V[k][j].
            for (int j = 0; j < p; ++j)
            {
                std::fill(X[j].begin(), X[j].end(), 0.0);
                for (int k = 0; k < p; ++k)
                {
                    const double s = V[k * p + j];
                    if (s == 0.0) continue;
                    for (int i = 0; i < N; ++i) X[j][i] += s * Y[k][i];
                }
            }
            theta = eval;

            // Convergence: relative change of the wanted eigenvalues.
            double maxrel = 0.0;
            const int want = std::min(p, std::max(1, nModes));
            for (int j = 0; j < want; ++j)
            {
                const double denom = std::abs(theta[j]) + 1e-30;
                maxrel = std::max(maxrel, std::abs(theta[j] - prevTheta[j]) / denom);
            }
            prevTheta = theta;
            if (outer > 2 && maxrel < 1e-6) break;
        }

        // Emit the lowest nModes modes, expanded to the full grid.
        const int emit = std::min(p, std::max(1, nModes));
        res.modes.resize(emit);
        for (int m = 0; m < emit; ++m)
        {
            EigenMode& mode = res.modes[m];
            mode.lambda = theta[m];
            mode.psi.assign(std::size_t(nx) * ny * nz, 0.0);
            // sign convention: largest-magnitude sample positive
            double amax = 0.0; int dmax = 0;
            for (int d = 0; d < N; ++d)
                if (std::abs(X[m][d]) > amax) { amax = std::abs(X[m][d]); dmax = d; }
            const double sgn = (X[m][dmax] < 0.0) ? -1.0 : 1.0;
            for (int d = 0; d < N; ++d)
                mode.psi[res.idx(ci[d], cj[d], ck[d])] = sgn * X[m][d];
        }
        return res;
    }

} // namespace waveguide
