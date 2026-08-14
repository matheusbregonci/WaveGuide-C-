#include "MicrostripXsec.hpp"

#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <algorithm>

namespace waveguide
{
    namespace
    {
        constexpr double kEps0 = 8.8541878128e-12;
        constexpr double kC0   = 299792458.0;

        // Solve div(eps * grad phi) = 0 on an nz-by-ny cell grid (uniform square
        // cell dx) with a finite-volume 5-point stencil. The strip cells are held
        // at phi = 1; the outer box (ground at the bottom, walls on the sides and
        // top) is held at phi = 0. Returns C = eps0 * dx^2 * sum(er * |grad phi|^2)
        // -- the per-unit-length capacitance for a 1 V drive.
        //
        // Only the free (non-strip) cells are unknowns, so the system matrix is
        // symmetric positive-definite and plain conjugate gradient converges.
        double solveCap(int nz, int ny, double dx,
                        const std::vector<double>& er,       // relative permittivity per cell
                        const std::vector<std::uint8_t>& strip, // 1 = conductor (phi = 1)
                        int& itersOut)
        {
            const int N = nz * ny;
            auto id = [&](int i, int j){ return j * nz + i; };
            // Harmonic mean of two cell permittivities = the face conductance for
            // piecewise-constant media (series capacitance across the face).
            auto gf = [](double a, double b){ return 2.0 * a * b / (a + b + 1e-30); };

            // Enumerate the free cells.
            std::vector<int> freeOf(N, -1);
            std::vector<int> cellOf; cellOf.reserve(N);
            for (int c = 0; c < N; ++c)
                if (!strip[c]) { freeOf[c] = int(cellOf.size()); cellOf.push_back(c); }
            const int M = int(cellOf.size());
            if (M == 0) { itersOut = 0; return 0.0; }

            // Precompute the diagonal, the RHS (from the phi=1 strip neighbours),
            // and the free-cell adjacency (up to 4 neighbours with their g).
            std::vector<double> diag(M, 0.0), b(M, 0.0);
            std::vector<std::array<int,4>>    nbIdx(M, {-1,-1,-1,-1});
            std::vector<std::array<double,4>> nbG(M,   {0,0,0,0});
            for (int m = 0; m < M; ++m) {
                const int c = cellOf[m];
                const int i = c % nz, j = c / nz;
                const double ec = er[c];
                int slot = 0;
                const int di[4] = {+1,-1,0,0}, dj[4] = {0,0,+1,-1};
                for (int d = 0; d < 4; ++d) {
                    const int ii = i + di[d], jj = j + dj[d];
                    double en, g;
                    if (ii < 0 || ii >= nz || jj < 0 || jj >= ny) {
                        // ground / wall ghost at phi = 0: adds to diag, not to b.
                        en = ec; g = gf(ec, en); diag[m] += g;
                    } else {
                        const int nc = id(ii, jj); en = er[nc]; g = gf(ec, en);
                        diag[m] += g;
                        if (strip[nc]) b[m] += g * 1.0;        // phi = 1 -> RHS
                        else { nbIdx[m][slot] = freeOf[nc]; nbG[m][slot] = g; ++slot; }
                    }
                }
            }

            // A * x  (symmetric): diag*x - sum_neighbours g*x_neighbour.
            auto applyA = [&](const std::vector<double>& x, std::vector<double>& y){
                for (int m = 0; m < M; ++m) {
                    double s = diag[m] * x[m];
                    for (int k = 0; k < 4 && nbIdx[m][k] >= 0; ++k) s -= nbG[m][k] * x[nbIdx[m][k]];
                    y[m] = s;
                }
            };

            // Conjugate gradient.
            std::vector<double> x(M, 0.0), r(b), p(b), Ap(M, 0.0);
            auto dot = [&](const std::vector<double>& a, const std::vector<double>& c){
                double s = 0; for (int m = 0; m < M; ++m) s += a[m] * c[m]; return s; };
            double rs = dot(r, r);
            const double bnorm = std::sqrt(std::max(rs, 1e-300));
            int it = 0; const int maxIt = 20 * M + 500;
            for (; it < maxIt; ++it) {
                applyA(p, Ap);
                const double alpha = rs / std::max(dot(p, Ap), 1e-300);
                for (int m = 0; m < M; ++m) { x[m] += alpha * p[m]; r[m] -= alpha * Ap[m]; }
                const double rsNew = dot(r, r);
                if (std::sqrt(rsNew) <= 1e-8 * bnorm) { ++it; break; }
                const double beta = rsNew / rs;
                for (int m = 0; m < M; ++m) p[m] = r[m] + beta * p[m];
                rs = rsNew;
            }
            itersOut = it;

            // Reconstruct the full potential (free values, strip = 1, ghost = 0).
            std::vector<double> phi(N, 0.0);
            for (int m = 0; m < M; ++m) phi[cellOf[m]] = x[m];
            for (int c = 0; c < N; ++c) if (strip[c]) phi[c] = 1.0;

            // C = 2 * U / V^2 (V = 1). U = 1/2 eps0 sum er |grad phi|^2 dx^2, so
            // C = eps0 * dx^2 * sum er |grad phi|^2. Central differences, phi = 0
            // outside the box.
            auto at = [&](int i, int j){ return (i<0||i>=nz||j<0||j>=ny) ? 0.0 : phi[id(i,j)]; };
            double sum = 0.0;
            for (int j = 0; j < ny; ++j) for (int i = 0; i < nz; ++i) {
                const double gx = (at(i+1,j) - at(i-1,j)) / (2.0 * dx);
                const double gy = (at(i,j+1) - at(i,j-1)) / (2.0 * dx);
                sum += er[id(i,j)] * (gx*gx + gy*gy);
            }
            return kEps0 * (dx * dx) * sum;
        }
    } // namespace

    XsecResult solveMicrostripXsec(double W, double H, double epsr, int cellsPerH)
    {
        XsecResult R;
        if (W <= 0 || H <= 0 || cellsPerH < 4) return R;
        const double dx = H / cellsPerH;
        // Box margins: several substrate heights (and strip widths) away from the
        // conductor so the walls barely perturb the fringing field.
        const double sideM = std::max(6.0 * H, 3.0 * W);
        const double topM  = std::max(6.0 * H, 2.0 * W);
        const double Zdom = W + 2.0 * sideM;
        const double Ydom = H + topM;
        const int nz = std::max(8, int(Zdom / dx + 0.5));
        const int ny = std::max(8, int(Ydom / dx + 0.5));
        const int N = nz * ny;

        // Strip sits at the top of the substrate (row nearest y = H), spanning
        // |z| <= W/2. Ground is the bottom boundary (y = 0), i.e. below row 0.
        const int jStrip = std::clamp(int(H / dx + 0.5) - 1, 1, ny - 2);
        const double z0 = -Zdom * 0.5;   // z of cell i = z0 + (i+0.5)*dx

        std::vector<std::uint8_t> strip(N, 0);
        for (int i = 0; i < nz; ++i) {
            const double zc = z0 + (i + 0.5) * dx;
            if (std::fabs(zc) <= 0.5 * W) strip[jStrip * nz + i] = 1;
        }

        // Permittivity fields: dielectric (er below y=H) and air (er=1 everywhere).
        std::vector<double> erDiel(N, 1.0), erAir(N, 1.0);
        for (int j = 0; j < ny; ++j) {
            const double yc = (j + 0.5) * dx;
            if (yc < H) for (int i = 0; i < nz; ++i) erDiel[j * nz + i] = epsr;
        }

        int itd = 0, ita = 0;
        R.C    = solveCap(nz, ny, dx, erDiel, strip, itd);
        R.Cair = solveCap(nz, ny, dx, erAir,  strip, ita);
        if (R.C <= 0 || R.Cair <= 0) return R;
        R.eeff = R.C / R.Cair;
        R.Z0   = 1.0 / (kC0 * std::sqrt(R.C * R.Cair));
        R.nz = nz; R.ny = ny; R.iters = itd; R.ok = true;
        return R;
    }

    XsecResult hammerstadMicrostrip(double W, double H, double epsr)
    {
        XsecResult R;
        if (W <= 0 || H <= 0) return R;
        const double u = W / H;
        const double a = 1.0 + (1.0 / 12.0) / std::sqrt(1.0 + 12.0 / u); // (1+12H/W)^-1/2 helper
        double eeff = 0.5 * (epsr + 1.0) + 0.5 * (epsr - 1.0) / std::sqrt(1.0 + 12.0 / u);
        if (u < 1.0) eeff += 0.5 * (epsr - 1.0) * 0.04 * (1.0 - u) * (1.0 - u);
        double Z0;
        const double kEta = 376.730313668; // free-space impedance
        if (u <= 1.0)
            Z0 = (60.0 / std::sqrt(eeff)) * std::log(8.0 / u + u / 4.0);
        else
            Z0 = (kEta / std::sqrt(eeff)) / (u + 1.393 + 0.667 * std::log(u + 1.444));
        (void)a;
        R.eeff = eeff; R.Z0 = Z0; R.ok = true;
        return R;
    }
}
