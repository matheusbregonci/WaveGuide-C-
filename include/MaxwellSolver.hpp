#pragma once

#include "Geometry.hpp" // VoxelMask, Aabb

#include <cstdint>
#include <vector>

// Full vector Maxwell cavity eigensolver (Phase 4).
//
// Solves the vector curl-curl eigenproblem for a PEC cavity of arbitrary voxel
// shape on a Yee grid (E on edges):
//
//     ∇×∇×E = k² E   in Ω,     n̂×E = 0 on the walls
//
// The discrete curl-curl (CᵀC) has a large null space of gradient fields; a
// grad-div penalty (β DᵀD) lifts those spurious modes, and a final divergence
// filter keeps only the physical (solenoidal) modes. Eigenpairs are found
// matrix-free with the shared subspace solver. Unlike the scalar Helmholtz
// solver this returns the true vector E field (and therefore the TE modes too,
// e.g. the dominant rectangular-cavity TE101). Validated against the analytic
// cavity spectrum f_mnl = (c/2)·sqrt((m/a)²+(n/b)²+(l/d)²).

namespace waveguide
{
    struct MaxwellMode
    {
        double lambda = 0.0;                 // k^2 [1/m^2]
        double divRatio = 0.0;               // ||curlE||/(||curlE||+||divE||); ~1 physical
        std::vector<double> ex, ey, ez;      // cell-centered E field (nx*ny*nz)
    };

    struct MaxwellResult
    {
        int nx = 0, ny = 0, nz = 0;
        Aabb box{};
        double dx = 0, dy = 0, dz = 0;
        std::vector<std::uint8_t> solid;
        std::vector<MaxwellMode> modes;      // ascending lambda (physical only)

        std::size_t idx(int i, int j, int k) const
        {
            return (std::size_t(k) * ny + j) * nx + i;
        }
        double frequency(int mode) const;
    };

    class MaxwellSolver
    {
    public:
        // Solve the lowest `nModes` physical cavity modes.
        //   penalty : grad-div penalty weight (β)
        static MaxwellResult solveCavity(const VoxelMask& mask, int nModes,
                                         double penalty = 3.0,
                                         int cgIters = 600, int subIters = 60);

        // Diagnostic: Rayleigh quotient (k^2) of the analytic TE_{m0l} field
        // (Ey = sin(m*pi*x/a) sin(l*pi*z/d)) under the discrete curl-curl, plus
        // its divergence ratio. Used to validate the operator.
        static double debugTE(const VoxelMask& mask, int m, int l, double* divRatio);
    };

} // namespace waveguide
