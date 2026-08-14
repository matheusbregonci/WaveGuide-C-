#pragma once

#include "Geometry.hpp" // VoxelMask, Aabb

#include <vector>

// Self-contained numerical eigenmode solver for arbitrary voxel geometries.
//
// It solves the scalar Helmholtz eigenproblem on the solid region of a
// VoxelMask with Dirichlet walls:
//
//     -∇²ψ = k² ψ   in Ω,     ψ = 0 on ∂Ω
//
// discretized with the standard 7-point finite-difference stencil (FDFD).
// The lowest few (k², ψ) pairs are found matrix-free with subspace inverse
// iteration (a block of conjugate-gradient solves + Rayleigh–Ritz), so no
// external linear-algebra library is needed.
//
// Physical meaning: for a rectangular box these eigenvalues are exactly
// (mπ/a)² + (nπ/b)² + (lπ/d)² — the cavity resonance wavenumbers — so the
// solver reproduces the analytic result. For an arbitrary shape (e.g. a
// T-junction) it gives the scalar resonance modes; the resonant frequency of
// a mode is f = k·c / (2π). (A full vector Maxwell solver is a later step.)

namespace waveguide
{

    struct EigenMode
    {
        double lambda = 0.0;      // eigenvalue k^2 [1/m^2]
        std::vector<double> psi;  // scalar field over the whole grid (0 outside)
    };

    struct HelmholtzResult
    {
        int nx = 0, ny = 0, nz = 0;
        Aabb box{};                       // grid extent (meters)
        double dx = 0, dy = 0, dz = 0;    // cell sizes
        std::vector<std::uint8_t> solid;  // occupancy mask (nx*ny*nz)
        std::vector<EigenMode> modes;     // ascending lambda

        // Node center position of grid index (i,j,k), in meters.
        double px(int i) const { return box.xmin + (i + 0.5) * dx; }
        double py(int j) const { return box.ymin + (j + 0.5) * dy; }
        double pz(int k) const { return box.zmin + (k + 0.5) * dz; }
        std::size_t idx(int i, int j, int k) const
        {
            return (std::size_t(k) * ny + j) * nx + i;
        }
        // Resonant frequency of a mode [Hz], vacuum light speed.
        double frequency(int mode) const;
    };

    class HelmholtzSolver
    {
    public:
        // Solve the lowest `nModes` Dirichlet modes on the mask.
        //   cgIters : max conjugate-gradient iterations per inverse-iteration solve
        //   subIters: subspace (block inverse) iterations
        static HelmholtzResult solveDirichlet(const VoxelMask& mask,
                                              int nModes,
                                              int cgIters = 400,
                                              int subIters = 30);
    };

} // namespace waveguide
