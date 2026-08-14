#pragma once

#include "FieldSource.hpp"
#include "HelmholtzSolver.hpp"
#include "MaxwellSolver.hpp"

#include <cstdint>
#include <vector>

// A FieldSource backed by a numerically solved cavity mode.
//
// It can wrap either:
//   - a scalar Helmholtz mode (Phase 3): the visualized field is grad(psi),
//     the transverse-E analog of a TM mode; or
//   - a full vector Maxwell mode (Phase 4): the visualized field is E itself.
//
// Both are presented as a per-cell vector field (gx,gy,gz) sampled by trilinear
// interpolation. Model coordinates are u in [0,W]x[0,H]x[0,D] (the shape's
// bounding box shifted to the origin), matching the [0,size] convention the
// renderer already centers.

namespace waveguide
{

    class NumericalModel : public FieldSource
    {
    public:
        explicit NumericalModel(const HelmholtzResult& r); // scalar: field = grad(psi)
        explicit NumericalModel(const MaxwellResult& r);   // vector: field = E

        void setMode(int mode);
        int mode() const { return mode_; }
        int modeCount() const { return int(modes_.size()); }
        bool isVector() const { return vector_; }

        // ---- FieldSource ----
        std::array<double, 3> fieldVector(double x, double y, double z,
                                          double phase = 0.0) const override;
        std::pair<double, double> transverseField(double x, double y, double z,
                                                   double phase = 0.0) const override;
        std::vector<Particle> sampleGrid(int nx, int ny, int nz,
                                         bool cutawayOn = true,
                                         float minIntensity = 0.05f,
                                         double phase = 0.0) const override;
        Bounds bounds() const override { return {float(W_), float(H_), float(D_)}; }
        double peakField() const override
        { return modes_.empty() ? 1.0 : modes_[mode_].peak; }
        FieldKind fieldKind() const override { return FieldKind::Electric; }
        double cutoffWavenumber() const override;
        double resonantFrequency() const override;
        double epsilonRel() const override { return 1.0; }
        double muRel() const override { return 1.0; }
        bool inside(double x, double y, double z) const override;

    private:
        struct Mode { double lambda; std::vector<double> gx, gy, gz; double peak; };

        void initGrid(int nx, int ny, int nz, const Aabb& box,
                      std::vector<std::uint8_t> solid);
        void addPeak(Mode& m);
        double sampleField(const std::vector<double>& A,
                           double ux, double uy, double uz) const;
        std::array<double, 3> vecAt(double ux, double uy, double uz) const;
        std::size_t gidx(int i, int j, int k) const { return (std::size_t(k) * ny_ + j) * nx_ + i; }

        bool vector_ = false;
        int nx_ = 0, ny_ = 0, nz_ = 0;
        double dx_ = 0, dy_ = 0, dz_ = 0;
        double W_ = 0, H_ = 0, D_ = 0;
        std::vector<std::uint8_t> solid_;
        std::vector<Mode> modes_;
        int mode_ = 0;
    };

} // namespace waveguide
