#pragma once

#include "Types.hpp"
#include "FieldSource.hpp"

#include <array>
#include <complex>
#include <cstdint>
#include <utility>
#include <vector>

// C++ port of Cilindrico_model.py (WaveGuidesSimulator).
//
// Circular metallic waveguide: TE_nm and TM_nm modes.
//   - n : azimuthal order  (0..2)
//   - m : radial  order    (1..3)
// Uses tabulated zeros of J_n (TM) and J'_n (TE) for the cutoff:
//   TE: k_c = p'_nm / R
//   TM: k_c = p_nm  / R

namespace waveguide
{
    class CylindricalModel : public FieldSource
    {
    public:
        //   cavity : if true, the guide is shorted at z = 0 and z = length,
        //            producing a resonant TE_nml / TM_nml standing wave along
        //            z. `modeL` is the number of half-wave variations along the
        //            length. When false the model is the open waveguide.
        CylindricalModel(double radiusMM = 25.0,
                         double lengthMeters = 0.11,
                         double frequency = 12e9,
                         double epsilonRel = 1.0,
                         double muRel = 1.0,
                         int modeN = 1, // azimuthal
                         int modeM = 1, // radial
                         ModeType type = ModeType::TE,
                         FieldKind field = FieldKind::Electric,
                         bool cavity = false,
                         int modeL = 1);

        double peakField() const override { ensurePeak(); return cachedMaxE_; }
        FieldKind fieldKind() const override { return field_; }

        bool isCavity() const { return cavity_; }
        int modeL() const { return l_; }

        // Resonant frequency f_nml of the closed cavity (Hz). For the open
        // waveguide this returns the transverse cutoff frequency (l = 0).
        double resonantFrequency() const override;

        Bounds bounds() const override
        {
            return {float(2.0 * raio_), float(2.0 * raio_), float(length_)};
        }

        // Transverse coordinates are CENTRED on the axis, unlike the rectangular
        // guide's corner-origin convention. See FieldSource::domain().
        void domain(double& x0, double& x1, double& y0, double& y1,
                    double& z0, double& z1) const override
        {
            x0 = -raio_;  x1 = raio_;
            y0 = -raio_;  y1 = raio_;
            z0 = 0.0;     z1 = length_;
        }

        // Bounding box for the cylinder is [-R,R] x [-R,R] x [0,L]; cells
        // outside the disk (x² + y² > R²) are culled. Sphere-size and color
        // are still driven by intensity in the shader.
        std::vector<Particle> sampleGrid(int nx, int ny, int nz,
                                         bool cutawayOn = true,
                                         float minIntensity = 0.05f,
                                         double phase = 0.0) const override;

        double radius() const { return raio_; }
        double cutoffKc() const { return k_c_; }
        double cutoffWavenumber() const override { return k_c_; }
        double epsilonRel() const override { return epsilon_; }
        double muRel() const override { return mu_; }

        // Circular cross-section: inside the disk x^2 + y^2 <= R^2, 0 <= z <= L.
        bool inside(double x, double y, double z) const override
        {
            return x * x + y * y <= raio_ * raio_ && z >= 0.0 && z <= length_;
        }

        // Real transverse (Vx, Vy) in Cartesian of the chosen field.
        std::pair<double, double> transverseField(double x, double y, double z,
                                                  double phase = 0.0) const override;

        std::array<double, 3> fieldVector(double x, double y, double z,
                                          double phase = 0.0) const override;

    private:
        double magnitudeE(double x, double y, double z, double phase) const;

        // Fill cachedMaxE_ on first use (const + mutable cache).
        void ensurePeak() const;

        // Real Cartesian field vector of the CLOSED cavity at (x,y,z) and time
        // `phase`. Transverse components carry sin(lπz/d), the axial component
        // cos/sin(lπz/d) per the shorting-wall boundary conditions; E and H are
        // 90° apart in time. Used by the public accessors when cavity_ is set.
        std::array<double, 3> cavityFieldVec(double x, double y, double z,
                                             double phase) const;

        double raio_;
        double length_;
        double frequencia_;
        double epsilon_;
        double mu_;
        int n_;
        int m_;
        ModeType type_;
        FieldKind field_;
        double A_ = 1.0;
        double B_ = 1.0;
        bool cavity_ = false;
        int l_ = 1; // z half-wave index (cavity only)

        double omega_;
        double k_;
        double k_c_;
        std::complex<double> beta_;

        mutable double cachedMaxE_ = 0.0;
    };

} // namespace waveguide
