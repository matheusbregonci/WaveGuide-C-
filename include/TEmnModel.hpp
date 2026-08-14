#pragma once

#include "Types.hpp"
#include "FieldSource.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// C++ port of TEmn_model.py (WaveGuidesSimulator).
//
// Solves the TEmn mode fields of a rectangular metallic waveguide and
// produces particle clouds suitable for point-cloud visualization in the
// style of kavan010/Atoms:
//   - stochastic sampling of positions with density proportional to |E|^2
//   - per-particle color taken from a "fire" heatmap applied to |E|
//   - optional interactive cutaway to expose the interior of the guide
//
// Units:
//   - width (a), height (b), depth (L)     : meters
//   - frequency                              : Hertz
//   - relative permittivity / permeability   : dimensionless

namespace waveguide
{

    // Particle, ModeType, FieldKind and Bounds now live in Types.hpp.

    class TEmnModel : public FieldSource
    {
    public:
        // Defaults match the Python reference:
        //   width  = 22.86 mm, height = 10.16 mm
        //   f      = 12 GHz,   depth  = 0.11 m
        //   mu_r = eps_r = 1,  m = 1, n = 0   (TE10 mode)
        //   cavity : if true, the guide is shorted at z = 0 and z = depth,
        //            producing a resonant standing wave along z (TE_mnl /
        //            TM_mnl cavity mode). `modeL` is the number of half-wave
        //            variations along the depth. When false the model is the
        //            open propagating waveguide (unchanged behavior).
        //   powerW : time-average power the propagating mode carries down the
        //            guide, in watts. The mode amplitude A_ is solved from it so
        //            the returned E/H are in real V/m and A/m -- this is the same
        //            convention HFSS uses for a wave port (1 W incident by
        //            default), which is what makes the two comparable. Ignored
        //            for a cavity or an evanescent mode (see physicalUnits()).
        TEmnModel(double widthMM = 22.86,
                  double heightMM = 10.16,
                  double frequency = 12e9,
                  double epsilonRel = 1.0,
                  double muRel = 1.0,
                  int modeM = 3,
                  int modeN = 3,
                  ModeType type = ModeType::TE,
                  FieldKind field = FieldKind::Electric,
                  double depthMM = 300.0,
                  bool cavity = false,
                  int modeL = 1,
                  double powerW = 1.0);

        double peakField() const override { ensurePeak(); return cachedMaxE_; }
        FieldKind fieldKind() const override { return field_; }

        // True once A_ has been solved from powerW, i.e. the field values are
        // real V/m and A/m. False for a cavity (stores energy, carries no net
        // power) and below cutoff (an evanescent mode carries none either): the
        // amplitude then stays at the arbitrary A_ = 1 and the UI must say so.
        bool physicalUnits() const override { return unitsPhysical_; }

        // Mode amplitude actually in use (Hz amplitude for TE, Ez for TM).
        double amplitude() const { return A_; }
        double powerW() const { return powerW_; }

        bool isCavity() const { return cavity_; }
        int modeL() const { return l_; }

        // Resonant frequency f_mnl of the closed cavity (Hz). For the open
        // waveguide this returns the transverse cutoff frequency (l = 0).
        double resonantFrequency() const override;

        // Transverse cutoff wavenumber k_c [rad/m].
        double cutoffWavenumber() const override { return std::sqrt(k_c_sq_); }

        // Rectangular cross-section: inside the [0,a] x [0,b] x [0,depth] box.
        bool inside(double x, double y, double z) const override
        {
            return x >= 0.0 && x <= largura_ && y >= 0.0 && y <= altura_ &&
                   z >= 0.0 && z <= profundidade_;
        }

        // Stochastic particle cloud: rejection sampling with density ~ |E|^2.
        //   count       : how many particles to emit
        //   cutawayOn   : if true, removes the upper-far corner of the box so
        //                 the interior is visible (same idea as Atoms)
        //   seed        : RNG seed (0 = random)
        std::vector<Particle> sampleProbabilistic(int count,
                                                  bool cutawayOn = true,
                                                  uint32_t seed = 0) const;

        // Regular 3D grid sampling: evaluates |E| on a uniform lattice and
        // returns one particle per cell with normalized intensity in [0,1].
        // Cells below `minIntensity` are dropped so the empty background
        // stays empty. This is the method that actually reveals the mode
        // oscillation pattern — size-modulation in the shader makes nodes
        // shrink and antinodes swell.
        std::vector<Particle> sampleGrid(int nx,
                                         int ny,
                                         int nz,
                                         bool cutawayOn = true,
                                         float minIntensity = 0.05f,
                                         double phase = 0.0) const override;

        Bounds bounds() const override
        {
            return {float(largura_), float(altura_), float(profundidade_)};
        }

        // Physics accessors.
        double cutoffKcSquared() const
        {
            return k_c_sq_;
        }
        std::complex<double> beta() const
        {
            return beta_;
        }
        double omega() const
        {
            return omega_;
        }
        double epsilonRel() const override { return epsilon_; }
        double muRel() const override { return mu_; }

    private:
        // Parameters
        double largura_;
        double altura_;
        double profundidade_;
        double frequencia_;
        double epsilon_;
        double mu_;
        int m_;
        int n_;
        ModeType type_;
        FieldKind field_;
        double A_ = 1.0;
        bool cavity_ = false;
        int l_ = 1; // z half-wave index (cavity only)
        double powerW_ = 1.0;          // requested transported power
        bool unitsPhysical_ = false;   // A_ solved from powerW_ (=> SI units)

        // Derived physics
        double omega_;
        double k_;
        double k_c_sq_;             // (mπ/a)^2 + (nπ/b)^2
        std::complex<double> beta_; // sqrt(k^2 - k_c^2)

        mutable double cachedMaxE_ = 0.0;

        // Evaluate real-valued |E|(x,y,z,t) in SI units. The time axis
        // is carried as a unitless phase angle (radians) so the caller
        // can advance it at any visualization speed.
        double magnitudeE(double x, double y, double z, double phase = 0.0) const;

        // Solve A_ from powerW_ (see the definition for the derivation).
        void solveAmplitudeFromPower();

        // Fill cachedMaxE_ on first use. const + mutable cache so peakField()
        // can stay a plain accessor from the caller's point of view.
        void ensurePeak() const;

        // Real Cartesian field vector of the CLOSED cavity at (x,y,z) and
        // time `phase`. Transverse components carry sin(lπz/d), the axial
        // component carries cos/sin(lπz/d) per the shorting-wall boundary
        // conditions; E and H oscillate 90° apart in time. Chosen field
        // (E or H) follows field_. Used by all three public accessors when
        // cavity_ is set.
        std::array<double, 3> cavityFieldVec(double x, double y, double z,
                                             double phase) const;

    public:
        // Real transverse vector components (Vx, Vy) in Cartesian — either
        // (Ex,Ey) or (Hx,Hy) depending on `field_`. Used to draw the
        // cross-section vector plot.
        std::pair<double, double> transverseField(double x, double y, double z,
                                                  double phase = 0.0) const override;

        // Real 3D Cartesian components of the chosen field.
        std::array<double, 3> fieldVector(double x, double y, double z,
                                          double phase = 0.0) const override;
    };

} // namespace waveguide
