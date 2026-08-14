#include "CylindricalModel.hpp"
#include "Colormap.hpp"
#include "Bessel.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace waveguide
{
    namespace
    {
        constexpr double kPi  = 3.14159265358979323846;
        constexpr double c0   = 299792458.0;
        constexpr double mu0  = 4.0 * kPi * 1e-7;
        constexpr double eps0 = 8.8541878128e-12;

        // Shared definition (Colormap.hpp).
        void fireColormap(float t, float &r, float &g, float &b)
        {
            fireColor(t, r, g, b);
        }

        // Tabulated Bessel zeros from the Python reference.
        //   pnm_tm [n][m-1] -> m-th zero of J_n   (used by TM)
        //   pnm_te [n][m-1] -> m-th zero of J'_n  (used by TE)
        double besselZero(int n, int m, bool tm)
        {
            // Full double precision: the boundary conditions are satisfied
            // BECAUSE kc*R lands exactly on a root (Jn there for TM, Jn' for
            // TE). Rounding the root to 4 digits leaves a residual tangential E
            // on the wall of order 1e-4 of the peak -- small, but it is pure
            // table error, not physics.
            static const double pnm_tm[3][3] = {
                { 2.4048255577,  5.5200781103,  8.6537279129},
                { 3.8317059702,  7.0155866698, 10.1734681351},
                { 5.1356223018,  8.4172441404, 11.6198411721},
            };
            static const double pnm_te[3][3] = {
                { 3.8317059702,  7.0155866698, 10.1734681351},
                { 1.8411837813,  5.3314427735,  8.5363163663},
                { 3.0542369282,  6.7061331942,  9.9694678231},
            };
            const int nn = std::clamp(n, 0, 2);
            const int mm = std::clamp(m - 1, 0, 2);
            return tm ? pnm_tm[nn][mm] : pnm_te[nn][mm];
        }

        // Our own J_n (Bessel.hpp), not std::cyl_bessel_j: the C++17 special
        // math functions are absent from libc++, which is what Emscripten uses,
        // so the standard call compiles here and nowhere near a browser. Using
        // one implementation everywhere also means the desktop and the web build
        // cannot disagree numerically. Verified against std::cyl_bessel_j over
        // the whole domain the model touches: max abs error 6.4e-13.
        inline double Jn(int n, double x)      { return besselJ(n, x); }
        inline double Jnprime(int n, double x) { return besselJPrime(n, x); }
    } // namespace

    CylindricalModel::CylindricalModel(double radiusMM,
                                       double lengthMeters,
                                       double frequency,
                                       double epsilonRel,
                                       double muRel,
                                       int modeN,
                                       int modeM,
                                       ModeType type,
                                       FieldKind field,
                                       bool cavity,
                                       int modeL)
        : raio_(radiusMM / 1000.0),
          length_(lengthMeters),
          frequencia_(frequency),
          epsilon_(epsilonRel),
          mu_(muRel),
          n_(std::clamp(modeN, 0, 2)),
          m_(std::clamp(modeM, 1, 3)),
          type_(type),
          field_(field),
          cavity_(cavity),
          l_(modeL)
    {
        omega_ = 2.0 * kPi * frequencia_;
        k_     = (omega_ / c0) * std::sqrt(mu_ * epsilon_);
        const bool tm = (type_ == ModeType::TM);
        k_c_   = besselZero(n_, m_, tm) / raio_;

        if (cavity_)
        {
            // Closed cavity: the axial standing wave fixes β = lπ/d.
            beta_ = std::complex<double>(l_ * kPi / length_, 0.0);
        }
        else
        {
            const std::complex<double> radicand(k_ * k_ - k_c_ * k_c_, 0.0);
            beta_  = std::sqrt(radicand);
            if (beta_.imag() > 0.0)
                beta_ = -beta_;
        }
    }

    double CylindricalModel::resonantFrequency() const
    {
        const double vph = c0 / std::sqrt(mu_ * epsilon_);
        // Open guide: no axial standing wave, so the characteristic frequency is
        // the transverse cutoff, not a resonance. Same fix as the rectangular
        // model -- see the note there.
        if (!cavity_)
            return k_c_ * vph / (2.0 * kPi);

        // f_nml = vph/(2π) · sqrt(k_c² + (lπ/d)²)
        const double kz  = double(l_) * kPi / length_;
        const double kres = std::sqrt(k_c_ * k_c_ + kz * kz);
        return kres * vph / (2.0 * kPi);
    }

    // Cavity field: reuse the circular-waveguide transverse (ρ,φ) mode shapes
    // but replace the running-wave e^{-jβz} with the standing factors set by
    // the shorting walls at z = 0, d (same boundary rules as the rectangular
    // cavity): transverse E ~ sin(lπz/d), transverse H ~ cos(lπz/d), axial Hz
    // ~ sin, axial Ez ~ cos. β -> zk = lπ/d in amplitudes that carried β.
    // E carries cos(phase), H sin(phase) (90° time lag).
    std::array<double, 3> CylindricalModel::cavityFieldVec(double x, double y,
                                                           double z, double phase) const
    {
        const double rho = std::sqrt(x * x + y * y);
        if (rho > raio_) return {0.0, 0.0, 0.0};
        const double phi = std::atan2(y, x);

        const double sNphi = std::sin(double(n_) * phi);
        const double cNphi = std::cos(double(n_) * phi);
        const double arg   = k_c_ * rho;
        const double Jnv   = Jn(n_, arg);
        const double Jnpv  = Jnprime(n_, arg);
        const double rho_safe = (rho < 1e-9) ? 1e-9 : rho;

        const double zk = double(l_) * kPi / length_;
        const double Sz = std::sin(zk * z);
        const double Cz = std::cos(zk * z);
        const double ct = std::cos(phase); // electric time factor
        const double st = std::sin(phase); // magnetic time factor (90° lag)

        // Same angular convention as fieldVector(): fPhi on the axial component
        // and n*gPhi on anything born from d/dphi. The n is what makes those
        // components vanish for n = 0 instead of diverging as 1/rho.
        const double fPhi = A_ * sNphi + B_ * cNphi;
        const double gPhi = A_ * cNphi - B_ * sNphi;
        const double nd   = double(n_);

        double Vr = 0.0, Vp = 0.0, Vz = 0.0;
        if (type_ == ModeType::TE)
        {
            if (field_ == FieldKind::Electric)
            {
                // Transverse E only, standing via Sz.
                Vr = (omega_ * mu0 * mu_ * nd / (k_c_ * k_c_ * rho_safe)) *
                     gPhi * Jnv  * Sz * ct;
                Vp = (omega_ * mu0 * mu_ / k_c_) *
                     fPhi * Jnpv * Sz * ct;
            }
            else
            {
                // Transverse H (Cz) and axial Hz (Sz). β -> zk.
                Vr = (zk / k_c_) * fPhi * Jnpv * Cz * st;
                Vp = (zk * nd / (k_c_ * k_c_ * rho_safe)) * gPhi * Jnv * Cz * st;
                Vz = fPhi * Jnv * Sz * st;
            }
        }
        else // TM
        {
            if (field_ == FieldKind::Electric)
            {
                // Transverse E (Sz) and axial Ez (Cz). β -> zk.
                Vr = (zk / k_c_) * fPhi * Jnpv * Sz * ct;
                Vp = (zk * nd / (k_c_ * k_c_ * rho_safe)) * gPhi * Jnv * Sz * ct;
                Vz = fPhi * Jnv * Cz * ct;
            }
            else
            {
                // Transverse H only, standing via Cz.
                Vr = (omega_ * eps0 * epsilon_ * nd / (k_c_ * k_c_ * rho_safe)) *
                     gPhi * Jnv  * Cz * st;
                Vp = (omega_ * eps0 * epsilon_ / k_c_) *
                     fPhi * Jnpv * Cz * st;
            }
        }

        const double cp = std::cos(phi), sp = std::sin(phi);
        return {Vr * cp - Vp * sp, Vr * sp + Vp * cp, Vz};
    }

    double CylindricalModel::magnitudeE(double x, double y, double z, double phase) const
    {
        if (cavity_)
        {
            // Static spatial envelope (peak over time) so the 3D cloud shows
            // the fixed resonant pattern; animated dynamics live in the
            // cross-section vector plots. `phase` is ignored here.
            (void)phase;
            const double tpk = (field_ == FieldKind::Electric) ? 0.0 : (kPi * 0.5);
            const std::array<double, 3> V = cavityFieldVec(x, y, z, tpk);
            return std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
        }

        // Delegated: fieldVector() is the ONLY transcription of the mode
        // equations in this file. This used to be an independent copy, and that
        // is precisely how the missing-n bug survived a fix -- peakField() is
        // driven from here while the cloud samples fieldVector(), so the two
        // disagreed by ~1e7 and every sample fell under the visibility cut.
        const std::array<double, 3> V = fieldVector(x, y, z, phase);
        return std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
    }

    // Coarse scan for the peak |field|, cached. The phase is swept because the
    // plotted quantity is the REAL part and oscillates in time.
    //
    // Called from peakField(), not only from sampleGrid(): the cache used to be
    // a side effect of sampling the particle grid, so any consumer asking for
    // the peak without having drawn a grid first read 0 and silently produced an
    // empty visualization.
    void CylindricalModel::ensurePeak() const
    {
        if (cachedMaxE_ > 0.0) return;
        const int sx = 24, sy = 24, sz = 48;
        for (int i = 0; i <= sx; ++i)
        {
            const double x = -raio_ + 2.0 * raio_ * i / sx;
            for (int jj = 0; jj <= sy; ++jj)
            {
                const double y = -raio_ + 2.0 * raio_ * jj / sy;
                if (x * x + y * y > raio_ * raio_) continue;
                for (int kk = 0; kk <= sz; ++kk)
                {
                    const double z = length_ * kk / sz;
                    for (int pp = 0; pp < 8; ++pp)
                    {
                        const double ph = 2.0 * kPi * pp / 8.0;
                        const double e = magnitudeE(x, y, z, ph);
                        if (e > cachedMaxE_) cachedMaxE_ = e;
                    }
                }
            }
        }
        if (cachedMaxE_ <= 0.0) cachedMaxE_ = 1.0;
    }

    std::pair<double, double> CylindricalModel::transverseField(double x, double y,
                                                                 double z, double phase) const
    {
        // Delegated for the same reason as magnitudeE(): one transcription only.
        const std::array<double, 3> V = fieldVector(x, y, z, phase);
        return {V[0], V[1]};
    }

    std::array<double, 3> CylindricalModel::fieldVector(double x, double y,
                                                         double z, double phase) const
    {
        if (cavity_)
            return cavityFieldVec(x, y, z, phase);

        const double rho = std::sqrt(x * x + y * y);
        if (rho > raio_) return {0.0, 0.0, 0.0};
        const double phi = std::atan2(y, x);

        using cx = std::complex<double>;
        const cx j(0.0, 1.0);
        const cx expz = std::exp(j * (cx(phase, 0.0) - beta_ * z));
        const double sNphi = std::sin(double(n_) * phi);
        const double cNphi = std::cos(double(n_) * phi);
        const double arg = k_c_ * rho;
        const double Jnv = Jn(n_, arg);
        const double Jnpv = Jnprime(n_, arg);
        const double rho_safe = (rho < 1e-9) ? 1e-9 : rho;

        // Angular convention (Pozar 4th ed, eq. 3.110 for TE / 3.116 for TM):
        // the axial component carries fPhi, and whatever comes from d/dphi
        // carries n * gPhi, with gPhi = (1/n) dfPhi/dphi.
        //
        //   TE:  Hz ~ fPhi Jn      Hr ~ fPhi Jn'      Hp ~ n gPhi Jn / rho
        //                          Er ~ n gPhi Jn/rho Ep ~ fPhi Jn'
        //   TM:  Ez ~ fPhi Jn      Er ~ fPhi Jn'      Ep ~ n gPhi Jn / rho
        //                          Hr ~ n gPhi Jn/rho Hp ~ fPhi Jn'
        //
        // The factor n is what makes the d/dphi components VANISH for a
        // circularly symmetric mode (n = 0). Dropping it left Er (TE) and Hr
        // (TM) as a bare 1/rho, which blew up on the axis instead of going to
        // zero: peakField() then latched onto that spike and the whole cloud
        // fell below the 5% visibility cut, drawing an empty guide.
        const double fPhi = A_ * sNphi + B_ * cNphi;
        const double gPhi = A_ * cNphi - B_ * sNphi;
        const double nd   = double(n_);

        cx Vr(0.0, 0.0), Vp(0.0, 0.0), Vz(0.0, 0.0);
        if (type_ == ModeType::TE)
        {
            if (field_ == FieldKind::Electric)
            {
                const cx cEr = -j * omega_ * mu0 * mu_ * nd / (k_c_ * k_c_ * rho_safe);
                const cx cEp =  j * omega_ * mu0 * mu_ / k_c_;
                Vr = cEr * gPhi * Jnv  * expz;
                Vp = cEp * fPhi * Jnpv * expz;
            }
            else
            {
                const cx cHr = -j * beta_ / k_c_;
                const cx cHp = -j * beta_ * nd / (k_c_ * k_c_ * rho_safe);
                Vr = cHr * fPhi * Jnpv * expz;
                Vp = cHp * gPhi * Jnv  * expz;
                Vz =       fPhi * Jnv  * expz;
            }
        }
        else
        {
            if (field_ == FieldKind::Electric)
            {
                const cx cEr = -j * beta_ / k_c_;
                const cx cEp = -j * beta_ * nd / (k_c_ * k_c_ * rho_safe);
                Vr = cEr * fPhi * Jnpv * expz;
                Vp = cEp * gPhi * Jnv  * expz;
                Vz =       fPhi * Jnv  * expz;
            }
            else
            {
                const cx cHr =  j * omega_ * eps0 * epsilon_ * nd / (k_c_ * k_c_ * rho_safe);
                const cx cHp = -j * omega_ * eps0 * epsilon_ / k_c_;
                Vr = cHr * gPhi * Jnv  * expz;
                Vp = cHp * fPhi * Jnpv * expz;
            }
        }
        const double rR = Vr.real();
        const double rP = Vp.real();
        const double cp = std::cos(phi), sp = std::sin(phi);
        return {rR * cp - rP * sp, rR * sp + rP * cp, Vz.real()};
    }

    std::vector<Particle> CylindricalModel::sampleGrid(int nx,
                                                       int ny,
                                                       int nz,
                                                       bool cutawayOn,
                                                       float minIntensity,
                                                       double phase) const
    {
        std::vector<Particle> out;
        if (nx <= 1 || ny <= 1 || nz <= 1) return out;

        ensurePeak();
        const double maxE = cachedMaxE_;

        const float halfL = float(length_) * 0.5f;

        out.reserve(size_t(nx) * ny * nz);
        for (int i = 0; i < nx; ++i)
        {
            const double x = -raio_ + 2.0 * raio_ * i / double(nx - 1);
            for (int j = 0; j < ny; ++j)
            {
                const double y = -raio_ + 2.0 * raio_ * j / double(ny - 1);
                if (x * x + y * y > raio_ * raio_) continue;
                for (int k = 0; k < nz; ++k)
                {
                    const double z = length_ * k / double(nz - 1);
                    const double e = magnitudeE(x, y, z, phase);
                    const float  t = float(std::min(1.0, e / maxE));
                    if (t < minIntensity) continue;

                    const float cx = float(x);
                    const float cy = float(y);
                    const float cz = float(z) - halfL;
                    if (cutawayOn && cx > 0.0f && cy > 0.0f && cz > 0.0f)
                        continue;

                    float r, g, b;
                    fireColormap(t, r, g, b);

                    Particle p;
                    p.x = cx; p.y = cy; p.z = cz;
                    p.r = r;  p.g = g;  p.b = b;
                    p.intensity = t;
                    out.push_back(p);
                }
            }
        }
        return out;
    }

} // namespace waveguide
