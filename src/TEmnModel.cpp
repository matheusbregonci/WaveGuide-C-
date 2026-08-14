#include "TEmnModel.hpp"
#include "Colormap.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

namespace waveguide
{

    namespace
    {
        constexpr double kPi  = 3.14159265358979323846;
        constexpr double c0   = 299792458.0;
        constexpr double mu0  = 4.0 * kPi * 1e-7;
        constexpr double eps0 = 8.8541878128e-12;

        // Shared definition (Colormap.hpp). Was a private black->purple->white
        // ramp, which is not what any of the plots use any more.
        void fireColormap(float t, float &r, float &g, float &b)
        {
            fireColor(t, r, g, b);
        }
    } // namespace

    TEmnModel::TEmnModel(double widthMM,
                         double heightMM,
                         double frequency,
                         double epsilonRel,
                         double muRel,
                         int modeM,
                         int modeN,
                         ModeType type,
                         FieldKind field,
                         double depthMM,
                         bool cavity,
                         int modeL,
                         double powerW)
        : largura_(widthMM / 1000.0),
          altura_(heightMM / 1000.0),
          profundidade_(depthMM / 1000.0),
          frequencia_(frequency),
          epsilon_(epsilonRel),
          mu_(muRel),
          m_(modeM),
          n_(modeN),
          type_(type),
          field_(field),
          cavity_(cavity),
          l_(modeL),
          powerW_(powerW)
    {
        omega_ = 2.0 * kPi * frequencia_;
        k_ = (omega_ / c0) * std::sqrt(mu_ * epsilon_);

        const double kx = m_ * kPi / largura_;
        const double ky = n_ * kPi / altura_;
        k_c_sq_ = kx * kx + ky * ky;

        if (cavity_)
        {
            // Closed cavity: the axial standing wave fixes the propagation
            // constant at β = lπ/d regardless of the drive frequency.
            beta_ = std::complex<double>(l_ * kPi / profundidade_, 0.0);
        }
        else
        {
            const std::complex<double> radicand(k_ * k_ - k_c_sq_, 0.0);
            beta_ = std::sqrt(radicand);
            // Evanescent modes (f < fc): ensure field decays along +z.
            // std::sqrt of negative real gives +j*alpha; we need -j*alpha
            // so that e^{-j*beta*z} = e^{-alpha*z} (decaying, not growing).
            if (beta_.imag() > 0.0)
                beta_ = -beta_;
        }

        solveAmplitudeFromPower();
    }

    // Pin the mode amplitude A_ to the requested transported power, turning the
    // eigenfunction into a field in real V/m and A/m.
    //
    // A propagating mode carries
    //     P = (1/2) Re integral (E x H*).z dA = (1/(2 Zw)) integral (|Ex|^2+|Ey|^2) dA
    // because the transverse E and H of a single mode are related by the wave
    // impedance Zw everywhere (Ex = Zw Hy, Ey = -Zw Hx). Both transverse E
    // components are A times a product of sin/cos, so the area integral is just a
    // product of the elementary integrals below and P comes out proportional to
    // A^2 -- invert that for A.
    //
    // Cross-checks against Pozar for TE10 (m=1, n=0):
    //     P = |A|^2 omega mu beta a^3 b / (4 pi^2)                  (eq. 3.86)
    //
    // Only meaningful for a propagating waveguide mode. A cavity stores energy
    // instead of transporting it, and an evanescent mode transports none (beta is
    // purely imaginary, so E x H* is imaginary and Re{} vanishes). Both keep the
    // arbitrary A_ = 1 and report physicalUnits() == false.
    void TEmnModel::solveAmplitudeFromPower()
    {
        A_ = 1.0;
        unitsPhysical_ = false;

        if (cavity_ || !(powerW_ > 0.0))
            return;
        const double betaR = beta_.real();
        if (!(betaR > 0.0) || std::fabs(beta_.imag()) > 1e-9 * std::fabs(betaR))
            return;              // at or below cutoff: no real power flow
        if (!(k_c_sq_ > 0.0))
            return;              // TEM placeholder: no transverse variation

        const double a = largura_, b = altura_;
        const double mpi_a = double(m_) * kPi / a;
        const double npi_b = double(n_) * kPi / b;

        // Elementary integrals over the cross-section. cos^2 integrates to the
        // full side when the index is 0 (the integrand is 1), to half otherwise;
        // sin^2 vanishes when the index is 0.
        const double Cm = (m_ == 0) ? a : a * 0.5;
        const double Sm = (m_ == 0) ? 0.0 : a * 0.5;
        const double Cn = (n_ == 0) ? b : b * 0.5;
        const double Sn = (n_ == 0) ? 0.0 : b * 0.5;

        const double muAbs  = mu0 * mu_;
        const double epsAbs = eps0 * epsilon_;

        double Zw = 0.0, integ = 0.0;
        if (type_ == ModeType::TE)
        {
            // Ex ~ (omega mu npi_b / kc^2) cos(mx) sin(ny)   -> Cm * Sn
            // Ey ~ (omega mu mpi_a / kc^2) sin(mx) cos(ny)   -> Sm * Cn
            Zw = omega_ * muAbs / betaR;
            const double pref = omega_ * muAbs / k_c_sq_;
            integ = pref * pref * (npi_b * npi_b * Cm * Sn +
                                   mpi_a * mpi_a * Sm * Cn);
        }
        else
        {
            // Ex ~ (beta mpi_a / kc^2) cos(mx) sin(ny)       -> Cm * Sn
            // Ey ~ (beta npi_b / kc^2) sin(mx) cos(ny)       -> Sm * Cn
            Zw = betaR / (omega_ * epsAbs);
            const double pref = betaR / k_c_sq_;
            integ = pref * pref * (mpi_a * mpi_a * Cm * Sn +
                                   npi_b * npi_b * Sm * Cn);
        }

        if (!(integ > 0.0) || !(Zw > 0.0))
            return;

        A_ = std::sqrt(2.0 * Zw * powerW_ / integ);
        unitsPhysical_ = true;
        cachedMaxE_ = 0.0;   // amplitude changed: force the peak scan to re-run
    }

    double TEmnModel::resonantFrequency() const
    {
        const double vph = c0 / std::sqrt(mu_ * epsilon_);
        // An OPEN guide has no axial standing wave, so there is no l and no
        // resonance: the characteristic frequency is the transverse cutoff.
        // Adding (l pi/d)^2 unconditionally -- as this did -- inflated the
        // reported f_c by a few tenths of a percent for every open guide, and
        // the exported report printed that as the cutoff.
        if (!cavity_)
            return std::sqrt(k_c_sq_) * vph / (2.0 * kPi);

        // f_mnl = vph/(2π) · sqrt((mπ/a)² + (nπ/b)² + (lπ/d)²)
        const double kz  = double(l_) * kPi / profundidade_;
        const double kres = std::sqrt(k_c_sq_ + kz * kz);
        return kres * vph / (2.0 * kPi);
    }

    // Cavity field: reuse the waveguide transverse (x,y) mode shapes but
    // replace the running-wave e^{-jβz} with the standing factors demanded
    // by the shorting walls at z = 0, d. Boundary conditions:
    //   tangential E = 0 at z = 0,d  -> transverse E  ~ sin(lπz/d)
    //   normal   H = 0 at z = 0,d     -> axial Hz      ~ sin(lπz/d)
    //   transverse H is a maximum      -> transverse H  ~ cos(lπz/d)
    //   (TM dual: axial Ez ~ cos(lπz/d), transverse E ~ sin)
    // E and H are 90° apart in time, so E carries cos(phase) and H sin(phase).
    std::array<double, 3> TEmnModel::cavityFieldVec(double x, double y,
                                                    double z, double phase) const
    {
        const double mx = m_ * kPi * x / largura_;
        const double ny = n_ * kPi * y / altura_;
        const double cos_mx = std::cos(mx), sin_mx = std::sin(mx);
        const double cos_ny = std::cos(ny), sin_ny = std::sin(ny);

        const double zk = double(l_) * kPi / profundidade_;
        const double Sz = std::sin(zk * z); // vanishes on the end walls
        const double Cz = std::cos(zk * z);

        const double denom = (k_c_sq_ > 0.0) ? k_c_sq_ : 1.0;
        const double mpi_a = double(m_) * kPi / largura_;
        const double npi_b = double(n_) * kPi / altura_;

        const double ct = std::cos(phase); // electric time factor
        const double st = std::sin(phase); // magnetic time factor (90° lag)

        double Vx = 0.0, Vy = 0.0, Vz = 0.0;
        if (type_ == ModeType::TE)
        {
            if (field_ == FieldKind::Electric)
            {
                // Transverse E only (no Ez for TE). Standing in z via Sz.
                // Ex carries sin(ny): it is tangential to the y = 0 and y = b
                // walls and must vanish there.
                Vx =  omega_ * mu0 * mu_ * npi_b / denom * A_ * cos_mx * sin_ny * Sz * ct;
                Vy = -omega_ * mu0 * mu_ * mpi_a / denom * A_ * sin_mx * cos_ny * Sz * ct;
            }
            else
            {
                // Transverse H (Cz) and axial Hz (Sz). β -> zk in amplitude.
                Vx = zk * mpi_a / denom * A_ * sin_mx * cos_ny * Cz * st;
                Vy = zk * npi_b / denom * A_ * cos_mx * sin_ny * Cz * st;
                Vz = A_ * cos_mx * cos_ny * Sz * st;
            }
        }
        else // TM
        {
            if (field_ == FieldKind::Electric)
            {
                // Transverse E (Sz) and axial Ez (Cz). β -> zk in amplitude.
                Vx = zk * mpi_a / denom * A_ * cos_mx * sin_ny * Sz * ct;
                Vy = zk * npi_b / denom * A_ * sin_mx * cos_ny * Sz * ct;
                Vz = A_ * sin_mx * sin_ny * Cz * ct;
            }
            else
            {
                // Transverse H only (no Hz for TM). Standing in z via Cz.
                // Hy carries sin(ny), mirroring Hx's sin(mx): both follow from
                // H_t = (j omega eps / kc^2) z_hat x grad(Ez) with
                // Ez ~ sin(mx) sin(ny).
                Vx =  omega_ * eps0 * epsilon_ * npi_b / denom * A_ * sin_mx * cos_ny * Cz * st;
                Vy = -omega_ * eps0 * epsilon_ * mpi_a / denom * A_ * cos_mx * sin_ny * Cz * st;
            }
        }
        return {Vx, Vy, Vz};
    }

    // Instantaneous magnitude |E|(x,y,z) at t = 0 using the REAL part of
    // each component (matches TEmn_model.py's np.real convention). This
    // preserves the standing-wave pattern along z — you see the actual
    // oscillation of the propagating mode inside the guide.
    double TEmnModel::magnitudeE(double x, double y, double z, double phase) const
    {
        if (cavity_)
        {
            // Time-independent spatial envelope: evaluate the field at the
            // instant its own temporal factor peaks (E at cos=1, H at sin=1)
            // so the 3D cloud shows the fixed resonant pattern instead of
            // blinking as the standing wave breathes. `phase` is ignored here;
            // the animated dynamics live in the cross-section vector plots.
            (void)phase;
            const double tpk = (field_ == FieldKind::Electric) ? 0.0 : (kPi * 0.5);
            const std::array<double, 3> V = cavityFieldVec(x, y, z, tpk);
            return std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
        }

        // Propagating guide: same single transcription of the mode equations
        // that fieldVector() holds.
        const std::array<double, 3> V = fieldVector(x, y, z, phase);
        return std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
    }

    // Transverse pair, delegated to fieldVector so there is exactly ONE copy of
    // the mode equations in this file. The three accessors used to carry three
    // independent transcriptions of the same formulas, which is how the Ex/Hy
    // sin-vs-cos errors survived: fixing one copy left the other two wrong.
    // Coarse scan for the peak |field|, cached. The phase is swept too: the
    // plotted quantity is the REAL part, which oscillates in time, so a scan at
    // t = 0 alone would miss the maximum of any component in quadrature.
    //
    // Called from peakField() rather than only from sampleGrid(): the cache used
    // to be filled as a side effect of sampling the particle grid, so any
    // consumer that asked for the peak WITHOUT having drawn a grid first got 0
    // and silently produced an empty visualization.
    void TEmnModel::ensurePeak() const
    {
        if (cachedMaxE_ > 0.0) return;
        const int sx = 32, sy = 32, sz = 64;
        for (int i = 0; i <= sx; ++i)
        {
            const double x = largura_ * i / sx;
            for (int jj = 0; jj <= sy; ++jj)
            {
                const double y = altura_ * jj / sy;
                for (int kk = 0; kk <= sz; ++kk)
                {
                    const double z = profundidade_ * kk / sz;
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

    std::pair<double, double> TEmnModel::transverseField(double x, double y,
                                                          double z, double phase) const
    {
        const std::array<double, 3> V = fieldVector(x, y, z, phase);
        return {V[0], V[1]};
    }

    std::array<double, 3> TEmnModel::fieldVector(double x, double y,
                                                  double z, double phase) const
    {
        if (cavity_)
            return cavityFieldVec(x, y, z, phase);

        using cx = std::complex<double>;
        const cx j(0.0, 1.0);
        const double mx = m_ * kPi * x / largura_;
        const double ny = n_ * kPi * y / altura_;
        const double cos_mx = std::cos(mx), cos_ny = std::cos(ny);
        const double sin_mx = std::sin(mx), sin_ny = std::sin(ny);
        const cx expz = std::exp(j * (cx(phase, 0.0) - beta_ * z));
        const double denom = (k_c_sq_ > 0.0) ? k_c_sq_ : 1.0;
        const double mpi_a = double(m_) * kPi / largura_;
        const double npi_b = double(n_) * kPi / altura_;

        // Field equations follow Pozar (Microwave Engineering, 4th ed, sec 3.3).
        // Prefactors use absolute mu = mu0*mur and eps = eps0*epsr (SI units).
        // Every transverse component pairs sin against cos in each coordinate:
        // that is what makes the tangential E vanish on all four walls. A stray
        // cos where a sin belongs puts a tangential E maximum ON the conductor.
        cx Vx(0.0, 0.0), Vy(0.0, 0.0), Vz(0.0, 0.0);
        if (type_ == ModeType::TE)
        {
            // Hz = A cos(mx) cos(ny) e^{-j beta z}
            if (field_ == FieldKind::Electric)
            {
                // Ex =  j omega mu (npi/b) / kc^2 * A cos(mx) sin(ny)
                // Ey = -j omega mu (mpi/a) / kc^2 * A sin(mx) cos(ny)
                Vx =  j * omega_ * mu0 * mu_ * npi_b / denom * A_ * cos_mx * sin_ny * expz;
                Vy = -j * omega_ * mu0 * mu_ * mpi_a / denom * A_ * sin_mx * cos_ny * expz;
            }
            else
            {
                // Hx = j beta (mpi/a) / kc^2 * A sin(mx) cos(ny)
                // Hy = j beta (npi/b) / kc^2 * A cos(mx) sin(ny)
                Vx = j * beta_ * mpi_a / denom * A_ * sin_mx * cos_ny * expz;
                Vy = j * beta_ * npi_b / denom * A_ * cos_mx * sin_ny * expz;
                Vz = A_ * cos_mx * cos_ny * expz;
            }
        }
        else
        {
            // Ez = A sin(mx) sin(ny) e^{-j beta z}
            if (field_ == FieldKind::Electric)
            {
                // Ex = -j beta (mpi/a) / kc^2 * A cos(mx) sin(ny)
                // Ey = -j beta (npi/b) / kc^2 * A sin(mx) cos(ny)
                Vx = -j * beta_ * mpi_a / denom * A_ * cos_mx * sin_ny * expz;
                Vy = -j * beta_ * npi_b / denom * A_ * sin_mx * cos_ny * expz;
                Vz =  A_ * sin_mx * sin_ny * expz;
            }
            else
            {
                // Hx =  j omega eps (npi/b) / kc^2 * A sin(mx) cos(ny)
                // Hy = -j omega eps (mpi/a) / kc^2 * A cos(mx) sin(ny)
                // The j matters: it puts H in quadrature with Ez, as Maxwell
                // requires. Dropping it (a real prefactor) does not change the
                // shape of an H-only plot but is wrong in absolute phase, and
                // it broke E x H* in the power integral.
                Vx =  j * omega_ * eps0 * epsilon_ * npi_b / denom * A_ * sin_mx * cos_ny * expz;
                Vy = -j * omega_ * eps0 * epsilon_ * mpi_a / denom * A_ * cos_mx * sin_ny * expz;
            }
        }
        return {Vx.real(), Vy.real(), Vz.real()};
    }

    std::vector<Particle> TEmnModel::sampleProbabilistic(int count,
                                                         bool cutawayOn,
                                                         uint32_t seed) const
    {
        std::vector<Particle> out;
        if (count <= 0)
            return out;
        out.reserve(size_t(count));

        // --------- estimate max of |E|^2 over the full 3D box ---------
        // The real-part magnitude carries a sin(beta*z)/cos(beta*z) factor,
        // so a scan at z = 0 alone may miss the peak entirely (sin(0) = 0).
        // A 3D scan is cheap (tens of thousands of evaluations) and robust
        // for any TEmn mode at any guide length.
        double maxE2 = 0.0;
        const int scanXY = 48;
        const int scanZ = 96;
        for (int i = 0; i <= scanXY; ++i)
        {
            const double x = largura_ * i / scanXY;
            for (int jj = 0; jj <= scanXY; ++jj)
            {
                const double y = altura_ * jj / scanXY;
                for (int kk = 0; kk <= scanZ; ++kk)
                {
                    const double z = profundidade_ * kk / scanZ;
                    const double e = magnitudeE(x, y, z);
                    const double e2 = e * e;
                    if (e2 > maxE2)
                        maxE2 = e2;
                }
            }
        }
        if (maxE2 <= 0.0)
            maxE2 = 1.0;
        const double maxE2Safe = maxE2 * 1.05;

        // --------- RNG ---------
        std::mt19937 rng;
        if (seed == 0)
        {
            std::random_device rd;
            rng.seed(rd());
        }
        else
        {
            rng.seed(seed);
        }
        std::uniform_real_distribution<double> uni(0.0, 1.0);

        // --------- rejection sampling ---------
        const float halfW = float(largura_) * 0.5f;
        const float halfH = float(altura_) * 0.5f;
        const float halfL = float(profundidade_) * 0.5f;

        // Cutaway: remove the upper-far corner (x > 0 && y > 0 && z > 0 in
        // centered coordinates). This matches the Atoms-style quarter cut.
        const auto inCutaway = [](float cx, float cy, float cz)
        {
            return (cx > 0.0f) && (cy > 0.0f) && (cz > 0.0f);
        };

        // maxE over the scan, used to normalize colors to [0,1].
        const double maxE = std::sqrt(maxE2);

        int accepted = 0;
        int attempts = 0;
        const int attemptCap = count * 60; // safety cap

        while (accepted < count && attempts < attemptCap)
        {
            ++attempts;
            const double x = uni(rng) * largura_;
            const double y = uni(rng) * altura_;
            const double z = uni(rng) * profundidade_;

            const double e = magnitudeE(x, y, z);
            const double e2 = e * e;
            if (uni(rng) > (e2 / maxE2Safe))
                continue;

            const float cx = float(x) - halfW;
            const float cy = float(y) - halfH;
            const float cz = float(z) - halfL;
            if (cutawayOn && inCutaway(cx, cy, cz))
                continue;

            float r, g, b;
            fireColormap(float(e / maxE), r, g, b);

            Particle p;
            p.x = cx;
            p.y = cy;
            p.z = cz;
            p.r = r;
            p.g = g;
            p.b = b;
            p.intensity = float(e / maxE);
            out.push_back(p);
            ++accepted;
        }

        return out;
    }

    std::vector<Particle> TEmnModel::sampleGrid(int nx,
                                                int ny,
                                                int nz,
                                                bool cutawayOn,
                                                float minIntensity,
                                                double phase) const
    {
        std::vector<Particle> out;
        if (nx <= 1 || ny <= 1 || nz <= 1)
            return out;

        ensurePeak();
        const double maxE = cachedMaxE_;

        const float halfW = float(largura_) * 0.5f;
        const float halfH = float(altura_) * 0.5f;
        const float halfL = float(profundidade_) * 0.5f;

        out.reserve(size_t(nx) * ny * nz);
        for (int i = 0; i < nx; ++i)
        {
            const double x = largura_ * i / double(nx - 1);
            for (int j = 0; j < ny; ++j)
            {
                const double y = altura_ * j / double(ny - 1);
                for (int k = 0; k < nz; ++k)
                {
                    const double z = profundidade_ * k / double(nz - 1);
                    const double e = magnitudeE(x, y, z, phase);
                    const float t = float(std::min(1.0, e / maxE));
                    if (t < minIntensity)
                        continue;

                    const float cx = float(x) - halfW;
                    const float cy = float(y) - halfH;
                    const float cz = float(z) - halfL;
                    if (cutawayOn && cx > 0.0f && cy > 0.0f && cz > 0.0f)
                        continue;

                    float r, g, b;
                    fireColormap(t, r, g, b);

                    Particle p;
                    p.x = cx;
                    p.y = cy;
                    p.z = cz;
                    p.r = r;
                    p.g = g;
                    p.b = b;
                    p.intensity = t;
                    out.push_back(p);
                }
            }
        }
        return out;
    }

} // namespace waveguide
