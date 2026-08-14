#include "NumericalModel.hpp"
#include "Colormap.hpp"

#include <algorithm>
#include <cmath>

namespace waveguide
{
    namespace
    {
        constexpr double kC0 = 299792458.0;
        constexpr double kPi = 3.14159265358979323846;

        // Shared definition (Colormap.hpp).
        void fireColormap(float t, float& r, float& g, float& b)
        {
            fireColor(t, r, g, b);
        }
    } // namespace

    void NumericalModel::initGrid(int nx, int ny, int nz, const Aabb& box,
                                  std::vector<std::uint8_t> solid)
    {
        nx_ = nx; ny_ = ny; nz_ = nz;
        W_ = box.sizeX(); H_ = box.sizeY(); D_ = box.sizeZ();
        dx_ = W_ / std::max(1, nx); dy_ = H_ / std::max(1, ny); dz_ = D_ / std::max(1, nz);
        solid_ = std::move(solid);
    }

    void NumericalModel::addPeak(Mode& m)
    {
        double pk = 1e-30;
        for (std::size_t g = 0; g < m.gx.size(); ++g)
        {
            const double v = std::sqrt(m.gx[g]*m.gx[g] + m.gy[g]*m.gy[g] + m.gz[g]*m.gz[g]);
            if (v > pk) pk = v;
        }
        m.peak = pk;
    }

    NumericalModel::NumericalModel(const HelmholtzResult& r)
    {
        vector_ = false;
        initGrid(r.nx, r.ny, r.nz, r.box, r.solid);
        const int nx = nx_, ny = ny_, nz = nz_;
        const double ix = 1.0 / (2.0 * dx_), iy = 1.0 / (2.0 * dy_), iz = 1.0 / (2.0 * dz_);
        for (const EigenMode& em : r.modes)
        {
            Mode m; m.lambda = em.lambda;
            m.gx.assign(em.psi.size(), 0.0); m.gy = m.gx; m.gz = m.gx;
            auto P = [&](int i, int j, int k) -> double {
                if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz) return 0.0;
                return em.psi[gidx(i, j, k)]; };
            for (int k = 0; k < nz; ++k) for (int j = 0; j < ny; ++j) for (int i = 0; i < nx; ++i)
            {
                const std::size_t g = gidx(i, j, k);
                m.gx[g] = (P(i+1,j,k) - P(i-1,j,k)) * ix;
                m.gy[g] = (P(i,j+1,k) - P(i,j-1,k)) * iy;
                m.gz[g] = (P(i,j,k+1) - P(i,j,k-1)) * iz;
            }
            addPeak(m);
            modes_.push_back(std::move(m));
        }
        mode_ = 0;
    }

    NumericalModel::NumericalModel(const MaxwellResult& r)
    {
        vector_ = true;
        initGrid(r.nx, r.ny, r.nz, r.box, r.solid);
        for (const MaxwellMode& em : r.modes)
        {
            Mode m; m.lambda = em.lambda;
            m.gx = em.ex; m.gy = em.ey; m.gz = em.ez;
            addPeak(m);
            modes_.push_back(std::move(m));
        }
        mode_ = 0;
    }

    void NumericalModel::setMode(int mode)
    {
        if (modes_.empty()) { mode_ = 0; return; }
        mode_ = std::clamp(mode, 0, int(modes_.size()) - 1);
    }

    double NumericalModel::cutoffWavenumber() const
    {
        if (modes_.empty()) return 0.0;
        return std::sqrt(std::max(0.0, modes_[mode_].lambda));
    }
    double NumericalModel::resonantFrequency() const
    {
        return cutoffWavenumber() * kC0 / (2.0 * kPi);
    }

    bool NumericalModel::inside(double x, double y, double z) const
    {
        const int i = int(x / std::max(1e-12, dx_));
        const int j = int(y / std::max(1e-12, dy_));
        const int k = int(z / std::max(1e-12, dz_));
        if (i < 0 || i >= nx_ || j < 0 || j >= ny_ || k < 0 || k >= nz_) return false;
        return solid_[gidx(i, j, k)] != 0;
    }

    double NumericalModel::sampleField(const std::vector<double>& A,
                                       double ux, double uy, double uz) const
    {
        const int nx = nx_, ny = ny_, nz = nz_;
        double fx = ux / dx_ - 0.5, fy = uy / dy_ - 0.5, fz = uz / dz_ - 0.5;
        fx = std::clamp(fx, 0.0, double(nx - 1));
        fy = std::clamp(fy, 0.0, double(ny - 1));
        fz = std::clamp(fz, 0.0, double(nz - 1));
        const int i0 = std::min(int(fx), std::max(0, nx - 2));
        const int j0 = std::min(int(fy), std::max(0, ny - 2));
        const int k0 = std::min(int(fz), std::max(0, nz - 2));
        const int i1 = std::min(i0 + 1, nx - 1), j1 = std::min(j0 + 1, ny - 1), k1 = std::min(k0 + 1, nz - 1);
        const double tx = fx - i0, ty = fy - j0, tz = fz - k0;
        auto V = [&](int i, int j, int k) { return A[gidx(i, j, k)]; };
        const double c00 = V(i0,j0,k0)*(1-tx) + V(i1,j0,k0)*tx;
        const double c10 = V(i0,j1,k0)*(1-tx) + V(i1,j1,k0)*tx;
        const double c01 = V(i0,j0,k1)*(1-tx) + V(i1,j0,k1)*tx;
        const double c11 = V(i0,j1,k1)*(1-tx) + V(i1,j1,k1)*tx;
        const double c0 = c00*(1-ty) + c10*ty;
        const double c1 = c01*(1-ty) + c11*ty;
        return c0*(1-tz) + c1*tz;
    }

    std::array<double, 3> NumericalModel::vecAt(double ux, double uy, double uz) const
    {
        if (modes_.empty()) return {0.0, 0.0, 0.0};
        const Mode& m = modes_[mode_];
        return {sampleField(m.gx, ux, uy, uz),
                sampleField(m.gy, ux, uy, uz),
                sampleField(m.gz, ux, uy, uz)};
    }

    std::array<double, 3> NumericalModel::fieldVector(double x, double y, double z,
                                                      double phase) const
    {
        const std::array<double, 3> v = vecAt(x, y, z);
        const double ct = std::cos(phase);
        return {v[0] * ct, v[1] * ct, v[2] * ct};
    }

    std::pair<double, double> NumericalModel::transverseField(double x, double y, double z,
                                                              double phase) const
    {
        const std::array<double, 3> v = fieldVector(x, y, z, phase);
        return {v[0], v[1]};
    }

    std::vector<Particle> NumericalModel::sampleGrid(int nx, int ny, int nz,
                                                     bool cutawayOn,
                                                     float minIntensity,
                                                     double phase) const
    {
        (void)phase;
        std::vector<Particle> out;
        if (nx <= 1 || ny <= 1 || nz <= 1 || modes_.empty()) return out;
        const double peak = modes_[mode_].peak;
        const float halfW = float(W_) * 0.5f, halfH = float(H_) * 0.5f, halfD = float(D_) * 0.5f;
        out.reserve(std::size_t(nx) * ny * nz);
        for (int i = 0; i < nx; ++i)
        {
            const double ux = W_ * (i + 0.5) / nx;
            for (int j = 0; j < ny; ++j)
            {
                const double uy = H_ * (j + 0.5) / ny;
                for (int k = 0; k < nz; ++k)
                {
                    const double uz = D_ * (k + 0.5) / nz;
                    if (!inside(ux, uy, uz)) continue;
                    const std::array<double, 3> g = vecAt(ux, uy, uz);
                    const double mag = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
                    const float t = float(std::min(1.0, mag / peak));
                    if (t < minIntensity) continue;
                    const float cx = float(ux) - halfW, cy = float(uy) - halfH, cz = float(uz) - halfD;
                    if (cutawayOn && cx > 0.0f && cy > 0.0f && cz > 0.0f) continue;
                    float r, gg, bb; fireColormap(t, r, gg, bb);
                    Particle p; p.x = cx; p.y = cy; p.z = cz;
                    p.r = r; p.g = gg; p.b = bb; p.intensity = t;
                    out.push_back(p);
                }
            }
        }
        return out;
    }

} // namespace waveguide
