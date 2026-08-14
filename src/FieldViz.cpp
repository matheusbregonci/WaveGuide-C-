#include "FieldViz.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <random>

namespace waveguide {

namespace {

constexpr double kPi = 3.14159265358979323846;

// Project a field vector onto the cut plane. Which two components are "in
// plane" depends on the cut, and getting this wrong is silent, so it lives in
// exactly one place.
inline void inPlanePair(int plane, const std::array<double, 3>& V,
                        double& a, double& b)
{
    if (plane == 0)      { a = V[0]; b = V[1]; }   // XY
    else if (plane == 1) { a = V[2]; b = V[0]; }   // ZX
    else                 { a = V[2]; b = V[1]; }   // ZY
}

// Map plane coordinates (u, v) plus the slice position to model coordinates.
inline void toModel(int plane, double u, double v, double slice,
                    double& x, double& y, double& z)
{
    if (plane == 0)      { x = u; y = v; z = slice; }
    else if (plane == 1) { z = u; x = v; y = slice; }
    else                 { z = u; y = v; x = slice; }
}

} // namespace

// ---------------------------------------------------------------- cloud ----

CloudResult bakeCloud(const FieldSource& src, const CloudParams& p)
{
    CloudResult out;

    // Sample in the source's OWN frame. Assuming [0, extent] here fed [0, 2R]
    // coordinates to the cylinder's centred inside(), which kept only the
    // quarter disc near the corner and drew it offset from the mesh.
    double x0, x1, y0, y1, z0, z1;
    src.domain(x0, x1, y0, y1, z0, z1);
    const double W = x1 - x0, H = y1 - y0, D = z1 - z0;
    const double xc = 0.5 * (x0 + x1), yc = 0.5 * (y0 + y1), zc = 0.5 * (z0 + z1);

    const double peak = src.peakField();
    // Field not built yet: bail rather than burning the whole rejection budget
    // on an all-zero volume.
    if (peak <= 1e-20) return out;
    out.invPeak = float(1.0 / peak);

    const int K = std::max(1, p.harmonics);
    const double invK = 1.0 / K;
    const int target = std::max(1000, p.targetPoints);
    const int attemptCap = target * 25;

    std::mt19937 rng(p.seed);
    std::uniform_real_distribution<double> Ux(x0, x1), Uy(y0, y1), Uz(z0, z1);
    out.points.reserve(target);

    int attempts = 0;
    for (; attempts < attemptCap && int(out.points.size()) < target; ++attempts) {
        const double u = Ux(rng), v = Uy(rng), w = Uz(rng);
        if (!src.inside(u, v, w)) continue;
        const float cx = float(u - xc), cy = float(v - yc), cz = float(w - zc);
        if (p.cutaway && cx > 0.0f && cy > 0.0f && cz > 0.0f) continue;

        // Fourier terms of |field|^2 over the animation phase. Storing these
        // lets the GPU reconstruct the instantaneous intensity every frame
        // without the CPU resampling the field.
        double dc = 0, a2 = 0, b2 = 0;
        for (int kk = 0; kk < K; ++kk) {
            const double ph = 2.0 * kPi * kk / K;
            const std::array<double, 3> f = src.fieldVector(u, v, w, ph);
            const double m2 = f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
            dc += m2; a2 += m2 * std::cos(2 * ph); b2 += m2 * std::sin(2 * ph);
        }
        dc *= invK; a2 *= 2 * invK; b2 *= 2 * invK;
        const double env = std::sqrt(dc + std::hypot(a2, b2));
        if (env * out.invPeak < p.minIntensity) continue;

        CloudPoint cp;
        cp.x = cx; cp.y = cy; cp.z = cz;
        cp.dc = float(dc); cp.a2 = float(a2); cp.b2 = float(b2);
        out.points.push_back(cp);
    }

    // Accepted points fill accepted/attempts of the box, so their mean spacing
    // is ~cbrt(volume / attempts) -- which is what sets a sensible dot size.
    out.meanSpacing = float(std::cbrt(std::max(1e-30, W * H * D) /
                                      std::max(1, attempts)));
    return out;
}

// -------------------------------------------------------------- sections ----

SectionSpan sectionSpan(const FieldSource& src, int plane, double manualFrac)
{
    SectionSpan s;
    double x0, x1, y0, y1, z0, z1;
    src.domain(x0, x1, y0, y1, z0, z1);

    if (plane == 0)      { s.uMin = x0; s.uMax = x1; s.vMin = y0; s.vMax = y1;
                           s.wMin = z0; s.wMax = z1; }
    else if (plane == 1) { s.uMin = z0; s.uMax = z1; s.vMin = x0; s.vMax = x1;
                           s.wMin = y0; s.wMax = y1; }
    else                 { s.uMin = z0; s.uMax = z1; s.vMin = y0; s.vMax = y1;
                           s.wMin = x0; s.wMax = x1; }

    if (manualFrac >= 0.0) {
        s.slice = s.wMin + manualFrac * (s.wMax - s.wMin);
        return s;
    }

    // Antinode search: scan the perpendicular axis for the strongest in-plane
    // field, over a coarse grid and a few phases.
    const int NW = 17, CU = 12, CV = 8, NP = 6;
    const double Uspan = s.uMax - s.uMin, Vspan = s.vMax - s.vMin;
    double bestW = 0.5 * (s.wMin + s.wMax), bestMax = 1e-30;
    for (int iw = 0; iw < NW; ++iw) {
        const double w = s.wMin + (iw / double(NW - 1)) * (s.wMax - s.wMin);
        double wmax = 0.0;
        for (int pp = 0; pp < NP; ++pp) {
            const double ph = 2.0 * kPi * pp / NP;
            for (int iv = 0; iv < CV; ++iv)
            for (int iu = 0; iu < CU; ++iu) {
                const double u = s.uMin + (iu + 0.5) / CU * Uspan;
                const double v = s.vMin + (iv + 0.5) / CV * Vspan;
                double x, y, z; toModel(plane, u, v, w, x, y, z);
                if (!src.inside(x, y, z)) continue;
                double a, b; inPlanePair(plane, src.fieldVector(x, y, z, ph), a, b);
                wmax = std::max(wmax, std::sqrt(a*a + b*b));
            }
        }
        if (wmax > bestMax) { bestMax = wmax; bestW = w; }
    }
    s.slice = bestW;
    return s;
}

namespace {

// Shared sampler: field at plane coordinates, zero outside the guide.
inline std::array<double, 3> sampleAt(const FieldSource& src, int plane,
                                      const SectionSpan& s,
                                      double u, double v, double ph)
{
    double x, y, z; toModel(plane, u, v, s.slice, x, y, z);
    if (!src.inside(x, y, z)) return {0.0, 0.0, 0.0};
    return src.fieldVector(x, y, z, ph);
}

} // namespace

double sectionReference(const FieldSource& src, int plane, const SectionSpan& s,
                        int nu, int nv, int phaseSamples)
{
    const double Uspan = s.uMax - s.uMin, Vspan = s.vMax - s.vMin;
    double mx = 1e-30;
    for (int pp = 0; pp < std::max(1, phaseSamples); ++pp) {
        const double ph = 2.0 * kPi * pp / std::max(1, phaseSamples);
        for (int iv = 0; iv < nv; ++iv)
        for (int iu = 0; iu < nu; ++iu) {
            double a, b;
            inPlanePair(plane, sampleAt(src, plane, s,
                                        s.uMin + (iu + 0.5) / nu * Uspan,
                                        s.vMin + (iv + 0.5) / nv * Vspan, ph), a, b);
            mx = std::max(mx, std::sqrt(a*a + b*b));
        }
    }
    return mx;
}

double sectionInstantMax(const FieldSource& src, int plane, const SectionSpan& s,
                         double phase, int nu, int nv)
{
    const double Uspan = s.uMax - s.uMin, Vspan = s.vMax - s.vMin;
    double mx = 1e-30;
    for (int iv = 0; iv < nv; ++iv)
    for (int iu = 0; iu < nu; ++iu) {
        double a, b;
        inPlanePair(plane, sampleAt(src, plane, s,
                                    s.uMin + (iu + 0.5) / nu * Uspan,
                                    s.vMin + (iv + 0.5) / nv * Vspan, phase), a, b);
        mx = std::max(mx, std::hypot(a, b));
    }
    return mx;
}

std::vector<SectionArrow> sampleSection(const FieldSource& src, int plane,
                                        const SectionSpan& s, double phase,
                                        int nu, int nv,
                                        double refInPlane, double refTotal)
{
    std::vector<SectionArrow> out;
    out.reserve(size_t(nu) * size_t(nv));
    const double Uspan = s.uMax - s.uMin, Vspan = s.vMax - s.vMin;
    const double rIn = std::max(1e-30, refInPlane);
    const double rTot = std::max(1e-30, refTotal);

    for (int iv = 0; iv < nv; ++iv)
    for (int iu = 0; iu < nu; ++iu) {
        const double u = s.uMin + (iu + 0.5) / nu * Uspan;
        const double v = s.vMin + (iv + 0.5) / nv * Vspan;
        const std::array<double, 3> V = sampleAt(src, plane, s, u, v, phase);
        double a, b; inPlanePair(plane, V, a, b);
        const double mg = std::sqrt(a*a + b*b);
        if (mg / rIn < 0.02) continue;

        SectionArrow ar;
        ar.u = float(u); ar.v = float(v);
        ar.du = float(a / mg); ar.dv = float(b / mg);
        ar.inPlane = float(std::min(1.0, mg / rIn));
        ar.total = float(std::min(1.0,
            std::sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]) / rTot));
        out.push_back(ar);
    }
    return out;
}

std::vector<Streamline> traceSection(const FieldSource& src, int plane,
                                     const SectionSpan& s, double phase,
                                     int seedsU, int seedsV,
                                     double refInPlane, double refTotal,
                                     double halfArcFrac)
{
    std::vector<Streamline> out;
    const double Uspan = s.uMax - s.uMin, Vspan = s.vMax - s.vMin;
    const double diag = std::sqrt(Uspan*Uspan + Vspan*Vspan);
    const double h = diag / 300.0;
    const int maxSteps = std::max(6, int(halfArcFrac * diag / h));
    const double rIn = std::max(1e-30, refInPlane);
    const double rTot = std::max(1e-30, refTotal);

    // Cutoff measured against the field RIGHT NOW (see the header): a cut whose
    // whole pattern pulses through zero would otherwise blank twice per cycle.
    const double eps = 1e-4 * sectionInstantMax(src, plane, s, phase,
                                                std::max(4, seedsU),
                                                std::max(4, seedsV));

    auto vec = [&](double u, double v, double& a, double& b, double& m3) {
        const std::array<double, 3> V = sampleAt(src, plane, s, u, v, phase);
        inPlanePair(plane, V, a, b);
        m3 = std::sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
    };

    for (int j = 0; j < seedsV; ++j)
    for (int i = 0; i < seedsU; ++i) {
        const double su = s.uMin + (i + 0.5) / double(seedsU) * Uspan;
        const double sv = s.vMin + (j + 0.5) / double(seedsV) * Vspan;
        { double a, b, m3; vec(su, sv, a, b, m3);
          if (std::hypot(a, b) < eps) continue; }

        // Integrate both ways, then splice so the result runs along +field.
        std::vector<float> bu, bv, bi, bt;   // backward half
        std::vector<float> fu, fv, fi, ft;   // forward half
        for (int dir = -1; dir <= 1; dir += 2) {
            std::vector<float>& U = (dir < 0) ? bu : fu;
            std::vector<float>& V_ = (dir < 0) ? bv : fv;
            std::vector<float>& I = (dir < 0) ? bi : fi;
            std::vector<float>& T = (dir < 0) ? bt : ft;
            double u = su, v = sv;
            { double a, b, m3; vec(u, v, a, b, m3);
              U.push_back(float(u)); V_.push_back(float(v));
              I.push_back(float(std::min(1.0, std::hypot(a, b) / rIn)));
              T.push_back(float(std::min(1.0, m3 / rTot))); }
            for (int st = 0; st < maxSteps; ++st) {
                double a, b, m3; vec(u, v, a, b, m3);
                const double m = std::hypot(a, b);
                if (m < eps) break;
                double a2, b2, m3b;
                vec(u + dir*(a/m)*h*0.5, v + dir*(b/m)*h*0.5, a2, b2, m3b);
                const double m2 = std::hypot(a2, b2);
                if (m2 < eps) break;
                const double nu = u + dir*(a2/m2)*h;
                const double nv = v + dir*(b2/m2)*h;
                if (nu < s.uMin || nu > s.uMax || nv < s.vMin || nv > s.vMax) break;
                double xx, yy, zz; toModel(plane, nu, nv, s.slice, xx, yy, zz);
                if (!src.inside(xx, yy, zz)) break;
                U.push_back(float(nu)); V_.push_back(float(nv));
                I.push_back(float(std::min(1.0, m / rIn)));
                T.push_back(float(std::min(1.0, m3 / rTot)));
                u = nu; v = nv;
            }
        }

        Streamline sl;
        const size_t n = bu.size() + (fu.empty() ? 0 : fu.size() - 1);
        sl.u.reserve(n); sl.v.reserve(n); sl.inPlane.reserve(n); sl.total.reserve(n);
        for (size_t k = bu.size(); k-- > 0; ) {
            sl.u.push_back(bu[k]); sl.v.push_back(bv[k]);
            sl.inPlane.push_back(bi[k]); sl.total.push_back(bt[k]);
        }
        for (size_t k = 1; k < fu.size(); ++k) {
            sl.u.push_back(fu[k]); sl.v.push_back(fv[k]);
            sl.inPlane.push_back(fi[k]); sl.total.push_back(ft[k]);
        }
        if (sl.u.size() >= 2) out.push_back(std::move(sl));
    }
    return out;
}

} // namespace waveguide
