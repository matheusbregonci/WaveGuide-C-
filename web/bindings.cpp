// WebAssembly surface for the Waveguide/Cavity domain.
//
// The physics is NOT reimplemented here. This file only adapts the same
// TEmnModel / CylindricalModel / FieldViz that the desktop build uses, so the
// browser and the desktop cannot drift apart numerically -- which is the whole
// reason the visualization geometry was pulled out of main.cpp first.
//
// Everything that returns bulk data returns a copied Float32Array. A
// typed_memory_view would avoid the copy, but it aliases WASM linear memory,
// which is invalidated the moment the heap grows -- a use-after-free that shows
// up as silently corrupted geometry rather than a crash. The copies here are a
// few hundred kB per rebuild and only happen when parameters change.

#include <emscripten/bind.h>
#include <emscripten/val.h>

#include "CylindricalModel.hpp"
#include "FieldViz.hpp"
#include "TEmnModel.hpp"

#include <cmath>
#include <memory>
#include <vector>

using namespace emscripten;
using namespace waveguide;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kC0 = 299792458.0;

val toF32(const std::vector<float>& v)
{
    // .slice() with no arguments copies the view out of WASM memory.
    return val(typed_memory_view(v.size(), v.data()))
        .call<val>("slice");
}

} // namespace

// Plain-old-data configuration, mirroring the desktop UI controls.
struct GuideConfig {
    int    geometry  = 0;      // 0 = rectangular, 1 = cylindrical
    double widthMM   = 22.86;  // rectangular a
    double heightMM  = 10.16;  // rectangular b
    double radiusMM  = 23.83;  // cylindrical R
    double depthMM   = 300.0;  // guide length, both geometries
    double freqGHz   = 12.0;
    double epsR      = 1.0;
    double muR       = 1.0;
    int    modeM     = 1;
    int    modeN     = 0;
    int    modeL     = 1;      // axial index, cavity only
    int    modeType  = 0;      // 0 = TE, 1 = TM
    int    field     = 0;      // 0 = E, 1 = H
    int    structure = 0;      // 0 = waveguide, 1 = cavity
    double powerW    = 1.0;
};

class Guide {
public:
    // Rebuilding the model is what the desktop does on every parameter change
    // too: the constructor is where beta, the cutoff and the power-normalized
    // amplitude are derived, so there is no partially-updated state to get wrong.
    void configure(const GuideConfig& c)
    {
        cfg_ = c;
        const ModeType  mt = (c.modeType == 0) ? ModeType::TE : ModeType::TM;
        const FieldKind fk = (c.field == 0) ? FieldKind::Electric : FieldKind::Magnetic;
        const bool cav = (c.structure == 1);

        if (c.geometry == 0) {
            // TE allows either index to be zero but not both; TM needs both >= 1.
            int m = c.modeM, n = c.modeN;
            if (mt == ModeType::TM) { m = std::max(1, m); n = std::max(1, n); }
            else if (m == 0 && n == 0) m = 1;
            rect_ = std::make_unique<TEmnModel>(c.widthMM, c.heightMM, c.freqGHz * 1e9,
                                                c.epsR, c.muR, m, n, mt, fk,
                                                c.depthMM, cav, c.modeL, c.powerW);
            cyl_.reset();
        } else {
            // Circular: the radial index counts Bessel roots and starts at 1.
            cyl_ = std::make_unique<CylindricalModel>(c.radiusMM, c.depthMM / 1000.0,
                                                      c.freqGHz * 1e9, c.epsR, c.muR,
                                                      c.modeN, std::max(1, c.modeM),
                                                      mt, fk, cav, c.modeL);
            rect_.reset();
        }
        // The cached spans and references belong to the PREVIOUS model; keeping
        // them would silently draw the new mode against the old normalization.
        for (int i = 0; i < 3; ++i) { span_[i] = SectionSpan{}; ref_[i] = 1.0; }
    }

    // ---- scalars ----
    double peakField()  const { return src().peakField(); }
    bool   physicalUnits() const { return src().physicalUnits(); }
    double amplitude()  const { return rect_ ? rect_->amplitude() : 1.0; }
    double cutoffHz()   const { return src().resonantFrequency(); }
    double epsilonRel() const { return src().epsilonRel(); }
    double muRel()      const { return src().muRel(); }

    // Real part of the propagation constant; 0 below cutoff.
    double beta() const
    {
        const double v = kC0 / std::sqrt(src().epsilonRel() * src().muRel());
        const double k = 2.0 * kPi * cfg_.freqGHz * 1e9 / v;
        const double kc = src().cutoffWavenumber();
        const double d = k * k - kc * kc;
        return (d > 0.0) ? std::sqrt(d) : 0.0;
    }
    bool propagating() const { return beta() > 0.0; }

    // ---- bulk geometry ----

    // Flat [x, y, z, dc, a2, b2] per point, centred for display. dc/a2/b2 are
    // the Fourier terms of |field|^2 in the animation phase, so the caller
    // reconstructs the instantaneous intensity on the GPU without resampling.
    val cloud(int points, int harmonics, bool cutaway)
    {
        CloudParams p;
        p.targetPoints = points;
        p.harmonics = harmonics;
        p.cutaway = cutaway;
        const CloudResult r = bakeCloud(src(), p);
        std::vector<float> out;
        out.reserve(r.points.size() * 6);
        for (const CloudPoint& c : r.points) {
            out.push_back(c.x); out.push_back(c.y); out.push_back(c.z);
            out.push_back(c.dc); out.push_back(c.a2); out.push_back(c.b2);
        }
        lastInvPeak_ = r.invPeak;
        lastSpacing_ = r.meanSpacing;
        return toF32(out);
    }
    float cloudInvPeak()     const { return lastInvPeak_; }
    float cloudMeanSpacing() const { return lastSpacing_; }

    // Physical extent of a cut plane plus the slice it settled on. Returned as
    // an object because the front end needs all of it to set up axes.
    val sectionInfo(int plane, double manualFrac)
    {
        const SectionSpan s = sectionSpan(src(), plane, manualFrac);
        span_[plane] = s;
        ref_[plane] = sectionReference(src(), plane, s, nu(plane), nv(plane));
        val o = val::object();
        o.set("uMin", s.uMin); o.set("uMax", s.uMax);
        o.set("vMin", s.vMin); o.set("vMax", s.vMax);
        o.set("slice", s.slice);
        o.set("wMin", s.wMin); o.set("wMax", s.wMax);
        o.set("reference", ref_[plane]);
        return o;
    }

    // Flat [u, v, du, dv, inPlane, total] per arrow.
    val sectionArrows(int plane, double phase)
    {
        ensureSection(plane);
        const auto a = sampleSection(src(), plane, span_[plane], phase,
                                     nu(plane), nv(plane),
                                     ref_[plane], src().peakField());
        std::vector<float> out;
        out.reserve(a.size() * 6);
        for (const SectionArrow& s : a) {
            out.push_back(s.u); out.push_back(s.v);
            out.push_back(s.du); out.push_back(s.dv);
            out.push_back(s.inPlane); out.push_back(s.total);
        }
        return toF32(out);
    }

    // Streamlines as one flat vertex array plus the index where each line
    // starts, so the front end can issue one LINE_STRIP per range without
    // walking a nested structure across the JS boundary.
    val sectionLines(int plane, double phase, int seedsU, int seedsV)
    {
        ensureSection(plane);
        const auto lines = traceSection(src(), plane, span_[plane], phase,
                                        seedsU, seedsV,
                                        ref_[plane], src().peakField());
        std::vector<float> verts;
        std::vector<float> starts;
        for (const Streamline& L : lines) {
            starts.push_back(float(verts.size() / 4));
            for (size_t i = 0; i < L.u.size(); ++i) {
                verts.push_back(L.u[i]);       verts.push_back(L.v[i]);
                verts.push_back(L.inPlane[i]); verts.push_back(L.total[i]);
            }
        }
        starts.push_back(float(verts.size() / 4));   // sentinel end
        val o = val::object();
        o.set("verts", toF32(verts));      // 4 floats per vertex
        o.set("starts", toF32(starts));    // n+1 entries
        return o;
    }

    // Single-point probe, for readouts and tooltips.
    val fieldAt(double x, double y, double z, double phase) const
    {
        const std::array<double, 3> V = src().fieldVector(x, y, z, phase);
        val o = val::object();
        o.set("x", V[0]); o.set("y", V[1]); o.set("z", V[2]);
        o.set("mag", std::sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]));
        return o;
    }

    // Domain box in model coordinates -- the front end needs it to build the
    // outline and to centre the camera. Rectangular spans [0,a]x[0,b]x[0,d]
    // while the cylinder is centred on its axis, so this must be asked for
    // rather than derived from the extent.
    val domain() const
    {
        double x0, x1, y0, y1, z0, z1;
        src().domain(x0, x1, y0, y1, z0, z1);
        val o = val::object();
        o.set("x0", x0); o.set("x1", x1);
        o.set("y0", y0); o.set("y1", y1);
        o.set("z0", z0); o.set("z1", z1);
        return o;
    }

private:
    // configure() has to have run. Returning a reference to a null unique_ptr
    // would be undefined behaviour that shows up as garbage numbers rather than
    // a crash, so fall back to a default guide and keep the failure visible in
    // the values instead of corrupting memory.
    const FieldSource& src() const
    {
        if (!rect_ && !cyl_) {
            static const TEmnModel fallback;
            return fallback;
        }
        return rect_ ? static_cast<const FieldSource&>(*rect_)
                     : static_cast<const FieldSource&>(*cyl_);
    }

    // sectionArrows/sectionLines need the span and the reference that
    // sectionInfo computes. Rather than trust the caller to have asked in the
    // right order -- the exact trap that made peakField() return 0 for anyone
    // who had not drawn a particle grid first -- fill them on demand.
    void ensureSection(int plane)
    {
        if (span_[plane].uMax > span_[plane].uMin) return;
        span_[plane] = sectionSpan(src(), plane, -1.0);
        ref_[plane] = sectionReference(src(), plane, span_[plane],
                                       nu(plane), nv(plane));
    }
    // ZX/ZY run along the long z axis, so they get more columns.
    static int nu(int plane) { return (plane == 0) ? 18 : 44; }
    static int nv(int plane) { return (plane == 0) ? 18 : 14; }

    GuideConfig cfg_;
    std::unique_ptr<TEmnModel> rect_;
    std::unique_ptr<CylindricalModel> cyl_;
    SectionSpan span_[3];
    double ref_[3] = {1.0, 1.0, 1.0};
    float lastInvPeak_ = 1.0f;
    float lastSpacing_ = 1e-3f;
};

EMSCRIPTEN_BINDINGS(waveguide_module)
{
    value_object<GuideConfig>("GuideConfig")
        .field("geometry",  &GuideConfig::geometry)
        .field("widthMM",   &GuideConfig::widthMM)
        .field("heightMM",  &GuideConfig::heightMM)
        .field("radiusMM",  &GuideConfig::radiusMM)
        .field("depthMM",   &GuideConfig::depthMM)
        .field("freqGHz",   &GuideConfig::freqGHz)
        .field("epsR",      &GuideConfig::epsR)
        .field("muR",       &GuideConfig::muR)
        .field("modeM",     &GuideConfig::modeM)
        .field("modeN",     &GuideConfig::modeN)
        .field("modeL",     &GuideConfig::modeL)
        .field("modeType",  &GuideConfig::modeType)
        .field("field",     &GuideConfig::field)
        .field("structure", &GuideConfig::structure)
        .field("powerW",    &GuideConfig::powerW);

    class_<Guide>("Guide")
        .constructor<>()
        .function("configure",        &Guide::configure)
        .function("peakField",        &Guide::peakField)
        .function("physicalUnits",    &Guide::physicalUnits)
        .function("amplitude",        &Guide::amplitude)
        .function("cutoffHz",         &Guide::cutoffHz)
        .function("epsilonRel",       &Guide::epsilonRel)
        .function("muRel",            &Guide::muRel)
        .function("beta",             &Guide::beta)
        .function("propagating",      &Guide::propagating)
        .function("cloud",            &Guide::cloud)
        .function("cloudInvPeak",     &Guide::cloudInvPeak)
        .function("cloudMeanSpacing", &Guide::cloudMeanSpacing)
        .function("sectionInfo",      &Guide::sectionInfo)
        .function("sectionArrows",    &Guide::sectionArrows)
        .function("sectionLines",     &Guide::sectionLines)
        .function("fieldAt",          &Guide::fieldAt)
        .function("domain",           &Guide::domain);
}
