#pragma once

#include "FieldSource.hpp"
#include "Types.hpp"

#include <cstdint>
#include <vector>

// Field visualization GEOMETRY, with no renderer and no UI attached.
//
// Everything here turns a FieldSource into plain buffers -- point clouds,
// polylines, arrow lists -- expressed in the source's own physical coordinates.
// Mapping those to pixels, uploading them to a GPU or drawing them with ImGui is
// the caller's business.
//
// This split exists because the same computations have to serve three very
// different front ends: the desktop OpenGL view, the software rasterizer in the
// exporter, and (next) a WebAssembly build driving WebGL from TypeScript. They
// used to live inside main.cpp's drawing lambdas, captured by reference from UI
// globals, which made them impossible to reuse or to test without a window.
//
// Nothing in this header includes GL, ImGui, or any platform header.

namespace waveguide {

// ---------------------------------------------------------------- cloud ----

struct CloudParams {
    int      targetPoints   = 60000;   // accepted (VISIBLE) points aimed for
    int      harmonics      = 8;       // phase samples used to build the envelope
    bool     cutaway        = false;   // drop the +x/+y/+z octant
    float    minIntensity   = 0.05f;   // reject below this fraction of the peak
    uint32_t seed           = 0xC0FFEEu;
};

struct CloudResult {
    std::vector<CloudPoint> points;
    float invPeak     = 1.0f;   // 1 / peakField(), the shader's normalizer
    float meanSpacing = 1e-3f;  // physical spacing -> on-screen dot size
};

// Monte-Carlo sample of the field volume. Random positions rather than a
// lattice, so the cloud has no visible planes of points. The RNG is seeded
// deterministically: the same field must give the same cloud every frame, or
// the points would swim.
CloudResult bakeCloud(const FieldSource& src, const CloudParams& p);

// -------------------------------------------------------------- sections ----

// A cut plane: 0 = XY (vary x,y), 1 = ZX (vary z,x), 2 = ZY (vary z,y).
struct SectionSpan {
    double uMin = 0.0, uMax = 0.0;   // horizontal physical range
    double vMin = 0.0, vMax = 0.0;   // vertical physical range
    double slice = 0.0;              // perpendicular coordinate of the cut
    double wMin = 0.0, wMax = 0.0;   // valid range of that perpendicular axis
};

// Physical extent of a cut plane, and where to put it.
//
// With manualFrac < 0 the slice is placed at an ANTINODE: the plane is scanned
// along its perpendicular axis for the position with the strongest in-plane
// field. Placing it at the geometric centre instead is a coin flip that lands on
// a node about as often as not, blanking the whole plot.
SectionSpan sectionSpan(const FieldSource& src, int plane, double manualFrac = -1.0);

// Largest in-plane magnitude anywhere in the cut over a FULL cycle. Used to
// scale arrow length and streamline head size. It must be a max over phase as
// well as position: a per-frame maximum would rescale every frame and make the
// arrows breathe when the field does.
double sectionReference(const FieldSource& src, int plane, const SectionSpan& s,
                        int nu, int nv, int phaseSamples = 12);

// Largest in-plane magnitude at ONE instant. The weak-field cutoff for tracing
// has to be measured against this, not against the cycle maximum: a transverse
// cut is one fixed shape times sin(wt - bz), so it collapses to zero twice per
// cycle and a cutoff pinned to the cycle peak would blank every line as it
// passed through.
double sectionInstantMax(const FieldSource& src, int plane, const SectionSpan& s,
                         double phase, int nu, int nv);

// One arrow of the quiver plot. Position and direction are physical; the two
// scalars are already normalized and ready to drive length and colour.
//   inPlane = |V_t| / sectionReference   -> length
//   total   = |V|   / peakField          -> colour, matching the 3D view
struct SectionArrow {
    float u = 0, v = 0;      // cell centre, physical
    float du = 0, dv = 0;    // unit in-plane direction
    float inPlane = 0;
    float total = 0;
};

std::vector<SectionArrow> sampleSection(const FieldSource& src, int plane,
                                        const SectionSpan& s, double phase,
                                        int nu, int nv,
                                        double refInPlane, double refTotal);

// A traced field line. Vertices are physical (u,v); the parallel scalar arrays
// carry the same two normalizations the arrows use.
struct Streamline {
    std::vector<float> u, v;
    std::vector<float> inPlane;
    std::vector<float> total;
};

// Even-spaced field lines through the cut.
//
// Seeds sit on a FIXED lattice and every seed is used, and each line runs a
// fixed arclength. There is deliberately no occupancy grid: the classic
// evenly-spaced placement (reject a seed near an existing line, stop a line when
// it meets one) couples every line to the ones traced before it, so on a moving
// field one flipped test cascades and the lines wink in and out. Independent
// lines vary as smoothly as the field does.
//
// Lines are ordered along +field, so a caller drawing arrowheads can take the
// direction straight from consecutive vertices.
std::vector<Streamline> traceSection(const FieldSource& src, int plane,
                                     const SectionSpan& s, double phase,
                                     int seedsU, int seedsV,
                                     double refInPlane, double refTotal,
                                     double halfArcFrac = 0.16);

} // namespace waveguide
