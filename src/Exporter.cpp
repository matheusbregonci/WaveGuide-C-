#include "Exporter.hpp"
#include "Colormap.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_EASY_FONT_IMPLEMENTATION
#include "stb_easy_font.h"

#include "gif.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <glm/gtc/matrix_transform.hpp>

namespace waveguide {

namespace {

constexpr double kPi = 3.14159265358979323846;

// ---------- Software RGBA canvas ----------
struct Canvas {
    int w = 0, h = 0;
    std::vector<uint8_t> px; // RGBA, top-left origin

    void resize(int ww, int hh, uint8_t r = 0, uint8_t g = 0, uint8_t b = 0) {
        w = ww; h = hh;
        px.assign(size_t(w) * h * 4, 0);
        for (size_t i = 0; i < px.size(); i += 4) {
            px[i]=r; px[i+1]=g; px[i+2]=b; px[i+3]=255;
        }
    }

    inline void setPx(int x, int y, uint8_t r, uint8_t g, uint8_t b, uint8_t a=255) {
        if ((unsigned)x >= (unsigned)w || (unsigned)y >= (unsigned)h) return;
        uint8_t* p = &px[(y * w + x) * 4];
        if (a == 255) { p[0]=r; p[1]=g; p[2]=b; p[3]=255; return; }
        const float af = a / 255.0f;
        p[0] = uint8_t(p[0]*(1-af) + r*af);
        p[1] = uint8_t(p[1]*(1-af) + g*af);
        p[2] = uint8_t(p[2]*(1-af) + b*af);
        p[3] = 255;
    }

    void fillRect(int x0, int y0, int x1, int y1,
                  uint8_t r, uint8_t g, uint8_t b, uint8_t a=255) {
        if (x1 < x0) std::swap(x0, x1);
        if (y1 < y0) std::swap(y0, y1);
        x0 = std::max(0, x0); y0 = std::max(0, y0);
        x1 = std::min(w - 1, x1); y1 = std::min(h - 1, y1);
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x)
                setPx(x, y, r, g, b, a);
    }

    void rect(int x0, int y0, int x1, int y1,
              uint8_t r, uint8_t g, uint8_t b) {
        fillRect(x0, y0, x1, y0, r, g, b);
        fillRect(x0, y1, x1, y1, r, g, b);
        fillRect(x0, y0, x0, y1, r, g, b);
        fillRect(x1, y0, x1, y1, r, g, b);
    }

    // Wu-ish line: simple DDA with thickness via 3x3 stamp.
    void line(float x0, float y0, float x1, float y1,
              uint8_t r, uint8_t g, uint8_t b, float thickness=1.0f) {
        const float dx = x1 - x0, dy = y1 - y0;
        const int steps = int(std::ceil(std::max(std::fabs(dx), std::fabs(dy))));
        if (steps <= 0) { setPx(int(x0), int(y0), r, g, b); return; }
        const float sx = dx / steps, sy = dy / steps;
        const int half = int(std::max(0.0f, (thickness - 1.0f) * 0.5f));
        for (int i = 0; i <= steps; ++i) {
            const int xi = int(std::round(x0 + sx * i));
            const int yi = int(std::round(y0 + sy * i));
            for (int oy = -half; oy <= half; ++oy)
                for (int ox = -half; ox <= half; ++ox)
                    setPx(xi + ox, yi + oy, r, g, b);
        }
    }

    void circle(int cx, int cy, int rad,
                uint8_t r, uint8_t g, uint8_t b) {
        // Midpoint circle
        int x = rad, y = 0, err = 0;
        while (x >= y) {
            setPx(cx + x, cy + y, r, g, b);
            setPx(cx + y, cy + x, r, g, b);
            setPx(cx - y, cy + x, r, g, b);
            setPx(cx - x, cy + y, r, g, b);
            setPx(cx - x, cy - y, r, g, b);
            setPx(cx - y, cy - x, r, g, b);
            setPx(cx + y, cy - x, r, g, b);
            setPx(cx + x, cy - y, r, g, b);
            ++y;
            if (err <= 0) err += 2 * y + 1;
            if (err > 0)  { --x; err -= 2 * x + 1; }
        }
    }

    // Text via stb_easy_font: each glyph is expressed as small quads which
    // we rasterize as tiny filled rects (sufficient for UI labels). The glyphs
    // are laid out at the origin and then scaled about (x, y), so `scale` gives
    // a genuinely larger font rather than a stretched bitmap.
    void text(int x, int y, const char* s,
              uint8_t r, uint8_t g, uint8_t b, float scale = 1.0f) {
        static thread_local char buf[99999];
        unsigned char col[4] = {r, g, b, 255};
        const int nq = stb_easy_font_print(0.0f, 0.0f,
                                           const_cast<char*>(s),
                                           col, buf, int(sizeof(buf)));
        const float* v = reinterpret_cast<const float*>(buf);
        for (int q = 0; q < nq; ++q) {
            const float x0 = x + v[0] * scale, y0 = y + v[1] * scale;
            const float x1 = x + v[8] * scale, y1 = y + v[9] * scale;
            fillRect(int(x0), int(y0), int(x1), int(y1), r, g, b);
            v += 16; // 4 verts * 4 floats
        }
    }
    static int textW(const char* s, float scale = 1.0f) {
        return int(stb_easy_font_width(const_cast<char*>(s)) * scale + 0.5f);
    }
    static int textH(float scale = 1.0f) { return int(7.0f * scale + 0.5f); }
};

// Single shared definition (Colormap.hpp). This used to hold its own
// black->purple->white stops under a comment claiming they matched the 3D
// view; they had not for a long time, so every exported cross section, field
// line and colour bar was drawn on a different scale from the 3D cloud sitting
// next to it in the same document.
void fireMap(float t, float& r, float& g, float& b) { fireColor(t, r, g, b); }

// ---------- Canvas composition ----------
void blit(Canvas& dst, int dx, int dy, const Canvas& src) {
    for (int y = 0; y < src.h; ++y) {
        const int yy = dy + y;
        if (yy < 0 || yy >= dst.h) continue;
        for (int x = 0; x < src.w; ++x) {
            const int xx = dx + x;
            if (xx < 0 || xx >= dst.w) continue;
            const uint8_t* s = &src.px[(y * src.w + x) * 4];
            uint8_t* d = &dst.px[(yy * dst.w + xx) * 4];
            d[0]=s[0]; d[1]=s[1]; d[2]=s[2]; d[3]=255;
        }
    }
}

// ---------- Item drawers ----------

// Color scale bar (static)
void drawColorBar(Canvas& c, const ExportContext& ctx) {
    c.resize(220, 340, 15, 15, 20);
    const bool isE = (ctx.fieldKind == 0);
    const double peak = (ctx.geometry == 0) ? ctx.rect->peakField()
                                             : ctx.cyl ->peakField();
    const char* unit = isE ? "V/m" : "A/m";
    const char* sym  = isE ? "|E|" : "|H|";
    char hdr[48]; std::snprintf(hdr, sizeof(hdr), "%s (%s)", sym, unit);
    c.text(10, 10, hdr, 255, 255, 255);

    const int x0 = 16, y0 = 30;
    const int barW = 44, barH = 280;
    const int x1 = x0 + barW, y1 = y0 + barH;
    for (int y = 0; y < barH; ++y) {
        const float t = 1.0f - float(y) / float(barH - 1);
        float r, g, b; fireMap(t, r, g, b);
        c.fillRect(x0, y0 + y, x1, y0 + y,
                   uint8_t(r*255), uint8_t(g*255), uint8_t(b*255));
    }
    c.rect(x0, y0, x1, y1, 200, 200, 200);

    char l1[32], l2[32], l3[32];
    std::snprintf(l1, sizeof(l1), "%.3g %s", peak, unit);
    std::snprintf(l2, sizeof(l2), "%.3g %s", peak * 0.5, unit);
    std::snprintf(l3, sizeof(l3), "0 %s", unit);
    c.text(x1 + 6, y0 - 4, l1, 255, 255, 255);
    c.text(x1 + 6, y0 + barH/2 - 4, l2, 255, 255, 255);
    c.text(x1 + 6, y1 - 8, l3, 255, 255, 255);
}

// Width of the scale strip appended to every field figure.
constexpr int kBarPaneW = 150;

// Append the colour scale to the right of a finished figure. Every plot that
// maps magnitude to colour carries its own scale, so a figure lifted out of the
// document into a paper still says what its colours mean -- no cross-reference
// to a separate legend, and nothing to caption, since a colour ramp with
// numbers on it explains itself.
void attachColorBar(Canvas& fig, const ExportContext& ctx) {
    if (fig.w <= 0 || fig.h <= 0) return;

    // Match the figure's own background and pick a readable ink for it: the
    // cross sections are near-black, the 3D render is near-white.
    const uint8_t br = fig.px[0], bg_ = fig.px[1], bb_ = fig.px[2];
    const float lum = (0.299f * br + 0.587f * bg_ + 0.114f * bb_) / 255.0f;
    const uint8_t ink = (lum > 0.5f) ? 0 : 255;

    Canvas out;
    out.resize(fig.w + kBarPaneW, fig.h, br, bg_, bb_);
    blit(out, 0, 0, fig);

    const FieldSource* src = (ctx.geometry == 0)
        ? static_cast<const FieldSource*>(ctx.rect)
        : static_cast<const FieldSource*>(ctx.cyl);
    const double peak = src->peakField();
    const bool phys = src->physicalUnits();
    const char* unit = phys ? ((ctx.fieldKind == 0) ? "V/m" : "A/m") : "u.a.";

    const int barW = 34;
    const int x0 = fig.w + 16, x1 = x0 + barW;
    const int y0 = 22, y1 = fig.h - 22;
    const int barH = std::max(1, y1 - y0);
    for (int y = 0; y < barH; ++y) {
        const float t = 1.0f - float(y) / float(std::max(1, barH - 1));
        float r, g, b; fireMap(t, r, g, b);
        out.fillRect(x0, y0 + y, x1, y0 + y,
                     uint8_t(r*255), uint8_t(g*255), uint8_t(b*255));
    }
    out.rect(x0, y0, x1, y1, ink, ink, ink);

    const float ts = 1.5f;
    char buf[48];
    for (int k = 0; k <= 4; ++k) {
        const float f = k / 4.0f;                     // 0 at the bottom
        const int py = y1 - int(f * barH + 0.5f);
        out.line(float(x1), float(py), float(x1 + 5), float(py), ink, ink, ink, 1.0f);
        std::snprintf(buf, sizeof(buf), "%.3g", peak * f);
        out.text(x1 + 9, py - Canvas::textH(ts) / 2, buf, ink, ink, ink, ts);
    }
    // Unit once, at the top of the ramp; the numbers below are all in it.
    out.text(x0, y0 - 16, unit, ink, ink, ink, ts);

    fig = std::move(out);
}

// Spectrum (static)
void drawSpectrum(Canvas& c, const ExportContext& ctx) {
    c.resize(480, 260, 15, 15, 20);
    const double c0 = 299792458.0;
    double kc, er, mr;
    if (ctx.geometry == 0) {
        kc = std::sqrt(ctx.rect->cutoffKcSquared());
        er = ctx.rect->epsilonRel();
        mr = ctx.rect->muRel();
    } else {
        kc = ctx.cyl->cutoffKc();
        er = ctx.cyl->epsilonRel();
        mr = ctx.cyl->muRel();
    }
    const double vph = c0 / std::sqrt(er * mr);
    const double fc = kc * vph / (2.0 * kPi);
    const double fMax = (fc > 0) ? 2.0 * fc : 1.0;

    char line1[96];
    if (fc >= 1e9) std::snprintf(line1, sizeof(line1), "f_c = %.3f GHz", fc*1e-9);
    else           std::snprintf(line1, sizeof(line1), "f_c = %.3f MHz", fc*1e-6);
    c.text(10, 10, line1, 255, 255, 255);
    char line2[96];
    std::snprintf(line2, sizeof(line2), "f = %.3f GHz   v_ph = %.3e m/s",
                  ctx.freqHz * 1e-9, vph);
    c.text(10, 26, line2, 220, 220, 220);

    const int px0 = 20, py0 = 50, px1 = 460, py1 = 240;
    c.rect(px0, py0, px1, py1, 120, 120, 140);
    // axes labels
    c.text(px0, py1 + 4, "0", 200,200,200);
    char fmx[32];
    std::snprintf(fmx, sizeof(fmx), "%.2f GHz", fMax * 1e-9);
    c.text(px1 - 50, py1 + 4, fmx, 200,200,200);
    c.text(px0 - 10, py0 - 12, "beta/k_c  alpha/k_c", 200,200,200);

    // Plot beta & alpha curves normalized to [0,1.5]
    const int N = (px1 - px0);
    float prevBx = -1, prevBy = -1, prevAx = -1, prevAy = -1;
    for (int i = 0; i < N; ++i) {
        const double f = fMax * double(i) / double(N - 1);
        const double kk = 2.0 * kPi * f / vph;
        const double d = kk * kk - kc * kc;
        double bn = 0.0, an = 0.0;
        if (d >= 0) bn = std::sqrt(d) / (kc > 0 ? kc : 1.0);
        else        an = std::sqrt(-d) / (kc > 0 ? kc : 1.0);
        const float xi = float(px0 + i);
        const float yiB = float(py1 - std::min(1.5, bn) / 1.5 * (py1 - py0));
        const float yiA = float(py1 - std::min(1.5, an) / 1.5 * (py1 - py0));
        if (i > 0) {
            c.line(prevBx, prevBy, xi, yiB, 80, 220, 120, 1.5f);
            c.line(prevAx, prevAy, xi, yiA, 255, 140, 60, 1.5f);
        }
        prevBx = xi; prevBy = yiB;
        prevAx = xi; prevAy = yiA;
    }
    // Cursor at current frequency
    const double frac = std::min(1.0, ctx.freqHz / fMax);
    const int cx = int(px0 + frac * (px1 - px0));
    c.line(float(cx), float(py0), float(cx), float(py1), 255, 255, 120, 1.0f);

    c.text(px0 + 8, py0 + 6,  "green = beta/k_c (propagating)", 80, 220, 120);
    c.text(px0 + 8, py0 + 18, "orange= alpha/k_c (evanescent)", 255, 140, 60);
}

// Physical span of a cut plane in model coordinates. Kept in one place so the
// size planner and the drawing agree on the aspect ratio.
void sectionSpan(const ExportContext& ctx, int plane,
                 double& uMin, double& uMax, double& vMin, double& vMax) {
    const Bounds tb = (ctx.geometry == 0) ? ctx.rect->bounds() : ctx.cyl->bounds();
    const double rad = (ctx.geometry == 1) ? ctx.cyl->radius() : 0.0;
    if (plane == 0) {
        uMin = (ctx.geometry == 0) ? 0.0 : -rad;
        uMax = (ctx.geometry == 0) ? double(tb.width) : rad;
        vMin = (ctx.geometry == 0) ? 0.0 : -rad;
        vMax = (ctx.geometry == 0) ? double(tb.height) : rad;
    } else if (plane == 1) {
        uMin = 0.0; uMax = double(tb.depth);
        vMin = (ctx.geometry == 0) ? 0.0 : -rad;
        vMax = (ctx.geometry == 0) ? double(tb.width) : rad;
    } else {
        uMin = 0.0; uMax = double(tb.depth);
        vMin = (ctx.geometry == 0) ? 0.0 : -rad;
        vMax = (ctx.geometry == 0) ? double(tb.height) : rad;
    }
}

// Plot margins, in pixels, reserved around the data area for the axes.
constexpr int kMargL = 62, kMargR = 16, kMargT = 26, kMargB = 42;

// "Nice" tick step for a span: 1, 2, 5 x 10^n, aiming at ~`want` ticks.
double niceStep(double span, int want) {
    if (span <= 0.0) return 1.0;
    const double raw = span / std::max(1, want);
    const double mag = std::pow(10.0, std::floor(std::log10(raw)));
    const double n = raw / mag;
    double s = 1.0;
    if      (n < 1.5) s = 1.0;
    else if (n < 3.5) s = 2.0;
    else if (n < 7.5) s = 5.0;
    else              s = 10.0;
    return s * mag;
}

// Axis ruler in millimetres along the bottom (u) and left (v) edges of the data
// area. Ticks land on round values and ALWAYS include the coordinate origin
// when it falls inside the span -- z = 0 is the mouth of the guide and r = 0 is
// the cylinder axis, and both are the references a reader measures against.
void drawAxes(Canvas& c, int x0, int y0, int x1, int y1,
              double uMin, double uMax, double vMin, double vMax,
              const char* uLabel, const char* vLabel) {
    const uint8_t ar = 170, ag = 170, ab = 180;
    c.line(float(x0), float(y1), float(x1), float(y1), ar, ag, ab, 1.0f);
    c.line(float(x0), float(y0), float(x0), float(y1), ar, ag, ab, 1.0f);

    const double u0mm = uMin * 1e3, u1mm = uMax * 1e3;
    const double v0mm = vMin * 1e3, v1mm = vMax * 1e3;
    const double du = niceStep(u1mm - u0mm, 6), dv = niceStep(v1mm - v0mm, 4);
    char buf[32];

    for (double t = std::ceil(u0mm / du) * du; t <= u1mm + 1e-9; t += du) {
        const int px = x0 + int((t - u0mm) / (u1mm - u0mm) * (x1 - x0) + 0.5);
        const bool origin = std::fabs(t) < 1e-9;
        c.line(float(px), float(y1), float(px), float(y1 + (origin ? 8 : 5)),
               ar, ag, ab, 1.0f);
        std::snprintf(buf, sizeof(buf), "%g", t);
        c.text(px - 8, y1 + 10, buf, origin ? 255 : 200, origin ? 255 : 200,
               origin ? 255 : 210);
    }
    for (double t = std::ceil(v0mm / dv) * dv; t <= v1mm + 1e-9; t += dv) {
        const int py = y1 - int((t - v0mm) / (v1mm - v0mm) * (y1 - y0) + 0.5);
        const bool origin = std::fabs(t) < 1e-9;
        c.line(float(x0 - (origin ? 8 : 5)), float(py), float(x0), float(py),
               ar, ag, ab, 1.0f);
        std::snprintf(buf, sizeof(buf), "%g", t);
        c.text(x0 - 40, py - 4, buf, origin ? 255 : 200, origin ? 255 : 200,
               origin ? 255 : 210);
    }
    c.text((x0 + x1) / 2 - 24, y1 + 24, uLabel, 220, 220, 230);
    c.text(4, y0 - 14, vLabel, 220, 220, 230);
}

// Cross section plane sampler — shared draw helper for XY/ZX/ZY.
// plane: 0 XY (z = L/2), 1 ZX (y = mid), 2 ZY (x = mid)
// W/H are the FULL canvas size; the data area is what is left after the axis
// margins, so the caller controls how much of the figure the plot fills.
void drawSection(Canvas& c, const ExportContext& ctx, int plane, double phase,
                 int W, int H) {
    const Bounds tb = (ctx.geometry == 0) ? ctx.rect->bounds() : ctx.cyl->bounds();
    const double rad = (ctx.geometry == 1) ? ctx.cyl->radius() : 0.0;

    double uMin, uMax, vMin, vMax;
    sectionSpan(ctx, plane, uMin, uMax, vMin, vMax);
    const double Uspan = uMax - uMin;
    const double Vspan = vMax - vMin;

    c.resize(W, H, 15, 15, 20);

    const char* titles[3] = {
        "Corte XY (z = L/2)",
        "Corte ZX (y = meio)",
        "Corte ZY (x = meio)"
    };
    c.text(kMargL, 6, titles[plane], 255, 255, 255);

    const int x0 = kMargL, y0 = kMargT;
    const int x1 = W - kMargR, y1 = H - kMargB;
    // Outline
    if (plane == 0 && ctx.geometry == 1) {
        c.circle((x0+x1)/2, (y0+y1)/2, std::min(x1-x0, y1-y0)/2 - 1, 180,180,180);
    } else {
        c.rect(x0, y0, x1, y1, 180, 180, 180);
    }

    // Axis names follow the plane: the horizontal axis of ZX/ZY is the guide
    // axis z, whose zero sits at the mouth of the guide.
    const char* uL = (plane == 0) ? "x (mm)" : "z (mm)";
    const char* vL = (plane == 0) ? "y (mm)" : (plane == 1 ? "x (mm)" : "y (mm)");
    drawAxes(c, x0, y0, x1, y1, uMin, uMax, vMin, vMax, uL, vL);

    auto sampleAt = [&](double up, double vp, double ph) -> std::array<double,3> {
        double xp = 0, yp = 0, zp = 0;
        if (plane == 0) {
            xp = up; yp = vp; zp = double(tb.depth) * 0.5;
        } else if (plane == 1) {
            zp = up; xp = vp;
            yp = (ctx.geometry == 0) ? double(tb.height) * 0.5 : 0.0;
        } else {
            zp = up; yp = vp;
            xp = (ctx.geometry == 0) ? double(tb.width) * 0.5 : 0.0;
        }
        if (ctx.geometry == 1 && plane == 0 && xp*xp + yp*yp > rad*rad)
            return {0.0, 0.0, 0.0};
        return (ctx.geometry == 0)
            ? ctx.rect->fieldVector(xp, yp, zp, ph)
            : ctx.cyl ->fieldVector(xp, yp, zp, ph);
    };

    // Arrow grid follows the data area so a bigger figure gets more arrows
    // rather than the same 20x20 blown up into sparse giants. ~26 px per cell.
    const int NU = std::clamp((x1 - x0) / 26, 10, 60);
    const int NV = std::clamp((y1 - y0) / 26, 6,  40);

    // Phase-independent reference
    double ref = 1e-30;
    for (int pp = 0; pp < 12; ++pp) {
        const double ph = 2.0 * kPi * pp / 12.0;
        for (int iv = 0; iv < NV; ++iv)
        for (int iu = 0; iu < NU; ++iu) {
            const double au = (iu + 0.5) / double(NU);
            const double av = (iv + 0.5) / double(NV);
            const double up = uMin + au * Uspan;
            const double vp = vMin + av * Vspan;
            auto V = sampleAt(up, vp, ph);
            double a, b;
            if (plane == 0)      { a = V[0]; b = V[1]; }
            else if (plane == 1) { a = V[2]; b = V[0]; }
            else                 { a = V[2]; b = V[1]; }
            const double mg = std::sqrt(a*a + b*b);
            if (mg > ref) ref = mg;
        }
    }

    const float plotW = float(x1 - x0);
    const float plotH = float(y1 - y0);
    const float cellW = plotW / NU, cellH = plotH / NV;
    const float arrowMax = 0.9f * std::min(cellW, cellH);

    // Colour reference: the same global peak the 3D cloud and the colour bar
    // divide by, so a hue means one value across the whole document. `ref`
    // above stays the IN-PLANE maximum and only sets arrow/head LENGTH.
    const FieldSource* srcC = (ctx.geometry == 0)
        ? static_cast<const FieldSource*>(ctx.rect)
        : static_cast<const FieldSource*>(ctx.cyl);
    const double colorRef = std::max(1e-30, srcC->peakField());
    auto planeOf = [&](const std::array<double,3>& V, double& a, double& b) {
        if (plane == 0)      { a = V[0]; b = V[1]; }
        else if (plane == 1) { a = V[2]; b = V[0]; }
        else                 { a = V[2]; b = V[1]; }
    };
    auto toScreenUV = [&](double u, double v, float& sx, float& sy) {
        sx = x0 + float((u - uMin) / Uspan * plotW);
        sy = y0 + float((1.0 - (v - vMin) / Vspan) * plotH);
    };

    // ---- Streamline mode (mirrors the on-screen "Field lines" checkbox) ----
    // Same construction as the live plot: a fixed seed lattice and a bounded
    // arclength, so each line is an independent function of the field. The
    // animated dashes are dropped -- on a still figure they would read as
    // broken lines, and direction is already carried by the arrowheads.
    if (ctx.sectionStreamlines) {
        const double diag = std::sqrt(Uspan*Uspan + Vspan*Vspan);
        const double hStep = diag / 300.0;
        const double halfArc = 0.16 * diag;
        const int maxSteps = std::max(6, int(halfArc / hStep));

        double instMax = 1e-30;
        for (int iv = 0; iv < NV; ++iv)
        for (int iu = 0; iu < NU; ++iu) {
            double a, b;
            planeOf(sampleAt(uMin + (iu+0.5)/double(NU)*Uspan,
                             vMin + (iv+0.5)/double(NV)*Vspan, phase), a, b);
            instMax = std::max(instMax, std::hypot(a, b));
        }
        const double eps = 1e-4 * instMax;

        const float seedPx =
            std::clamp(46.0f / std::max(0.25f, ctx.sectionLineDensity), 16.0f, 160.0f);
        const int NSU = std::max(2, int(plotW / seedPx));
        const int NSV = std::max(2, int(plotH / seedPx));
        const float headEvery = 42.0f, headMax = 8.0f;

        for (int j = 0; j < NSV; ++j)
        for (int i = 0; i < NSU; ++i) {
            const double su = uMin + (i + 0.5)/double(NSU) * Uspan;
            const double sv = vMin + (j + 0.5)/double(NSV) * Vspan;
            if (ctx.geometry == 1 && plane == 0 && su*su + sv*sv > rad*rad) continue;
            { double a, b; planeOf(sampleAt(su, sv, phase), a, b);
              if (std::hypot(a, b) < eps) continue; }

            for (int dir = -1; dir <= 1; dir += 2) {
                double u = su, v = sv;
                float px, py; toScreenUV(u, v, px, py);
                float sinceHead = headEvery * 0.5f;
                for (int s = 0; s < maxSteps; ++s) {
                    std::array<double,3> V = sampleAt(u, v, phase);
                    double a, b; planeOf(V, a, b);
                    const double m = std::hypot(a, b);
                    if (m < eps) break;
                    double a2, b2;
                    planeOf(sampleAt(u + dir*(a/m)*hStep*0.5,
                                     v + dir*(b/m)*hStep*0.5, phase), a2, b2);
                    const double m2 = std::hypot(a2, b2);
                    if (m2 < eps) break;
                    const double nu = u + dir*(a2/m2)*hStep;
                    const double nv = v + dir*(b2/m2)*hStep;
                    if (nu < uMin || nu > uMax || nv < vMin || nv > vMax) break;
                    if (ctx.geometry == 1 && plane == 0 && nu*nu + nv*nv > rad*rad) break;
                    float qx, qy; toScreenUV(nu, nv, qx, qy);

                    const double m3 = std::sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
                    float rr, gg, bb;
                    fireMap(float(std::min(1.0, m3 / colorRef)), rr, gg, bb);
                    const uint8_t R = uint8_t(rr*255), G = uint8_t(gg*255), B = uint8_t(bb*255);
                    c.line(px, py, qx, qy, R, G, B, 1.6f);

                    // Heads point along +field, so trace the backward half with
                    // the arrow reversed or it would claim the wrong direction.
                    const float ddx = (dir > 0) ? (qx - px) : (px - qx);
                    const float ddy = (dir > 0) ? (qy - py) : (py - qy);
                    const float seg = std::sqrt(ddx*ddx + ddy*ddy);
                    sinceHead += seg;
                    if (sinceHead >= headEvery && seg > 1e-3f) {
                        sinceHead = 0.0f;
                        const float hx = (dir > 0) ? qx : px, hy = (dir > 0) ? qy : py;
                        const float ux = ddx/seg, uy = ddy/seg;
                        const float nx = -uy, ny = ux;
                        const float hl = headMax *
                            (0.30f + 0.70f * std::sqrt(float(std::min(1.0, m / ref))));
                        c.line(hx, hy, hx - ux*hl + nx*hl*0.55f,
                                       hy - uy*hl + ny*hl*0.55f, R, G, B, 1.6f);
                        c.line(hx, hy, hx - ux*hl - nx*hl*0.55f,
                                       hy - uy*hl - ny*hl*0.55f, R, G, B, 1.6f);
                    }
                    px = qx; py = qy; u = nu; v = nv;
                }
            }
        }
        return;
    }

    for (int iv = 0; iv < NV; ++iv) {
        for (int iu = 0; iu < NU; ++iu) {
            const double au = (iu + 0.5) / double(NU);
            const double av = (iv + 0.5) / double(NV);
            const double up = uMin + au * Uspan;
            const double vp = vMin + av * Vspan;
            auto V = sampleAt(up, vp, phase);
            double a, b;
            if (plane == 0)      { a = V[0]; b = V[1]; }
            else if (plane == 1) { a = V[2]; b = V[0]; }
            else                 { a = V[2]; b = V[1]; }
            const double mg = std::sqrt(a*a + b*b);
            const double nrm = mg / ref;
            if (nrm < 0.02) continue;
            const float cx = x0 + float(au * plotW);
            const float cy = y0 + float((1.0 - av) * plotH);
            const float scale = float(std::sqrt(nrm) / nrm) * arrowMax;
            const float dx = float(a / ref) * scale;
            const float dy = -float(b / ref) * scale;
            const float ax = cx - dx * 0.5f, ay = cy - dy * 0.5f;
            const float bx = cx + dx * 0.5f, by = cy + dy * 0.5f;
            // Length from the in-plane pair (readability), colour from the FULL
            // |field| on the global scale, exactly as the live cut does -- so a
            // hue matches the 3D figure and the colour bar in the same report.
            const double m3 = std::sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
            float rr, gg, bb;
            fireMap(float(std::min(1.0, m3 / colorRef)), rr, gg, bb);
            const uint8_t R = uint8_t(rr*255), G = uint8_t(gg*255), B = uint8_t(bb*255);
            c.line(ax, ay, bx, by, R, G, B, 1.5f);
            // arrowhead
            const float len = std::sqrt(dx*dx + dy*dy);
            if (len > 1e-3f) {
                const float ux = dx/len, uy = dy/len;
                const float px = -uy, py = ux;
                const float h = 3.0f;
                c.line(bx, by, bx - ux*h + px*h*0.5f, by - uy*h + py*h*0.5f, R,G,B, 1.0f);
                c.line(bx, by, bx - ux*h - px*h*0.5f, by - uy*h - py*h*0.5f, R,G,B, 1.0f);
            }
        }
    }
}

// ---------- 3D field lines ----------
// The live view draws these through ImGui's draw list, which the exporter has
// no access to, so the trace is reproduced here against the software canvas.
// Seeding is a fixed lattice with a bounded arclength (the same choice made for
// the on-screen cross sections): each line depends only on the field, so the
// figure is reproducible and does not depend on tracing order.
void drawFieldLines3D(Canvas& c, const ExportContext& ctx,
                      const glm::mat4& view, const glm::mat4& proj,
                      double phase, float density)
{
    const FieldSource* src = (ctx.geometry == 0)
        ? static_cast<const FieldSource*>(ctx.rect)
        : static_cast<const FieldSource*>(ctx.cyl);
    double x0, x1, y0, y1, z0, z1;
    src->domain(x0, x1, y0, y1, z0, z1);
    const double xc = 0.5*(x0+x1), yc = 0.5*(y0+y1), zc = 0.5*(z0+z1);
    const double peak = src->peakField();
    if (peak <= 0.0) return;

    const glm::vec4 vp(0.0f, 0.0f, float(c.w), float(c.h));
    auto project = [&](double x, double y, double z, float& sx, float& sy) -> bool {
        const glm::vec3 w(float(x - xc), float(y - yc), float(z - zc));
        const glm::vec3 s = glm::project(w, view, proj, vp);
        if (s.z < 0.0f || s.z > 1.0f) return false;
        sx = s.x; sy = float(c.h) - s.y;    // canvas y grows downward
        return true;
    };
    auto fieldAt = [&](double x, double y, double z) {
        return (ctx.geometry == 0) ? ctx.rect->fieldVector(x, y, z, phase)
                                   : ctx.cyl ->fieldVector(x, y, z, phase);
    };

    const double dx = x1-x0, dy = y1-y0, dz = z1-z0;
    const double diag = std::sqrt(dx*dx + dy*dy + dz*dz);
    const double h = diag / 260.0;
    const double halfArc = 0.16 * diag;
    const int maxSteps = std::max(6, int(halfArc / h));
    const double eps = 1e-4 * peak;

    const float d = std::clamp(density, 0.3f, 2.5f);
    const int NX = std::max(2, int(7 * d)), NY = std::max(2, int(5 * d)),
              NZ = std::max(3, int(18 * d));

    for (int iz = 0; iz < NZ; ++iz)
    for (int iy = 0; iy < NY; ++iy)
    for (int ix = 0; ix < NX; ++ix) {
        const double sx = x0 + (ix + 0.5) / NX * dx;
        const double sy = y0 + (iy + 0.5) / NY * dy;
        const double sz = z0 + (iz + 0.5) / NZ * dz;
        if (!src->inside(sx, sy, sz)) continue;
        for (int dir = -1; dir <= 1; dir += 2) {
            double x = sx, y = sy, z = sz;
            float px = 0, py = 0;
            bool have = project(x, y, z, px, py);
            for (int s = 0; s < maxSteps; ++s) {
                auto F = fieldAt(x, y, z);
                double m = std::sqrt(F[0]*F[0] + F[1]*F[1] + F[2]*F[2]);
                if (m < eps) break;
                auto F2 = fieldAt(x + dir*F[0]/m*h*0.5,
                                  y + dir*F[1]/m*h*0.5,
                                  z + dir*F[2]/m*h*0.5);
                double m2 = std::sqrt(F2[0]*F2[0] + F2[1]*F2[1] + F2[2]*F2[2]);
                if (m2 < eps) break;
                const double nx = x + dir*F2[0]/m2*h;
                const double ny = y + dir*F2[1]/m2*h;
                const double nz = z + dir*F2[2]/m2*h;
                if (!src->inside(nx, ny, nz)) break;
                float qx = 0, qy = 0;
                const bool ok = project(nx, ny, nz, qx, qy);
                if (have && ok) {
                    float r, g, b;
                    fireMap(float(std::min(1.0, m / peak)), r, g, b);
                    c.line(px, py, qx, qy,
                           uint8_t(r*255), uint8_t(g*255), uint8_t(b*255), 1.3f);
                }
                px = qx; py = qy; have = ok;
                x = nx; y = ny; z = nz;
            }
        }
    }
}

// X/Y/Z orientation gizmo, bottom-left. Uses only the view ROTATION, so the
// arms show the axis directions with honest foreshortening but a fixed on-screen
// size, independent of how far the camera ended up.
void drawGizmo3D(Canvas& c, const glm::mat4& view) {
    const glm::vec3 right(view[0][0], view[1][0], view[2][0]);
    const glm::vec3 up   (view[0][1], view[1][1], view[2][1]);
    const float L = 46.0f;
    const float ox = 58.0f, oy = float(c.h) - 52.0f;
    const glm::vec3 axes[3] = {{1,0,0},{0,1,0},{0,0,1}};
    // Darkened from the on-screen palette: these sit on the near-white
    // background of the 3D render, where a light green or sky blue washes out.
    const uint8_t col[3][3] = {{200,30,30},{20,140,50},{30,80,190}};
    const char* nm[3] = {"X","Y","Z"};
    const float ts = 1.8f;
    for (int i = 0; i < 3; ++i) {
        // Canvas y grows downward, hence the sign flip on the up component.
        const float dx = glm::dot(axes[i], right) * L;
        const float dy = -glm::dot(axes[i], up) * L;
        c.line(ox, oy, ox + dx, oy + dy, col[i][0], col[i][1], col[i][2], 2.5f);
        int tx = int(ox + dx * 1.22f) - Canvas::textW(nm[i], ts) / 2;
        int ty = int(oy + dy * 1.22f) - Canvas::textH(ts) / 2;
        tx = std::clamp(tx, 2, std::max(2, c.w - Canvas::textW(nm[i], ts) - 2));
        ty = std::clamp(ty, 2, std::max(2, c.h - Canvas::textH(ts) - 2));
        c.text(tx, ty, nm[i], col[i][0], col[i][1], col[i][2], ts);
    }
}

// Edge dimensions, the same ones the live view annotates. Each label sits at
// the projected midpoint of one representative edge of the bounding box, with
// the edge itself highlighted so the number is unambiguous about WHICH span it
// measures.
void drawDimensions3D(Canvas& c, const ExportContext& ctx,
                      const glm::mat4& view, const glm::mat4& proj) {
    const Bounds tb = (ctx.geometry == 0) ? ctx.rect->bounds() : ctx.cyl->bounds();
    const float hw = tb.width * 0.5f, hh = tb.height * 0.5f, hd = tb.depth * 0.5f;
    const glm::vec4 vp(0.0f, 0.0f, float(c.w), float(c.h));
    auto toScreen = [&](const glm::vec3& w, float& sx, float& sy) -> bool {
        const glm::vec3 s = glm::project(w, view, proj, vp);
        if (s.z < 0.0f || s.z > 1.0f) return false;
        sx = s.x; sy = float(c.h) - s.y;
        return true;
    };

    // Which edge carries which label. The export camera is fixed (yaw 35,
    // pitch 25) with the eye at +x/+y/+z, so on screen +x runs right, +z runs
    // LEFT and -z is the far ("back") side. Each span therefore sits on the
    // back-facing edge that is not hidden by the body:
    //   a -> top-back edge      (y = +hh, z = -hd), label above
    //   b -> right-back edge    (x = +hw, z = -hd), label to the right
    //   L -> bottom-right edge  (x = +hw, y = -hh), label below
    // align: -1 puts the text to the right of the anchor, 0 centres it.
    struct Dim { glm::vec3 a, b; std::string text; int ox, oy, align; };
    std::vector<Dim> dims;
    char buf[64];
    if (ctx.geometry == 0) {
        std::snprintf(buf, sizeof(buf), "a = %.2f mm", tb.width * 1e3);
        dims.push_back({{-hw,  hh, -hd}, { hw,  hh, -hd}, buf,   0, -30,  0});
        std::snprintf(buf, sizeof(buf), "b = %.2f mm", tb.height * 1e3);
        dims.push_back({{ hw, -hh, -hd}, { hw,  hh, -hd}, buf,  16,  -6, -1});
        std::snprintf(buf, sizeof(buf), "L = %.1f mm", tb.depth * 1e3);
        dims.push_back({{ hw, -hh, -hd}, { hw, -hh,  hd}, buf,   0,  16,  0});
    } else {
        std::snprintf(buf, sizeof(buf), "R = %.2f mm", ctx.cyl->radius() * 1e3);
        dims.push_back({{0.0f, 0.0f, -hd}, { hw, 0.0f, -hd}, buf,  16, -6, -1});
        std::snprintf(buf, sizeof(buf), "L = %.1f mm", tb.depth * 1e3);
        dims.push_back({{ hw, -hh, -hd}, { hw, -hh,  hd}, buf,   0, 16,  0});
    }

    // The 3D view clears to a near-white background, so the old pale-yellow
    // leader lines and cream text were nearly invisible on it.
    const uint8_t lr = 25, lg = 25, lb = 30;
    const float ts = 2.2f;                    // font scale for the labels

    for (const Dim& d : dims) {
        float ax, ay, bx, by;
        if (!toScreen(d.a, ax, ay) || !toScreen(d.b, bx, by)) continue;
        c.line(ax, ay, bx, by, lr, lg, lb, 2.0f);
        // End caps make the measured span read as a dimension, not an edge.
        for (int e = 0; e < 2; ++e) {
            const float px = e ? bx : ax, py = e ? by : ay;
            c.line(px - 4, py - 4, px + 4, py + 4, lr, lg, lb, 1.5f);
            c.line(px - 4, py + 4, px + 4, py - 4, lr, lg, lb, 1.5f);
        }
        const int tw = Canvas::textW(d.text.c_str(), ts);
        const int th = Canvas::textH(ts);
        int tx = int(0.5f * (ax + bx)) + d.ox - (d.align == 0 ? tw / 2 : 0);
        int ty = int(0.5f * (ay + by)) + d.oy;
        // Never let a label run off the page: the figure has to contain the
        // number that explains it, and a label clipped by the PNG edge is worse
        // than one nudged a few pixels inward.
        tx = std::clamp(tx, 3, std::max(3, c.w - tw - 3));
        ty = std::clamp(ty, 3, std::max(3, c.h - th - 3));
        // Opaque plate behind the text. Needed because the clamp above can push
        // a label back over the guide when there is not enough margin for it
        // (a 2.2x "b = 10.16 mm" is ~158 px wide, more than the border), and a
        // number sitting on the field cloud is unreadable without it.
        c.fillRect(tx - 5, ty - 4, tx + tw + 5, ty + th + 4, 255, 255, 255, 235);
        c.text(tx, ty, d.text.c_str(), 0, 0, 0, ts);
    }
}

// One 3D figure, honouring the "Modo 3D" the user picked on screen.
//   0 = intensity cloud   1 = field lines   2 = both
// For lines-only the offscreen pass still runs, with an empty particle set, so
// the guide outline and floor grid stay as the spatial reference.
void renderScene3D(Canvas& sub, const ExportContext& ctx, int w, int h,
                   const glm::mat4& view, const glm::mat4& proj,
                   float sphereRadius, double phase, float lineDensity)
{
    (void)sphereRadius;
    Renderer::OffscreenOpts oo;
    oo.showGrid     = true;
    // A cylindrical guide gets its own translucent shell; the wireframe box is
    // the rectangular outline and was being drawn for BOTH before.
    oo.showBox      = (ctx.geometry == 0);
    oo.showCylinder = (ctx.geometry == 1);
    // Intensity comes from the same baked cloud the UI displays, animated by
    // the same uPhase, normalized by the same 1/peak and shaded by the same
    // fire palette as the colour bar. The old path drew instanced solid
    // spheres coloured by the model's internal black->purple->white map, which
    // matched neither the screen nor the legend.
    oo.useCloud        = (ctx.view3D != 1);
    oo.cloudPhase      = float(phase);
    oo.cloudInvPeak    = ctx.cloudInvPeak;
    oo.cloudOpaque     = ctx.cloudOpaque;
    // Same expression as the live view: dot size follows the mean spacing of
    // the accepted samples, scaled into pixels through the projection.
    oo.cloudPointScale = (ctx.cloudMeanSpacing * 1.5f) * proj[1][1] * float(h) * 0.5f;

    std::vector<uint8_t> rgba;
    ctx.renderer->renderOffscreen(w, h, view, proj, oo, rgba);
    sub.w = w; sub.h = h; sub.px = std::move(rgba);
    if (ctx.view3D != 0)
        drawFieldLines3D(sub, ctx, view, proj, phase, lineDensity);
    // Overlays last, so neither the cloud nor the field lines can bury the
    // scale reference or the orientation.
    drawDimensions3D(sub, ctx, view, proj);
    drawGizmo3D(sub, view);
}

// ---------- Item plan ----------
struct Item {
    enum Kind { Scene3D, XY, ZX, ZY, Spectrum, ColorBar } kind;
    int w = 0, h = 0;    // per-item canvas size
    int animated = 0;    // 1 if its content depends on phase
};

std::vector<Item> planItems(const ExportContext& ctx, const ExportOptions& o) {
    std::vector<Item> its;
    if (o.incScene3D)  its.push_back({Item::Scene3D, 0, 0, 1});
    if (o.incXY)       its.push_back({Item::XY,      0, 0, 1});
    if (o.incZX)       its.push_back({Item::ZX,      0, 0, 1});
    if (o.incZY)       its.push_back({Item::ZY,      0, 0, 1});
    if (o.incSpectrum) its.push_back({Item::Spectrum,0, 0, 0});
    // No standalone ColorBar item: the scale is now attached to the side of
    // every field figure instead of standing alone as its own plate.

    // The scale strip is appended to every field figure afterwards, so reserve
    // its width here or the finished figure would overrun the requested box.
    const int barPane = o.incColorBar ? kBarPaneW : 0;
    const int BW = std::max(320, o.imageW) - barPane;
    const int BH = std::max(240, o.imageH);

    // Fit the DATA AREA to the box while keeping the plane's true aspect ratio,
    // then add the axis margins back. Whichever of width/height binds first is
    // driven to the box edge, so a long ZX cut spans the full width and a square
    // XY cut spans the full height -- each uses as much of the figure as its own
    // proportions allow, instead of a fixed 520 px that wasted most of both.
    auto sectSize = [&](int plane, int& w, int& h) {
        double uMin, uMax, vMin, vMax;
        sectionSpan(ctx, plane, uMin, uMax, vMin, vMax);
        const double Us = std::max(1e-12, uMax - uMin);
        const double Vs = std::max(1e-12, vMax - vMin);
        const double availW = BW - (kMargL + kMargR);
        const double availH = BH - (kMargT + kMargB);
        const double scale = std::min(availW / Us, availH / Vs);
        w = int(Us * scale + 0.5) + kMargL + kMargR;
        h = int(Vs * scale + 0.5) + kMargT + kMargB;
        // A very long guide would otherwise leave a data strip a few px tall.
        const int minData = 150;
        if (h - (kMargT + kMargB) < minData) h = minData + kMargT + kMargB;
    };

    for (auto& it : its) {
        switch (it.kind) {
            case Item::Scene3D:  it.w = BW; it.h = BH; break;
            case Item::XY:       sectSize(0, it.w, it.h); break;
            case Item::ZX:       sectSize(1, it.w, it.h); break;
            case Item::ZY:       sectSize(2, it.w, it.h); break;
            case Item::Spectrum: it.w = std::min(BW, 900);
                                 it.h = std::min(BH, 420); break;
            default:             it.w = std::min(BW, 260);
                                 it.h = std::min(BH, 420); break;
        }
    }
    return its;
}

// Render all items into one stacked canvas for the given phase.
void composeFrame(Canvas& out,
                  const ExportContext& ctx,
                  const std::vector<Item>& items,
                  const std::vector<int>& yOffsets,
                  int pageW, int pageH,
                  double phase,
                  const glm::mat4& view3d,
                  const glm::mat4& proj3d,
                  float sphereRadius3d,
                  const std::vector<Particle>& particles3d,
                  bool showFloor3d,
                  bool attachScale)
{
    (void)particles3d;
    out.resize(pageW, pageH, 8, 8, 12);

    for (size_t i = 0; i < items.size(); ++i) {
        const Item& it = items[i];
        if (yOffsets[i] < 0) continue;
        const int dy = yOffsets[i];
        Canvas sub;
        if (it.kind == Item::Scene3D) {
            (void)showFloor3d;
            renderScene3D(sub, ctx, it.w, it.h, view3d, proj3d,
                          sphereRadius3d, phase, ctx.lineDensity);
        }
        else if (it.kind == Item::XY)       drawSection(sub, ctx, 0, phase, it.w, it.h);
        else if (it.kind == Item::ZX)       drawSection(sub, ctx, 1, phase, it.w, it.h);
        else if (it.kind == Item::ZY)       drawSection(sub, ctx, 2, phase, it.w, it.h);
        else if (it.kind == Item::Spectrum)   drawSpectrum(sub, ctx);
        else                                  drawColorBar(sub, ctx);
        // Field plots carry their own scale; the spectrum is a frequency plot
        // with no field colours, so it gets none.
        if (attachScale && it.kind != Item::Spectrum && it.kind != Item::ColorBar)
            attachColorBar(sub, ctx);
        const int dx = (pageW - sub.w) / 2;
        blit(out, dx, dy, sub);
    }
}

// ==================== LaTeX project ====================

std::string num(double v, int prec = 3) {
    char b[64];
    std::snprintf(b, sizeof(b), "%.*g", prec, v);
    return std::string(b);
}

const char* kModeTypeName[2] = {"TE", "TM"};

// "TE$_{11}$", or "TE$_{110}$" for a cavity (the axial index only exists there).
std::string modeLabel(const ExportContext& ctx) {
    std::ostringstream s;
    s << kModeTypeName[ctx.modeType == 0 ? 0 : 1] << "$_{"
      << ctx.modeM << ctx.modeN;
    if (ctx.structure == 1) s << ctx.modeL;
    s << "}$";
    return s.str();
}

struct Derived {
    double fcHz = 0.0;       // cutoff (guide) or resonance (cavity)
    double kc = 0.0;
    double beta = 0.0;       // real part; 0 below cutoff
    double alpha = 0.0;      // decay constant when evanescent
    double lambdaG = 0.0;    // guide wavelength (m), 0 below cutoff
    double Zw = 0.0;         // wave impedance (ohm), 0 below cutoff
    bool   propagating = false;
    double vPhaseMedium = 0.0;
};

Derived derive(const ExportContext& ctx) {
    Derived d;
    const double c0 = 299792458.0, mu0 = 4.0 * kPi * 1e-7, eps0 = 8.8541878128e-12;
    const FieldSource* src = (ctx.geometry == 0)
        ? static_cast<const FieldSource*>(ctx.rect)
        : static_cast<const FieldSource*>(ctx.cyl);
    const double er = src->epsilonRel(), mr = src->muRel();
    d.vPhaseMedium = c0 / std::sqrt(er * mr);
    d.kc = src->cutoffWavenumber();
    d.fcHz = src->resonantFrequency();
    const double k = 2.0 * kPi * ctx.freqHz / d.vPhaseMedium;
    const double disc = k * k - d.kc * d.kc;
    if (disc > 0.0) {
        d.propagating = true;
        d.beta = std::sqrt(disc);
        d.lambdaG = 2.0 * kPi / d.beta;
        const double omega = 2.0 * kPi * ctx.freqHz;
        d.Zw = (ctx.modeType == 0) ? omega * mu0 * mr / d.beta
                                   : d.beta / (omega * eps0 * er);
    } else {
        d.alpha = std::sqrt(-disc);
    }
    return d;
}

// Every mode of the active geometry whose cutoff sits below fMax, sorted.
// This is the table that tells you whether the guide is single-mode at the
// operating point -- the thing that decides whether a measured or simulated
// field can be compared against one analytic mode at all.
struct ModeRow { double fcHz; std::string name; bool propagates; };

std::vector<ModeRow> modeSpectrum(const ExportContext& ctx, double fMax) {
    std::vector<ModeRow> rows;
    const double c0 = 299792458.0;
    const FieldSource* src = (ctx.geometry == 0)
        ? static_cast<const FieldSource*>(ctx.rect)
        : static_cast<const FieldSource*>(ctx.cyl);
    const double v = c0 / std::sqrt(src->epsilonRel() * src->muRel());
    if (ctx.geometry == 0) {
        const Bounds b = ctx.rect->bounds();
        const double a = b.width, h = b.height;
        for (int m = 0; m <= 4; ++m)
        for (int n = 0; n <= 4; ++n) {
            if (m == 0 && n == 0) continue;
            const double kc = std::sqrt(std::pow(m * kPi / a, 2.0) +
                                        std::pow(n * kPi / h, 2.0));
            const double fc = kc * v / (2.0 * kPi);
            if (fc > fMax) continue;
            std::ostringstream s; s << "TE$_{" << m << n << "}$";
            rows.push_back({fc, s.str(), fc < ctx.freqHz});
            if (m >= 1 && n >= 1) {
                std::ostringstream t; t << "TM$_{" << m << n << "}$";
                rows.push_back({fc, t.str(), fc < ctx.freqHz});
            }
        }
    } else {
        // Same tabulated Bessel zeros the model uses.
        static const double tm[3][3] = {
            { 2.4048255577,  5.5200781103,  8.6537279129},
            { 3.8317059702,  7.0155866698, 10.1734681351},
            { 5.1356223018,  8.4172441404, 11.6198411721}};
        static const double te[3][3] = {
            { 3.8317059702,  7.0155866698, 10.1734681351},
            { 1.8411837813,  5.3314427735,  8.5363163663},
            { 3.0542369282,  6.7061331942,  9.9694678231}};
        const double R = ctx.cyl->radius();
        for (int n = 0; n < 3; ++n)
        for (int m = 0; m < 3; ++m) {
            const double fcTE = (te[n][m] / R) * v / (2.0 * kPi);
            const double fcTM = (tm[n][m] / R) * v / (2.0 * kPi);
            if (fcTE <= fMax) {
                std::ostringstream s; s << "TE$_{" << n << (m + 1) << "}$";
                rows.push_back({fcTE, s.str(), fcTE < ctx.freqHz});
            }
            if (fcTM <= fMax) {
                std::ostringstream s; s << "TM$_{" << n << (m + 1) << "}$";
                rows.push_back({fcTM, s.str(), fcTM < ctx.freqHz});
            }
        }
    }
    std::sort(rows.begin(), rows.end(),
              [](const ModeRow& a, const ModeRow& b) { return a.fcHz < b.fcHz; });
    return rows;
}

struct FigureRef { std::string file, caption, label; };

std::string buildLatex(const ExportContext& ctx, const ExportOptions& o,
                       const std::vector<FigureRef>& figs)
{
    const Derived d = derive(ctx);
    const FieldSource* src = (ctx.geometry == 0)
        ? static_cast<const FieldSource*>(ctx.rect)
        : static_cast<const FieldSource*>(ctx.cyl);
    const bool cav   = (ctx.structure == 1);
    const bool physU = src->physicalUnits();
    const char* fldSym  = (ctx.fieldKind == 0) ? "|E|" : "|H|";
    const std::string unit = physU ? ((ctx.fieldKind == 0) ? "V/m" : "A/m")
                                   : "a.u.";
    const Bounds bb = src->bounds();

    std::ostringstream t;
    t << "% Generated by WaveguideTEmn. Build with: pdflatex main.tex\n"
      << "\\documentclass[11pt,a4paper]{article}\n"
      << "\\usepackage[utf8]{inputenc}\n"
      << "\\usepackage[T1]{fontenc}\n"
      << "\\usepackage{graphicx}\n"
      << "\\usepackage{booktabs}\n"
      << "\\usepackage{amsmath}\n"
      << "\\usepackage[margin=2.5cm]{geometry}\n"
      << "\\usepackage{caption}\n"
      << "\\graphicspath{{pictures/}}\n\n";

    std::string title = o.title;
    if (title.empty()) {
        std::ostringstream ts;
        ts << (cav ? "Cavidade " : "Guia de onda ")
           << (ctx.geometry == 0 ? "retangular" : "cil\\'indrica")
           << " --- modo " << modeLabel(ctx);
        title = ts.str();
    }
    t << "\\title{" << title << "}\n";
    t << "\\author{" << (o.author.empty() ? std::string("") : o.author) << "}\n";
    t << "\\date{\\today}\n\n\\begin{document}\n\\maketitle\n\n";

    // ---- 1. Geometry and medium ----
    t << "\\section{Configura\\c{c}\\~ao da simula\\c{c}\\~ao}\n\n"
      << "\\begin{table}[h]\n\\centering\n"
      << "\\caption{Geometria e meio de preenchimento.}\n"
      << "\\label{tab:geometria}\n"
      << "\\begin{tabular}{ll}\n\\toprule\nPar\\^ametro & Valor \\\\\n\\midrule\n";
    t << "Geometria & " << (ctx.geometry == 0 ? "Retangular" : "Cil\\'indrica") << " \\\\\n";
    t << "Estrutura & " << (cav ? "Cavidade (curto em ambas as extremidades)"
                                : "Guia aberta (onda viajante)") << " \\\\\n";
    if (ctx.geometry == 0) {
        t << "Largura $a$ & " << num(bb.width * 1e3) << " mm \\\\\n"
          << "Altura $b$ & "  << num(bb.height * 1e3) << " mm \\\\\n"
          << "Comprimento $d$ & " << num(bb.depth * 1e3) << " mm \\\\\n";
    } else {
        t << "Raio $R$ & " << num(ctx.cyl->radius() * 1e3) << " mm \\\\\n"
          << "Comprimento $d$ & " << num(bb.depth * 1e3) << " mm \\\\\n";
    }
    t << "Meio & " << ctx.mediumName << " \\\\\n"
      << "$\\varepsilon_r$ & " << num(src->epsilonRel()) << " \\\\\n"
      << "$\\mu_r$ & " << num(src->muRel()) << " \\\\\n"
      << "Velocidade de fase no meio & " << num(d.vPhaseMedium / 1e8, 4)
      << "$\\times 10^{8}$ m/s \\\\\n"
      << "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n";

    // ---- 2. Mode and operating point ----
    t << "\\begin{table}[h]\n\\centering\n"
      << "\\caption{Modo e ponto de opera\\c{c}\\~ao.}\n"
      << "\\label{tab:modo}\n"
      << "\\begin{tabular}{ll}\n\\toprule\nPar\\^ametro & Valor \\\\\n\\midrule\n";
    t << "Modo & " << modeLabel(ctx) << " \\\\\n"
      << "Campo apresentado & " << (ctx.fieldKind == 0 ? "El\\'etrico $\\mathbf{E}$"
                                                       : "Magn\\'etico $\\mathbf{H}$")
      << " \\\\\n"
      << "Frequ\\^encia de opera\\c{c}\\~ao & " << num(ctx.freqHz / 1e9, 5) << " GHz \\\\\n"
      << (cav ? "Frequ\\^encia de resson\\^ancia $f_{res}$ & "
              : "Frequ\\^encia de corte $f_c$ & ")
      << num(d.fcHz / 1e9, 5) << " GHz \\\\\n"
      << "$k_c$ & " << num(d.kc, 5) << " rad/m \\\\\n";
    if (!cav) {
        t << "Regime & " << (d.propagating ? "Propagante" : "Evanescente (abaixo do corte)")
          << " \\\\\n";
        if (d.propagating) {
            t << "$\\beta$ & " << num(d.beta, 5) << " rad/m \\\\\n"
              << "Comprimento de onda guiado $\\lambda_g$ & "
              << num(d.lambdaG * 1e3, 4) << " mm \\\\\n"
              << "Imped\\^ancia de onda $Z_{" << kModeTypeName[ctx.modeType == 0 ? 0 : 1]
              << "}$ & " << num(d.Zw, 5) << " $\\Omega$ \\\\\n";
        } else {
            t << "Constante de atenua\\c{c}\\~ao $\\alpha$ & " << num(d.alpha, 5)
              << " Np/m \\\\\n"
              << "Profundidade de penetra\\c{c}\\~ao $1/\\alpha$ & "
              << num(1e3 / std::max(1e-30, d.alpha), 4) << " mm \\\\\n";
        }
    }
    t << "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n";

    // ---- 3. Amplitude scale ----
    t << "\\begin{table}[h]\n\\centering\n"
      << "\\caption{Escala de amplitude. "
      << (physU
          ? "A amplitude modal foi resolvida a partir da pot\\^encia transportada, "
            "de modo que os valores est\\~ao em unidades SI (mesma conven\\c{c}\\~ao de "
            "uma \\textit{wave port} do HFSS, cujo padr\\~ao \\'e 1 W incidente)."
          : "A amplitude n\\~ao p\\^ode ser fixada em pot\\^encia (cavidade ou modo "
            "evanescente n\\~ao transportam pot\\^encia l\\'iquida), portanto a escala "
            "\\'e arbitr\\'aria.")
      << "}\n\\label{tab:amplitude}\n"
      << "\\begin{tabular}{ll}\n\\toprule\nPar\\^ametro & Valor \\\\\n\\midrule\n";
    if (physU) {
        t << "Pot\\^encia transportada & " << num(ctx.powerW, 4) << " W \\\\\n";
        if (ctx.geometry == 0)
            t << "Amplitude modal $A$ & " << num(ctx.rect->amplitude(), 5) << " \\\\\n";
    }
    t << "Unidades & " << (physU ? "F\\'isicas (SI)" : "Arbitr\\'arias (u.a.)") << " \\\\\n"
      << "Pico de $" << fldSym << "$ & " << num(src->peakField(), 5) << " " << unit
      << " \\\\\n"
      << "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n";

    // ---- 4. Mode spectrum ----
    {
        const auto rows = modeSpectrum(ctx, ctx.freqHz * 1.6);
        int nProp = 0;
        for (const auto& r : rows) if (r.propagates) ++nProp;
        t << "\\begin{table}[h]\n\\centering\n"
          << "\\caption{Modos com corte at\\'e $1{,}6\\,f$. Em " << num(ctx.freqHz / 1e9, 5)
          << " GHz, " << nProp << " modo(s) propagam"
          << (nProp > 1
              ? ". A guia \\'e multimodo neste ponto: uma solu\\c{c}\\~ao num\\'erica "
                "excitada por porta conter\\'a uma superposi\\c{c}\\~ao, e n\\~ao este "
                "modo isolado."
              : ". A guia \\'e monomodo neste ponto.")
          << "}\n\\label{tab:espectro}\n"
          << "\\begin{tabular}{lrl}\n\\toprule\nModo & $f_c$ (GHz) & Regime \\\\\n\\midrule\n";
        for (const auto& r : rows)
            t << r.name << " & " << num(r.fcHz / 1e9, 4) << " & "
              << (r.propagates ? "propaga" : "evanescente") << " \\\\\n";
        t << "\\bottomrule\n\\end{tabular}\n\\end{table}\n\n\\clearpage\n\n";
    }

    // ---- 5. Figures ----
    t << "\\section{Resultados}\n\n";
    for (const auto& f : figs) {
        t << "\\begin{figure}[htbp]\n\\centering\n"
          << "\\includegraphics[width=\\linewidth]{" << f.file << "}\n"
          << "\\caption{" << f.caption << "}\n"
          << "\\label{fig:" << f.label << "}\n"
          << "\\end{figure}\n\n";
    }

    t << "\\end{document}\n";
    return t.str();
}

} // namespace

bool runExport(const ExportContext& ctx, const ExportOptions& options)
{
    if (!ctx.renderer) return false;

    auto items = planItems(ctx, options);
    if (items.empty()) return false;

    // If the color bar will be overlaid on the 3D scene, skip its own row.
    bool hasScene3D = false;
    for (auto& it : items) if (it.kind == Item::Scene3D) { hasScene3D = true; break; }

    // Layout: stacked vertically with 12 px padding
    const int pad = 12;
    int pageW = 0, pageH = pad;
    std::vector<int> yOff(items.size(), 0);
    for (size_t i = 0; i < items.size(); ++i) {
        if (hasScene3D && items[i].kind == Item::ColorBar) { yOff[i] = -1; continue; }
        yOff[i] = pageH;
        pageH += items[i].h + pad;
        // Field items grow by the scale strip appended to their right.
        const int wEff = items[i].w +
            ((options.incColorBar && items[i].kind != Item::Spectrum &&
              items[i].kind != Item::ColorBar) ? kBarPaneW : 0);
        if (wEff > pageW) pageW = wEff;
    }
    pageW += pad * 2;

    // Default 3D camera: orbit at (yaw=35, pitch=25)
    const Bounds tb = (ctx.geometry == 0) ? ctx.rect->bounds() : ctx.cyl->bounds();
    const float yawR = glm::radians(35.0f), pitR = glm::radians(25.0f);
    const glm::vec3 target(0.0f, 0.0f, 0.0f);
    const glm::vec3 eyeDir(std::cos(pitR) * std::sin(yawR),
                           std::sin(pitR),
                           std::cos(pitR) * std::cos(yawR));
    // Aspect must follow the ACTUAL scene size, which is now driven by the
    // requested figure box. A fixed 720x480 here stretched the geometry as soon
    // as the figure stopped being 3:2.
    int sceneW = options.imageW, sceneH = options.imageH;
    for (const Item& it : items)
        if (it.kind == Item::Scene3D) { sceneW = it.w; sceneH = it.h; break; }
    const float fovY = glm::radians(45.0f);
    const float aspect = float(sceneW) / float(std::max(1, sceneH));
    const glm::mat4 proj3d = glm::perspective(fovY, aspect, 0.001f, 1000.0f);

    // Frame the guide tightly AND centre it. Two separate problems:
    //
    //  - zoom: backing off a fixed 2.2 x maxDim ignores the other two extents
    //    and the figure's aspect ratio, so a long thin guide became a small
    //    diagonal streak in a mostly empty frame.
    //  - centring: the box is centred on the origin, but the PROJECTED
    //    silhouette of a box seen obliquely is not symmetric about the centre
    //    of projection. Aiming the camera at the box centre therefore leaves
    //    the drawing visibly off-centre; the silhouette is what has to be
    //    centred, not the object.
    //
    // Both are solved together: project the eight corners, shift the camera so
    // the silhouette's bounding box is centred, rescale the distance so it
    // fills `fill` of the half-viewport, and repeat. The two corrections
    // interact (moving the camera changes the silhouette), so it is iterated
    // rather than solved in one shot; a handful of rounds is plenty.
    const glm::vec3 fwd = -eyeDir;                    // camera looks this way
    const glm::vec3 right = glm::normalize(glm::cross(fwd, glm::vec3(0, 1, 0)));
    const glm::vec3 upv = glm::cross(right, fwd);
    const float tanY = std::tan(fovY * 0.5f);
    const float tanX = tanY * aspect;

    glm::vec3 corners[8];
    {
        const float hw = tb.width * 0.5f, hh = tb.height * 0.5f, hd = tb.depth * 0.5f;
        for (int c = 0; c < 8; ++c)
            corners[c] = glm::vec3((c & 1) ? hw : -hw,
                                   (c & 2) ? hh : -hh,
                                   (c & 4) ? hd : -hd);
    }

    // Corner coordinates in the camera's fixed basis, independent of where the
    // camera ends up: eyeDir is perpendicular to `right` and `upv`, so sliding
    // the camera sideways by (cx, cy) shifts every corner by exactly (-cx, -cy)
    // and leaves depth alone.
    float cxs[8], cys[8], czs[8];
    for (int i = 0; i < 8; ++i) {
        cxs[i] = glm::dot(corners[i], right);
        cys[i] = glm::dot(corners[i], upv);
        czs[i] = glm::dot(corners[i], fwd);
    }
    // Border left around the box. It is not just cosmetic: the dimension labels
    // are drawn OUTSIDE the outline (above the top edge, right of the right
    // edge, below the bottom one), so the frame has to reserve room for them or
    // they get pushed back inward by the anti-clipping clamp and end up sitting
    // on the guide. ~7 % per side is enough for the 2.2x labels.
    const float fill = 0.86f;
    const float tX = tanX * fill, tY = tanY * fill;

    // At a distance d, corner i sits at depth (czs[i] + d) and the containment
    // condition |x - cx| <= tX*depth becomes an INTERVAL for cx. The box fits
    // iff the intersection of the eight intervals is non-empty, in x and in y.
    // Both bounds move monotonically with d, so the smallest workable distance
    // is a clean bisection -- and the midpoint of the feasible interval is the
    // centring, so framing and centring fall out of the same solve.
    auto bounds1D = [&](float d, const float* v, float t,
                        float& lo, float& hi) {
        lo = -1e30f; hi = 1e30f;
        for (int i = 0; i < 8; ++i) {
            const float depth = czs[i] + d;
            lo = std::max(lo, v[i] - t * depth);
            hi = std::min(hi, v[i] + t * depth);
        }
    };
    float minDepth = 0.0f;
    for (int i = 0; i < 8; ++i) minDepth = std::max(minDepth, -czs[i]);
    auto feasible = [&](float d) {
        if (d <= minDepth) return false;      // a corner would be behind the eye
        float lx, hx, ly, hy;
        bounds1D(d, cxs, tX, lx, hx);
        bounds1D(d, cys, tY, ly, hy);
        return lx <= hx && ly <= hy;
    };

    float lo = minDepth + 1e-6f;
    float hi = 8.0f * std::max({tb.width, tb.height, tb.depth}) + minDepth;
    while (!feasible(hi) && hi < 1e6f) hi *= 2.0f;
    for (int i = 0; i < 60; ++i) {
        const float mid = 0.5f * (lo + hi);
        if (feasible(mid)) hi = mid; else lo = mid;
    }
    const float dist = hi;

    // Centre each axis on the SILHOUETTE, not on the feasible interval. They
    // coincide only on the axis that binds; on the slack axis the interval is
    // wide and its midpoint leaves the drawing visibly off-centre (the box sits
    // at different depths, so equal clearance in world units is not equal
    // clearance on screen). Solve instead for the offset that makes the
    // projected extent symmetric: ndc_i(c) = (v_i - c)/(t*(z_i+dist)) falls
    // monotonically with c, so max+min does too, and a bisection nails it.
    auto centreAxis = [&](const float* v, float t, float lo_, float hi_) {
        auto skew = [&](float c) {
            float mn = 1e30f, mx = -1e30f;
            for (int i = 0; i < 8; ++i) {
                const float n = (v[i] - c) / (t * (czs[i] + dist));
                mn = std::min(mn, n); mx = std::max(mx, n);
            }
            return mx + mn;
        };
        for (int i = 0; i < 50; ++i) {
            const float mid = 0.5f * (lo_ + hi_);
            if (skew(mid) > 0.0f) lo_ = mid; else hi_ = mid;
        }
        return 0.5f * (lo_ + hi_);
    };
    float lx, hx2, ly, hy2;
    bounds1D(dist, cxs, tX, lx, hx2);
    bounds1D(dist, cys, tY, ly, hy2);
    const glm::vec3 target3d = right * centreAxis(cxs, tX, lx, hx2)
                             + upv   * centreAxis(cys, tY, ly, hy2);
    const glm::vec3 eye = target3d + dist * eyeDir;
    const glm::mat4 view3d = glm::lookAt(eye, target3d, glm::vec3(0, 1, 0));
    const float sphereRadius = (tb.depth / 110.0f) * 0.9f;

    // ---------------- LaTeX project ----------------
    // Each selected item becomes its OWN png under pictures/ and its own float
    // in the document, which is what a paper or a book chapter needs: figures
    // that can be placed, sized and referenced independently. The stacked
    // single-image export below stays available for a quick look.
    if (options.latexProject) {
        namespace fs = std::filesystem;
        const fs::path root = fs::path(options.filename);
        const fs::path pics = root / "pictures";
        std::error_code ec;
        fs::create_directories(pics, ec);
        if (ec) return false;

        std::vector<FigureRef> figs;
        const bool multiPhase = options.phasesDeg.size() > 1;

        auto itemInfo = [&](Item::Kind k) -> std::pair<std::string, std::string> {
            switch (k) {
                case Item::Scene3D:  return {"cena3d",   "Distribui\\c{c}\\~ao de intensidade em 3D"};
                case Item::XY:       return {"corte_xy", "Corte transversal XY"};
                case Item::ZX:       return {"corte_zx", "Corte longitudinal ZX"};
                case Item::ZY:       return {"corte_zy", "Corte longitudinal ZY"};
                case Item::Spectrum: return {"espectro", "Espectro e frequ\\^encia de corte"};
                default:             return {"escala_cores", "Escala de cores"};
            }
        };

        for (const Item& it : items) {
            const auto info = itemInfo(it.kind);
            // Static items are written once no matter how many phases were asked
            // for -- repeating an identical figure would only pad the document.
            const std::vector<double> phases =
                it.animated ? options.phasesDeg : std::vector<double>{0.0};

            for (size_t pi = 0; pi < phases.size(); ++pi) {
                const double phase = phases[pi] * kPi / 180.0;
                Canvas sub;
                if (it.kind == Item::Scene3D)
                    renderScene3D(sub, ctx, it.w, it.h, view3d, proj3d,
                                  sphereRadius, phase, ctx.lineDensity);
                else if (it.kind == Item::XY)       drawSection(sub, ctx, 0, phase, it.w, it.h);
                else if (it.kind == Item::ZX)       drawSection(sub, ctx, 1, phase, it.w, it.h);
                else if (it.kind == Item::ZY)       drawSection(sub, ctx, 2, phase, it.w, it.h);
                else if (it.kind == Item::Spectrum) drawSpectrum(sub, ctx);
                else                                drawColorBar(sub, ctx);
                // Every field figure carries its own scale, so a figure pulled
                // into a paper on its own still explains its colours.
                if (options.incColorBar && it.kind != Item::Spectrum &&
                    it.kind != Item::ColorBar)
                    attachColorBar(sub, ctx);

                std::string stem = info.first;
                std::string cap  = info.second;
                if (it.animated && multiPhase) {
                    char sfx[32];
                    std::snprintf(sfx, sizeof(sfx), "_fase%03d", int(phases[pi] + 0.5));
                    stem += sfx;
                    cap  += ", fase " + num(phases[pi], 4) + "$^\\circ$";
                }
                const std::string file = stem + ".png";
                if (!stbi_write_png((pics / file).string().c_str(),
                                    sub.w, sub.h, 4, sub.px.data(), sub.w * 4))
                    return false;
                figs.push_back({file, cap + ".", stem});
            }
        }

        // The animation itself, when asked for. LaTeX cannot embed a GIF, so it
        // ships beside the source rather than inside the document.
        if (options.isGif) {
            const int frames = std::max(1, options.frames);
            GifWriter gw{};
            const std::string gifPath =
                (root / (fs::path(options.filename).filename().string() + ".gif")).string();
            if (GifBegin(&gw, gifPath.c_str(), uint32_t(pageW), uint32_t(pageH),
                         uint32_t(options.delayCs))) {
                Canvas gpage;
                for (int f = 0; f < frames; ++f) {
                    const double ph = 2.0 * kPi * double(f) / double(frames);
                    composeFrame(gpage, ctx, items, yOff, pageW, pageH,
                                 ph, view3d, proj3d, sphereRadius, {}, true,
                                 options.incColorBar);
                    GifWriteFrame(&gw, gpage.px.data(), uint32_t(pageW),
                                  uint32_t(pageH), uint32_t(options.delayCs));
                }
                GifEnd(&gw);
            }
        }

        std::ofstream tex((root / "main.tex").string(), std::ios::binary);
        if (!tex) return false;
        tex << buildLatex(ctx, options, figs);
        return tex.good();
    }

    Canvas page;
    if (options.isGif) {
        const int frames = std::max(1, options.frames);
        GifWriter gw{};
        const std::string fname = options.filename + ".gif";
        if (!GifBegin(&gw, fname.c_str(),
                      uint32_t(pageW), uint32_t(pageH),
                      uint32_t(options.delayCs)))
            return false;
        for (int f = 0; f < frames; ++f) {
            const double phase = 2.0 * kPi * double(f) / double(frames);
            composeFrame(page, ctx, items, yOff, pageW, pageH,
                         phase, view3d, proj3d, sphereRadius,
                         {}, true, options.incColorBar);
            GifWriteFrame(&gw, page.px.data(),
                          uint32_t(pageW), uint32_t(pageH),
                          uint32_t(options.delayCs));
        }
        GifEnd(&gw);
        return true;
    } else {
        composeFrame(page, ctx, items, yOff, pageW, pageH,
                     0.0, view3d, proj3d, sphereRadius, {}, true,
                     options.incColorBar);
        const std::string fname = options.filename + ".png";
        return stbi_write_png(fname.c_str(), pageW, pageH, 4,
                              page.px.data(), pageW * 4) != 0;
    }
}

} // namespace waveguide
