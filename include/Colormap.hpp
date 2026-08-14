#pragma once

// THE field colormap. Blue (low) -> cyan -> green -> amber -> red -> dark red.
//
// Chosen for a light background: a black-to-white ramp loses its low end on the
// pale viewport and its high end on white paper. Every plot that maps a field
// magnitude to colour must use this and nothing else -- the 3D cloud, the cross
// sections, the field lines, the colour bar and the exported figures all share
// one scale, so a hue means the same value wherever it appears.
//
// It used to be copied by hand into eleven places; five of them still held an
// older black->purple->white ramp, which is why an exported cross section and
// the colour bar printed beside it disagreed with the 3D view.
//
// ONE copy cannot include this header: shaders/cloud.frag, because it is GLSL.
// Change the stops here and mirror them there.

namespace waveguide {

inline void fireColor(float t, float& r, float& g, float& b)
{
    if (t < 0.0f) t = 0.0f;
    if (t > 1.0f) t = 1.0f;
    static const float stops[6][3] = {
        {0.20f, 0.45f, 0.95f},   // blue      (low)
        {0.10f, 0.72f, 0.85f},   // cyan
        {0.20f, 0.72f, 0.25f},   // green
        {0.98f, 0.70f, 0.10f},   // amber
        {0.92f, 0.25f, 0.10f},   // red
        {0.55f, 0.00f, 0.08f},   // dark red  (high)
    };
    const float s = t * 5.0f;
    const int   i = (s < 4.0f) ? int(s) : 4;
    const float l = s - float(i);
    r = stops[i][0] + l * (stops[i + 1][0] - stops[i][0]);
    g = stops[i][1] + l * (stops[i + 1][1] - stops[i][1]);
    b = stops[i][2] + l * (stops[i + 1][2] - stops[i][2]);
}

} // namespace waveguide
