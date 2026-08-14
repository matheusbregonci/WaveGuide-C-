#pragma once

// Shared value types used across the field models, the renderer and the
// (future) numerical solvers. Kept in one small header so a new FieldSource
// implementation can depend on these without pulling in a concrete model.

namespace waveguide
{

    // One particle ready to be uploaded to the GPU.
    // Positions are CENTERED on the guide (so the model sits around the origin)
    // to make camera framing and the cutaway easy. `intensity` is the
    // normalized |field| in [0,1], used in the shader to scale the sphere size
    // so nodes shrink to nothing while antinodes swell.
    struct Particle
    {
        float x, y, z;   // centered position in meters
        float r, g, b;   // flat color (fire heatmap)
        float intensity; // normalized magnitude in [0,1]
    };

    enum class ModeType
    {
        TE,
        TM
    };

    enum class FieldKind
    {
        Electric, // |E| in V/m
        Magnetic  // |H| in A/m
    };

    // Axis-aligned extent of a model, in meters. For a cylinder, width and
    // height are the bounding diameter (2R).
    struct Bounds
    {
        float width;  // a  (meters)
        float height; // b  (meters)
        float depth;  // L  (meters)
    };

    // A baked cloud point: centered position + the Fourier decomposition of the
    // instantaneous |field|^2 in the animation phase: |field|^2(phase) = dc +
    // a2*cos(2*phase) + b2*sin(2*phase). The shader reconstructs the intensity
    // per frame from a phase uniform, so the cloud is sampled ONCE (on change)
    // and animated on the GPU (no per-frame CPU resampling / re-upload).
    struct CloudPoint
    {
        float x, y, z;   // centered position in meters
        float dc, a2, b2;
    };

} // namespace waveguide
