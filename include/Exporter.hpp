#pragma once

#include "CylindricalModel.hpp"
#include "TEmnModel.hpp"
#include "Renderer.hpp"

#include <glm/glm.hpp>
#include <string>
#include <vector>

namespace waveguide {

struct ExportOptions
{
    bool        isGif        = true;   // GIF animated; otherwise PNG static
    bool        incScene3D   = true;
    bool        incXY        = true;
    bool        incZX        = false;
    bool        incZY        = false;
    bool        incSpectrum  = true;
    bool        incColorBar  = true;
    int         frames       = 48;     // GIF frame count
    int         delayCs      = 4;      // GIF frame delay in centiseconds
    std::string filename     = "waveguide_export";

    // Emit a self-contained LaTeX project instead of one stacked image:
    //
    //   <filename>/
    //     main.tex          document with the configuration tables and figures
    //     pictures/*.png    one file per selected item, each its own figure
    //     <filename>.gif    the animation too, when isGif (LaTeX cannot embed
    //                       it, but it is useful next to the paper source)
    //
    // Static figures are what a paper or a textbook needs, so every item is
    // also written as a PNG even in GIF mode.
    bool        latexProject = false;
    // Phases (in degrees) at which each animated item is exported. One entry
    // gives a single figure per item; several give a phase series, which is how
    // a travelling mode is usually shown on paper.
    std::vector<double> phasesDeg = {0.0};
    std::string title        = "";     // document title; auto-built when empty
    std::string author       = "";

    // Target size of each exported figure. Items are sized to use as much of
    // this box as their true aspect ratio allows: a long ZX cut gets the full
    // width and stays short, a square XY cut is limited by the height instead.
    int         imageW       = 1200;
    int         imageH       = 900;
};

struct ExportContext
{
    TEmnModel*        rect      = nullptr;
    CylindricalModel* cyl       = nullptr;
    Renderer*         renderer  = nullptr;
    int               geometry  = 0;   // 0 = rect, 1 = cyl
    int               fieldKind = 0;   // 0 = E, 1 = H
    int               modeType  = 0;
    int               modeM     = 1;
    int               modeN     = 1;
    int               modeL     = 1;   // axial index (cavity only)
    int               structure = 0;   // 0 = waveguide, 1 = cavity
    double            freqHz    = 12e9;

    // Reported in the LaTeX tables. Everything else the report needs (bounds,
    // cutoff, peak field, amplitude, units) is read straight off the model, so
    // the tables cannot drift from what was actually simulated.
    std::string       mediumName = "Vacuum / air";
    double            powerW     = 1.0;  // requested transported power

    // Mirrors the on-screen "Modo 3D" choice so the report shows what you were
    // actually looking at: 0 = intensity cloud, 1 = field lines, 2 = both.
    int               view3D      = 0;
    float             lineDensity = 1.0f;  // 3D field-line density

    // State of the baked point-sprite cloud, so the exported 3D figure is the
    // same draw call the UI makes rather than the old instanced-sphere path
    // (different geometry AND a different colormap from the colour bar).
    float             cloudInvPeak     = 1.0f;
    float             cloudMeanSpacing = 1e-3f;
    bool              cloudOpaque      = true;

    // Cross sections: streamlines instead of the arrow grid, mirroring the
    // "Field lines (streamlines)" checkbox and its density slider.
    bool              sectionStreamlines = false;
    float             sectionLineDensity = 1.0f;
};

// Runs a blocking export pass. Returns true on success and writes the file
// next to the working directory using `options.filename` + ".png"/".gif".
bool runExport(const ExportContext& ctx, const ExportOptions& options);

} // namespace waveguide
