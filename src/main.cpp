// Waveguide TEmn field visualizer (Atoms-style).
// ------------------------------------------------
// Stochastic particle cloud of |E| inside a rectangular metallic waveguide
// running a TEmn mode. Density is proportional to |E|^2, color is the fire
// heatmap applied to |E|, spheres are instanced low-poly meshes.

#include "Camera.hpp"
#include "CylindricalModel.hpp"
#include "Exporter.hpp"
#include "FdtdSim.hpp"
#include "MicrostripSim.hpp"
#include "MicrostripXsec.hpp"
#include "Geometry.hpp"
#include "HelmholtzSolver.hpp"
#include "MaxwellSolver.hpp"
#include "NumericalModel.hpp"
#include "Renderer.hpp"
#include "Colormap.hpp"
#include "FieldViz.hpp"
#include "TEmnModel.hpp"
#include "gl_loader.h"

#include <glm/ext/matrix_projection.hpp>

#include <GLFW/glfw3.h>

#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace waveguide;

namespace
{
    // Directory that holds the running executable (independent of the current
    // working directory). Used to locate the shaders folder so the app runs
    // no matter where it is launched from — double-click, shortcut, or any cwd.
    std::string executableDir()
    {
#ifdef _WIN32
        char buf[MAX_PATH];
        const DWORD n = GetModuleFileNameA(nullptr, buf, MAX_PATH);
        const std::string p(buf, (n > 0) ? n : 0);
        const size_t slash = p.find_last_of("\\/");
        return (slash == std::string::npos) ? std::string(".") : p.substr(0, slash);
#else
        return ".";
#endif
    }

    // Pick the first candidate directory that actually contains the shaders.
    // Prefer the folder next to the executable (where CMake copies them), then
    // fall back to cwd-relative paths for dev runs from the project root.
    std::string resolveShaderDir()
    {
        const std::string exeDir = executableDir();
        const std::string candidates[] = {
            exeDir + "/shaders",
            "shaders",
            exeDir + "/../shaders",
            "../shaders",
        };
        for (const std::string& c : candidates)
        {
            std::ifstream probe(c + "/point.vert");
            if (probe.good())
                return c;
        }
        return exeDir + "/shaders"; // best guess; init() will report if missing
    }

    OrbitCamera gCamera;
    bool gMouseDown = false;
    double gLastX = 0.0, gLastY = 0.0;

    // Port-configuration mode: left-click a bounding-box face to cycle its role.
    bool gPortConfig = false;
    bool gPortClickPending = false;
    double gPortClickX = 0.0, gPortClickY = 0.0;
    double gPressX = 0.0, gPressY = 0.0;

    void glfwErrorCallback(int code, const char *desc)
    {
        std::fprintf(stderr, "GLFW error %d: %s\n", code, desc);
    }

    void framebufferSizeCallback(GLFWwindow *, int w, int h)
    {
        glViewport(0, 0, w, h);
    }

    void mouseButtonCallback(GLFWwindow *win, int button, int action, int /*mods*/)
    {
        if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse)
        {
            gMouseDown = false;
            return;
        }
        if (button == GLFW_MOUSE_BUTTON_LEFT)
        {
            gMouseDown = (action == GLFW_PRESS);
            if (gMouseDown) {
                glfwGetCursorPos(win, &gLastX, &gLastY);
                gPressX = gLastX; gPressY = gLastY;
            } else if (gPortConfig) {
                // Release: if it was a click (little movement), request a face pick.
                double cx, cy; glfwGetCursorPos(win, &cx, &cy);
                if (std::abs(cx - gPressX) < 6.0 && std::abs(cy - gPressY) < 6.0) {
                    gPortClickPending = true; gPortClickX = cx; gPortClickY = cy;
                }
            }
        }
    }

    void cursorPosCallback(GLFWwindow *, double x, double y)
    {
        if (!gMouseDown)
            return;
        if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse)
            return;
        const double dx = x - gLastX;
        const double dy = y - gLastY;
        gLastX = x;
        gLastY = y;
        gCamera.orbit(float(dx), float(dy));
    }

    void scrollCallback(GLFWwindow *, double /*xoff*/, double yoff)
    {
        if (ImGui::GetCurrentContext() && ImGui::GetIO().WantCaptureMouse)
            return;
        gCamera.zoom(float(yoff));
    }

    // Hot-swappable visualization params.
    bool gCutawayOn = true;
    // Only the Z count survives, and only to size the (now unused) sphere
    // instances. The X/Y counts went with the particle-restore call that used
    // to run after an export; nothing samples on that grid any more.
    int gGridNz = 110;
    // Field cloud is now a Monte-Carlo point set: points drawn uniformly at
    // random in the bounding box (deterministic seed, so positions are stable
    // across frames), culled by inside()/cutaway/threshold. gCloudSamples is the
    // target number of VISIBLE points -- rejection sampling keeps drawing until it
    // has that many, so a thin field region (e.g. a microstrip in a big air box)
    // still fills in instead of coming out sparse. Resolution = this count.
    int   gCloudSamples = 60000;
    float gCloudMeanSpacing = 1e-3f;   // mean spacing of accepted points -> dot size
    // Opaque cloud: near dots occlude far ones, so a pixel shows the field at
    // ONE point instead of the sum along the whole line of sight.
    bool  gCloudOpaque = true;
    float gPhaseSpeed = 2.5f; // radians per second
    bool gAnimate = true;
    // Cross-section overlay: false = arrow (quiver) field, true = streamlines
    // (integral field lines traced from the vector field).
    bool gFieldLines = false;
    // Animated "current" along the (frozen-geometry) streamlines: dashes scroll
    // in the field direction, faster where the field is stronger. gFlowPhase is
    // a free-running clock advanced every frame; gFlowSpeed is a UI multiplier.
    float  gFlowSpeed = 1.0f;
    double gFlowPhase = 0.0;
    // Cross-section streamline density. Seeds sit on a fixed lattice measured in
    // SCREEN pixels, so a wide window would otherwise take proportionally more
    // lines and turn into a thicket. 1.0 = one seed per ~46 px.
    float  gSecLineDensity = 1.0f;
    // 3D field-line density multiplier: scales the seed + spacing grids (more
    // lines = more per-frame tracing cost).
    float  gFieldLineDensity = 1.0f;
    // Colour the microstrip copper by the field intensity on its surface (the
    // "skin" / surface-current view) instead of a flat copper colour.
    bool   gCopperSkin = false;

    // Slice planes shown in the 3D view. When gManualSlice is false the slice
    // is auto-placed at an antinode; when true, the sliders below drive it and
    // the cross-section plots follow. gSliceZ/Y/X are fractions in [0,1] of the
    // corresponding dimension (z for XY, y for ZX, x for ZY).
    bool  gShowSlicePlanes = true;
    bool  gManualSlice     = false;
    float gSliceZ = 0.5f, gSliceY = 0.5f, gSliceX = 0.5f;

    // Main 3D view mode: 0 = particle cloud, 1 = 3D field lines (streamlines
    // traced through the volume), 2 = cloud + field lines.
    int gView3D = 0;


    // Pre-configured rectangular (WR) waveguides. fmin/fmax are the recommended
    // single-mode operating band (GHz) from the standard tables; the last entry
    // is the user-editable "custom" size. Cylindrical presets can be added the
    // same way once that table is available.
    struct RectPreset { const char* label; float a_mm, b_mm, fmin, fmax; };
    const RectPreset kRectPresets[] = {
        {"WR-42 (K)",     10.70f,  4.30f, 18.00f, 26.50f},
        {"WR-62 (Ku)",    15.80f,  7.90f, 12.40f, 18.00f},
        {"WR-90 (X)",     22.86f, 10.16f,  8.20f, 12.40f},
        {"WR-112 (W)",    28.50f, 12.62f,  7.05f, 10.00f},
        {"WR-137 (C)",    34.85f, 15.80f,  5.85f,  8.20f},
        {"Personalizada", 22.86f, 10.16f,  0.00f,  0.00f},
    };
    const int kNumRectPresets = 6;
    const int kCustomRectIndex = 5;

    // Pre-configured circular waveguides (radius-based); last entry is custom.
    struct CylPreset { const char* label; float R_mm, fmin, fmax; };
    const CylPreset kCylPresets[] = {
        {"Guia 1 (X)",   23.83f,  8.50f, 11.60f},
        {"Guia 2 (Ku)",  15.08f, 13.40f, 18.00f},
        {"Guia 3 (K)",   10.06f, 20.00f, 24.50f},
        {"Guia 4 (Ka)",   6.35f, 33.00f, 38.50f},
        {"Guia 5 (Q)",    5.56f, 38.50f, 43.00f},
        {"Personalizada", 23.00f, 0.00f,  0.00f},
    };
    const int kNumCylPresets = 6;
    const int kCustomCylIndex = 5;

    // One step of the CSG geometry builder (Phase 1/2). Parameters are in mm.
    //   box:      p = {cx, cy, cz, sx, sy, sz}
    //   cylinder: p = {cx, cy, z0, z1, radius, -}
    struct BuildStep
    {
        int   type = 0; // 0 = box, 1 = cylinder (along z)
        int   op   = 0; // 0 = add (union), 1 = subtract
        float p[6] = {0.f, 0.f, 100.f, 22.86f, 10.16f, 200.f};
    };

    // ---- Persist CSG shapes to shapes/<name>.csg (one line per step) ----
    const char* kShapeDir = "shapes";

    void saveShape(const std::string& name, const std::vector<BuildStep>& steps)
    {
        std::error_code ec;
        std::filesystem::create_directories(kShapeDir, ec);
        std::ofstream f(std::string(kShapeDir) + "/" + name + ".csg");
        if (!f) return;
        f << "# waveguide CSG shape (type op cx cy cz sx sy sz | mm)\n";
        for (const BuildStep& s : steps)
        {
            f << s.type << ' ' << s.op;
            for (int i = 0; i < 6; ++i) f << ' ' << s.p[i];
            f << '\n';
        }
    }

    bool loadShape(const std::string& name, std::vector<BuildStep>& steps)
    {
        std::ifstream f(std::string(kShapeDir) + "/" + name + ".csg");
        if (!f) return false;
        std::vector<BuildStep> out;
        std::string line;
        while (std::getline(f, line))
        {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream is(line);
            BuildStep s;
            if (is >> s.type >> s.op >> s.p[0] >> s.p[1] >> s.p[2] >> s.p[3] >> s.p[4] >> s.p[5])
                out.push_back(s);
        }
        if (out.empty()) return false;
        steps = out;
        return true;
    }

    std::vector<std::string> listShapes()
    {
        std::vector<std::string> names;
        std::error_code ec;
        if (!std::filesystem::exists(kShapeDir, ec)) return names;
        for (const auto& e : std::filesystem::directory_iterator(kShapeDir, ec))
            if (e.path().extension() == ".csg") names.push_back(e.path().stem().string());
        std::sort(names.begin(), names.end());
        return names;
    }

    // Detect the real port openings: connected patches of solid on each
    // bounding-box face (an arm end), not the whole face. Transverse extents are
    // in model coords [0,size]. Default roles: input at the z-min opening.
    std::vector<FdtdPort> detectPorts(const Geometry& g, int res)
    {
        std::vector<FdtdPort> ports;
        const Aabb bb = g.bounds();
        const double sx = bb.sizeX(), sy = bb.sizeY(), sz = bb.sizeZ();
        const double smax = std::max({sx, sy, sz, 1e-9});
        const int nx = std::max(4, int(res*sx/smax)), ny = std::max(4, int(res*sy/smax)), nz = std::max(4, int(res*sz/smax));
        const VoxelMask m = g.voxelize(nx, ny, nz, bb);
        auto solid = [&](int i, int j, int k) {
            if (i<0||i>=nx||j<0||j>=ny||k<0||k>=nz) return false;
            return m.occ[(std::size_t(k)*ny+j)*nx+i] != 0; };
        for (int axis = 0; axis < 3; ++axis) for (int side = 0; side < 2; ++side)
        {
            const int dim = (axis==0)?nx:(axis==1)?ny:nz;
            const int plane = (side==0)?0:dim-1;
            int tu, tv; double du, dv;
            if (axis==0){ tu=ny; tv=nz; du=sy/ny; dv=sz/nz; }
            else if (axis==1){ tu=nx; tv=nz; du=sx/nx; dv=sz/nz; }
            else { tu=nx; tv=ny; du=sx/nx; dv=sy/ny; }
            auto solAt = [&](int a, int b) {
                int i,j,k;
                if (axis==0){ i=plane; j=a; k=b; } else if (axis==1){ i=a; j=plane; k=b; } else { i=a; j=b; k=plane; }
                return solid(i,j,k); };
            std::vector<int> lab(std::size_t(tu)*tv, -1);
            for (int b0 = 0; b0 < tv; ++b0) for (int a0 = 0; a0 < tu; ++a0)
            {
                if (!solAt(a0,b0) || lab[b0*tu+a0]!=-1) continue;
                std::vector<std::pair<int,int>> stack{{a0,b0}};
                lab[b0*tu+a0]=1;
                double umin=1e30,umax=-1e30,vmin=1e30,vmax=-1e30;
                while (!stack.empty()) {
                    const int ca=stack.back().first, cb=stack.back().second; stack.pop_back();
                    umin=std::min(umin,ca*du); umax=std::max(umax,(ca+1)*du);
                    vmin=std::min(vmin,cb*dv); vmax=std::max(vmax,(cb+1)*dv);
                    const int dd[4][2]={{1,0},{-1,0},{0,1},{0,-1}};
                    for (auto& d : dd) { const int na=ca+d[0], nb=cb+d[1];
                        if (na<0||na>=tu||nb<0||nb>=tv) continue;
                        if (solAt(na,nb) && lab[nb*tu+na]==-1){ lab[nb*tu+na]=1; stack.push_back({na,nb}); } }
                }
                FdtdPort p; p.axis=axis; p.side=side; p.role=0;
                p.uMin=umin; p.uMax=umax; p.vMin=vmin; p.vMax=vmax;
                ports.push_back(p);
            }
        }
        int inIdx = -1;
        for (std::size_t i = 0; i < ports.size(); ++i) if (ports[i].axis==2 && ports[i].side==0){ inIdx=int(i); break; }
        if (inIdx<0 && !ports.empty()) inIdx=0;
        for (std::size_t i = 0; i < ports.size(); ++i) ports[i].role = (int(i)==inIdx)?1:2;
        return ports;
    }

    struct Resample
    {
        bool needed = false;
        TEmnModel *model = nullptr;
        Renderer *rend = nullptr;
    } gResample;

    void keyCallback(GLFWwindow *win, int key, int /*sc*/, int action, int /*mods*/)
    {
        if (action != GLFW_PRESS && action != GLFW_REPEAT)
            return;

        if (key == GLFW_KEY_ESCAPE)
        {
            glfwSetWindowShouldClose(win, GLFW_TRUE);
            return;
        }
        if (key == GLFW_KEY_C)
        {
            gCutawayOn = !gCutawayOn;
            gResample.needed = true;
        }
        if (key == GLFW_KEY_EQUAL || key == GLFW_KEY_KP_ADD)   // denser Monte-Carlo cloud
            gCloudSamples = std::min(int(gCloudSamples * 1.3f) + 1000, 500000);
        if (key == GLFW_KEY_MINUS || key == GLFW_KEY_KP_SUBTRACT)
            gCloudSamples = std::max(int(gCloudSamples / 1.3f), 5000);
        if (key == GLFW_KEY_SPACE)
            gAnimate = !gAnimate;
        if (key == GLFW_KEY_LEFT_BRACKET)
            gPhaseSpeed = std::max(gPhaseSpeed - 0.5f, 0.0f);
        if (key == GLFW_KEY_RIGHT_BRACKET)
            gPhaseSpeed = std::min(gPhaseSpeed + 0.5f, 30.0f);
    }
} // namespace

int main()
{
    glfwSetErrorCallback(glfwErrorCallback);
    if (!glfwInit())
    {
        std::fprintf(stderr, "Failed to init GLFW\n");
        return EXIT_FAILURE;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    const int winW = 1280;
    const int winH = 800;
    GLFWwindow *window = glfwCreateWindow(winW, winH,
                                          "Waveguide TEmn | |E| Probabilistic Cloud", nullptr, nullptr);
    if (!window)
    {
        std::fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return EXIT_FAILURE;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glfwSetFramebufferSizeCallback(window, framebufferSizeCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetCursorPosCallback(window, cursorPosCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetKeyCallback(window, keyCallback);

    if (!loadOpenGLFunctions(reinterpret_cast<GLLoadProc>(glfwGetProcAddress)))
    {
        std::fprintf(stderr, "Failed to load OpenGL function pointers\n");
        glfwTerminate();
        return EXIT_FAILURE;
    }

    std::printf("OpenGL version: %s\n", glGetString(GL_VERSION));
    std::printf("GLSL   version: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
    std::printf("Renderer      : %s\n", glGetString(GL_RENDERER));

    // ---------------- ImGui init ----------------
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsLight();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330 core");
    ImGui::GetIO().ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // dockable panels

    // ---------------- physics: live-editable params ----------------
    float  uiWidthMM   = 22.86f;
    float  uiHeightMM  = 10.16f;
    float  uiFreqGHz   = 12.0f;
    float  uiSweepMinGHz = 8.0f;   // FDTD S-parameter sweep band
    float  uiSweepMaxGHz = 16.0f;
    float  uiCWfreqGHz   = 10.0f;  // single-frequency CW animation

    // Filling medium (relative permittivity/permeability). The analytic models
    // are fully eps_r/mu_r aware: the medium slows the wave (v = c/sqrt(eps*mu)),
    // lowering the cutoff / resonant frequency and rescaling the fields.
    int    uiDielectric = 0;   // index into the preset table below
    float  uiEpsR       = 1.0f;
    float  uiMuR        = 1.0f;
    // Label of the selected medium, mirrored out of the preset table so the
    // exporter can name it in the report without duplicating the table.
    std::string gMediumName = "Vacuum / air";
    int    uiM         = 1;
    int    uiN         = 0;  // default TE10: the dominant, above-cutoff mode
    int    uiL         = 1;  // axial half-wave index (cavity only)
    int    uiModeType  = 0; // 0 = TE, 1 = TM
    int    uiGeometry  = 0; // 0 = Rectangular, 1 = Cylindrical
    int    uiStructure = 0; // 0 = Waveguide (open), 1 = Cavity (shorted ends)
    int    uiField     = 0; // 0 = Electric (V/m), 1 = Magnetic (A/m)
    // Time-average power the mode carries down the guide. The analytic models
    // are eigenfunctions with a free amplitude; pinning it to a power is what
    // turns the plotted numbers into real V/m and A/m. 1 W matches the HFSS
    // wave-port default, so the two tools become directly comparable.
    float  uiPowerW    = 1.0f;
    float  uiRadiusMM  = 23.83f; // matches the default cylindrical preset (Guia 1)
    float  uiDepthMM   = 300.0f; // rectangular depth d (also the shown length)
    float  uiLengthMM  = 300.0f; // cylindrical length d
    int    uiRectPreset = 2;     // default WR-90 (X) — matches the rect dims above
    int    uiCylPreset  = 0;     // default Guia 1 (X)

    auto asField = [&]() { return uiField == 0 ? FieldKind::Electric : FieldKind::Magnetic; };
    auto asMode  = [&]() { return uiModeType == 0 ? ModeType::TE : ModeType::TM; };
    auto isCavity = [&]() { return uiStructure == 1; };

    TEmnModel rectModel(uiWidthMM, uiHeightMM, uiFreqGHz * 1e9,
                        uiEpsR, uiMuR, uiM, uiN, asMode(), asField(),
                        uiDepthMM, isCavity(), uiL, double(uiPowerW));
    CylindricalModel cylModel(uiRadiusMM, uiLengthMM / 1000.0, uiFreqGHz * 1e9,
                              uiEpsR, uiMuR, uiN, uiM, asMode(), asField(),
                              isCavity(), uiL);

    // Phase 3: numerically solved field over an arbitrary CSG shape. When
    // present and enabled, it becomes the active FieldSource.
    std::unique_ptr<NumericalModel> numModel;
    bool gUseNumerical = false;
    int  gNumMode = 0;
    int  gSolveRes = 32;    // grid cells on the longest axis
    int  gSolverType = 0;   // 0 = scalar Helmholtz, 1 = vector Maxwell

    // Exactly one simulation domain is live at a time. gDomain is the single
    // authority over what renders; switchDomain() (below) enforces exclusivity so
    // a closed panel can never leave its simulation on screen.
    enum class SimDomain { Waveguide, Geometry, Microstrip };
    SimDomain gDomain = SimDomain::Waveguide;

    // Phase 5: time-domain FDTD — a pulse/CW source launched into the geometry,
    // advanced each frame so the wave visibly propagates through the circuit.
    std::unique_ptr<FdtdSim> gFdtd;
    bool gUseFdtd = false;       // FDTD is the active/displayed source
    bool gFdtdStepping = false;  // time is advancing (frozen once the sweep converges)
    bool gSweepDone = false;     // pulse rang down -> S-curve is final & stable
    int  gFdtdSteps = 12;   // time steps advanced per frame

    // VNA-style CW sweep: one steady-state measurement per frequency, filling the
    // S-curve left to right. gS21db/gS11db hold the measured points (NaN = not yet
    // measured); gCWindex is the frequency currently being measured.
    bool gCWsweep = false;       // a CW sweep is in progress
    bool gCWhold  = false;       // continuous CW at a single fixed frequency (animation)
    int  gCWindex = 0;
    int  gCWnf = 61;             // frequency points in the sweep (user-adjustable)
    int  gCWsteps = 80;          // FDTD steps advanced per frame during the sweep
    int  gCWmaxBlocks = 150;     // give up steady-state search after this many windows
    std::vector<double> gCWfreqs;
    std::vector<float>  gS21db, gS11db;

    // Persistent snapshot of the broadband-pulse S-curve, so it can be overlaid
    // against the VNA sweep for comparison even after the sim is rebuilt.
    bool gPulseActive = false;   // last run was the pulse method (keep snapshotting)
    std::vector<double> gPulseFreqs;
    std::vector<float>  gPulseS21db, gPulseS11db;
    std::vector<FdtdPort> gPorts; // real openings of the solid, with roles

    // ---- Microstrip (open FDTD): draw the trace, run the quasi-TEM wave -------
    std::unique_ptr<MicrostripSim> gMicro;
    // THRU de-embedding reference: same stack + source + ports as gMicro but a
    // straight feed-width line spanning the whole domain (filter body removed) --
    // an ideal pass-through. The sweep divides its complex S-params out, so the
    // plotted S11/S21 are referenced past the access lines (to an equal-length
    // uniform line). Rebuilt with gMicro; stepped in lock-step during a sweep.
    std::unique_ptr<MicrostripSim> gMsThru;
    bool  gMsDeembed = true;       // divide the sweep by the THRU reference
    // Passivity check: max over the band of |S11|^2 + |S21|^2 from the RAW
    // (self-consistent, non-de-embedded) S-params. A passive 2-port must have
    // this <= 1 (=1 lossless, <1 lossy); > 1 means energy is being created.
    // gMsThruPassivity is the same for the straight THRU line (which is passive by
    // construction) -- if IT exceeds 1 too, the S-param EXTRACTION is the culprit,
    // not the field.
    float gMsRawPassivity = 0.0f;
    float gMsThruPassivity = 0.0f;
    bool  gUseMicro = false;       // microstrip sim is the active/displayed source
    bool  gMicroStepping = false;
    float uiMsHsubMM   = 1.6f;     // substrate thickness
    float uiMsEpsSub   = 4.4f;     // substrate permittivity (FR-4)
    float uiMsStripWmm = 3.0f;     // trace width (straight preset)
    float uiMsLenMM    = 40.0f;    // trace length (straight preset)
    float uiMsAirMM    = 6.0f;     // air height above the trace
    int   uiMsSubCells = 4;        // cells across the substrate (sets dx)
    float uiMsFcGHz    = 5.0f;     // pulse center frequency
    float uiMsSrcFrac  = 0.0f;     // launch-port at the input extremity (0% length)
    float uiMsSenseFrac= 1.0f;     // sense-port at the output extremity (100% length)
    float uiMsFilterScale = 1.0f;  // low-pass preset: mm per drawing unit (1.0 = the
                                   // paper's millimetre dimensions).
    float uiMsLengthFactor = 1.4f; // Duenas "lengthening factor": stretches the
                                   // propagation-direction (x) lengths only, to
                                   // correct the too-fast numerical phase velocity
                                   // of the coarse mesh (dips move DOWN in freq).
                                   // Leaves the transverse widths -> impedances fixed.
    bool  showMicrostrip = true;   // Microstrip lives in its own window

    // ---- 2D microstrip trace designer (top view, x-z plane) ----
    // A layout of copper patches (rectangles + circles, mm) that union into an
    // arbitrary trace. x = propagation (left = input, right = output), z = width.
    // 0=rect (w x d), 1=ring/disk. For a ring: r = centerline (mean) radius and
    // w = trace width, so outer = r + w/2 and inner = r - w/2 (thick unused now).
    struct MsShape { int type=0; float x=0, z=0, w=4, d=1, r=1, thick=0; };
    std::vector<MsShape> gMsDesign;
    int    gMsDesignSel = -1;
    bool   gMsDesignOpen = false;      // dedicated 2D designer window open
    float  gMsZoom = 1.0f;             // canvas zoom (x auto-fit basis)
    ImVec2 gMsPan = ImVec2(0.0f,0.0f); // canvas pan offset (screen px)
    int    gMsDragHandle = -1;         // -1 none/move, >=0 resize handle, -2 pan

    // 2D cross-section (quasi-static) solver result: Z0 / eps_eff of a straight
    // line, validated against Hammerstad. Phase 1 of the frequency-domain solver.
    XsecResult gXsec;         // FD (finite-difference) result
    XsecResult gXsecHam;      // Hammerstad closed form (comparison)
    int  gXsecCells = 32;     // grid cells across the substrate height

    // Microstrip frequency sweep: retunes fc across a band, letting the field
    // settle at each point, and records Z0 / eps_eff / L' / C' vs frequency.
    bool  gMsSweep = false, gMsSweepDone = false, gMsWindowOpen = false;
    int   gMsSweepIdx = 0, gMsSweepStepsDone = 0;
    int   gMsSweepSettle = 1500;    // steps to reach steady state before measuring
    int   gMsSweepTarget = 3000;    // total steps per frequency point
    int   gMsSweepChunk  = 250;     // steps advanced per frame during the sweep
    float uiMsSweepMinGHz = 1.0f, uiMsSweepMaxGHz = 15.0f;
    int   uiMsSweepPts = 15;
    std::vector<float> gMsFreqGHz, gMsZ0, gMsEeff, gMsLp, gMsCp, gMsS11, gMsS21;

    // ---- Developer mode (Ctrl+D): extra "plots aleatorios" -------------------
    bool  gDevMode    = false;
    bool  gShowDipole = false;     // rotating magnetic dipole (pulsar) plot
    double gDipoleT   = 0.0;       // animation time (P units)
    float uiDipAlphaDeg = 30.0f;   // obliquity angle alpha
    float uiDipPeriod   = 4.0f;    // rotation period P
    float uiDipMag      = 1.0f;    // dipole magnitude m
    bool  uiDipAnimate  = true;
    float uiDipSpeed    = 1.0f;
    bool  uiDipFieldLines = true;
    int   uiDipFieldAz  = 16;      // field lines per shell (azimuths)
    int   uiDipShellsN  = 6;       // number of nested field-line shells
    bool  uiDipCones    = true;    // pulsar radiation cones at the poles
    float uiDipConeLen  = 1.6f;    // cone length

    // ---- Thru calibration ----------------------------------------------------
    // Our raw S normalizes by the source DRIVE waveform, which sits at an
    // arbitrary level (a scale factor K), so |S| can float above 0 dB. A thru
    // calibration runs the SAME source into a straight matched line, measures the
    // launched incident wave, and stores H(f) = X_incident / X_source. Dividing
    // every raw measurement by H(f) removes K: a straight thru then reads 0 dB
    // and passive structures read <= 0 dB, exactly like a real VNA.
    bool gHaveCal = false;   // reference measured
    bool gApplyCal = false;  // apply it to the plotted S
    std::vector<double> gCalFreqs, gCalHre, gCalHim; // H(f) over the reference band

    // Interpolate the complex calibration H at a frequency (nearest-linear).
    auto calH = [&](double fHz, double& hr, double& hi) {
        hr = 1.0; hi = 0.0;
        const int n = int(gCalFreqs.size());
        if (n < 2) return;
        if (fHz <= gCalFreqs.front()) { hr = gCalHre.front(); hi = gCalHim.front(); return; }
        if (fHz >= gCalFreqs.back())  { hr = gCalHre.back();  hi = gCalHim.back();  return; }
        int lo = 0, hiIdx = n - 1;
        while (hiIdx - lo > 1) { int mid = (lo + hiIdx) / 2; (gCalFreqs[mid] <= fHz ? lo : hiIdx) = mid; }
        const double t = (fHz - gCalFreqs[lo]) / (gCalFreqs[hiIdx] - gCalFreqs[lo] + 1e-30);
        hr = gCalHre[lo] + t * (gCalHre[hiIdx] - gCalHre[lo]);
        hi = gCalHim[lo] + t * (gCalHim[hiIdx] - gCalHim[lo]);
    };

    // Turn raw port/source spectra into calibrated S21/S11 in dB. Xout, Xin, Xsrc
    // are complex; when calibration is active they are divided by H(f).
    auto calS = [&](double orr, double oii, double irr, double iii,
                    double cr,  double ci,  double fHz,
                    float& s21dB, float& s11dB) {
        const double cd = cr*cr + ci*ci + 1e-30;
        double r21r = (orr*cr + oii*ci)/cd, r21i = (oii*cr - orr*ci)/cd; // Xout/Xsrc
        double rinr = (irr*cr + iii*ci)/cd, rini = (iii*cr - irr*ci)/cd; // Xin /Xsrc
        double hr = 1.0, hi = 0.0;
        if (gHaveCal && gApplyCal) calH(fHz, hr, hi);
        const double hd = hr*hr + hi*hi + 1e-300;
        const double s21r = (r21r*hr + r21i*hi)/hd, s21i = (r21i*hr - r21r*hi)/hd;
        const double inr  = (rinr*hr + rini*hi)/hd, ini  = (rini*hr - rinr*hi)/hd;
        const double s11r = inr - 1.0, s11i = ini;   // S11 = Xin/Xinc - 1
        s21dB = float(20.0*std::log10(std::max(std::hypot(s21r,s21i), 1e-12)));
        s11dB = float(20.0*std::log10(std::max(std::hypot(s11r,s11i), 1e-12)));
    };

    // The single seam between the visualization and the physics: everything
    // below samples the field through this FieldSource* instead of naming a
    // concrete model, so the numerical solver drops in here.
    auto active = [&]() -> FieldSource* {
        if (gUseMicro && gMicro) return gMicro.get();
        if (gUseFdtd && gFdtd) return gFdtd.get();
        if (gUseNumerical && numModel) return numModel.get();
        return uiGeometry == 0 ? static_cast<FieldSource*>(&rectModel)
                               : static_cast<FieldSource*>(&cylModel);
    };
    auto activeBounds = [&]() { return active()->bounds(); };
    auto activeSample = [&](int nx, int ny, int nz, bool co, float mi, double ph) {
        return active()->sampleGrid(nx, ny, nz, co, mi, ph);
    };
    auto activePeak = [&]() { return active()->peakField(); };

    // ---- Geometry builder (Phase 2): a CSG shape assembled from primitives.
    // Purely a design/preview tool for now — the field is still the analytic
    // model; the numerical solver over this shape is Phase 3. Seeded with a
    // T-junction so enabling the builder shows something immediately.
    bool gBuilderOn = false;
    std::vector<BuildStep> gSteps = {
        {0, 0, {0.f, 0.f,     100.f, 22.86f, 10.16f, 200.f}},   // main guide (z)
        {0, 0, {0.f, 30.08f,  100.f, 22.86f, 50.00f, 22.86f}},  // stub (+y)
    };
    auto buildGeo = [&]() -> Geometry {
        Geometry g;
        for (const BuildStep& s : gSteps) {
            const float* p = s.p;
            Geometry prim = (s.type == 0)
                ? Geometry::box(p[0]/1e3, p[1]/1e3, p[2]/1e3,
                                p[3]/1e3, p[4]/1e3, p[5]/1e3)
                : Geometry::cylinderZ(p[0]/1e3, p[1]/1e3, p[2]/1e3, p[3]/1e3, p[4]/1e3);
            g = (s.op == 0) ? g.unite(prim) : g.subtract(prim);
        }
        return g;
    };

    // The cloud is now a baked point-sprite buffer (updateCloud), so the sphere
    // instance buffer starts empty.
    const std::string shaderDir = resolveShaderDir();
    Renderer renderer;
    if (!renderer.init(shaderDir, {}, activeBounds()))
    {
        std::fprintf(stderr, "Renderer init failed\n");
        glfwTerminate();
        return EXIT_FAILURE;
    }

    // ---- baked field cloud: sampled once (on change), animated on the GPU ----
    float gCloudInvPeak = 1.0f;
    bool  gCloudDirty = true;
    // The sampling itself now lives in FieldViz (no renderer, no UI), so the
    // exporter and the coming WebAssembly build can produce the same cloud.
    // This lambda is only the glue: UI state in, GPU buffer out.
    auto bakeCloud = [&](int K) {
        CloudParams cp;
        cp.targetPoints = gCloudSamples;
        cp.harmonics    = K;
        cp.cutaway      = gCutawayOn;
        // Qualified: the local lambda would otherwise shadow the free function.
        const CloudResult cr = waveguide::bakeCloud(*active(), cp);
        gCloudInvPeak     = cr.invPeak;
        gCloudMeanSpacing = cr.meanSpacing;
        renderer.updateCloud(cr.points);
    };

    gResample.model = &rectModel;
    gResample.rend = &renderer;

    // Build initial cylinder mesh for cylindrical geometry.
    {
        Bounds cb = cylModel.bounds();
        renderer.updateCylinder(cb.width * 0.5f, cb.depth);
    }

    // Frame the waveguide (centered on origin in the new model).
    Bounds b = activeBounds();
    gCamera.setTarget({0.0f, 0.0f, 0.0f});
    const float maxDim = std::max({b.width, b.height, b.depth});
    gCamera.setDistance(maxDim * 2.2f);

    // Sphere radius in world units — tuned to look like small glowing dots.
    // On a regular grid the spacing along z is ~depth/nz. Make the max
    // sphere radius roughly half that so antinodes look full and nodes
    // disappear (radius scales with intensity in the vertex shader).
    float sphereRadius = (b.depth / float(gGridNz)) * 0.9f;

    std::printf(
        "Controls:\n"
        "  Mouse drag   - orbit\n"
        "  Scroll       - zoom\n"
        "  C            - toggle cutaway\n"
        "  +/-          - more/fewer particles\n"
        "  ESC          - quit\n");

    glEnable(GL_DEPTH_TEST);

    double prevTime = glfwGetTime();
    double phase = 0.0;
    bool showXY = true, showZX = false, showZY = false;
    bool showFloor    = true;
    bool showAxis     = true;
    bool showControls = true;
    bool showColorBar = true;
    bool showSpectrum = false;   // Spectrum/cutoff window starts closed (open via View menu)

    // Switch the live simulation domain. Clears every other domain's active /
    // stepping flags first, so leaving microstrip (or geometry) always removes
    // its object, overlay and dimension labels from the scene, then reframes the
    // camera on whatever is now active.
    auto switchDomain = [&](SimDomain d) {
        gDomain = d;
        gUseMicro = false;     gMicroStepping = false;
        gUseFdtd  = false;     gFdtdStepping  = false;
        gUseNumerical = false;
        gBuilderOn = false;
        gPortConfig = false;   // geometry port-config overlay must not leak out
        gCWsweep = false; gCWhold = false; gPulseActive = false; gSweepDone = false;
        switch (d) {
            case SimDomain::Waveguide:  showControls   = true; break;
            case SimDomain::Geometry:   gBuilderOn     = true; break;
            case SimDomain::Microstrip: showMicrostrip = true;
                if (gMicro) { gUseMicro = true; gMicroStepping = true; }
                break;
        }
        const Bounds nb = activeBounds();
        gCamera.setTarget({0.0f, 0.0f, 0.0f});
        gCamera.setDistance(std::max({nb.width, nb.height, nb.depth}) * 2.2f);
        gCloudDirty = true;
    };

    // Field-visualization controls shared by every domain. The cross sections
    // operate on the FieldSource interface, so they work identically for the
    // analytic waveguide, a built geometry and a microstrip. Rendered in the
    // Simulacao panel below whichever domain is active (an "inherited" block).
    auto drawVizControls = [&]() {
        ImGui::SeparatorText("Plot 3D");
        // 0 = intensity only (cloud), 1 = field lines only, 2 = both.
        const char* view3DItems[] = { "Intensidade (nuvem)", "Linhas de campo",
                                      "Intensidade + linhas" };
        ImGui::Combo("Modo 3D", &gView3D, view3DItems, 3);
        // Monte-Carlo cloud resolution = target number of visible points.
        ImGui::SliderInt("Resolucao da nuvem (MC)", &gCloudSamples, 5000, 400000, "%d pontos",
                         ImGuiSliderFlags_Logarithmic);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(
            "Pontos sorteados aleatoriamente no volume (Monte Carlo), nao numa\n"
            "grade -- evita os 'planos' de pontos. E o numero de pontos VISIVEIS\n"
            "buscado (a nuvem preenche mesmo campos finos). Mais pontos = nuvem\n"
            "mais densa e menores. Tambem ajustavel com as teclas + / -.");
        ImGui::Checkbox("Nuvem opaca", &gCloudOpaque);
        if (ImGui::IsItemHovered()) ImGui::SetTooltip(
            "Ligado: os pontos da frente escondem os de tras, e a nuvem tem\n"
            "superficie -- o que voce ve e o campo NAQUELA face.\n"
            "Desligado: nuvem translucida, cada ponto soma com tudo atras dele.\n"
            "Mostra o interior, mas a intensidade lida na tela e um acumulado ao\n"
            "longo da linha de visada, nao o valor local. Use 'Cutaway' para\n"
            "abrir a nuvem opaca e ver o miolo sem perder essa leitura.");
        if (gView3D != 0)   // field lines active -> expose their density
            ImGui::SliderFloat("Densidade de linhas", &gFieldLineDensity, 0.4f, 2.5f, "%.2f");
        ImGui::SeparatorText("Cross sections");
        ImGui::Checkbox("XY (transverse)", &showXY);
        ImGui::Checkbox("ZX plane",        &showZX);
        ImGui::Checkbox("ZY plane",        &showZY);
        ImGui::Checkbox("Field lines (streamlines)", &gFieldLines);
        if (gFieldLines)
        {
            ImGui::SliderFloat("Flow speed", &gFlowSpeed, 0.0f, 4.0f, "%.2f");
            ImGui::SliderFloat("Densidade (secoes)", &gSecLineDensity, 0.3f, 2.5f, "%.2f");
        }
        ImGui::Checkbox("Show slice planes (3D)", &gShowSlicePlanes);
        ImGui::Checkbox("Manual slice position",  &gManualSlice);
        if (gManualSlice) {
            if (showXY) ImGui::SliderFloat("XY slice z", &gSliceZ, 0.0f, 1.0f, "%.2f");
            if (showZX) ImGui::SliderFloat("ZX slice y", &gSliceY, 0.0f, 1.0f, "%.2f");
            if (showZY) ImGui::SliderFloat("ZY slice x", &gSliceX, 0.0f, 1.0f, "%.2f");
        }
    };

    // The active domain only has a field to sample when its own source is live.
    // Without this, active() falls back to the analytic waveguide, so an empty
    // Microstrip/Geometry tab would show the cavity field (cloud + cross sections).
    auto hasFieldForDomain = [&]() {
        return (gDomain == SimDomain::Waveguide) ||
               (gDomain == SimDomain::Geometry   && (gUseNumerical || gUseFdtd)) ||
               (gDomain == SimDomain::Microstrip && gUseMicro && gMicro);
    };

    while (!glfwWindowShouldClose(window))
    {
        const double now = glfwGetTime();
        const double dt = now - prevTime;
        prevTime = now;
        if (gAnimate)
            phase += double(gPhaseSpeed) * dt;
        // Flow-dash clock runs independently of the field animation so the
        // "current" keeps flowing even when the field itself is paused.
        gFlowPhase += dt;

        // Developer "plots aleatorios": advance the rotating-dipole clock.
        if (gDevMode && gShowDipole && uiDipAnimate)
            gDipoleT += dt * double(uiDipSpeed);

        // Advance the time-domain simulation each frame.
        if (gUseFdtd && gFdtd && gFdtdStepping && gAnimate) {
            if (gCWhold) {
                // Continuous CW at one fixed frequency: keep stepping forever so
                // the wave visibly travels through the geometry (no freeze).
                gFdtd->step(gFdtdSteps);
                gFdtd->syncField();
                gCloudDirty = true;
            }
            else if (gCWsweep) {
                // VNA-style CW sweep: hold each frequency until the field reaches
                // steady state (lock-in windows agree) or we hit the block cap,
                // then record S(f_k), store the point, and step to f_{k+1}.
                gFdtd->step(gCWsteps);
                gFdtd->syncField();
                gCloudDirty = true; // cloud shows the CW field at f_k
                if (gFdtd->cwSteady() || gFdtd->cwBlocks() > gCWmaxBlocks) {
                    double cr, ci; gFdtd->cwSourceAmp(cr, ci);
                    // S21 -> first output port; S11 at the input port.
                    int outP = -1, inP = -1;
                    for (int i = 0; i < gFdtd->portCount(); ++i) {
                        if (gFdtd->port(i).role == 2 && outP < 0) outP = i;
                        if (gFdtd->port(i).role == 1 && inP  < 0) inP  = i;
                    }
                    if (gCWindex >= 0 && gCWindex < gCWnf) {
                        double orr = 0, oii = 0, irr = 0, iii = 0;
                        if (outP >= 0) gFdtd->cwPortAmp(outP, orr, oii);
                        if (inP  >= 0) gFdtd->cwPortAmp(inP,  irr, iii);
                        float s21, s11;
                        calS(orr, oii, irr, iii, cr, ci, gCWfreqs[gCWindex], s21, s11);
                        if (outP >= 0) gS21db[gCWindex] = s21;
                        if (inP  >= 0) gS11db[gCWindex] = s11;
                    }
                    ++gCWindex;
                    if (gCWindex >= gCWnf) {
                        gCWsweep = false; gFdtdStepping = false; gSweepDone = true;
                    } else {
                        gFdtd->setCW(gCWfreqs[gCWindex]);
                    }
                }
            } else {
                // Broadband pulse: run until the field rings down, then freeze so
                // the running DFT (and the S-curve) no longer changes.
                gFdtd->step(gFdtdSteps);
                gFdtd->syncField();
                if (gFdtd->pulsePast() && gFdtd->decayRatio() < 1e-3) {
                    gFdtdStepping = false;
                    gSweepDone = true;
                    gCloudDirty = true; // one final bake of the settled field
                }
            }
        }

        // Advance the open microstrip FDTD (quasi-TEM wave traveling down the
        // trace). Continuous stepping so the pulse visibly propagates.
        if (gUseMicro && gMicro) {
            // Honor the global E/H selector, even while paused (re-project the
            // stored fields and re-bake so the switch is immediate).
            if (gMicro->fieldKind() != asField()) {
                gMicro->setDisplayField(asField());
                gMicro->syncField();
                gCloudDirty = true;
            }
            if (gMsSweep) {
                // Broadband pulse: one run fills the whole band. Advance the field,
                // refresh the S-curve from the running DFT (it refines as the pulse
                // rings down), and finish once the transient has essentially left.
                gMicro->step(gMsSweepChunk);
                // The THRU (a uniform line) rings down far sooner than a resonant
                // filter; once its transient is gone its DFT is converged, so stop
                // stepping it -- saves ~half the sweep work on high-Q designs.
                if (gMsThru && !(gMsThru->pulsePast() && gMsThru->pulseDecay() < 1e-2))
                    gMsThru->step(gMsSweepChunk);
                gMsSweepStepsDone += gMsSweepChunk;
                // Refresh the field view only every few chunks during the sweep:
                // syncField and the Monte-Carlo re-bake are each a full-grid pass
                // and otherwise dominate the per-frame overhead (the S-curve the
                // user is watching still refines every chunk from the running DFT).
                static int msViz = 0;
                if (++msViz >= 6) { msViz = 0; gMicro->syncField(); gCloudDirty = true; }
                const int nf = gMicro->pulseFreqCount();
                float rawPassivity = 0.0f, thruPassivity = 0.0f;
                for (int m = 0; m < nf && m < int(gMsS11.size()); ++m) {
                    // Reference Z0 = the THRU line's feed impedance at this freq
                    // (robust even when the filter reflects almost everything).
                    const double z0ref = gMsThru ? gMsThru->pulseZ0(m) : 0.0;
                    double s11r, s11i, s21r, s21i;
                    gMicro->pulseS(m, s11r, s11i, s21r, s21i, z0ref);
                    // Passivity from the RAW S-params (before de-embed): |S11|^2+|S21|^2.
                    rawPassivity = std::max(rawPassivity,
                        float((s11r*s11r+s11i*s11i) + (s21r*s21r+s21i*s21i)));
                    std::complex<double> S11(s11r, s11i), S21(s21r, s21i);
                    // THRU S-params (its own feed Z0) -- used for the passivity
                    // sanity line and the transmission de-embed below.
                    std::complex<double> T21(0.0, 0.0);
                    if (gMsThru) {
                        double t11r, t11i, t21r, t21i;
                        gMsThru->pulseS(m, t11r, t11i, t21r, t21i);
                        thruPassivity = std::max(thruPassivity,
                            float((t11r*t11r+t11i*t11i) + (t21r*t21r+t21i*t21i)));
                        T21 = std::complex<double>(t21r, t21i);
                    }
                    // THRU de-embed (transmission response normalisation): divide
                    // S21 by the equal-length uniform-line reference. S11 stays RAW
                    // (dividing reflection by the through response is not a valid
                    // 1-port de-embed and pushes |S11| above 0 dB).
                    if (gMsDeembed && gMsThru && std::abs(T21) > 1e-9) S21 /= T21;
                    gMsS11[m] = float(20.0*std::log10(std::max(std::abs(S11), 1e-6)));
                    gMsS21[m] = float(20.0*std::log10(std::max(std::abs(S21), 1e-6)));
                }
                gMsRawPassivity = rawPassivity; gMsThruPassivity = thruPassivity;
                // Stop when the pulse has passed and the field rung down, with a
                // hard step cap as a safety net for a slow-decaying resonance.
                if ((gMicro->pulsePast() && gMicro->pulseDecay() < 1e-2) ||
                    gMsSweepStepsDone > 40000) {
                    gMsSweep = false; gMsSweepDone = true; gMicroStepping = false;
                    gMicro->syncField(); gCloudDirty = true;   // final settled view
                }
            } else if (gMicroStepping && gAnimate) {
                gMicro->step(gFdtdSteps);
                gMicro->syncField();
                gCloudDirty = true;
            }
        }

        // Snapshot the pulse S-curve each frame while the pulse method is active,
        // so it persists for comparison against a later VNA sweep. (Frozen once
        // the user starts a CW sweep, since gPulseActive goes false.)
        if (gPulseActive && gFdtd && gFdtd->numFreqs() > 0) {
            const int nf = gFdtd->numFreqs();
            int outP = -1, inP = -1;
            for (int i = 0; i < gFdtd->portCount(); ++i) {
                if (gFdtd->port(i).role == 2 && outP < 0) outP = i;
                if (gFdtd->port(i).role == 1 && inP  < 0) inP  = i;
            }
            gPulseFreqs.resize(nf);
            gPulseS21db.assign(nf, std::numeric_limits<float>::quiet_NaN());
            gPulseS11db.assign(nf, std::numeric_limits<float>::quiet_NaN());
            for (int m = 0; m < nf; ++m) {
                gPulseFreqs[m] = gFdtd->freqHz(m);
                double cr, ci; gFdtd->sourceSpectrum(m, cr, ci);
                double orr = 0, oii = 0, irr = 0, iii = 0;
                if (outP >= 0) gFdtd->portSpectrum(outP, m, orr, oii);
                if (inP  >= 0) gFdtd->portSpectrum(inP,  m, irr, iii);
                float s21, s11;
                calS(orr, oii, irr, iii, cr, ci, gFdtd->freqHz(m), s21, s11);
                if (outP >= 0) gPulseS21db[m] = s21;
                if (inP  >= 0) gPulseS11db[m] = s11;
            }
        }

        // ---- ImGui frame ----
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Full-viewport dock space with a transparent central node: panels dock to
        // the edges while the 3D scene shows through the empty center.
        ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport(),
                                     ImGuiDockNodeFlags_PassthruCentralNode);

        // Ctrl+D toggles developer mode (unless typing in a text field).
        if (!ImGui::GetIO().WantTextInput && ImGui::GetIO().KeyCtrl &&
            ImGui::IsKeyPressed(ImGuiKey_D, false))
            gDevMode = !gDevMode;
        // C toggles the cartesian floor grid.
        if (!ImGui::GetIO().WantTextInput && !ImGui::GetIO().KeyCtrl &&
            ImGui::IsKeyPressed(ImGuiKey_C, false))
            showFloor = !showFloor;

        // When the developer plot is active it takes over the viewport: the whole
        // waveguide visualization (3D field, cross-sections, colorbar, overlays)
        // is suppressed so only the dev scene shows.
        const bool devDipole = gDevMode && gShowDipole;

        // ---- Top menu bar (View) ----
        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("View")) {
                ImGui::MenuItem("Simulation floor",   nullptr, &showFloor);
                ImGui::MenuItem("Coordinate axis",    nullptr, &showAxis);
                ImGui::Separator();
                ImGui::MenuItem("Color scale",        nullptr, &showColorBar);
                ImGui::Separator();
                if (ImGui::BeginMenu("Windows")) {
                    ImGui::MenuItem("XY cross section", nullptr, &showXY);
                    ImGui::MenuItem("ZX cross section", nullptr, &showZX);
                    ImGui::MenuItem("ZY cross section", nullptr, &showZY);
                    ImGui::MenuItem("Spectrum / cutoff", nullptr, &showSpectrum);
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            // Developer tools (enabled with Ctrl+D).
            if (gDevMode) {
                ImGui::Separator();
                ImGui::TextColored(ImVec4(0.9f,0.5f,0.2f,1.0f), "[DEV]");
                if (ImGui::Button("plots aleatorios")) {
                    gShowDipole = !gShowDipole;
                    if (gShowDipole) { gCamera.setTarget({0.0f,0.0f,0.0f}); gCamera.setDistance(4.2f); }
                }
            }
            ImGui::EndMainMenuBar();
        }

        // ---- Domain selector (tabs): picks the one live simulation ----
        if (!devDipole) {
            // Auto-fit the panel to its content, but never taller than the app
            // window (WorkSize.y already excludes the menu bar). Past that cap the
            // window clamps and AlwaysAutoResize hands over to a vertical scrollbar.
            const ImGuiViewport* vp = ImGui::GetMainViewport();
            const float maxPanelH = vp->WorkSize.y - 20.0f;
            ImGui::SetNextWindowPos(ImVec2(10, 26), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSizeConstraints(ImVec2(340, 0),
                                                ImVec2(10000.0f, maxPanelH));
            ImGui::Begin("Simulacao", nullptr,
                         ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_AlwaysAutoResize);
            if (ImGui::BeginTabBar("##domainTabs")) {
                auto domainTab = [&](const char* name, SimDomain d) {
                    if (ImGui::BeginTabItem(name)) {
                        if (gDomain != d) switchDomain(d);
                        ImGui::EndTabItem();
                    }
                };
                domainTab("Waveguide / Cavidade", SimDomain::Waveguide);
                domainTab("Geometria",            SimDomain::Geometry);
                domainTab("Microstrip",           SimDomain::Microstrip);
                ImGui::EndTabBar();
            }
            ImGui::End();
        }

        // ---- [DEV] Rotating magnetic dipole control window ----
        if (gDevMode && gShowDipole) {
            const double PI = 3.14159265358979;
            ImGui::SetNextWindowSize(ImVec2(350, 0), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowPos(ImVec2(20, 60), ImGuiCond_FirstUseEver);
            ImGui::Begin("plots aleatorios - Dipolo magnetico girante");
            ImGui::TextWrapped("Estrela compacta: o momento de dipolo magnetico m esta "
                "desalinhado do eixo de rotacao z pela obliquidade alpha, e gira com periodo P. "
                "So a parte transversal girante (m_x, m_y) radia.");
            ImGui::Separator();
            ImGui::SliderFloat("alpha (graus)", &uiDipAlphaDeg, 0.0f, 90.0f, "%.1f");
            ImGui::SliderFloat("periodo P", &uiDipPeriod, 0.5f, 20.0f, "%.2f");
            ImGui::SliderFloat("|m|", &uiDipMag, 0.1f, 2.0f, "%.2f");
            ImGui::Checkbox("Animar", &uiDipAnimate); ImGui::SameLine();
            ImGui::Checkbox("Linhas de campo", &uiDipFieldLines);
            if (uiDipFieldLines) {
                ImGui::SetNextItemWidth(120); ImGui::SliderInt("linhas/casca", &uiDipFieldAz, 4, 40);
                ImGui::SameLine(); ImGui::SetNextItemWidth(110); ImGui::SliderInt("cascas", &uiDipShellsN, 1, 10);
            }
            ImGui::Checkbox("Cones de radiacao", &uiDipCones); ImGui::SameLine();
            ImGui::SetNextItemWidth(120); ImGui::SliderFloat("compr. cone", &uiDipConeLen, 0.5f, 3.0f, "%.1f");
            ImGui::TextDisabled("abertura do cone ~ |m| (%.0f graus)",
                                std::clamp(double(uiDipMag)*0.35, 0.05, 0.9) * 180.0 / 3.14159265);
            ImGui::SetNextItemWidth(140); ImGui::SliderFloat("velocidade", &uiDipSpeed, 0.1f, 5.0f, "%.1f");

            const double A = uiDipAlphaDeg*PI/180.0;
            const double w = 2.0*PI/std::max(1e-3f, uiDipPeriod);
            const double phi = w*gDipoleT;
            const double m = uiDipMag;
            ImGui::Separator();
            ImGui::Text("m(t) = (%.3f, %.3f, %.3f)",
                        m*std::sin(A)*std::cos(phi), m*std::sin(A)*std::sin(phi), m*std::cos(A));
            ImGui::TextDisabled("m_z = m cos(alpha) = %.3f (estacionario, nao radia)", m*std::cos(A));
            const double mdd2 = 16.0*std::pow(PI,4)*m*m*std::sin(A)*std::sin(A)/std::pow(std::max(1e-3f,uiDipPeriod),4);
            ImGui::Text("|m''|^2 = 16*pi^4*m^2*sin^2(a)/P^4 = %.4g", mdd2);
            ImGui::TextDisabled("(constante no tempo -> potencia radiada ~ constante)");

            static float sx[128], sy[128];
            for (int i=0;i<128;++i){ const double t=uiDipPeriod*i/127.0, ph=w*t;
                sx[i]=float(m*std::sin(A)*std::cos(ph)); sy[i]=float(m*std::sin(A)*std::sin(ph)); }
            const float amp = float(std::max(1e-3, m*std::sin(A)));
            ImGui::PlotLines("m_x(t)", sx, 128, 0, nullptr, -amp, amp, ImVec2(0,45));
            ImGui::PlotLines("m_y(t)", sy, 128, 0, nullptr, -amp, amp, ImVec2(0,45));
            ImGui::TextDisabled("Vista 3D: esfera=estrela, eixo z, cone de precessao, "
                                "vetor m(t) e linhas de campo girantes.");
            ImGui::End();
        }

        bool rebuild = false;
        // relock: a parameter that changes the mode's characteristic frequency
        // was edited, so we re-snap the frequency to f_res (cavity) or f_c
        // (waveguide). The frequency field itself sets rebuild without relock
        // (a manual detune), so a typed/dragged frequency is preserved.
        bool relock = false;
        if (gDomain == SimDomain::Waveguide && !devDipole) {
        // Append into the domain-selector panel: the parametrization for the
        // active tab renders right below the tab bar, in the same window.
        ImGui::Begin("Simulacao");
        ImGui::SeparatorText("Waveguide / Cavidade");

        // Slider + typed-input row: drag the slider OR type a value and press
        // Enter to apply. Typed values are clamped to the [lo,hi] range.
        auto rowF = [&](const char* label, float* v, float lo, float hi,
                        const char* fmt) -> bool {
            bool ch = false;
            ImGui::PushID(label);
            ImGui::SetNextItemWidth(140.0f);
            if (ImGui::SliderFloat("##sl", v, lo, hi, fmt)) ch = true;
            ImGui::SameLine(0.0f, 6.0f);
            ImGui::SetNextItemWidth(68.0f);
            if (ImGui::InputFloat("##in", v, 0.0f, 0.0f, fmt,
                                  ImGuiInputTextFlags_EnterReturnsTrue)) {
                *v = std::clamp(*v, lo, hi); ch = true;
            }
            ImGui::SameLine(0.0f, 6.0f); ImGui::TextUnformatted(label);
            ImGui::PopID();
            return ch;
        };
        auto rowI = [&](const char* label, int* v, int lo, int hi) -> bool {
            bool ch = false;
            ImGui::PushID(label);
            ImGui::SetNextItemWidth(140.0f);
            if (ImGui::SliderInt("##sl", v, lo, hi)) ch = true;
            ImGui::SameLine(0.0f, 6.0f);
            ImGui::SetNextItemWidth(68.0f);
            if (ImGui::InputInt("##in", v, 0, 0, ImGuiInputTextFlags_EnterReturnsTrue)) {
                *v = std::clamp(*v, lo, hi); ch = true;
            }
            ImGui::SameLine(0.0f, 6.0f); ImGui::TextUnformatted(label);
            ImGui::PopID();
            return ch;
        };
        // Frequency row: slider + input like rowF, but the track/grab is colored
        // by band and a zone bar is drawn under it. red = below cutoff or above
        // the recommended range; yellow = cutoff .. recommended min; green =
        // inside the recommended band.
        auto freqRow = [&](double fc, double fmin, double fmax) -> bool {
            const float lo = 1.0f, hi = 60.0f;
            auto zoneOf = [&](double f) -> int {
                if (f < fc) return 0; if (f < fmin) return 1;
                if (f <= fmax) return 2; return 0; };
            const ImVec4 red(0.86f,0.28f,0.28f,1), yel(0.92f,0.82f,0.20f,1),
                         grn(0.32f,0.76f,0.38f,1);
            const int z = zoneOf(uiFreqGHz);
            const ImVec4 zc = (z==0)? red : (z==1)? yel : grn;
            auto sh = [](ImVec4 c, float k){ return ImVec4(c.x*k,c.y*k,c.z*k,1.0f); };
            bool ch = false;
            ImGui::PushID("freqrow");
            ImGui::PushStyleColor(ImGuiCol_SliderGrab,       zc);
            ImGui::PushStyleColor(ImGuiCol_SliderGrabActive, sh(zc,1.15f));
            ImGui::PushStyleColor(ImGuiCol_FrameBg,          sh(zc,0.30f));
            ImGui::PushStyleColor(ImGuiCol_FrameBgHovered,   sh(zc,0.45f));
            ImGui::PushStyleColor(ImGuiCol_FrameBgActive,    sh(zc,0.40f));
            ImGui::SetNextItemWidth(140.0f);
            if (ImGui::SliderFloat("##fsl", &uiFreqGHz, lo, hi, "%.2f")) ch = true;
            ImGui::PopStyleColor(5);
            ImGui::SameLine(0.0f, 6.0f);
            ImGui::SetNextItemWidth(68.0f);
            if (ImGui::InputFloat("##fin", &uiFreqGHz, 0.0f, 0.0f, "%.2f",
                                  ImGuiInputTextFlags_EnterReturnsTrue)) {
                uiFreqGHz = std::clamp(uiFreqGHz, lo, hi); ch = true;
            }
            ImGui::SameLine(0.0f, 6.0f); ImGui::TextUnformatted("Frequency (GHz)");
            // Zone bar aligned under the slider.
            ImGui::Dummy(ImVec2(140.0f, 6.0f));
            const ImVec2 bmin = ImGui::GetItemRectMin(), bmax = ImGui::GetItemRectMax();
            ImDrawList* dl = ImGui::GetWindowDrawList();
            auto xOf = [&](double f){ double t=(f-lo)/(hi-lo); t=std::clamp(t,0.0,1.0);
                                      return bmin.x + float(t)*(bmax.x-bmin.x); };
            const ImU32 cR=ImGui::GetColorU32(red), cY=ImGui::GetColorU32(yel),
                        cG=ImGui::GetColorU32(grn);
            const float xfc=xOf(fc), xmn=xOf(std::max(fc,fmin)),
                        xmx=xOf(std::max(std::max(fc,fmin),fmax));
            dl->AddRectFilled(ImVec2(bmin.x,bmin.y), ImVec2(xfc,bmax.y), cR);
            if (xmn>xfc) dl->AddRectFilled(ImVec2(xfc,bmin.y), ImVec2(xmn,bmax.y), cY);
            if (xmx>xmn) dl->AddRectFilled(ImVec2(xmn,bmin.y), ImVec2(xmx,bmax.y), cG);
            dl->AddRectFilled(ImVec2(xmx,bmin.y), ImVec2(bmax.x,bmax.y), cR);
            const float xf=xOf(uiFreqGHz);
            dl->AddLine(ImVec2(xf,bmin.y-1.0f), ImVec2(xf,bmax.y+1.0f), IM_COL32(20,20,20,255), 1.6f);
            ImGui::PopID();
            return ch;
        };

        ImGui::SliderFloat("Speed (rad/s)", &gPhaseSpeed, 0.0f, 20.0f, "%.2f");
        ImGui::Checkbox("Animate", &gAnimate);
        ImGui::Separator();
        const char* geomItems[] = { "Rectangular", "Cylindrical" };
        if (ImGui::Combo("Geometry", &uiGeometry, geomItems, 2)) { rebuild = true; relock = true; }
        const char* structItems[] = { "Waveguide", "Cavity" };
        if (ImGui::Combo("Structure", &uiStructure, structItems, 2)) { rebuild = true; relock = true; }
        const char* modeItems[] = { "TE", "TM" };
        if (ImGui::Combo("Mode type", &uiModeType, modeItems, 2)) { rebuild = true; relock = true; }
        const char* fieldItems[] = { "Electric (V/m)", "Magnetic (A/m)" };
        if (ImGui::Combo("Field",     &uiField,    fieldItems, 2)) rebuild = true;

        // Filling medium: a table of common microwave dielectrics plus a Custom
        // entry that exposes eps_r / mu_r directly. Changing it relocks the
        // frequency to the new cutoff/resonance (the medium slows the wave).
        {
            struct Diel { const char* name; float er; float mr; };
            static const Diel kDiel[] = {
                {"Vacuum / air",        1.00f, 1.0f},
                {"PTFE (Teflon) 2.1",   2.10f, 1.0f},
                {"Polyethylene 2.25",   2.25f, 1.0f},
                {"Rexolite 2.53",       2.53f, 1.0f},
                {"Fused silica 3.8",    3.80f, 1.0f},
                {"FR-4 4.4",            4.40f, 1.0f},
                {"Alumina 9.8",         9.80f, 1.0f},
                {"Water ~78",          78.00f, 1.0f},
                {"Custom",              0.00f, 0.0f},
            };
            constexpr int nDiel = int(sizeof(kDiel) / sizeof(kDiel[0]));
            constexpr int customIdx = nDiel - 1;
            const char* names[nDiel];
            for (int i = 0; i < nDiel; ++i) names[i] = kDiel[i].name;
            gMediumName = names[std::clamp(uiDielectric, 0, nDiel - 1)];
            // Changing the medium relocks the frequency ONLY for a cavity, where
            // the mode is defined at its resonance f_mnl (which the medium does
            // shift). For an open waveguide the drive frequency is set by the
            // SOURCE and does not physically retune with the fill, so it is kept
            // fixed; only the cutoff and beta (guide wavelength) are recomputed.
            if (ImGui::Combo("Medium", &uiDielectric, names, nDiel)) {
                if (uiDielectric != customIdx) { uiEpsR = kDiel[uiDielectric].er; uiMuR = kDiel[uiDielectric].mr; }
                rebuild = true; relock = isCavity();
            }
            if (uiDielectric == customIdx) {
                ImGui::SetNextItemWidth(90);
                if (ImGui::InputFloat("eps_r", &uiEpsR, 0.0f, 0.0f, "%.3f")) {
                    if (uiEpsR < 1e-3f) uiEpsR = 1e-3f; rebuild = true; relock = isCavity(); }
                ImGui::SameLine(); ImGui::SetNextItemWidth(90);
                if (ImGui::InputFloat("mu_r", &uiMuR, 0.0f, 0.0f, "%.3f")) {
                    if (uiMuR < 1e-3f) uiMuR = 1e-3f; rebuild = true; relock = isCavity(); }
            } else {
                ImGui::TextDisabled("eps_r = %.2f, mu_r = %.2f  (wave slowed x%.2f)",
                                    uiEpsR, uiMuR, std::sqrt(uiEpsR * uiMuR));
            }
        }

        // Rectangular waveguide preset (WR types + custom). Picking a named
        // guide sets its a/b and drops the drive frequency into the recommended
        // band; "Personalizada" exposes the width/height sliders.
        if (uiGeometry == 0) {
            const char* presetLabels[kNumRectPresets];
            for (int i = 0; i < kNumRectPresets; ++i) presetLabels[i] = kRectPresets[i].label;
            if (ImGui::Combo("Waveguide", &uiRectPreset, presetLabels, kNumRectPresets)) {
                if (uiRectPreset != kCustomRectIndex) {
                    uiWidthMM  = kRectPresets[uiRectPreset].a_mm;
                    uiHeightMM = kRectPresets[uiRectPreset].b_mm;
                    if (!isCavity())
                        uiFreqGHz = std::clamp(0.5f * (kRectPresets[uiRectPreset].fmin +
                                                       kRectPresets[uiRectPreset].fmax),
                                               1.0f, 60.0f);
                    rebuild = true;
                    if (isCavity()) relock = true;
                }
            }
        } else {
            const char* presetLabels[kNumCylPresets];
            for (int i = 0; i < kNumCylPresets; ++i) presetLabels[i] = kCylPresets[i].label;
            if (ImGui::Combo("Guia", &uiCylPreset, presetLabels, kNumCylPresets)) {
                if (uiCylPreset != kCustomCylIndex) {
                    uiRadiusMM = kCylPresets[uiCylPreset].R_mm;
                    if (!isCavity())
                        uiFreqGHz = std::clamp(0.5f * (kCylPresets[uiCylPreset].fmin +
                                                       kCylPresets[uiCylPreset].fmax),
                                               1.0f, 60.0f);
                    rebuild = true;
                    if (isCavity()) relock = true;
                }
            }
        }

        if (uiGeometry == 0) {
            // Both indices start at 0: a TE mode may have EITHER index zero
            // (TE10, TE20 ... and equally TE01, TE02 ...), it just cannot have
            // both. M was pinned to >=1 here, which made the whole TE0n family
            // unreachable for no physical reason.
            if (rowI("Mode M", &uiM, 0, 5)) { rebuild = true; relock = true; }
            if (rowI("Mode N", &uiN, 0, 5)) { rebuild = true; relock = true; }
            // Say when the requested pair is not a mode, instead of silently
            // substituting one behind the user's back.
            if (uiModeType == 1 && (uiM < 1 || uiN < 1))
                ImGui::TextColored(ImVec4(0.98f,0.70f,0.10f,1.0f),
                    "TM exige M>=1 e N>=1 (Ez ~ sin.sin some) -> usando TM%d%d",
                    std::max(1, uiM), std::max(1, uiN));
            else if (uiModeType == 0 && uiM == 0 && uiN == 0)
                ImGui::TextColored(ImVec4(0.98f,0.70f,0.10f,1.0f),
                    "TE00 nao existe (kc=0, todo campo se anula) -> usando TE10");
        } else {
            // Cylindrical: n=azimuthal (0..2), m=radial (1..3)
            if (rowI("Mode n (azimuthal)", &uiN, 0, 2)) { rebuild = true; relock = true; }
            if (rowI("Mode m (radial)",    &uiM, 1, 3)) { rebuild = true; relock = true; }
        }
        // Axial index L (half-wave variations along the depth) — cavity only.
        if (isCavity()) {
            if (rowI("Mode L (axial)", &uiL, 0, 5)) { rebuild = true; relock = true; }
        }

        // Cutoff of the current mode and the recommended band (from the preset,
        // or derived as ~1.25..1.9 x f_c) for the coloured frequency slider.
        double fcGHz;
        {
            const double c0 = 299792458.0, kPiC = 3.14159265358979323846;
            const double kc = active()->cutoffWavenumber();
            const double vph = c0 / std::sqrt(active()->epsilonRel() * active()->muRel());
            fcGHz = kc * vph / (2.0 * kPiC) * 1e-9;
        }
        double recMin, recMax;
        if (uiGeometry == 0 && uiRectPreset != kCustomRectIndex) {
            recMin = kRectPresets[uiRectPreset].fmin;
            recMax = kRectPresets[uiRectPreset].fmax;
        } else if (uiGeometry == 1 && uiCylPreset != kCustomCylIndex) {
            recMin = kCylPresets[uiCylPreset].fmin;
            recMax = kCylPresets[uiCylPreset].fmax;
        } else {
            recMin = 1.25 * fcGHz; recMax = 1.90 * fcGHz;
        }

        ImGui::TextDisabled(isCavity()
            ? "cavity: f snaps to resonance on M/N/L change"
            : "green = recommended band, yellow = above cutoff, red = out of band");
        if (freqRow(fcGHz, recMin, recMax)) rebuild = true;

        if (uiGeometry == 0) {
            if (uiRectPreset == kCustomRectIndex) {
                if (rowF("Width a (mm)",  &uiWidthMM,  5.0f, 60.0f, "%.2f")) { rebuild = true; relock = true; }
                if (rowF("Height b (mm)", &uiHeightMM, 2.0f, 40.0f, "%.2f")) { rebuild = true; relock = true; }
            } else {
                ImGui::Text("a = %.2f mm    b = %.2f mm", uiWidthMM, uiHeightMM);
            }
            if (rowF("Depth d (mm)",  &uiDepthMM,  5.0f, 400.0f, "%.1f")) { rebuild = true; relock = true; }
        } else {
            if (uiCylPreset == kCustomCylIndex) {
                if (rowF("Radius (mm)", &uiRadiusMM, 2.0f, 80.0f, "%.2f")) { rebuild = true; relock = true; }
            } else {
                ImGui::Text("R = %.2f mm", uiRadiusMM);
            }
            if (rowF("Length d (mm)", &uiLengthMM, 5.0f, 400.0f, "%.1f")) { rebuild = true; relock = true; }
        }
        // ---- Amplitude scale ----
        // Without this the model is a bare eigenfunction with A = 1 and the
        // "V/m" / "A/m" on the colour bar is a lie. Solving A from the
        // transported power is the same normalization HFSS applies at a wave
        // port (1 W incident by default), so the numbers become comparable.
        ImGui::SeparatorText("Escala de amplitude");
        if (uiGeometry == 0) {
            if (rowF("Potencia (W)", &uiPowerW, 0.001f, 100.0f, "%.3f")) rebuild = true;
            if (rectModel.physicalUnits()) {
                ImGui::TextDisabled("A = %.4g -> campo em %s reais",
                                    rectModel.amplitude(),
                                    (uiField == 0) ? "V/m" : "A/m");
            } else {
                ImGui::TextColored(ImVec4(0.98f,0.70f,0.10f,1.0f),
                    isCavity() ? "Cavidade: nao transporta potencia."
                               : "Abaixo do corte: nao transporta potencia.");
                ImGui::TextDisabled("Amplitude arbitraria (A = 1), escala em u.a.");
            }
        } else {
            ImGui::TextDisabled("Guia cilindrica: amplitude arbitraria (u.a.).");
        }

        ImGui::Separator();
        ImGui::Checkbox("Cutaway", &gCutawayOn);
        ImGui::Text("Nuvem: %d amostras (Monte Carlo)", gCloudSamples);
        // Cross-section / field-viz controls are shared across all domains and
        // rendered once, below, via drawVizControls() (inherited block).
        ImGui::End();
        }

        // ---- Geometry builder window (Phase 2) ----
        if (gDomain == SimDomain::Geometry) {
            ImGui::Begin("Simulacao");   // append below the tab bar
            ImGui::SeparatorText("Geometria construida");
            ImGui::TextDisabled("Assemble a CSG shape (preview; the numerical");
            ImGui::TextDisabled("solver over it is Phase 3). Units: mm.");
            ImGui::Separator();

            int removeIdx = -1;
            for (size_t i = 0; i < gSteps.size(); ++i) {
                BuildStep& s = gSteps[i];
                ImGui::PushID(int(i));
                ImGui::Text("#%d", int(i) + 1); ImGui::SameLine();
                const char* types[] = { "Box", "Cylinder" };
                ImGui::SetNextItemWidth(90); ImGui::Combo("##t", &s.type, types, 2);
                ImGui::SameLine();
                const char* ops[] = { "Add", "Subtract" };
                ImGui::SetNextItemWidth(90); ImGui::Combo("##o", &s.op, ops, 2);
                ImGui::SameLine();
                if (ImGui::SmallButton("remove")) removeIdx = int(i);
                if (s.type == 0) {
                    ImGui::InputFloat3("center xyz", &s.p[0], "%.2f");
                    ImGui::InputFloat3("size xyz",   &s.p[3], "%.2f");
                } else {
                    ImGui::InputFloat2("axis xy",  &s.p[0], "%.2f");
                    ImGui::InputFloat("z0",        &s.p[2], 0, 0, "%.2f");
                    ImGui::InputFloat("z1",        &s.p[3], 0, 0, "%.2f");
                    ImGui::InputFloat("radius",    &s.p[4], 0, 0, "%.2f");
                }
                ImGui::Separator();
                ImGui::PopID();
            }
            if (removeIdx >= 0) gSteps.erase(gSteps.begin() + removeIdx);

            if (ImGui::Button("+ Box"))
                gSteps.push_back({0, 0, {0.f, 0.f, 100.f, 20.f, 10.f, 50.f}});
            ImGui::SameLine();
            if (ImGui::Button("+ Cylinder"))
                gSteps.push_back({1, 0, {0.f, 0.f, 0.f, 100.f, 10.f, 0.f}});
            ImGui::SameLine();
            if (ImGui::Button("Clear")) gSteps.clear();
            if (ImGui::Button("Load T-junction")) {
                gSteps = {
                    {0, 0, {0.f, 0.f,    100.f, 22.86f, 10.16f, 200.f}},
                    {0, 0, {0.f, 30.08f, 100.f, 22.86f, 50.00f, 22.86f}},
                };
            }

            // ---- Save / load shapes ----
            ImGui::Separator();
            static char shapeName[96] = "myshape";
            ImGui::SetNextItemWidth(180);
            ImGui::InputText("name", shapeName, sizeof(shapeName));
            ImGui::SameLine();
            if (ImGui::Button("Save") && shapeName[0] != '\0')
                saveShape(shapeName, gSteps);
            ImGui::SameLine();
            if (ImGui::Button("Load") && shapeName[0] != '\0')
                loadShape(shapeName, gSteps);
            const std::vector<std::string> saved = listShapes();
            if (!saved.empty()) {
                ImGui::TextDisabled("saved shapes (click to load):");
                ImGui::BeginChild("savedShapes", ImVec2(0, 84), true);
                for (const std::string& nm : saved) {
                    if (ImGui::Selectable(nm.c_str())) {
                        loadShape(nm, gSteps);
                        std::snprintf(shapeName, sizeof(shapeName), "%s", nm.c_str());
                    }
                }
                ImGui::EndChild();
            }

            ImGui::Separator();
            const Geometry g = buildGeo();
            const Aabb bb = g.bounds();
            const VoxelMask vm = g.voxelize(24, 24, 48);
            ImGui::Text("AABB: %.1f x %.1f x %.1f mm",
                        bb.sizeX() * 1e3, bb.sizeY() * 1e3, bb.sizeZ() * 1e3);
            ImGui::Text("solid voxels: %zu / %d",
                        vm.solidCount(), vm.nx * vm.ny * vm.nz);

            // ---- Numerical solve (Phase 3 scalar / Phase 4 vector) ----
            ImGui::Separator();
            ImGui::TextDisabled("Numerical eigenmode solver");
            const char* solverItems[] = { "Scalar (Helmholtz)", "Vector (Maxwell)" };
            ImGui::SetNextItemWidth(180);
            ImGui::Combo("Solver", &gSolverType, solverItems, 2);
            ImGui::SetNextItemWidth(140);
            ImGui::SliderInt("Grid res", &gSolveRes, 12, 64);
            if (ImGui::Button("Solve modes")) {
                // Grid proportional to the AABB, longest axis = gSolveRes.
                const double sx = bb.sizeX(), sy = bb.sizeY(), sz = bb.sizeZ();
                const double smax = std::max({sx, sy, sz, 1e-9});
                const int rnx = std::max(6, int(gSolveRes * sx / smax));
                const int rny = std::max(6, int(gSolveRes * sy / smax));
                const int rnz = std::max(6, int(gSolveRes * sz / smax));
                const VoxelMask solveMask = g.voxelize(rnx, rny, rnz, bb);
                if (gSolverType == 0) {
                    HelmholtzResult hr = HelmholtzSolver::solveDirichlet(solveMask, 6, 500, 40);
                    if (!hr.modes.empty()) numModel = std::make_unique<NumericalModel>(hr);
                } else {
                    MaxwellResult mr = MaxwellSolver::solveCavity(solveMask, 6, 1.0, 600, 50);
                    if (!mr.modes.empty()) numModel = std::make_unique<NumericalModel>(mr);
                }
                if (numModel) {
                    gNumMode = 0;
                    gUseNumerical = true;
                    // The numerical field uses [0,size] model coords like the
                    // rectangular model; force that convention so the section /
                    // 3D-view layout (which keys on uiGeometry) lines up.
                    uiGeometry = 0;
                    const Bounds nb = numModel->bounds();
                    gCamera.setTarget({0.0f, 0.0f, 0.0f});
                    gCamera.setDistance(std::max({nb.width, nb.height, nb.depth}) * 2.2f);
                }
            }
            if (numModel) {
                ImGui::SameLine();
                if (ImGui::Button(gUseNumerical ? "Back to analytic" : "Show numerical"))
                    gUseNumerical = !gUseNumerical;
                if (ImGui::SliderInt("Mode", &gNumMode, 0, numModel->modeCount() - 1))
                    numModel->setMode(gNumMode);
                ImGui::Text("f_res = %.3f GHz  (mode %d/%d)",
                            numModel->resonantFrequency() * 1e-9,
                            gNumMode + 1, numModel->modeCount());
            }

            // ---- Time-domain FDTD (Phase 5): wave propagating through it ----
            ImGui::Separator();
            ImGui::TextDisabled("Time-domain (FDTD) — wave through the circuit");
            ImGui::Checkbox("Configure ports (click an opening)", &gPortConfig);
            if (gPortConfig) {
                ImGui::SameLine();
                if (ImGui::Button("Detect ports")) gPorts = detectPorts(buildGeo(), 40);
                ImGui::TextDisabled("green=input, red=output, gray=off; click an opening to cycle");
            }
            // Sweep band: the broadband pulse covers [fmin,fmax]; the port DFTs
            // give S-parameters over exactly this range (choose it around the
            // resonance/cutoff you want to see the dip of).
            ImGui::SetNextItemWidth(90); ImGui::InputFloat("f_min (GHz)", &uiSweepMinGHz, 0.0f, 0.0f, "%.2f");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90); ImGui::InputFloat("f_max (GHz)", &uiSweepMaxGHz, 0.0f, 0.0f, "%.2f");
            if (uiSweepMinGHz < 0.1f) uiSweepMinGHz = 0.1f;
            if (uiSweepMaxGHz <= uiSweepMinGHz) uiSweepMaxGHz = uiSweepMinGHz + 1.0f;
            // Build a fresh FdtdSim from the current geometry/ports/medium.
            auto makeFdtd = [&]() {
                const double sx = bb.sizeX(), sy = bb.sizeY(), sz = bb.sizeZ();
                const double smax = std::max({sx, sy, sz, 1e-9});
                const int rnx = std::max(6, int(gSolveRes * sx / smax));
                const int rny = std::max(6, int(gSolveRes * sy / smax));
                const int rnz = std::max(6, int(gSolveRes * sz / smax));
                const VoxelMask fm = g.voxelize(rnx, rny, rnz, bb);
                if (gPorts.empty()) gPorts = detectPorts(g, 40);
                gFdtd = std::make_unique<FdtdSim>(fm, double(uiSweepMinGHz) * 1e9,
                                                  double(uiSweepMaxGHz) * 1e9, gPorts,
                                                  double(uiEpsR), double(uiMuR));
                gUseFdtd = true; gUseNumerical = false; uiGeometry = 0;
                const Bounds nb = gFdtd->bounds();
                gCamera.setTarget({0.0f, 0.0f, 0.0f});
                gCamera.setDistance(std::max({nb.width, nb.height, nb.depth}) * 2.2f);
            };
            if (ImGui::Button("Run FDTD")) {
                makeFdtd();
                gCWsweep = false; gCWhold = false; gFdtdStepping = true; gSweepDone = false;
                gPulseActive = true;             // snapshot the pulse curve each frame
                gPulseFreqs.clear(); gPulseS21db.clear(); gPulseS11db.clear();
            }
            ImGui::SameLine();
            // VNA-style: measure each frequency to steady state, filling the
            // S-curve left to right (each point final, no whole-curve wobble).
            if (ImGui::Button("Run CW sweep (VNA)")) {
                makeFdtd();
                gPulseActive = false;            // freeze any pulse snapshot for comparison
                gCWhold = false;
                if (gCWnf < 2) gCWnf = 2;
                gCWfreqs.resize(gCWnf);
                for (int m = 0; m < gCWnf; ++m)
                    gCWfreqs[m] = (double(uiSweepMinGHz) + (double(uiSweepMaxGHz) - double(uiSweepMinGHz)) * m / double(gCWnf - 1)) * 1e9;
                gS21db.assign(gCWnf, std::numeric_limits<float>::quiet_NaN());
                gS11db.assign(gCWnf, std::numeric_limits<float>::quiet_NaN());
                gCWindex = 0; gCWsweep = true; gFdtdStepping = true; gSweepDone = false;
                gFdtd->setCW(gCWfreqs[0]);
            }
            ImGui::SetNextItemWidth(140);
            ImGui::SliderInt("sweep points", &gCWnf, 11, 301);

            // Single-frequency CW: drive continuously at one frequency and watch
            // the wave travel through the built geometry (retunable live).
            ImGui::SetNextItemWidth(90);
            if (ImGui::InputFloat("CW freq (GHz)", &uiCWfreqGHz, 0.0f, 0.0f, "%.3f")) {
                if (uiCWfreqGHz < 0.1f) uiCWfreqGHz = 0.1f;
                if (gCWhold && gFdtd) gFdtd->setCW(double(uiCWfreqGHz) * 1e9);
            }
            ImGui::SameLine();
            if (ImGui::Button("Run CW @ freq")) {
                makeFdtd();
                gPulseActive = false; gCWsweep = false; gCWhold = true;
                gFdtdStepping = true; gSweepDone = false;
                gFdtd->setCW(double(uiCWfreqGHz) * 1e9);
            }

            // ---- Thru calibration: measure the incident wave on a straight,
            // matched line driven by the SAME source, so a thru reads 0 dB. -----
            ImGui::Separator();
            if (ImGui::Button("Calibrate reference (thru)")) {
                const FdtdPort* inp = nullptr;
                for (const FdtdPort& p : gPorts) if (p.role == 1) { inp = &p; break; }
                if (!inp) {
                    if (gPorts.empty()) gPorts = detectPorts(g, 40);
                    for (const FdtdPort& p : gPorts) if (p.role == 1) { inp = &p; break; }
                }
                if (inp) {
                    // Straight box: full length along the input axis, port
                    // cross-section transverse, its own bounds start at 0.
                    const int a = inp->axis;
                    const double wu = std::max(1e-4, inp->uMax - inp->uMin);
                    const double wv = std::max(1e-4, inp->vMax - inp->vMin);
                    const double La = (a==0)?bb.sizeX() : (a==1)?bb.sizeY() : bb.sizeZ();
                    double sx, sy, sz;
                    if (a == 2)      { sx = wu; sy = wv; sz = La; }
                    else if (a == 0) { sy = wu; sz = wv; sx = La; }
                    else             { sx = wu; sz = wv; sy = La; }
                    Geometry refG = Geometry::box(sx*0.5, sy*0.5, sz*0.5, sx, sy, sz);
                    const double smax = std::max({sx, sy, sz, 1e-9});
                    const int rnx = std::max(6, int(gSolveRes * sx / smax));
                    const int rny = std::max(6, int(gSolveRes * sy / smax));
                    const int rnz = std::max(6, int(gSolveRes * sz / smax));
                    const VoxelMask refMask = refG.voxelize(rnx, rny, rnz, refG.bounds());
                    std::vector<FdtdPort> refPorts;
                    FdtdPort pi; pi.axis = a; pi.side = inp->side; pi.role = 1;
                    pi.uMin = 0; pi.uMax = wu; pi.vMin = 0; pi.vMax = wv;
                    FdtdPort po = pi; po.side = 1 - inp->side; po.role = 2;
                    refPorts.push_back(pi); refPorts.push_back(po);
                    FdtdSim ref(refMask, double(uiSweepMinGHz)*1e9, double(uiSweepMaxGHz)*1e9,
                                refPorts, double(uiEpsR), double(uiMuR));
                    for (int it = 0; it < 5000; ++it) {   // run the pulse to ring-down
                        ref.step(20); ref.syncField();
                        if (ref.pulsePast() && ref.decayRatio() < 1e-3) break;
                    }
                    const int nf = ref.numFreqs();
                    gCalFreqs.resize(nf); gCalHre.resize(nf); gCalHim.resize(nf);
                    for (int m = 0; m < nf; ++m) {
                        double xir, xii; ref.portSpectrum(0, m, xir, xii); // input = incident
                        double xsr, xsi; ref.sourceSpectrum(m, xsr, xsi);
                        const double sd = xsr*xsr + xsi*xsi + 1e-30;
                        gCalFreqs[m] = ref.freqHz(m);
                        gCalHre[m] = (xir*xsr + xii*xsi)/sd;   // H = Xinc / Xsrc
                        gCalHim[m] = (xii*xsr - xir*xsi)/sd;
                    }
                    gHaveCal = true; gApplyCal = true;
                }
            }
            if (gHaveCal) {
                ImGui::SameLine();
                ImGui::Checkbox("Apply cal", &gApplyCal);
                ImGui::SameLine();
                ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1), "reference OK");
            } else {
                ImGui::SameLine();
                ImGui::TextDisabled("(raw S is uncalibrated; may exceed 0 dB)");
            }

            if (gFdtd) {
                ImGui::SameLine();
                if (ImGui::Button(gFdtdStepping ? "Pause" : "Resume")) gFdtdStepping = !gFdtdStepping;
                ImGui::SameLine();
                if (ImGui::Button("Reset")) {
                    gFdtd->reset(); gFdtdStepping = true; gSweepDone = false;
                    if (gCWsweep) {  // restart the sweep from the first point
                        gCWindex = 0;
                        gS21db.assign(gCWnf, std::numeric_limits<float>::quiet_NaN());
                        gS11db.assign(gCWnf, std::numeric_limits<float>::quiet_NaN());
                        if (!gCWfreqs.empty()) gFdtd->setCW(gCWfreqs[0]);
                    } else if (gCWhold) {          // restart the CW drive cleanly
                        gFdtd->setCW(double(uiCWfreqGHz) * 1e9);
                    }
                }
                if (gCWhold) {
                    ImGui::SetNextItemWidth(140);
                    ImGui::SliderInt("steps/frame", &gFdtdSteps, 1, 40);
                    ImGui::Text("t = %.3f ns  (step %ld)", gFdtd->simTime() * 1e9, gFdtd->stepCount());
                    ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                                       "Continuous CW @ %.3f GHz (edit above to retune).", uiCWfreqGHz);
                } else if (gCWsweep || (gSweepDone && !gCWfreqs.empty() && !gS21db.empty())) {
                    // VNA sweep status.
                    ImGui::SetNextItemWidth(120);
                    ImGui::SliderInt("sweep speed", &gCWsteps, 20, 240);
                    const double fnow = (gCWindex < int(gCWfreqs.size())) ? gCWfreqs[gCWindex] : gCWfreqs.back();
                    if (gCWsweep)
                        ImGui::TextDisabled("VNA sweep: point %d/%d  (%.2f GHz, %d windows, res %.1f%%)",
                                            gCWindex + 1, gCWnf, fnow * 1e-9,
                                            gFdtd->cwBlocks(), gFdtd->cwResidual() * 100.0);
                    else
                        ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                                           "VNA sweep complete (%d points).", gCWnf);
                } else {
                    ImGui::SetNextItemWidth(140);
                    ImGui::SliderInt("steps/frame", &gFdtdSteps, 1, 40);
                    ImGui::Text("t = %.3f ns  (step %ld)", gFdtd->simTime() * 1e9, gFdtd->stepCount());
                    const double dr = gFdtd->decayRatio();
                    const double drDb = 20.0 * std::log10(std::max(dr, 1e-9));
                    if (gSweepDone)
                        ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                                           "Pulse sweep converged (field %.0f dB).", drDb);
                    else if (gFdtdStepping)
                        ImGui::TextDisabled("Pulse %.1f-%.1f GHz... field %.0f dB (waiting for ring-down)",
                                            uiSweepMinGHz, uiSweepMaxGHz, drDb);
                    else
                        ImGui::TextDisabled("Paused.");
                }
            }
            ImGui::End();   // end "Geometry builder" window
        }

        // ---- Microstrip: its own window (split from the builder/waveguide) ----
        // A ground plane + dielectric substrate + copper trace in air, with
        // absorbing (CPML) outer boundaries. A quasi-TEM wave is launched under
        // the trace and travels down the line; the cloud shows the field
        // (including the fringing field) around the strip.
        if (gDomain == SimDomain::Microstrip && !devDipole) {
            ImGui::Begin("Simulacao");   // append below the tab bar
            ImGui::SeparatorText("Microstrip");
            ImGui::TextDisabled("Microstrip (open) — pre-drawn straight line");
            ImGui::SetNextItemWidth(80); ImGui::InputFloat("W strip(mm)", &uiMsStripWmm, 0,0,"%.3f");
            ImGui::SameLine(); ImGui::SetNextItemWidth(80); ImGui::InputFloat("length(mm)", &uiMsLenMM, 0,0,"%.1f");
            ImGui::SetNextItemWidth(80); ImGui::InputFloat("h_sub(mm)", &uiMsHsubMM, 0,0,"%.3f");
            ImGui::SameLine(); ImGui::SetNextItemWidth(70); ImGui::InputFloat("eps_sub", &uiMsEpsSub, 0,0,"%.2f");
            ImGui::SetNextItemWidth(80); ImGui::InputFloat("air(mm)", &uiMsAirMM, 0,0,"%.2f");
            ImGui::SameLine(); ImGui::SetNextItemWidth(70); ImGui::InputFloat("fc(GHz)", &uiMsFcGHz, 0,0,"%.2f");
            ImGui::SetNextItemWidth(140); ImGui::SliderInt("sub cells", &uiMsSubCells, 2, 40);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                "Celulas atraves da altura do substrato (define o dx fino).\n"
                "Custo cresce ~dx^-4; o grid tem um teto que engrossa dx de\n"
                "volta se ficar fino demais. Ctrl+clique para digitar o valor.");
            // Field type shown in 3D / cross-sections. The running sim re-projects
            // immediately (it already keys on this selector), even while paused.
            const char* msFieldItems[] = { "Electric |E|", "Magnetic |H|" };
            ImGui::SetNextItemWidth(140); ImGui::Combo("Field##ms", &uiField, msFieldItems, 2);
            ImGui::Checkbox("Field skin on copper", &gCopperSkin);
            if (gMicro) {
                ImGui::SameLine();
                bool conf = gMicro->conformal();
                if (ImGui::Checkbox("conformal PEC", &conf)) {
                    gMicro->setConformal(conf);
                    if (gMsThru) gMsThru->setConformal(conf);
                }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                    "Dey-Mittra: pondera o laco de Faraday pela fracao da celula\n"
                    "cortada pela borda do cobre, em vez de escadinha (staircase).\n"
                    "Melhora anel/curvas. Desligue para comparar A/B.");
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Colour the conductor by the intensity of the\n"
                                  "selected field (E or H) at its surface, instead\n"
                                  "of a flat copper colour (surface-current view).");

            // ---- 2D cross-section solver (Phase 1: quasi-static Z0 / eps_eff) ----
            ImGui::SetNextItemWidth(140); ImGui::SliderInt("xsec res (cells/H)", &gXsecCells, 8, 96);
            if (ImGui::Button("Solve Z0 / eps_eff (2D)")) {
                const double mm = 1e-3;
                gXsec    = solveMicrostripXsec(uiMsStripWmm*mm, uiMsHsubMM*mm, std::max(1.0f,uiMsEpsSub), gXsecCells);
                gXsecHam = hammerstadMicrostrip(uiMsStripWmm*mm, uiMsHsubMM*mm, std::max(1.0f,uiMsEpsSub));
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Quasi-static finite-difference solve of the line\n"
                                  "cross-section (W strip / h_sub / eps_sub above):\n"
                                  "capacitance with and without the dielectric ->\n"
                                  "eps_eff and Z0. Compared to the Hammerstad formula.");
            if (gXsec.ok) {
                ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                    "FD:        eps_eff = %.3f   Z0 = %.1f ohm", gXsec.eeff, gXsec.Z0);
                ImGui::TextColored(ImVec4(0.6f,0.78f,1.0f,1),
                    "Hammerstad: eps_eff = %.3f   Z0 = %.1f ohm", gXsecHam.eeff, gXsecHam.Z0);
                ImGui::TextDisabled("grid %dx%d, CG %d iters. eps_eff err %.1f%%, Z0 err %.1f%%",
                    gXsec.nz, gXsec.ny, gXsec.iters,
                    100.0*(gXsec.eeff-gXsecHam.eeff)/std::max(1e-9,gXsecHam.eeff),
                    100.0*(gXsec.Z0-gXsecHam.Z0)/std::max(1e-9,gXsecHam.Z0));
            }

            // Build & launch a microstrip from a set of trace footprints. Rects
            // are {x0,x1, w0,w1} mm; circles are {cx,cz,r} mm round patches (the
            // 2D designer). x = propagation, z = transverse width.
            auto runMicrostrip = [&](const std::vector<std::array<double,4>>& traces,
                                     const std::vector<std::array<double,4>>& circles
                                         = std::vector<std::array<double,4>>{}) {
                if (traces.empty() && circles.empty()) return;
                // Bounding box over every rect AND circle (x = propagation, z = width).
                double tx0=1e30, tx1=-1e30, tw0=1e30, tw1=-1e30;
                for (const auto& t : traces) {
                    tx0=std::min(tx0,t[0]); tx1=std::max(tx1,t[1]);
                    tw0=std::min(tw0,t[2]); tw1=std::max(tw1,t[3]);
                }
                for (const auto& c : circles) {
                    tx0=std::min(tx0,c[0]-c[2]); tx1=std::max(tx1,c[0]+c[2]);
                    tw0=std::min(tw0,c[1]-c[2]); tw1=std::max(tw1,c[1]+c[2]);
                }
                // Leftmost / rightmost trace extent -> where the input and output
                // leads + source attach (works for rects and circles alike).
                double inX=1e30, outX=-1e30, inZ0=tw0, inZ1=tw1, outZ0=tw0, outZ1=tw1;
                for (const auto& t : traces) {
                    if (t[0] < inX)  { inX = t[0];  inZ0 = t[2]; inZ1 = t[3]; }
                    if (t[1] > outX) { outX = t[1]; outZ0 = t[2]; outZ1 = t[3]; }
                }
                for (const auto& c : circles) {
                    if (c[0]-c[2] < inX)  { inX = c[0]-c[2];  inZ0 = c[1]-c[2]; inZ1 = c[1]+c[2]; }
                    if (c[0]+c[2] > outX) { outX = c[0]+c[2]; outZ0 = c[1]-c[2]; outZ1 = c[1]+c[2]; }
                }
                using P = MicrostripSim::Prim;
                const double mm=1e-3;
                const double hSub=uiMsHsubMM*mm, air=std::max(uiMsAirMM*mm, 1.0*mm), epsr=std::max(1.0f,uiMsEpsSub);
                double dx = hSub / std::max(2, uiMsSubCells);
                // Margin scales with the substrate (fringing ~ a few H) with a
                // cell-count floor for the CPML, instead of a fixed 6 mm that
                // dwarfs a thin-substrate design (e.g. the eps_r=10, H=0.2 mm LPF).
                const double marX=std::max(12.0*dx, 4.0*hSub), marW=std::max(12.0*dx, 4.0*hSub);
                const double groundT=std::max(dx, 0.4*mm), stripT=dx;
                double domX=(tx1-tx0)*mm + 2*marX;        // propagation (x)
                double domZ=(tw1-tw0)*mm + 2*marW;        // width (z)
                double domY=groundT+hSub+air;             // vertical (y)
                // Cap the (uniform-estimate) cell count. The graded Y mesh uses
                // fewer cells than this estimate, so the real grid stays well under.
                while ((domX/dx)*(domY/dx)*(domZ/dx) > 12.0e6) dx*=1.25; // cap grid
                const double shiftX=marX - tx0*mm, shiftW=marW - tw0*mm;
                const double yt0=groundT+hSub, yt1=groundT+hSub+stripT; // trace layer (y)
                std::vector<P> prims;
                { P g0; g0.kind=0; g0.xmin=0; g0.xmax=domX; g0.zmin=0; g0.zmax=domZ; g0.ymin=0; g0.ymax=groundT; g0.mat=MicrostripSim::Pec; prims.push_back(g0); }
                { P s0; s0.kind=0; s0.xmin=0; s0.xmax=domX; s0.zmin=0; s0.zmax=domZ; s0.ymin=groundT; s0.ymax=yt0; s0.mat=MicrostripSim::Dielectric; s0.epsr=epsr; prims.push_back(s0); }
                for (const auto& t : traces) { P tr; tr.kind=0;
                    tr.xmin=t[0]*mm+shiftX; tr.xmax=t[1]*mm+shiftX;
                    tr.zmin=t[2]*mm+shiftW; tr.zmax=t[3]*mm+shiftW;
                    tr.ymin=yt0; tr.ymax=yt1; tr.mat=MicrostripSim::Pec; prims.push_back(tr); }
                for (const auto& c : circles) { P d; d.kind=2;   // round trace patch (ring if c[3]>0)
                    d.cx=c[0]*mm+shiftX; d.cz=c[1]*mm+shiftW; d.radius=c[2]*mm; d.rinner=std::max(0.0,c[3])*mm;
                    d.ymin=yt0; d.ymax=yt1; d.mat=MicrostripSim::Pec; prims.push_back(d); }
                // Lead-in / lead-out: run the trace into the x-CPML at both ends
                // so the guided wave is absorbed (an approximately matched 2-port
                // for S-parameters) instead of reflecting off an open end.
                { P li; li.kind=0; li.xmin=0; li.xmax=marX+dx; li.ymin=yt0; li.ymax=yt1;
                  li.zmin=inZ0*mm+shiftW; li.zmax=inZ1*mm+shiftW;
                  li.mat=MicrostripSim::Pec; prims.push_back(li); }
                { P lo; lo.kind=0; lo.xmin=domX-marX-dx; lo.xmax=domX; lo.ymin=yt0; lo.ymax=yt1;
                  lo.zmin=outZ0*mm+shiftW; lo.zmax=outZ1*mm+shiftW;
                  lo.mat=MicrostripSim::Pec; prims.push_back(lo); }
                // Source = modal sheet under the input trace's footprint (its width
                // in z), spanning the substrate height in y. The launch port
                // position along x is chosen by the user (clamped past the margin).
                const double srcX=std::clamp(double(uiMsSrcFrac)*domX, marX+dx, domX-marX-dx);
                const double srcZ0=inZ0*mm + shiftW, srcZ1=inZ1*mm + shiftW;
                gMicro = std::make_unique<MicrostripSim>(domX,domY,domZ,dx, prims,
                            srcX, srcZ0, srcZ1, groundT, yt0, double(uiMsFcGHz)*1e9);
                // Ports at the copper extremities (the feed ends just inside the
                // PML), as far toward x=0 / x=domX as valid. Placing them here --
                // not a few cells in -- is what makes the 2-port S-params correct
                // (matches the physical feed-end port of the reference geometry).
                gMicro->setSourcePlaneX(0.0);      // input feed end (x -> 0)
                gMicro->setProbePlaneX(domX);      // output feed end (x -> domX)
                // THRU reference for de-embedding: identical stack + source + ports,
                // but a single straight feed-width trace across the whole domain in
                // place of the filter. Its sense plane is forced to gMicro's actual
                // sense x so both share the same reference planes.
                {
                    std::vector<P> tp;
                    { P g0; g0.kind=0; g0.xmin=0; g0.xmax=domX; g0.zmin=0; g0.zmax=domZ; g0.ymin=0; g0.ymax=groundT; g0.mat=MicrostripSim::Pec; tp.push_back(g0); }
                    { P s0; s0.kind=0; s0.xmin=0; s0.xmax=domX; s0.zmin=0; s0.zmax=domZ; s0.ymin=groundT; s0.ymax=yt0; s0.mat=MicrostripSim::Dielectric; s0.epsr=epsr; tp.push_back(s0); }
                    { P tr; tr.kind=0; tr.xmin=0; tr.xmax=domX; tr.ymin=yt0; tr.ymax=yt1;
                      tr.zmin=inZ0*mm+shiftW; tr.zmax=inZ1*mm+shiftW; tr.mat=MicrostripSim::Pec; tp.push_back(tr); }
                    gMsThru = std::make_unique<MicrostripSim>(domX,domY,domZ,dx, tp,
                                srcX, srcZ0, srcZ1, groundT, yt0, double(uiMsFcGHz)*1e9);
                    gMsThru->setSourcePlaneX(0.0);          // same extremity ports as gMicro
                    gMsThru->setProbePlaneX(domX);
                }
                gUseMicro=true; gMicroStepping=true; gUseFdtd=false; gUseNumerical=false;
                gCWsweep=false; gCWhold=false; gPulseActive=false; uiGeometry=0;
                const Bounds nb=gMicro->bounds();
                gCamera.setTarget({0,0,0});
                gCamera.setDistance(std::max({nb.width,nb.height,nb.depth})*2.2f);
            };

            if (ImGui::Button("Run straight microstrip")) {
                const double L=std::max(2.0f,uiMsLenMM), W=std::max(0.1f,uiMsStripWmm);
                runMicrostrip({ {-L/2, L/2, -W/2, W/2} });
            }
            ImGui::SameLine();
            if (ImGui::Button("From builder boxes")) {
                std::vector<std::array<double,4>> traces;
                for (const BuildStep& s : gSteps) { const float* p = s.p;
                    if (s.type==0) traces.push_back({p[0]-p[3]/2, p[0]+p[3]/2, p[1]-p[4]/2, p[1]+p[4]/2});
                    else           traces.push_back({p[0]-p[3],   p[0]+p[3],   p[1]-p[3],   p[1]+p[3]});
                }
                if (!traces.empty()) runMicrostrip(traces);
            }
            // Stepped-impedance low-pass filter (ref: Estrutura-do-Filtro-Passa-
            // Baixa.webp): a thin high-Z line carrying three wide low-Z shunt
            // stubs, symmetric about the line. Villegas/Duenas ISMOT'09: mm on an
            // eps_r=10, H=0.2 mm substrate (the button sets those), 1-15 GHz.
            ImGui::SetNextItemWidth(90);
            ImGui::InputFloat("filter scale (mm/un)", &uiMsFilterScale, 0, 0, "%.2f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("mm per drawing unit for the low-pass preset\n"
                                  "(1.0 = the paper's millimetre dimensions).\n"
                                  "Bigger -> the wide stubs hit a half-wave\n"
                                  "resonance inside the band (spurious notch).");
            ImGui::SetNextItemWidth(90);
            ImGui::InputFloat("length factor", &uiMsLengthFactor, 0, 0, "%.2f");
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Duenas lengthening factor: stretches the x\n"
                                  "(propagation) lengths only, to correct the\n"
                                  "too-fast numerical velocity of the coarse mesh.\n"
                                  "1.0 = raw geometry; ~1.4 pulls the dips down to\n"
                                  "the paper's frequencies. Widths/impedances unchanged.");
            if (ImGui::Button("Low-pass filter (stubs)")) {
                // Villegas/Duenas ISMOT'09 stepped-impedance LPF: the Fig. 2
                // dimensions are millimetres on an eps_r=10, H=0.2 mm substrate,
                // characterised 1-15 GHz. Set that substrate + band so the preset
                // reproduces the paper (thin feed line ~57 ohm -> ~50 ohm system).
                uiMsEpsSub      = 10.0f;
                uiMsHsubMM      = 0.2f;
                uiMsAirMM       = 1.0f;
                uiMsSubCells    = 6;
                uiMsFilterScale = 1.0f;     // drawing units are millimetres
                uiMsSweepMinGHz = 1.0f;
                uiMsSweepMaxGHz = 15.0f;
                uiMsFcGHz       = 8.0f;
                const double s  = std::max(0.2f, uiMsFilterScale);  // mm per drawing unit
                // Lengthening factor (Duenas): stretch propagation-direction lengths
                // to correct the coarse-mesh velocity; sx applies to x only, so the
                // transverse widths wl/ws (hence the section impedances) stay exact.
                const double Lf = std::max(0.5f, uiMsLengthFactor);
                const double sx = s * Lf;
                const double wl = 0.1435 * s * 0.5;     // half thin-line width (z, no Lf)
                const double ws = 0.99895 * s;          // half stub width (z, no Lf)
                // Short plain feed at each end so the launch plane (port 1) and
                // sense plane (port 2) sit on plain line clear of the stubs.
                const double feed = 1.0 * s;            // feed at each end (mm)
                const double b0  = feed;                // filter body start
                const double x1  = b0 + 0.4498 * sx;    // stub 1 start
                const double x2  = x1 + 0.400  * sx;     // stub 1 end
                const double x3  = x2 + 0.900  * sx;     // stub 2 start
                const double x4  = x3 + 0.500  * sx;     // stub 2 end
                const double x5  = x4 + 0.900  * sx;     // stub 3 start
                const double x6  = x5 + 0.400  * sx;     // stub 3 end
                const double bEnd = x6 + 0.4498 * sx;    // filter body end
                const double xEnd = bEnd + feed;         // output feed end
                // Ports at the extremities: source at the input end, sense at the
                // output end, straddling the whole filter (layoutPorts keeps each
                // reference plane on its own feed line).
                uiMsSrcFrac   = 0.0f;
                uiMsSenseFrac = 1.0f;
                runMicrostrip({
                    { 0.0, xEnd, -wl, wl },  // central thin line: feed + body + feed
                    { x1, x2, -ws, ws },     // wide shunt stub 1
                    { x3, x4, -ws, ws },     // wide shunt stub 2
                    { x5, x6, -ws, ws },     // wide shunt stub 3
                    { bEnd, xEnd, -wl, wl },  // output feed segment (keeps back() thin)
                });
            }
            if (ImGui::IsItemHovered())
                ImGui::SetTooltip("Filtro passa-baixa de impedancia escalonada:\n"
                                  "linha fina (alta-Z) + 3 stubs largos (baixa-Z).");

            // ---- 2D trace designer: dedicated, zoomable window -----------------
            if (ImGui::Button("Open 2D trace designer (top view)")) gMsDesignOpen = true;
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                "Zoomable canvas: drag shapes to move, drag edge/radius handles to\n"
                "resize, read the coupling gap live. Ring radius R is the centerline\n"
                "(outer = R+w/2, inner = R-w/2).");
            if (gMsDesignOpen) {
                ImGui::SetNextWindowSize(ImVec2(560, 580), ImGuiCond_FirstUseEver);
                ImGui::Begin("2D trace designer", &gMsDesignOpen, ImGuiWindowFlags_NoScrollWithMouse);
                ImGuiIO& io = ImGui::GetIO();
                // ---- Ring-resonator presets ------------------------------------
                // A microstrip ring couples to two collinear feed lines across a
                // gap (0 deg / 180 deg). It resonates when the mean circumference
                // 2*pi*R = n*lambda_g, so transmission (S21) peaks at f_n. Tune R,
                // gap and substrate afterwards; smaller gap = stronger coupling.
                auto ringSubstrate = [&]{
                    uiMsEpsSub=4.4f; uiMsHsubMM=0.8f; uiMsAirMM=4.0f; uiMsSubCells=8;
                    uiMsFcGHz=4.0f; uiMsSweepMinGHz=1.0f; uiMsSweepMaxGHz=10.0f;
                    uiMsFilterScale=1.0f; uiMsLengthFactor=1.0f;
                };
                auto makeRing = [&](float R, float w, float feedW, float gap){
                    gMsDesign.clear(); gMsDesignSel=-1;
                    const float ro=R+w*0.5f, feedLen=5.0f;
                    gMsDesign.push_back({1, 0.f,0.f, w,0.f, R, 0.f});             // ring: centerline R, width w
                    const float fx = ro + gap + feedLen*0.5f;
                    gMsDesign.push_back({0, -fx,0.f, feedLen, feedW, 0.f,0.f});   // input feed (left)
                    gMsDesign.push_back({0, +fx,0.f, feedLen, feedW, 0.f,0.f});   // output feed (right)
                };
                auto makeSquareRing = [&](float S, float w, float feedW, float gap){
                    gMsDesign.clear(); gMsDesignSel=-1;
                    const float h2=S*0.5f, outer=h2+w*0.5f, len=S+w, feedLen=5.0f;
                    gMsDesign.push_back({0, 0.f,+h2, len, w, 0.f,0.f});           // top bar
                    gMsDesign.push_back({0, 0.f,-h2, len, w, 0.f,0.f});           // bottom bar
                    gMsDesign.push_back({0, -h2,0.f, w, len, 0.f,0.f});           // left bar
                    gMsDesign.push_back({0, +h2,0.f, w, len, 0.f,0.f});           // right bar
                    const float fx = outer + gap + feedLen*0.5f;
                    gMsDesign.push_back({0, -fx,0.f, feedLen, feedW, 0.f,0.f});   // input feed
                    gMsDesign.push_back({0, +fx,0.f, feedLen, feedW, 0.f,0.f});   // output feed
                };
                ImGui::TextDisabled("Ring resonators (FR-4; tune R / gap / substrate):");
                if (ImGui::Button("Gap-coupled ring##rr"))     { ringSubstrate(); makeRing(6.0f, 1.5f, 1.5f, 0.3f); }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Circular ring, tight 0.3 mm gap (strong coupling,\nclear transmission peaks). R=6 mm, w=1.5 mm.");
                ImGui::SameLine();
                if (ImGui::Button("Loose ring##rr"))           { ringSubstrate(); makeRing(6.0f, 1.5f, 1.5f, 0.8f); }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Circular ring, wide 0.8 mm gap (weak coupling,\nhigher-Q sharper peaks). R=6 mm, w=1.5 mm.");
                ImGui::SameLine();
                if (ImGui::Button("Square loop##rr"))          { ringSubstrate(); makeSquareRing(9.5f, 1.5f, 1.5f, 0.3f); }
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Square-loop resonator (4 bars) - avoids circle\nstaircasing. Side=9.5 mm, w=1.5 mm, gap=0.3 mm.");
                ImGui::Separator();
                if (ImGui::Button("+ Rect##dsg"))   { gMsDesign.push_back({0, 0,0, 4,1, 0,0}); gMsDesignSel=int(gMsDesign.size())-1; }
                ImGui::SameLine();
                if (ImGui::Button("+ Ring##dsg"))   { gMsDesign.push_back({1, 0,0, 1,0, 4,0}); gMsDesignSel=int(gMsDesign.size())-1; }
                ImGui::SameLine();
                if (ImGui::Button("Clear##dsg"))    { gMsDesign.clear(); gMsDesignSel=-1; }
                ImGui::SameLine();
                if (ImGui::Button("Run design##dsg") && !gMsDesign.empty()) {
                    std::vector<std::array<double,4>> rects;
                    std::vector<std::array<double,4>> circs;
                    for (const MsShape& sh : gMsDesign) {
                        if (sh.type==0) rects.push_back({ sh.x-sh.w*0.5, sh.x+sh.w*0.5, sh.z-sh.d*0.5, sh.z+sh.d*0.5 });
                        else {
                            const double outer = std::max(0.1, double(sh.r)+sh.w*0.5);
                            const double inner = double(sh.r)-sh.w*0.5;   // ring if >0, else solid disk
                            circs.push_back({ sh.x, sh.z, outer, inner>0.05 ? inner : 0.0 });
                        }
                    }
                    runMicrostrip(rects, circs);
                }
                ImGui::SameLine();
                if (ImGui::Button("Fit view##dsg")) { gMsZoom=1.0f; gMsPan=ImVec2(0.0f,0.0f); }
                ImGui::SameLine(); ImGui::SetNextItemWidth(120);
                ImGui::SliderFloat("zoom##dsg", &gMsZoom, 0.1f, 40.0f, "%.2fx", ImGuiSliderFlags_Logarithmic);

                // Canvas: top view of the x-z plane. Auto-fit is the zoom=1 basis;
                // the mouse wheel zooms about the cursor, middle/empty drag pans.
                const ImVec2 avail = ImGui::GetContentRegionAvail();
                const ImVec2 cSize(std::max(300.0f, avail.x), std::max(240.0f, avail.y - 140.0f));
                const ImVec2 cA = ImGui::GetCursorScreenPos();
                ImGui::InvisibleButton("##dsgCanvas", cSize,
                    ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonMiddle);
                const bool hov = ImGui::IsItemHovered(), act = ImGui::IsItemActive();
                ImDrawList* dl = ImGui::GetWindowDrawList();
                const ImVec2 cB(cA.x+cSize.x, cA.y+cSize.y);
                dl->PushClipRect(cA, cB, true);
                dl->AddRectFilled(cA, cB, IM_COL32(25,25,30,255));
                dl->AddRect(cA, cB, IM_COL32(90,90,100,255));
                // World bounds = auto-fit basis (zoom=1). Zoom/pan applied on top.
                float wx0=-6, wx1=6, wz0=-4, wz1=4;
                if (!gMsDesign.empty()) {
                    wx0=1e9f; wx1=-1e9f; wz0=1e9f; wz1=-1e9f;
                    for (const MsShape& sh : gMsDesign) {
                        const float hx=(sh.type==0)?sh.w*0.5f:(sh.r+sh.w*0.5f);
                        const float hz=(sh.type==0)?sh.d*0.5f:(sh.r+sh.w*0.5f);
                        wx0=std::min(wx0,sh.x-hx); wx1=std::max(wx1,sh.x+hx);
                        wz0=std::min(wz0,sh.z-hz); wz1=std::max(wz1,sh.z+hz);
                    }
                    const float px=(wx1-wx0)*0.15f+0.5f, pz=(wz1-wz0)*0.15f+0.5f;
                    wx0-=px; wx1+=px; wz0-=pz; wz1+=pz;
                }
                const float wW=std::max(1e-3f,wx1-wx0), wH=std::max(1e-3f,wz1-wz0);
                const float sc0=std::min(cSize.x/wW, cSize.y/wH);
                const float wcx=(wx0+wx1)*0.5f, wcz=(wz0+wz1)*0.5f;
                const ImVec2 cc(cA.x+cSize.x*0.5f, cA.y+cSize.y*0.5f);
                auto s2w=[&](float sx,float sy){ const float s=sc0*gMsZoom;
                    return ImVec2((sx-cc.x-gMsPan.x)/s+wcx, (sy-cc.y-gMsPan.y)/s+wcz); };
                if (hov && io.MouseWheel!=0.0f) {                        // wheel zoom about cursor
                    const ImVec2 mwB = s2w(io.MousePos.x, io.MousePos.y);
                    gMsZoom = std::clamp<float>(gMsZoom * std::pow(1.15f, io.MouseWheel), 0.05f, 200.0f);
                    const float s=sc0*gMsZoom;
                    const ImVec2 sp(cc.x+gMsPan.x+(mwB.x-wcx)*s, cc.y+gMsPan.y+(mwB.y-wcz)*s);
                    gMsPan.x += io.MousePos.x - sp.x; gMsPan.y += io.MousePos.y - sp.y;
                }
                const float sc=sc0*gMsZoom;
                auto w2s=[&](float wx,float wz){ return ImVec2(cc.x+gMsPan.x+(wx-wcx)*sc, cc.y+gMsPan.y+(wz-wcz)*sc); };
                {   // grid (adaptive mm pitch) + z=0 centerline
                    const ImVec2 tl=s2w(cA.x,cA.y), br=s2w(cB.x,cB.y);
                    float gs=0.5f; while (gs*sc<40.0f) gs*=2.0f; while (gs*sc>140.0f) gs*=0.5f;
                    for (float gx=std::ceil(tl.x/gs)*gs; gx<=br.x; gx+=gs)
                        dl->AddLine(w2s(gx,tl.y), w2s(gx,br.y), IM_COL32(44,44,54,255));
                    for (float gz=std::ceil(tl.y/gs)*gs; gz<=br.y; gz+=gs)
                        dl->AddLine(w2s(tl.x,gz), w2s(br.x,gz), IM_COL32(44,44,54,255));
                    dl->AddLine(w2s(tl.x,0.0f), w2s(br.x,0.0f), IM_COL32(85,85,100,255)); // z=0
                }
                for (int i=0;i<int(gMsDesign.size());++i) {
                    const MsShape& sh=gMsDesign[i]; const bool sel=(i==gMsDesignSel);
                    const ImU32 fill=sel?IM_COL32(210,140,60,190):IM_COL32(180,120,50,140);
                    const ImU32 bd  =sel?IM_COL32(255,200,120,255):IM_COL32(210,150,80,220);
                    if (sh.type==0) {
                        ImVec2 a=w2s(sh.x-sh.w*0.5f,sh.z-sh.d*0.5f), b=w2s(sh.x+sh.w*0.5f,sh.z+sh.d*0.5f);
                        dl->AddRectFilled(a,b,fill); dl->AddRect(a,b,bd);
                    } else {
                        const ImVec2 ct=w2s(sh.x,sh.z);
                        const float outer=sh.r+sh.w*0.5f, inner=sh.r-sh.w*0.5f;
                        if (inner>0.02f) {                                 // ring: stroke the annulus at the centerline
                            dl->AddCircle(ct, sh.r*sc, fill, 64, std::max(1.0f,(outer-inner)*sc));
                            dl->AddCircle(ct, outer*sc, bd, 64);
                            dl->AddCircle(ct, inner*sc, bd, 64);
                        } else {
                            dl->AddCircleFilled(ct, outer*sc, fill, 64);
                            dl->AddCircle(ct, outer*sc, bd, 64);
                        }
                    }
                }
                // Resize handles for the selected shape (edges/corners for a rect,
                // outer/inner radius for a ring). Collected for hit-testing below.
                struct Hnd { ImVec2 p; int code; };
                std::vector<Hnd> handles;
                if (gMsDesignSel>=0 && gMsDesignSel<int(gMsDesign.size())) {
                    const MsShape& sh=gMsDesign[gMsDesignSel];
                    if (sh.type==0) {
                        const float hw=sh.w*0.5f, hd=sh.d*0.5f;
                        handles={{w2s(sh.x+hw,sh.z),0},{w2s(sh.x-hw,sh.z),1},
                                 {w2s(sh.x,sh.z+hd),2},{w2s(sh.x,sh.z-hd),3},
                                 {w2s(sh.x+hw,sh.z+hd),4},{w2s(sh.x+hw,sh.z-hd),5},
                                 {w2s(sh.x-hw,sh.z+hd),6},{w2s(sh.x-hw,sh.z-hd),7}};
                    } else {
                        handles={{w2s(sh.x+(sh.r+sh.w*0.5f),sh.z),10},     // outer edge (east)
                                 {w2s(sh.x,sh.z-(sh.r-sh.w*0.5f)),11}};    // inner edge (north)
                    }
                    for (const Hnd& h : handles) {
                        dl->AddCircleFilled(h.p,4.5f,IM_COL32(255,240,180,255));
                        dl->AddCircle(h.p,4.5f,IM_COL32(25,25,30,255));
                    }
                }
                // Live coupling gap: selected shape vs its nearest neighbour, edge
                // to edge along the line of centres, with a dimension line + label.
                auto reachRect=[&](const MsShape& s,float ux,float uz){
                    const float hw=s.w*0.5f, hd=s.d*0.5f, ax=std::fabs(ux), az=std::fabs(uz);
                    return std::min((ax>1e-6f)?hw/ax:1e9f, (az>1e-6f)?hd/az:1e9f);
                };
                if (gMsDesignSel>=0 && gMsDesignSel<int(gMsDesign.size()) && gMsDesign.size()>=2) {
                    const MsShape& S=gMsDesign[gMsDesignSel];
                    int best=-1; float bestGap=1e9f, bUx=0,bUz=0,brS=0,bCd=0,brT=0;
                    for (int j=0;j<int(gMsDesign.size());++j){ if(j==gMsDesignSel) continue;
                        const MsShape& T=gMsDesign[j];
                        const float dx=T.x-S.x, dz=T.z-S.z, cd=std::hypot(dx,dz);
                        if (cd<1e-6f) continue; const float ux=dx/cd, uz=dz/cd;
                        const float rS=(S.type==0)?reachRect(S,ux,uz):(S.r+S.w*0.5f);
                        const float rT=(T.type==0)?reachRect(T,-ux,-uz):(T.r+T.w*0.5f);
                        const float gap=cd-rS-rT;
                        if (gap<bestGap){ bestGap=gap; best=j; bUx=ux;bUz=uz;brS=rS;brT=rT;bCd=cd; }
                    }
                    if (best>=0) {
                        const ImVec2 pA=w2s(S.x+bUx*brS, S.z+bUz*brS);
                        const ImVec2 pB=w2s(S.x+bUx*(bCd-brT), S.z+bUz*(bCd-brT));
                        const ImU32 gc = bestGap>=0? IM_COL32(120,220,255,255):IM_COL32(255,120,120,255);
                        dl->AddLine(pA,pB,gc,1.6f);
                        ImVec2 pp(-(pB.y-pA.y),(pB.x-pA.x)); const float pl=std::hypot(pp.x,pp.y);
                        if (pl>1e-3f){ pp.x=pp.x/pl*4.0f; pp.y=pp.y/pl*4.0f;
                            dl->AddLine(ImVec2(pA.x-pp.x,pA.y-pp.y),ImVec2(pA.x+pp.x,pA.y+pp.y),gc,1.6f);
                            dl->AddLine(ImVec2(pB.x-pp.x,pB.y-pp.y),ImVec2(pB.x+pp.x,pB.y+pp.y),gc,1.6f); }
                        char gl[40]; std::snprintf(gl,sizeof(gl), bestGap>=0?"gap %.3f mm":"overlap %.3f mm", std::fabs(bestGap));
                        dl->AddText(ImVec2((pA.x+pB.x)*0.5f+6,(pA.y+pB.y)*0.5f-14), gc, gl);
                    }
                }
                // Dimension label on the selected shape.
                if (gMsDesignSel>=0 && gMsDesignSel<int(gMsDesign.size())) {
                    const MsShape& sh=gMsDesign[gMsDesignSel]; char lb[72];
                    if (sh.type==0) std::snprintf(lb,sizeof(lb),"%.2f x %.2f mm", sh.w, sh.d);
                    else std::snprintf(lb,sizeof(lb),"R %.2f  w %.2f  (out %.2f / in %.2f)", sh.r, sh.w, sh.r+sh.w*0.5f, sh.r-sh.w*0.5f);
                    dl->AddText(w2s(sh.x,sh.z), IM_COL32(240,240,245,255), lb);
                }
                dl->PopClipRect();

                // ---- Interaction ---------------------------------------------
                const ImVec2 mw = s2w(io.MousePos.x, io.MousePos.y);
                if (act && ImGui::IsMouseClicked(0)) {
                    gMsDragHandle=-1;
                    for (const Hnd& h : handles)                       // 1) grab a handle?
                        if (std::fabs(h.p.x-io.MousePos.x)<7.0f && std::fabs(h.p.y-io.MousePos.y)<7.0f) { gMsDragHandle=h.code; break; }
                    if (gMsDragHandle<0) {                             // 2) select a shape body, else pan
                        int hit=-1;
                        for (int i=int(gMsDesign.size())-1;i>=0;--i){ const MsShape& sh=gMsDesign[i];
                            const bool in=(sh.type==0)?(std::fabs(mw.x-sh.x)<=sh.w*0.5f && std::fabs(mw.y-sh.z)<=sh.d*0.5f)
                                                       :(std::hypot(mw.x-sh.x,mw.y-sh.z)<=sh.r+sh.w*0.5f);
                            if (in){ hit=i; break; } }
                        gMsDesignSel=hit; gMsDragHandle=(hit>=0)?-1:-2;
                    }
                }
                if (act && ImGui::IsMouseDragging(0)) {
                    if (gMsDragHandle==-2) { gMsPan.x+=io.MouseDelta.x; gMsPan.y+=io.MouseDelta.y; }   // pan
                    else if (gMsDesignSel>=0 && gMsDesignSel<int(gMsDesign.size())) {
                        MsShape& sh=gMsDesign[gMsDesignSel];
                        if (gMsDragHandle<0) { sh.x+=io.MouseDelta.x/sc; sh.z+=io.MouseDelta.y/sc; }   // move
                        else if (sh.type==0) {
                            const float L=sh.x-sh.w*0.5f, R=sh.x+sh.w*0.5f, T=sh.z-sh.d*0.5f, B=sh.z+sh.d*0.5f;
                            auto setX=[&](float l,float r){ l=std::min(l,r-0.05f); sh.w=r-l; sh.x=(l+r)*0.5f; };
                            auto setZ=[&](float t,float b){ t=std::min(t,b-0.05f); sh.d=b-t; sh.z=(t+b)*0.5f; };
                            switch(gMsDragHandle){
                                case 0: setX(L,mw.x); break;            case 1: setX(mw.x,R); break;
                                case 2: setZ(T,mw.y); break;            case 3: setZ(mw.y,B); break;
                                case 4: setX(L,mw.x); setZ(T,mw.y); break;
                                case 5: setX(L,mw.x); setZ(mw.y,B); break;
                                case 6: setX(mw.x,R); setZ(T,mw.y); break;
                                case 7: setX(mw.x,R); setZ(mw.y,B); break;
                            }
                        } else {
                            const float radial=std::hypot(mw.x-sh.x, mw.y-sh.z);
                            if (gMsDragHandle==10){ const float inner=sh.r-sh.w*0.5f, outer=std::max(radial,inner+0.05f);
                                sh.r=(outer+inner)*0.5f; sh.w=outer-inner; }
                            else if (gMsDragHandle==11){ const float outer=sh.r+sh.w*0.5f, inner=std::clamp(radial,0.0f,outer-0.05f);
                                sh.r=(outer+inner)*0.5f; sh.w=outer-inner; }
                        }
                    }
                }
                if (ImGui::IsMouseReleased(0)) gMsDragHandle=-1;
                if (hov && ImGui::IsMouseDragging(2)) { gMsPan.x+=io.MouseDelta.x; gMsPan.y+=io.MouseDelta.y; } // middle-drag pan

                // ---- Selected-shape numeric editor ---------------------------
                if (gMsDesignSel>=0 && gMsDesignSel<int(gMsDesign.size())) {
                    MsShape& sh=gMsDesign[gMsDesignSel];
                    ImGui::Text("Shape %d (%s) - mm", gMsDesignSel+1, sh.type==0?"rect":"ring");
                    ImGui::SetNextItemWidth(80); ImGui::InputFloat("x##dsg",&sh.x,0,0,"%.3f"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(80); ImGui::InputFloat("z##dsg",&sh.z,0,0,"%.3f");
                    if (sh.type==0) {
                        ImGui::SetNextItemWidth(80); ImGui::InputFloat("w (x)##dsg",&sh.w,0,0,"%.3f"); ImGui::SameLine();
                        ImGui::SetNextItemWidth(80); ImGui::InputFloat("d (z)##dsg",&sh.d,0,0,"%.3f");
                    } else {
                        ImGui::SetNextItemWidth(80); ImGui::InputFloat("R center##dsg",&sh.r,0,0,"%.3f"); ImGui::SameLine();
                        ImGui::SetNextItemWidth(80); ImGui::InputFloat("width##dsg",&sh.w,0,0,"%.3f");
                        ImGui::TextDisabled("outer = R+w/2 = %.3f, inner = R-w/2 = %.3f mm", sh.r+sh.w*0.5f, sh.r-sh.w*0.5f);
                    }
                    if (ImGui::Button("Delete##dsg")) { gMsDesign.erase(gMsDesign.begin()+gMsDesignSel); gMsDesignSel=-1; }
                } else ImGui::TextDisabled("Click a shape to select; drag to move, drag handles to resize.");
                ImGui::TextDisabled("Wheel = zoom, middle/empty drag = pan. x=propagation (left=in), z=width. Overlaps connect.");

                // ---- Save / load designs to mstrips/<name>.mstrip ------------
                // Mirrors the geometry tab's shape save/load: one line per shape
                // (type x z w d r thick | mm). Persists the drawn trace only, not
                // the substrate / sweep settings.
                auto saveMsDesign = [&](const std::string& name){
                    std::error_code ec; std::filesystem::create_directories("mstrips", ec);
                    std::ofstream f(std::string("mstrips/") + name + ".mstrip");
                    if (!f) return;
                    f << "# microstrip 2D design (type x z w d r thick | mm; type 0=rect 1=ring)\n";
                    for (const MsShape& s : gMsDesign)
                        f << s.type << ' ' << s.x << ' ' << s.z << ' ' << s.w << ' '
                          << s.d << ' ' << s.r << ' ' << s.thick << '\n';
                };
                auto loadMsDesign = [&](const std::string& name)->bool{
                    std::ifstream f(std::string("mstrips/") + name + ".mstrip");
                    if (!f) return false;
                    std::vector<MsShape> out; std::string line;
                    while (std::getline(f, line)) {
                        if (line.empty() || line[0]=='#') continue;
                        std::istringstream is(line); MsShape s;
                        if (is >> s.type >> s.x >> s.z >> s.w >> s.d >> s.r >> s.thick) out.push_back(s);
                    }
                    if (out.empty()) return false;
                    gMsDesign = out; gMsDesignSel = -1; return true;
                };
                auto listMsDesigns = [&](){
                    std::vector<std::string> names; std::error_code ec;
                    if (!std::filesystem::exists("mstrips", ec)) return names;
                    for (const auto& e : std::filesystem::directory_iterator("mstrips", ec))
                        if (e.path().extension()==".mstrip") names.push_back(e.path().stem().string());
                    std::sort(names.begin(), names.end()); return names;
                };
                ImGui::Separator();
                static char msName[96] = "mydesign";
                ImGui::SetNextItemWidth(160); ImGui::InputText("name##msdsg", msName, sizeof(msName));
                ImGui::SameLine();
                if (ImGui::Button("Save##msdsg") && msName[0]!='\0' && !gMsDesign.empty()) saveMsDesign(msName);
                ImGui::SameLine();
                if (ImGui::Button("Load##msdsg") && msName[0]!='\0') loadMsDesign(msName);
                const std::vector<std::string> savedMs = listMsDesigns();
                if (!savedMs.empty()) {
                    ImGui::TextDisabled("saved designs (click to load):");
                    ImGui::BeginChild("savedMsDesigns", ImVec2(0, 74), true);
                    for (const std::string& nm : savedMs)
                        if (ImGui::Selectable(nm.c_str())) {
                            loadMsDesign(nm); std::snprintf(msName, sizeof(msName), "%s", nm.c_str());
                        }
                    ImGui::EndChild();
                }
                ImGui::End();
            }
            if (gMicro) {
                ImGui::SameLine();
                if (ImGui::Button(gMicroStepping ? "Pause##ms" : "Resume##ms")) gMicroStepping=!gMicroStepping;
                ImGui::SameLine();
                if (ImGui::Button("Reset##ms")) { gMicro->reset(); gMicroStepping=true; }
                if (!gUseMicro) { ImGui::SameLine(); if (ImGui::Button("Show##ms")) gUseMicro=true; }
                ImGui::SetNextItemWidth(140); ImGui::SliderInt("steps/frame##ms", &gFdtdSteps, 1, 40);
                // Live retune of the running CW animation (ignored mid-sweep). The
                // field restarts and re-fills at the new frequency.
                ImGui::SetNextItemWidth(140);
                if (ImGui::SliderFloat("anim freq (GHz)##ms", &uiMsFcGHz, 0.5f, 20.0f, "%.2f") && !gMsSweep) {
                    gMicro->setFrequency(double(uiMsFcGHz) * 1e9);
                    gMicroStepping = true;
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Retunes the live CW animation frequency.\n"
                                      "The field relaunches at the new frequency.");
                // Exact-frequency entry (type a value + Enter). Not clamped to the
                // slider range, so you can go below 0.5 or above 20 GHz.
                ImGui::SameLine(); ImGui::SetNextItemWidth(90);
                if (ImGui::InputFloat("= GHz##msfc", &uiMsFcGHz, 0.0f, 0.0f, "%.4f",
                                      ImGuiInputTextFlags_EnterReturnsTrue) && !gMsSweep) {
                    uiMsFcGHz = std::clamp(uiMsFcGHz, 0.001f, 200.0f);
                    gMicro->setFrequency(double(uiMsFcGHz) * 1e9);
                    gMicroStepping = true;
                }
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("Digite a frequencia exata (GHz) e tecle Enter.\n"
                                      "Aceita valores fora da faixa do slider.");

                // Ports are fixed at the extremities (source at the input end,
                // sense at the output end); shown as spheres in 3D.
                ImGui::TextDisabled("ports: source x=%.1f mm (in), sense x=%.1f mm (out)",
                                    gMicro->sourcePlaneX()*1e3, gMicro->probePlaneX()*1e3);

                ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                                   "Microstrip: t=%.3f ns (step %ld)", gMicro->simTime()*1e9, gMicro->stepCount());
                ImGui::TextDisabled("Plot: %s (use the Field selector to switch)",
                                    gMicro->fieldKind()==FieldKind::Magnetic ? "|H| (A/m)" : "|E| (V/m)");

                // Transmission-line probe: launched current and per-length line
                // parameters extracted from the field (L' from the magnetic
                // energy, C' from the electric energy). These are FDTD estimates
                // that converge as the mesh ("sub cells") is refined.
                if (gMicro->stripCurrent() > 0.0) {
                    ImGui::Separator();
                    ImGui::Text("port  I = %.4g mA    V = %.4g mV",
                                gMicro->stripCurrent()*1e3, gMicro->portVoltage()*1e3);
                    ImGui::Text("L' = %.4g nH/m    C' = %.4g pF/m",
                                gMicro->inductancePerLength()*1e9, gMicro->capacitancePerLength()*1e12);
                    ImGui::SameLine();
                    ImGui::TextColored(ImVec4(0.6f,0.78f,1.0f,1), "   Z0 = %.1f ohm", gMicro->lineImpedance());
                    if (ImGui::IsItemHovered())
                        ImGui::SetTooltip("I = closed-loop integral of H around the trace\n"
                                          "V = -integral of E through the substrate\n"
                                          "L' = 2*Wm'/I^2,  C' = 2*We'/V^2,  Z0 = sqrt(L'/C')\n"
                                          "Estimates; raise 'sub cells' for better accuracy.");
                } else {
                    ImGui::TextDisabled("(L'/C'/Z0: run a few hundred steps to populate)");
                }

                // ---- Frequency sweep: Z0 and eps_eff vs frequency ----
                ImGui::Separator();
                ImGui::TextDisabled("Frequency sweep -> S11 (reflection), S21 (transmission)");
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip(
                        "One broadband Gaussian pulse excites the whole band; a running\n"
                        "DFT reads the V/I phasors at port 1 (launch) and port 2 (sense)\n"
                        "at every frequency at once -- ~10-40x faster than a CW-per-point\n"
                        "sweep. TEM split V+-=(V+-Z0 I)/2 gives S11=V1-/V1+, S21=V2+/V1+.\n"
                        "The trace runs into the CPML at both ends (~matched loads).\n"
                        "Qualitative: an ~-10 dB S11 / ~-2 dB S21 floor on a plain line\n"
                        "comes from the CPML mode match + coarse mesh; raise 'sub cells'.");
                ImGui::SetNextItemWidth(70); ImGui::InputFloat("f0(GHz)", &uiMsSweepMinGHz, 0,0,"%.1f");
                ImGui::SameLine(); ImGui::SetNextItemWidth(70); ImGui::InputFloat("f1(GHz)", &uiMsSweepMaxGHz, 0,0,"%.1f");
                ImGui::SameLine(); ImGui::SetNextItemWidth(90); ImGui::SliderInt("points", &uiMsSweepPts, 3, 40);
                ImGui::SameLine(); ImGui::Checkbox("de-embed##ms", &gMsDeembed);
                if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                    "THRU normalisation: divide the sweep by a straight feed-width line\n"
                    "of the same length (filter removed), cancelling the access-line\n"
                    "phase and loss. S11/S21 are then referenced to an equal-length\n"
                    "uniform line (0 dB S21 = as good as a plain through line).");
                if (!gMsSweep) {
                    if (ImGui::Button("Run sweep##ms")) {
                        const int n = std::clamp(uiMsSweepPts, 2, 60);
                        const float f0 = std::min(uiMsSweepMinGHz, uiMsSweepMaxGHz);
                        const float f1 = std::max(uiMsSweepMinGHz, uiMsSweepMaxGHz);
                        gMsFreqGHz.resize(n); gMsZ0.assign(n,0); gMsEeff.assign(n,0);
                        gMsLp.assign(n,0); gMsCp.assign(n,0); gMsS11.assign(n,0); gMsS21.assign(n,0);
                        for (int i=0;i<n;i++) gMsFreqGHz[i] = f0 + (f1-f0)*i/float(std::max(1,n-1));
                        // Broadband pulse: one run covers the whole band instead of
                        // a separate CW steady state per frequency (~10-40x faster).
                        std::vector<double> fHz(n);
                        for (int i=0;i<n;i++) fHz[i] = double(gMsFreqGHz[i]) * 1e9;
                        gMicro->startPulseSweep(fHz);
                        if (gMsThru) gMsThru->startPulseSweep(fHz);  // shadow run for de-embedding
                        gMsSweepIdx=n;             // all points refine together
                        gMsSweepStepsDone=0; gMsSweep=true; gMsSweepDone=false; gMsWindowOpen=false;
                        gMicroStepping=true;
                    }
                } else {
                    if (ImGui::Button("Stop sweep##ms")) { gMsSweep=false; gMicroStepping=false; }
                    ImGui::SameLine();
                    const double dec = gMicro ? gMicro->pulseDecay() : 1.0;
                    ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                        "pulse ring-down: field %.0f dB  (%d steps)",
                        20.0*std::log10(std::max(dec, 1e-6)), gMsSweepStepsDone);
                }
                const int shown = gMsSweep ? gMsSweepIdx : int(gMsZ0.size());
                if (shown > 0) {
                    char ov[64];
                    // ---- S11/S21 plot (dB vs GHz), styled like the FDTD geometry
                    // sweep: dB axis + gridlines, frequency axis in GHz, a marker at
                    // the frequency being swept, and auto S21-peak / S11-dip markers.
                    float dmin=1e30f, dmax=-1e30f;
                    for (int i=0;i<shown;i++)
                        for (float v : {gMsS11[i], gMsS21[i]}) {
                            if (v!=v) continue; const float c=std::clamp(v,-100.0f,10.0f);
                            dmin=std::min(dmin,c); dmax=std::max(dmax,c); }
                    float dbLo, dbHi;
                    if (dmax<dmin){ dbLo=-40.0f; dbHi=5.0f; }
                    else { dbLo=dmin-5.0f; dbHi=dmax+5.0f; if (dbHi-dbLo<1.0f) dbHi=dbLo+1.0f; }
                    double fLo=gMsFreqGHz.front(), fHi=gMsFreqGHz.back();
                    if (fHi<=fLo) fHi=fLo+1.0;

                    const ImU32 colOut=IM_COL32(90,220,120,255);   // S21 (green)
                    const ImU32 colRefl=IM_COL32(240,150,80,255);  // S11 (orange)
                    ImDrawList* dl=ImGui::GetWindowDrawList();
                    const ImVec2 org=ImGui::GetCursorScreenPos();
                    const float pw=ImGui::GetContentRegionAvail().x-8.0f;
                    const float ph=190.0f;
                    const ImVec2 p0=org, p1=ImVec2(org.x+pw, org.y+ph);
                    dl->AddRectFilled(p0,p1,IM_COL32(20,20,24,255));
                    dl->AddRect(p0,p1,IM_COL32(120,120,120,255));
                    auto yOf=[&](double db){ return p1.y-float(std::clamp(db,double(dbLo),double(dbHi))-dbLo)/(dbHi-dbLo)*ph; };
                    auto xOfF=[&](double f){ return p0.x+float((f-fLo)/(fHi-fLo))*pw; };
                    auto isNum=[](float v){ return v==v; };
                    {   // dB gridlines
                        const double span=std::max(double(dbHi-dbLo),1e-9), rough=span/5.0;
                        const double mag=std::pow(10.0,std::floor(std::log10(rough))), norm=rough/mag;
                        const double step=(norm<1.5?1.0:norm<3.0?2.0:norm<7.0?5.0:10.0)*mag;
                        for (double d=std::ceil(dbLo/step)*step; d<=dbHi+step*1e-6; d+=step){
                            const float y=yOf(d);
                            dl->AddLine(ImVec2(p0.x,y),ImVec2(p1.x,y),IM_COL32(60,60,66,255));
                            char lb[12]; std::snprintf(lb,sizeof lb,"%g",d);
                            dl->AddText(ImVec2(p0.x+2,y-12),IM_COL32(140,140,140,255),lb); }
                    }
                    dl->AddText(ImVec2(p0.x+2,p0.y+2),IM_COL32(170,170,175,255),"|S| (dB)");
                    {   // frequency gridlines (GHz)
                        const double span=std::max(fHi-fLo,1e-9), rough=span/6.0;
                        const double mag=std::pow(10.0,std::floor(std::log10(rough))), norm=rough/mag;
                        const double step=(norm<1.5?1.0:norm<3.0?2.0:norm<7.0?5.0:10.0)*mag;
                        const int dec = step>=1.0?1:step>=0.1?2:3;
                        char fmt[8]; std::snprintf(fmt,sizeof fmt,"%%.%df",dec);
                        for (double f=std::ceil(fLo/step)*step; f<=fHi+step*1e-6; f+=step){
                            const float x=xOfF(f);
                            dl->AddLine(ImVec2(x,p0.y),ImVec2(x,p1.y),IM_COL32(55,55,62,255),1.0f);
                            char lb[16]; std::snprintf(lb,sizeof lb,fmt,f);
                            const ImVec2 ts=ImGui::CalcTextSize(lb);
                            dl->AddText(ImVec2(x-ts.x*0.5f,p1.y+3.0f),IM_COL32(150,150,155,255),lb); }
                        const char* xlab="Frequency (GHz)";
                        const ImVec2 xs=ImGui::CalcTextSize(xlab);
                        dl->AddText(ImVec2(p0.x+(pw-xs.x)*0.5f,p1.y+17.0f),IM_COL32(170,170,175,255),xlab);
                    }
                    auto plotF=[&](const std::vector<float>& arr, ImU32 col){
                        ImVec2 prev; bool have=false;
                        for (int m=0;m<shown;m++){
                            if (!isNum(arr[m])){ have=false; continue; }
                            const ImVec2 pt(xOfF(gMsFreqGHz[m]), yOf(arr[m]));
                            if (have) dl->AddLine(prev,pt,col,2.0f);
                            dl->AddCircleFilled(pt,2.0f,col); prev=pt; have=true; } };
                    plotF(gMsS21, colOut);
                    plotF(gMsS11, colRefl);
                    if (gMsSweep && gMsSweepIdx<int(gMsFreqGHz.size())){   // current sweep freq
                        const float x=xOfF(gMsFreqGHz[gMsSweepIdx]);
                        dl->AddLine(ImVec2(x,p0.y),ImVec2(x,p1.y),IM_COL32(230,230,120,180),1.5f); }
                    auto extreme=[&](const std::vector<float>& arr, bool wantMax, double& fO, float& dO){
                        bool f=false; fO=0; dO=0;
                        for (int m=0;m<shown;m++){ if(!isNum(arr[m])) continue;
                            if(!f||(wantMax?arr[m]>dO:arr[m]<dO)){ dO=arr[m]; fO=gMsFreqGHz[m]; f=true; } }
                        return f; };
                    double fPk,fDip; float dbPk,dbDip;
                    const bool hasPk=extreme(gMsS21,true,fPk,dbPk);
                    const bool hasDip=extreme(gMsS11,false,fDip,dbDip);
                    auto marker=[&](double f, float db, ImU32 col, const char* tag, float ty){
                        const float x=xOfF(f);
                        dl->AddLine(ImVec2(x,p0.y),ImVec2(x,p1.y),col,1.0f);
                        dl->AddCircleFilled(ImVec2(x,yOf(db)),3.5f,col);
                        char t[48]; std::snprintf(t,sizeof t,"%s %.3f GHz",tag,f);
                        float tx=x+4; if (tx>p1.x-90) tx=x-92;
                        dl->AddText(ImVec2(tx,p0.y+ty),col,t); };
                    if (hasPk)  marker(fPk, dbPk, colOut,  "S21^", 4.0f);
                    if (hasDip) marker(fDip,dbDip,colRefl, "S11v", 20.0f);
                    ImGui::Dummy(ImVec2(pw, ph+36.0f));
                    if (hasDip) ImGui::TextColored(ImVec4(0.94f,0.59f,0.31f,1),"S11 dip: %.3f GHz (%.1f dB)",fDip,dbDip);
                    if (hasPk)  ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),"S21 peak: %.3f GHz (%.1f dB)",fPk,dbPk);
                    ImGui::TextDisabled("last: S11 %.1f dB, S21 %.1f dB @ %.1f GHz",
                        gMsS11[shown-1], gMsS21[shown-1], gMsFreqGHz[std::max(0,shown-1)]);
                    // Passivity: max |S11|^2+|S21|^2 (raw). <= 1 = energy conserved.
                    if (gMsRawPassivity > 0.0f) {
                        const bool ok = gMsRawPassivity <= 1.03f;   // ~3% numerical slack
                        ImGui::TextColored(ok ? ImVec4(0.55f,0.8f,0.95f,1) : ImVec4(1.0f,0.45f,0.35f,1),
                            "passividade (filtro) max |S11|^2+|S21|^2 = %.3f  %s",
                            gMsRawPassivity, ok ? "OK (<=1)" : "> 1: nao-passivo!");
                    }
                    if (gMsThruPassivity > 0.0f) {
                        const bool tok = gMsThruPassivity <= 1.03f;
                        ImGui::TextColored(tok ? ImVec4(0.55f,0.8f,0.95f,1) : ImVec4(1.0f,0.6f,0.3f,1),
                            "passividade (THRU reta) = %.3f  %s",
                            gMsThruPassivity, tok ? "OK -> extracao consistente" : "> 1 -> ARTEFATO de extracao (nao o campo)");
                    }

                    // Z0(f)/eps_eff(f) come from the CW lock-in's L'/C' extraction;
                    // the broadband pulse sweep leaves them zero, so only plot them
                    // when real data is present.
                    float z0min=1e30f,z0max=-1e30f,emin=1e30f,emax=-1e30f;
                    for (int i=0;i<shown;i++){ z0min=std::min(z0min,gMsZ0[i]); z0max=std::max(z0max,gMsZ0[i]);
                                               emin=std::min(emin,gMsEeff[i]); emax=std::max(emax,gMsEeff[i]); }
                    if (z0max > 1.0f) {
                        snprintf(ov,sizeof ov,"Z0 %.0f-%.0f ohm",z0min,z0max);
                        ImGui::PlotLines("Z0(f)##ms", gMsZ0.data(), shown, 0, ov, z0min*0.9f, z0max*1.1f, ImVec2(320,50));
                        snprintf(ov,sizeof ov,"eps_eff %.2f-%.2f",emin,emax);
                        ImGui::PlotLines("eeff(f)##ms", gMsEeff.data(), shown, 0, ov, emin*0.9f, emax*1.1f, ImVec2(320,50));
                    }
                    ImGui::TextDisabled("x: %.1f -> %.1f GHz", gMsFreqGHz.front(), gMsFreqGHz[std::max(0,shown-1)]);
                }
            }
            ImGui::End();
        }

        // ---- Shared viz controls (cross sections): appended to the Simulacao
        // panel for every domain, so they are available in all simulations ----
        if (!devDipole) {
            ImGui::Begin("Simulacao");
            drawVizControls();
            ImGui::End();
        }

        // ---- Port signals window (FDTD): waveform + relative power per port ----
        if (gUseFdtd && gFdtd && !devDipole) {
            ImGui::SetNextWindowSize(ImVec2(380, 0), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowPos(ImVec2(20, 420), ImGuiCond_FirstUseEver);
            ImGui::Begin("Port signals (FDTD)");
            const char* axName[3] = {"x","y","z"};
            const char* roleName[3] = {"off","IN","OUT"};
            float gmax = 1e-12f;
            for (int i = 0; i < gFdtd->portCount(); ++i)
                for (float v : gFdtd->portHistory(i)) gmax = std::max(gmax, std::abs(v));
            if (gFdtd->portCount() == 0)
                ImGui::TextDisabled("No ports configured (mark openings as IN/OUT).");
            for (int i = 0; i < gFdtd->portCount(); ++i) {
                const FdtdPort& p = gFdtd->port(i);
                const std::vector<float>& h = gFdtd->portHistory(i);
                char label[80];
                std::snprintf(label, sizeof(label), "P%d %s%s %s   rms %.3f",
                              i+1, axName[p.axis], p.side?"+":"-", roleName[p.role], gFdtd->portRms(i));
                ImGui::PushID(i);
                ImGui::PlotLines("##sig", h.data(), int(h.size()), gFdtd->portHistoryPos(i),
                                 label, -gmax, gmax, ImVec2(0, 58));
                ImGui::PopID();
            }
            if (gFdtd->portCount() > 0)
                ImGui::TextDisabled("IN = source; OUT = transmitted. rms ~ relative power.");
            ImGui::End();
        }

        // ---- S-parameters window: |S| in dB vs frequency --------------------
        // Two sources feed this plot:
        //  - VNA CW sweep (gCWfreqs set): one steady-state point per frequency,
        //    filling left to right; each point is final.
        //  - Broadband pulse (fallback): S(f)=X_port(f)/X_src(f) from the running
        //    port DFTs (whole curve refines together as the pulse rings down).
        // Both are 20*log10|S|, so resonances appear as dips.
        {
            const bool haveCW    = !gCWfreqs.empty();
            const bool havePulse = !gPulseFreqs.empty();
            if (gUseFdtd && (haveCW || havePulse) && !devDipole) {
            ImGui::SetNextWindowSize(ImVec2(480, 320), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowPos(ImVec2(410, 420), ImGuiCond_FirstUseEver);
            ImGui::Begin("S-parameters (dB)");

            // Auto-scale the dB axis to the data: (max dB)+5 on top, (min dB)-5 at
            // the bottom. Values are clamped to a sane window first so numerical
            // zeros / band-edge artifacts don't blow up the range.
            float dmin = 1e30f, dmax = -1e30f;
            auto scanRange = [&](const std::vector<float>& a){
                for (float v : a) { if (v != v) continue;
                    const float c = std::clamp(v, -100.0f, 10.0f);
                    dmin = std::min(dmin, c); dmax = std::max(dmax, c); } };
            scanRange(gS21db); scanRange(gS11db); scanRange(gPulseS21db); scanRange(gPulseS11db);
            float dbLo, dbHi;
            if (dmax < dmin) { dbLo = -40.0f; dbHi = 5.0f; }         // no data yet
            else { dbLo = dmin - 5.0f; dbHi = dmax + 5.0f;
                   if (dbHi - dbLo < 1.0f) dbHi = dbLo + 1.0f; }

            const ImU32 colOut  = IM_COL32(90, 220, 120, 255);   // S21 (green)
            const ImU32 colRefl = IM_COL32(240, 150, 80, 255);   // S11 (orange)
            const ImU32 colOutF  = IM_COL32(90, 220, 120, 120);  // faint = pulse
            const ImU32 colReflF = IM_COL32(240, 150, 80, 120);

            // Frequency axis: union of whatever data exists.
            double fLo = 1e30, fHi = -1e30;
            if (haveCW)    { fLo = std::min(fLo, gCWfreqs.front());    fHi = std::max(fHi, gCWfreqs.back()); }
            if (havePulse) { fLo = std::min(fLo, gPulseFreqs.front()); fHi = std::max(fHi, gPulseFreqs.back()); }
            if (fHi <= fLo) { fLo = double(uiSweepMinGHz)*1e9; fHi = double(uiSweepMaxGHz)*1e9; }

            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImVec2 org = ImGui::GetCursorScreenPos();
            const float pw = ImGui::GetContentRegionAvail().x - 8.0f;
            const float ph = 200.0f;
            const ImVec2 p0 = org, p1 = ImVec2(org.x + pw, org.y + ph);
            dl->AddRectFilled(p0, p1, IM_COL32(20,20,24,255));
            dl->AddRect(p0, p1, IM_COL32(120,120,120,255));
            auto yOf   = [&](double db){ return p1.y - float(std::clamp(db,double(dbLo),double(dbHi)) - dbLo)/(dbHi-dbLo)*ph; };
            auto xOfF  = [&](double f ){ return p0.x + float((f - fLo)/(fHi - fLo))*pw; };
            auto isNum = [](float v){ return v == v; };

            // Horizontal dB gridlines at "nice" steps across the auto-scaled range.
            {
                const double span = std::max(double(dbHi - dbLo), 1e-9);
                const double rough = span / 5.0;
                const double mag = std::pow(10.0, std::floor(std::log10(rough)));
                const double norm = rough / mag;
                const double step = (norm < 1.5 ? 1.0 : norm < 3.0 ? 2.0 : norm < 7.0 ? 5.0 : 10.0) * mag;
                for (double d = std::ceil(dbLo / step) * step; d <= dbHi + step * 1e-6; d += step) {
                    const float y = yOf(d);
                    dl->AddLine(ImVec2(p0.x, y), ImVec2(p1.x, y), IM_COL32(60,60,66,255));
                    char lb[12]; std::snprintf(lb, sizeof(lb), "%g", d);
                    dl->AddText(ImVec2(p0.x + 2, y - 12), IM_COL32(140,140,140,255), lb);
                }
            }
            // Y-axis label (top-left of the plot).
            dl->AddText(ImVec2(p0.x + 2, p0.y + 2), IM_COL32(170,170,175,255), "|S| (dB)");

            // Vertical frequency gridlines at "nice" GHz steps, each labeled with
            // its frequency just below the axis.
            {
                const double span = std::max(fHi - fLo, 1e-30);
                const double rough = span / 6.0;
                const double mag = std::pow(10.0, std::floor(std::log10(rough)));
                const double norm = rough / mag;
                const double step = (norm < 1.5 ? 1.0 : norm < 3.0 ? 2.0 : norm < 7.0 ? 5.0 : 10.0) * mag;
                const double stepGHz = step * 1e-9;
                const int dec = stepGHz >= 1.0 ? 1 : stepGHz >= 0.1 ? 2 : 3;
                char fmt[8]; std::snprintf(fmt, sizeof(fmt), "%%.%df", dec);
                for (double f = std::ceil(fLo / step) * step; f <= fHi + step * 1e-6; f += step) {
                    const float x = xOfF(f);
                    dl->AddLine(ImVec2(x, p0.y), ImVec2(x, p1.y), IM_COL32(55,55,62,255), 1.0f);
                    char lb[16]; std::snprintf(lb, sizeof(lb), fmt, f * 1e-9);
                    const ImVec2 ts = ImGui::CalcTextSize(lb);
                    dl->AddText(ImVec2(x - ts.x * 0.5f, p1.y + 3.0f), IM_COL32(150,150,155,255), lb);
                }
                // X-axis label, centered below the frequency numbers.
                const char* xlab = "Frequency (GHz)";
                const ImVec2 xs = ImGui::CalcTextSize(xlab);
                dl->AddText(ImVec2(p0.x + (pw - xs.x) * 0.5f, p1.y + 17.0f),
                            IM_COL32(170,170,175,255), xlab);
            }

            // Draw one dB-vs-frequency trace (skips NaN gaps).
            auto plotF = [&](const std::vector<double>& fs, const std::vector<float>& arr,
                             ImU32 col, float thick, bool dots) {
                ImVec2 prev; bool have = false;
                const int n = int(std::min(fs.size(), arr.size()));
                for (int m = 0; m < n; ++m) {
                    if (!isNum(arr[m])) { have = false; continue; }
                    const ImVec2 pt(xOfF(fs[m]), yOf(arr[m]));
                    if (have) dl->AddLine(prev, pt, col, thick);
                    if (dots) dl->AddCircleFilled(pt, 2.0f, col);
                    prev = pt; have = true;
                }
            };
            // Pulse first (faint, thin), VNA on top (bright, with dots).
            if (havePulse) { plotF(gPulseFreqs, gPulseS21db, colOutF, 1.5f, false);
                             plotF(gPulseFreqs, gPulseS11db, colReflF, 1.5f, false); }
            if (haveCW)    { plotF(gCWfreqs, gS21db, colOut, 2.0f, true);
                             plotF(gCWfreqs, gS11db, colRefl, 2.0f, true); }

            // Marker at the frequency currently being measured (VNA).
            if (gCWsweep && gCWindex < int(gCWfreqs.size())) {
                const float x = xOfF(gCWfreqs[gCWindex]);
                dl->AddLine(ImVec2(x, p0.y), ImVec2(x, p1.y), IM_COL32(230,230,120,180), 1.5f);
            }

            // ---- Automatic resonance markers: S21 peak + S11 dip -------------
            // Use the VNA data if present (final), else the pulse snapshot.
            const std::vector<double>& mf  = haveCW ? gCWfreqs   : gPulseFreqs;
            const std::vector<float>&  m21 = haveCW ? gS21db     : gPulseS21db;
            const std::vector<float>&  m11 = haveCW ? gS11db     : gPulseS11db;
            auto extreme = [&](const std::vector<float>& arr, bool wantMax,
                               double& fOut, float& dbOut) {
                bool found = false; fOut = 0; dbOut = 0;
                const int n = int(std::min(mf.size(), arr.size()));
                for (int m = 0; m < n; ++m) {
                    if (!isNum(arr[m])) continue;
                    if (!found || (wantMax ? arr[m] > dbOut : arr[m] < dbOut)) {
                        dbOut = arr[m]; fOut = mf[m]; found = true;
                    }
                }
                return found;
            };
            double fPk, fDip; float dbPk, dbDip;
            const bool hasPk  = extreme(m21, true,  fPk,  dbPk);   // S21 transmission peak
            const bool hasDip = extreme(m11, false, fDip, dbDip);  // S11 reflection dip
            auto marker = [&](double f, float db, ImU32 col, const char* tag, float ty) {
                const float x = xOfF(f);
                dl->AddLine(ImVec2(x, p0.y), ImVec2(x, p1.y), col, 1.0f);
                dl->AddCircleFilled(ImVec2(x, yOf(db)), 3.5f, col);
                char t[48]; std::snprintf(t, sizeof(t), "%s %.3f GHz", tag, f*1e-9);
                float tx = x + 4; if (tx > p1.x - 90) tx = x - 92;
                dl->AddText(ImVec2(tx, p0.y + ty), col, t);
            };
            if (hasPk)  marker(fPk,  dbPk,  colOut,  "S21^", 4.0f);
            if (hasDip) marker(fDip, dbDip, colRefl, "S11v", 20.0f);

            // Reserve room below the plot for the frequency numbers + axis label.
            ImGui::Dummy(ImVec2(pw, ph + 36.0f));

            // Readout of the auto-detected resonance features.
            if (hasDip)
                ImGui::TextColored(ImVec4(0.94f,0.59f,0.31f,1),
                                   "S11 dip: %.3f GHz (%.1f dB)", fDip*1e-9, dbDip);
            if (hasPk)
                ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                                   "S21 peak: %.3f GHz (%.1f dB)", fPk*1e-9, dbPk);

            if (havePulse && haveCW)
                ImGui::TextDisabled("bright+dots = VNA (CW), faint = FDTD (pulse)");
            else if (haveCW)
                ImGui::TextDisabled(gCWsweep ? "VNA sweep filling left->right..."
                                             : "VNA sweep complete. Run FDTD to overlay the pulse result.");
            else
                ImGui::TextDisabled("FDTD pulse. Run CW sweep to overlay the VNA result.");
            if (gHaveCal && gApplyCal)
                ImGui::TextColored(ImVec4(0.35f,0.86f,0.47f,1),
                                   "Thru-calibrated (a matched thru = 0 dB; passive <= 0 dB).");
            else
                ImGui::TextDisabled("Uncalibrated: |S| level is arbitrary and may exceed 0 dB. "
                                    "Click 'Calibrate reference (thru)'.");
            ImGui::End();
            }
        }

        // Rebuild the model immediately (before any panel reads it) so the
        // whole frame is consistent with the current UI parameters.
        if (rebuild) {
            const ModeType  mt = asMode();
            const FieldKind fk = asField();
            const bool cav = isCavity();

            auto buildRect = [&](double fGHz) {
                // Valid rectangular index pairs:
                //   TE - either index may be 0, but not both. TE00 would give
                //        kc = 0, and every transverse component carries a
                //        factor m or n, so the whole field vanishes.
                //   TM - both must be >= 1: Ez ~ sin(m pi x/a) sin(n pi y/b) is
                //        identically zero if either index is 0.
                // Only N was guarded before, which happened to be safe only
                // because the M slider could not reach 0 at all.
                int effM = uiM, effN = uiN;
                if (mt == ModeType::TM) {
                    effM = std::max(1, effM);
                    effN = std::max(1, effN);
                } else if (effM == 0 && effN == 0) {
                    effM = 1;                    // TE00 -> dominant TE10
                }
                rectModel = TEmnModel(uiWidthMM, uiHeightMM, fGHz * 1e9,
                                      uiEpsR, uiMuR, effM, effN, mt, fk,
                                      uiDepthMM, cav, uiL, double(uiPowerW));
            };
            auto buildCyl = [&](double fGHz) {
                cylModel = CylindricalModel(uiRadiusMM, uiLengthMM / 1000.0,
                                            fGHz * 1e9, uiEpsR, uiMuR,
                                            uiN, uiM, mt, fk, cav, uiL);
            };

            if (uiGeometry == 0) buildRect(double(uiFreqGHz));
            else                 buildCyl(double(uiFreqGHz));

            // Snap the drive frequency to the mode's characteristic frequency
            // whenever the geometry/mode/dimensions changed (relock): the
            // resonant frequency f_mnl for a closed cavity, or the cutoff
            // frequency f_c for an open waveguide. A manual frequency edit sets
            // rebuild WITHOUT relock, so a detuned value is preserved.
            if (relock) {
                double fHz;
                if (cav) {
                    fHz = active()->resonantFrequency();
                } else {
                    const double c0 = 299792458.0;
                    const double kPi = 3.14159265358979323846;
                    const double kc = active()->cutoffWavenumber();
                    const double vph = c0 / std::sqrt(active()->epsilonRel() *
                                                      active()->muRel());
                    // Snap to 0.2 GHz ABOVE cutoff so the mode is actually
                    // propagating (at exactly f_c, beta = 0 and the field is
                    // uniform along z).
                    fHz = kc * vph / (2.0 * kPi) + 0.2e9;
                }
                uiFreqGHz = std::clamp(float(fHz * 1e-9), 1.0f, 60.0f);
                if (uiGeometry == 0) buildRect(double(uiFreqGHz));
                else                 buildCyl(double(uiFreqGHz));
            }

            const Bounds nb = activeBounds();
            renderer.updateBounds(nb);
            if (uiGeometry == 1)
                renderer.updateCylinder(nb.width * 0.5f, nb.depth);
            sphereRadius = (nb.depth / float(gGridNz)) * 0.9f;
            const float md = std::max({nb.width, nb.height, nb.depth});
            gCamera.setDistance(md * 2.2f);
        }

        // ---- Color scale window ----
        if (showColorBar && !devDipole) {
            int fbW2, fbH2;
            glfwGetFramebufferSize(window, &fbW2, &fbH2);
            const float barW = 40.0f, barH = 280.0f;
            const float winW = 150.0f, winH = barH + 80.0f;
            ImGui::SetNextWindowPos(ImVec2(float(fbW2) - winW - 10.0f, 10.0f),
                                    ImGuiCond_Always);
            ImGui::SetNextWindowSize(ImVec2(winW, winH), ImGuiCond_Always);
            ImGui::Begin("Color scale", nullptr,
                         ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                         ImGuiWindowFlags_NoCollapse);

            // Only claim SI units when the source actually pinned its amplitude
            // to something physical. An unnormalized eigenmode is "u.a."
            // (unidades arbitrarias) -- printing A/m there is how our numbers
            // ended up ~21x off an HFSS plot that WAS power-normalized.
            const bool physUnits = active()->physicalUnits();
            const char* unit = physUnits ? ((uiField == 0) ? "V/m" : "A/m") : "u.a.";
            const char* sym  = (uiField == 0) ? "|E|" : "|H|";
            ImGui::Text("%s (%s)", sym, unit);

            // Draw a vertical gradient bar using the fire colormap stops.
            const ImVec2 p0 = ImGui::GetCursorScreenPos();
            const ImVec2 p1 = ImVec2(p0.x + barW, p0.y + barH);
            ImDrawList* dl = ImGui::GetWindowDrawList();

            const int steps = 64;
            auto fireAt = [&](float t, float& r, float& g, float& b) {
                fireColor(t, r, g, b);
            };
            for (int s = 0; s < steps; ++s) {
                const float t0 = 1.0f - float(s)     / float(steps);
                const float t1 = 1.0f - float(s + 1) / float(steps);
                float r0,g0,b0, r1,g1,b1;
                fireAt(t0, r0,g0,b0);
                fireAt(t1, r1,g1,b1);
                const ImU32 c0 = ImGui::ColorConvertFloat4ToU32(ImVec4(r0,g0,b0,1.0f));
                const ImU32 c1 = ImGui::ColorConvertFloat4ToU32(ImVec4(r1,g1,b1,1.0f));
                const float y0 = p0.y + barH * float(s)     / float(steps);
                const float y1 = p0.y + barH * float(s + 1) / float(steps);
                dl->AddRectFilledMultiColor(ImVec2(p0.x, y0), ImVec2(p1.x, y1),
                                            c0, c0, c1, c1);
            }
            dl->AddRect(p0, p1, IM_COL32(200,200,200,255));

            // Tick labels: peak, half, 0
            const double peak = activePeak();
            const float labelX = p1.x + 6.0f;
            auto fmtVal = [&](double v) {
                char buf[32];
                std::snprintf(buf, sizeof(buf), "%.3g %s", v, unit);
                return std::string(buf);
            };
            const ImU32 tickCol = IM_COL32(0, 0, 0, 255);   // black: readable on the light panel
            dl->AddText(ImVec2(labelX, p0.y - 4.0f),          tickCol,
                        fmtVal(peak).c_str());
            dl->AddText(ImVec2(labelX, p0.y + barH * 0.5f - 7.0f), tickCol,
                        fmtVal(peak * 0.5).c_str());
            dl->AddText(ImVec2(labelX, p1.y - 12.0f),         tickCol,
                        fmtVal(0.0).c_str());

            ImGui::Dummy(ImVec2(barW, barH + 10.0f));
            ImGui::End();
        }

        // Per-plane slice state, filled by drawSection and read by the 3D
        // slice-plane overlay below. secSlicePosC is the centered perpendicular
        // coordinate (z for XY, y for ZX, x for ZY); secShown marks open planes.
        double secSlicePosC[3] = {0.0, 0.0, 0.0};
        bool   secShown[3]     = {false, false, false};

        // ---- Cross-section vector windows (XY / ZX / ZY) ----
        // Available in every domain, but only when that domain has a live field.
        if (!devDipole && hasFieldForDomain()) {
            const Bounds tb = activeBounds();
            // Plane id: 0 = XY (vary x,y at z=L/2), 1 = ZX (vary z,x at y=mid),
            //           2 = ZY (vary z,y at x=mid)
            auto drawSection = [&](int plane, bool* pOpen, const char* title,
                                   float defPosX, float defPosY)
            {
                if (!*pOpen) return;

                // Geometry, slice placement and the magnitude reference now come
                // from FieldViz; what is left here is UI and drawing. Cached per
                // plane because the antinode scan and the reference sweep are the
                // expensive parts and only change with the mode.
                const int NU = (plane == 0) ? 18 : 44;   // ZX/ZY run along z
                const int NV = (plane == 0) ? 18 : 14;

                static SectionSpan spanC[3];
                static double ref[3] = {0,0,0};
                static int refKey[3][8] = {{-1,-1,-1,-1,-1,-1,-1,-1},
                                           {-1,-1,-1,-1,-1,-1,-1,-1},
                                           {-1,-1,-1,-1,-1,-1,-1,-1}};
                const float manualFrac = (plane == 0) ? gSliceZ
                                       : (plane == 1) ? gSliceY : gSliceX;
                const int sliceKey = gManualSlice ? int(manualFrac * 400.0f + 0.5f) : -1;
                const int key[8] = {uiGeometry, uiModeType, uiField, uiM, uiN,
                                    uiL, uiStructure, sliceKey};
                bool stale = ref[plane] <= 0.0;
                for (int kk = 0; kk < 8; ++kk)
                    if (refKey[plane][kk] != key[kk]) { stale = true; break; }
                if (stale) {
                    spanC[plane] = sectionSpan(*active(), plane,
                                               gManualSlice ? double(manualFrac) : -1.0);
                    ref[plane] = sectionReference(*active(), plane, spanC[plane], NU, NV);
                    for (int kk = 0; kk < 8; ++kk) refKey[plane][kk] = key[kk];
                }
                const SectionSpan span = spanC[plane];
                const double uMin = span.uMin, uMax = span.uMax;
                const double vMin = span.vMin, vMax = span.vMax;
                const double wMin = span.wMin, wMax = span.wMax;
                const double sliceNow = span.slice;
                const double maxMag = ref[plane];
                const double Uspan = uMax - uMin;
                const double Vspan = vMax - vMin;

                const float plotW = 320.0f;
                const float plotH = plotW * float(Vspan / Uspan);
                const float winW2 = plotW + 20.0f;
                const float winH2 = plotH + 60.0f;
                ImGui::SetNextWindowPos(ImVec2(defPosX, defPosY), ImGuiCond_FirstUseEver);
                ImGui::SetNextWindowSize(ImVec2(winW2, winH2), ImGuiCond_FirstUseEver);
                if (!ImGui::Begin(title, pOpen,
                                  ImGuiWindowFlags_NoCollapse)) {
                    ImGui::End();
                    return;
                }

                const double rad = cylModel.radius();

                // Two references, deliberately different:
                //   maxMag   - max of the IN-PLANE pair over this slice. Drives
                //              arrow LENGTH and streamline flow speed, i.e. how
                //              readable the drawing is. Purely cosmetic.
                //   colorRef - peakField() over the whole volume, all three
                //              components. This is the exact denominator the 3D
                //              cloud uses (uInvPeak in cloud.vert) and the one the
                //              "Color scale" legend labels in V/m or A/m, so a hue
                //              here means the same |field| it means in 3D.
                // Colouring by maxMag made the cut disagree with the 3D view by a
                // large factor whenever the out-of-plane component matters. TE11
                // magnetic is the worst case: Hz ~ cos(pi x/a)cos(pi y/b) is the
                // LARGEST component near cutoff, peaks in the corners where the
                // transverse pair is weakest, and is invisible in the (Hx,Hy) cut.
                const double sectionPeak = activePeak();
                const double colorRef = (sectionPeak > 1e-30) ? sectionPeak
                                                              : std::max(maxMag, 1e-30);
                // Publish the (centered) slice position so the 3D view can draw
                // the matching plane.
                secSlicePosC[plane] = sliceNow - 0.5 * (wMin + wMax);
                secShown[plane]     = true;

                // Title showing the actual slice position and mode. The arrows are
                // the in-plane pair; the colour is the full |field|, hence both
                // symbols in the label.
                const char* tsym = (uiField == 0) ? "E" : "H";
                const double frac = (wMax > wMin) ? (sliceNow - wMin)/(wMax - wMin) : 0.5;
                if (plane == 0)      ImGui::Text("%s_t  |  z-slice %.2f d  |  cor = |%s|", tsym, frac, tsym);
                else if (plane == 1) ImGui::Text("(%sz,%sx)  |  y-slice %.2f  |  cor = |%s|", tsym, tsym, frac, tsym);
                else                 ImGui::Text("(%sz,%sy)  |  x-slice %.2f  |  cor = |%s|", tsym, tsym, frac, tsym);

                const ImVec2 avail = ImGui::GetContentRegionAvail();
                const float W2 = std::max(60.0f, avail.x);
                const float H2 = std::max(60.0f, W2 * float(Vspan / Uspan));
                const ImVec2 p0 = ImGui::GetCursorScreenPos();
                const ImVec2 p1 = ImVec2(p0.x + W2, p0.y + H2);
                ImDrawList* dl2 = ImGui::GetWindowDrawList();
                dl2->AddRectFilled(p0, p1, IM_COL32(15, 15, 20, 255));

                // Outline
                if (plane == 0 && uiGeometry == 1) {
                    const ImVec2 c((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);
                    dl2->AddCircle(c, std::min(W2, H2) * 0.5f - 1.0f,
                                   IM_COL32(180,180,180,255), 64);
                } else {
                    dl2->AddRect(p0, p1, IM_COL32(180,180,180,255));
                }

                // Even after auto-placing the slice at an antinode, a component
                // can be identically zero across the plane (e.g. TE with L=0, or
                // the in-plane E of TE10 in the ZX cut). There is genuinely no
                // field to show, so draw the empty frame rather than dividing by
                // a ~zero reference (which would make the arrows explode).
                if (sectionPeak > 1e-12 && maxMag <= 1e-4 * sectionPeak) {
                    ImGui::Dummy(ImVec2(W2, H2));
                    ImGui::End();
                    return;
                }

                // ---- Field lines: even-spaced streamlines with animated flow ----
                // Geometry is traced at a FIXED reference phase (refPhase) so the
                // lines stay put instead of flickering; dashes then scroll along
                // them in the field direction, faster where the field is stronger,
                // to read as a "current". Two occupancy grids decouple the seed
                // spacing (coarse, d_sep) from the stop distance (fine, d_test <
                // d_sep) so lines run long and complete rather than ending in stubs.
                if (gFieldLines) {
                    // Geometry comes from FieldViz::traceSection -- see its
                    // header for why the lines are independent rather than
                    // greedily even-spaced. Everything below is presentation:
                    // mapping physical coordinates into the plot rect, the
                    // scrolling dashes, and the arrowheads.
                    auto toScreen = [&](float u, float v) {
                        return ImVec2(p0.x + float((u - uMin) / Uspan * W2),
                                      p0.y + float((1.0 - (v - vMin) / Vspan) * H2));
                    };
                    auto fireCol = [&](double t, float alpha) -> ImU32 {
                        float r, g, b; fireColor(float(t), r, g, b);
                        return ImGui::ColorConvertFloat4ToU32(ImVec4(r, g, b, alpha));
                    };

                    // Seed lattice in FIXED screen spacing, so the line count
                    // does not jump when the window is resized.
                    const float seedPx =
                        std::clamp(46.0f / std::max(0.25f, gSecLineDensity), 16.0f, 160.0f);
                    const int NSU = std::max(2, int(W2 / seedPx));
                    const int NSV = std::max(2, int(H2 / seedPx));
                    const std::vector<Streamline> lines =
                        traceSection(*active(), plane, span, phase,
                                     NSU, NSV, maxMag, colorRef);

                    // Flow scroll offset; dash period/floor in flow units.
                    const double tOffset = gFlowPhase * double(gFlowSpeed) * 2.5;
                    const float period = 30.0f, duty = 0.5f, Mfloor = 0.12f;
                    // Arrowheads every headEvery px of screen arclength, sized by
                    // the local in-plane magnitude. The dashes only imply motion;
                    // the heads state the SIGN, which matters because the field
                    // reverses every half cycle and a dashed line alone looks
                    // identical either way.
                    const float headEvery = 42.0f, headMax = 8.0f;

                    for (const Streamline& sl : lines) {
                        double flow = 0.0;
                        float sinceHead = headEvery * 0.5f;
                        for (size_t i = 0; i + 1 < sl.u.size(); ++i) {
                            const ImVec2 A = toScreen(sl.u[i],   sl.v[i]);
                            const ImVec2 B = toScreen(sl.u[i+1], sl.v[i+1]);
                            const float dx = B.x - A.x, dy = B.y - A.y;
                            const float segLen = std::sqrt(dx*dx + dy*dy);
                            const float M = std::max(Mfloor,
                                0.5f * (sl.inPlane[i] + sl.inPlane[i+1]));
                            const float fmid = float(flow) + 0.5f * segLen / M;
                            float q = fmid / period - float(tOffset);
                            q -= std::floor(q);
                            const bool on = q < duty;
                            const ImU32 col = on ? fireCol(sl.total[i], 1.0f)
                                                 : fireCol(sl.total[i], 0.22f);
                            dl2->AddLine(A, B, col, on ? 1.9f : 1.0f);
                            flow += segLen / M;

                            sinceHead += segLen;
                            if (sinceHead >= headEvery && segLen > 1e-3f) {
                                sinceHead = 0.0f;
                                const float ux = dx/segLen, uy = dy/segLen;
                                const float nx = -uy, ny = ux;
                                const float hl = headMax *
                                    (0.30f + 0.70f *
                                     std::sqrt(std::clamp(sl.inPlane[i], 0.0f, 1.0f)));
                                const ImU32 hc = fireCol(sl.total[i], 1.0f);
                                dl2->AddLine(B, ImVec2(B.x - ux*hl + nx*hl*0.55f,
                                                       B.y - uy*hl + ny*hl*0.55f), hc, 1.7f);
                                dl2->AddLine(B, ImVec2(B.x - ux*hl - nx*hl*0.55f,
                                                       B.y - uy*hl - ny*hl*0.55f), hc, 1.7f);
                            }
                        }
                    }

                    ImGui::Dummy(ImVec2(W2, H2));
                    ImGui::End();
                    return;
                }

                const float cellW = W2 / float(NU);
                const float cellH = H2 / float(NV);
                const float arrowMax = 0.9f * std::min(cellW, cellH);


                // Arrow list from FieldViz: positions and directions physical,
                // the two scalars already normalized (inPlane -> length,
                // total -> colour on the 3D/legend scale). Only the mapping to
                // the plot rect and the ImGui calls remain here.
                const std::vector<SectionArrow> arrows =
                    sampleSection(*active(), plane, span, phase, NU, NV,
                                  maxMag, colorRef);
                for (const SectionArrow& ar : arrows) {
                    const float cx = p0.x + float((ar.u - uMin) / Uspan * W2);
                    const float cy = p0.y + float((1.0 - (ar.v - vMin) / Vspan) * H2);
                    // sqrt keeps small arrows visible without letting big ones
                    // overrun their cell.
                    const float alen = std::sqrt(ar.inPlane) * arrowMax;
                    const float dx =  ar.du * alen;
                    const float dy = -ar.dv * alen;   // screen y grows downward
                    const float ax = cx - dx * 0.5f, ay = cy - dy * 0.5f;
                    const float bx = cx + dx * 0.5f, by = cy + dy * 0.5f;

                    float rr, gg, bb;
                    fireColor(ar.total, rr, gg, bb);
                    const ImU32 col = ImGui::ColorConvertFloat4ToU32(
                        ImVec4(rr, gg, bb, 1.0f));
                    dl2->AddLine(ImVec2(ax, ay), ImVec2(bx, by), col, 1.5f);

                    const float hlen = 3.0f;
                    const float len = std::sqrt(dx*dx + dy*dy);
                    if (len > 1e-3f) {
                        const float ux = dx / len, uy = dy / len;
                        const float px_ = -uy, py_ = ux;
                        dl2->AddLine(ImVec2(bx, by),
                            ImVec2(bx - ux*hlen + px_*hlen*0.5f,
                                   by - uy*hlen + py_*hlen*0.5f), col, 1.5f);
                        dl2->AddLine(ImVec2(bx, by),
                            ImVec2(bx - ux*hlen - px_*hlen*0.5f,
                                   by - uy*hlen - py_*hlen*0.5f), col, 1.5f);
                    }
                }

                ImGui::Dummy(ImVec2(W2, H2));
                ImGui::End();
            };

            int fbW3, fbH3;
            glfwGetFramebufferSize(window, &fbW3, &fbH3);
            const float rightX = float(fbW3) - 360.0f;
            drawSection(0, &showXY, "XY cross section", rightX, 380.0f);
            drawSection(1, &showZX, "ZX cross section", 20.0f, float(fbH3) - 220.0f);
            drawSection(2, &showZY, "ZY cross section", 20.0f, float(fbH3) - 120.0f);
        }
        // ---- Cutoff frequency & spectrum window ----
        if (showSpectrum && !devDipole) {
            const double c0 = 299792458.0;
            const double kc = active()->cutoffWavenumber();
            const double er = active()->epsilonRel();
            const double mr = active()->muRel();
            const double vph = c0 / std::sqrt(er * mr);
            const double fc  = kc * vph / (2.0 * 3.14159265358979323846);

            ImGui::SetNextWindowPos(ImVec2(10.0f, 360.0f), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(360.0f, 240.0f), ImGuiCond_FirstUseEver);
            ImGui::Begin("Spectrum / cutoff", nullptr, ImGuiWindowFlags_NoCollapse);

            auto fmtHz = [](double f) {
                char buf[32];
                if (f >= 1e9)      std::snprintf(buf, sizeof(buf), "%.3f GHz", f * 1e-9);
                else if (f >= 1e6) std::snprintf(buf, sizeof(buf), "%.3f MHz", f * 1e-6);
                else               std::snprintf(buf, sizeof(buf), "%.3f Hz",  f);
                return std::string(buf);
            };
            ImGui::Text("f_c = %s", fmtHz(fc).c_str());
            const double fNow = double(uiFreqGHz) * 1e9;
            ImGui::Text("f   = %s  (%s cutoff)", fmtHz(fNow).c_str(),
                        (fNow >= fc ? "above" : "below"));
            if (isCavity()) {
                const double fres = active()->resonantFrequency();
                ImGui::TextColored(ImVec4(0.6f, 0.9f, 1.0f, 1.0f),
                                   "f_res (%s%d%d%d) = %s",
                                   (uiModeType == 0 ? "TE" : "TM"),
                                   uiM, uiN, uiL, fmtHz(fres).c_str());
                const double detune = fNow - fres;
                ImGui::Text("detuning = %+.3f GHz", detune * 1e-9);
            }

            // Sample β(f) for f in [0, 2 f_c]. Below cutoff β is imaginary
            // (mode is evanescent) so we plot 0 there and plot α on a
            // secondary curve. Here we just show |β|/kc and α/kc side by side.
            const int N = 200;
            static float betaArr[200], alphaArr[200];
            const double fMax = (fc > 0.0) ? 2.0 * fc : 1.0;
            for (int i = 0; i < N; ++i) {
                const double f = fMax * double(i) / double(N - 1);
                const double kk = 2.0 * 3.14159265358979323846 * f / vph;
                const double d  = kk * kk - kc * kc;
                if (d >= 0.0) {
                    betaArr[i]  = float(std::sqrt(d) / (kc > 0.0 ? kc : 1.0));
                    alphaArr[i] = 0.0f;
                } else {
                    betaArr[i]  = 0.0f;
                    alphaArr[i] = float(std::sqrt(-d) / (kc > 0.0 ? kc : 1.0));
                }
            }
            char ov1[64]; std::snprintf(ov1, sizeof(ov1), "β/k_c (propagating)");
            char ov2[64]; std::snprintf(ov2, sizeof(ov2), "α/k_c (evanescent)");
            ImGui::PlotLines("##beta", betaArr, N, 0, ov1, 0.0f, 1.5f,
                             ImVec2(0, 70));
            ImGui::PlotLines("##alpha", alphaArr, N, 0, ov2, 0.0f, 1.5f,
                             ImVec2(0, 70));

            // Cursor at current f
            const float fracNow = float(std::min(1.0, fNow / fMax));
            ImGui::Text("f/f_max = %.2f   (f_max = %s)", fracNow,
                        fmtHz(fMax).c_str());
            ImGui::End();
        }

        // (Model rebuild was moved up to just after the controls window so the
        // color scale, cross-section and spectrum panels all use the freshly
        // rebuilt model this same frame — otherwise the section normalization
        // lags one model version behind and the arrows misbehave.)

        // Design mode = builder open with no numerical field yet: only the CSG
        // wireframe preview shows. Once solved (gUseNumerical) or FDTD, the field
        // is shown. The cloud is drawn from the baked point-sprite buffer.
        const bool designMode = gBuilderOn && !gUseNumerical && !gUseFdtd && !gUseMicro;
        const bool showCloud = !designMode && gView3D != 1 && !devDipole; // dev plot takes over

        // Keep the source's field peak fresh (cheap coarse scan, cached) so the
        // sections / streamlines / colorbar normalize correctly.
        activeSample(4, 4, 4, gCutawayOn, 0.05f, phase);

        // Re-bake the cloud only when something relevant changes (mode, grid,
        // cutaway, view, source). FDTD is time-domain so it re-bakes every frame.
        {
            static int pSamp=-1,pView=-1,pMode=-1;
            static bool pCut=false,pNum=false,pFdtd=false,pBuild=false; static FieldSource* pSrc=nullptr;
            FieldSource* s = active();
            if (rebuild || pSamp!=gCloudSamples||pCut!=gCutawayOn||
                pView!=gView3D||pMode!=gNumMode||pNum!=gUseNumerical||pFdtd!=gUseFdtd||
                pBuild!=gBuilderOn||pSrc!=s)
                gCloudDirty = true;
            pSamp=gCloudSamples;pCut=gCutawayOn;pView=gView3D;pMode=gNumMode;
            pNum=gUseNumerical;pFdtd=gUseFdtd;pBuild=gBuilderOn;pSrc=s;
        }
        const bool domainHasField = hasFieldForDomain();

        if (showCloud && domainHasField &&
            // During a microstrip sweep the cloud is refreshed on a throttle (via
            // gCloudDirty) rather than every frame, so a whole-band sweep is not
            // paying for a Monte-Carlo re-bake on each step chunk.
            (gCloudDirty || (gUseFdtd && gFdtdStepping) || (gUseMicro && gMicroStepping && !gMsSweep))) {
            // Time-domain sources (FDTD, microstrip) hold a real instantaneous
            // field, so one sample suffices; analytic phasor models need several
            // phases to bake the standing-wave envelope.
            bakeCloud((gUseFdtd || gUseMicro) ? 1 : 8);
            gCloudDirty = false;
        }

        int fbW, fbH;
        glfwGetFramebufferSize(window, &fbW, &fbH);
        const float aspect = (fbH > 0) ? float(fbW) / float(fbH) : 1.0f;

        glClearColor(0.90f, 0.91f, 0.93f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const glm::mat4 view = gCamera.viewMatrix();
        const glm::mat4 proj = gCamera.projectionMatrix(aspect);

        const bool analyticShape = gDomain == SimDomain::Waveguide &&
                                   !gUseNumerical && !gUseFdtd && !gUseMicro && !devDipole;
        renderer.draw(view, proj, sphereRadius, showFloor && !devDipole,
                      uiGeometry == 0 && analyticShape,
                      uiGeometry == 1 && analyticShape);

        // ---- Microstrip structure as SOLID 3D geometry (real depth/occlusion) ----
        // Built from the material primitives into triangle meshes with normals, so
        // the copper genuinely occludes the field cloud and the port spheres.
        // Rebuilt when the sim changes, or every frame while the field-skin colour
        // (|Js| on the copper) is live.
        if (gUseMicro && gMicro) {
            static const MicrostripSim* sPrevMs = nullptr;
            static bool sPrevSkin = false;
            if (sPrevMs != gMicro.get() || sPrevSkin != gCopperSkin || gCopperSkin) {
                sPrevMs = gMicro.get(); sPrevSkin = gCopperSkin;
                const Bounds sb = gMicro->bounds();
                const glm::vec3 sc(sb.width*0.5f, sb.height*0.5f, sb.depth*0.5f);
                const double kPiL = 3.14159265358979323846;
                std::vector<float> op, tr; op.reserve(8192); tr.reserve(2048);

                // Peak |Js| over the copper, for the field-skin normalisation. The
                // scan must be fine enough to actually catch the peak, otherwise the
                // colours saturate (everything red) or wash out.
                double jsMax = 1e-30;
                if (gCopperSkin) for (const MicrostripSim::Prim& p : gMicro->prims()) {
                    if (p.mat != MicrostripSim::Pec) continue;
                    if (!((p.kind==0 && p.ymin>1e-9) || p.kind==2)) continue;
                    const double ys = p.ymin - 1.5*(p.ymax - p.ymin);
                    if (p.kind==0) { for (int ix=0;ix<=64;++ix) for (int iz=0;iz<=16;++iz)
                        jsMax = std::max(jsMax, gMicro->surfaceCurrent(
                            p.xmin+(p.xmax-p.xmin)*ix/64.0, ys, p.zmin+(p.zmax-p.zmin)*iz/16.0)); }
                    else { const double R=p.radius; for (int ix=0;ix<=64;++ix) for (int iz=0;iz<=64;++iz) {
                        const double x=p.cx-R+2*R*ix/64.0, z=p.cz-R+2*R*iz/64.0;
                        const double d2=(x-p.cx)*(x-p.cx)+(z-p.cz)*(z-p.cz);
                        if (d2>R*R || d2<p.rinner*p.rinner) continue;
                        jsMax = std::max(jsMax, gMicro->surfaceCurrent(x,ys,z)); } }
                }
                double skinYs = 0.0;   // y at which |Js| is sampled for the current prim

                auto fireRGB = [](double t, float& r, float& g, float& b){
                    // `!(t >= 0)` also catches NaN: int(NaN) below would index the
                    // stop table out of bounds and take the whole app down.
                    if (!(t >= 0.0)) t = 0.0;
                    if (t > 1.0) t = 1.0;
                    fireColor(float(t), r, g, b);
                };
                // |Js| colour sampled just BELOW the copper (tangential E/H is forced
                // to zero ON a PEC face, so sampling on it would always read ~0).
                auto skinColAt = [&](double x, double z, const glm::vec3& base) -> glm::vec3 {
                    if (!gCopperSkin || jsMax <= 1e-30) return base;
                    float r,g,b; fireRGB(gMicro->surfaceCurrent(x, skinYs, z)/jsMax, r,g,b);
                    return glm::vec3(r,g,b);
                };
                auto emitTriC = [&](std::vector<float>& out, glm::vec3 a, glm::vec3 b, glm::vec3 c,
                                    glm::vec3 n, glm::vec3 col){
                    const glm::vec3 v[3] = {a,b,c};
                    for (int i=0;i<3;++i){
                        out.push_back(v[i].x - sc.x); out.push_back(v[i].y - sc.y); out.push_back(v[i].z - sc.z);
                        out.push_back(n.x); out.push_back(n.y); out.push_back(n.z);
                        out.push_back(col.r); out.push_back(col.g); out.push_back(col.b);
                    }
                };
                // Quad subdivided NU x NV; with skin on, each sub-cell gets one flat
                // colour sampled at its centre (same look the 2D overlay had).
                auto emitQuad = [&](std::vector<float>& out, glm::vec3 p00, glm::vec3 p10,
                                    glm::vec3 p11, glm::vec3 p01, glm::vec3 n,
                                    glm::vec3 base, bool skin, int NU, int NV){
                    auto L=[](const glm::vec3& a, const glm::vec3& b, float t){ return a + (b-a)*t; };
                    auto P=[&](float u,float v){ return L(L(p00,p10,u), L(p01,p11,u), v); };
                    for (int i=0;i<NU;++i) for (int j=0;j<NV;++j){
                        const float u0=float(i)/NU, u1=float(i+1)/NU, v0=float(j)/NV, v1=float(j+1)/NV;
                        const glm::vec3 a=P(u0,v0), b=P(u1,v0), c=P(u1,v1), d=P(u0,v1);
                        const glm::vec3 m=P(0.5f*(u0+u1), 0.5f*(v0+v1));
                        const glm::vec3 col = skin ? skinColAt(m.x, m.z, base) : base;
                        emitTriC(out,a,b,c,n,col); emitTriC(out,a,c,d,n,col);
                    }
                };
                auto emitBox = [&](std::vector<float>& out, double x0,double x1,double y0,double y1,
                                   double z0,double z1, glm::vec3 base, bool skin){
                    const glm::vec3 a000(x0,y0,z0), a100(x1,y0,z0), a110(x1,y1,z0), a010(x0,y1,z0);
                    const glm::vec3 a001(x0,y0,z1), a101(x1,y0,z1), a111(x1,y1,z1), a011(x0,y1,z1);
                    // Only the top face needs the fine grid (that is the visible skin);
                    // the sides of a trace are sub-cell slivers.
                    int NX = 1, NZ = 1;
                    if (skin) {
                        NX = std::clamp(int(220.0*(x1-x0)/std::max(1e-9,double(sb.width))), 4, 220);
                        NZ = std::clamp(int(90.0 *(z1-z0)/std::max(1e-9,double(sb.depth))), 2, 90);
                    }
                    emitQuad(out, a010,a110,a111,a011, glm::vec3(0,1,0),  base, skin, NX, NZ); // top
                    emitQuad(out, a000,a100,a101,a001, glm::vec3(0,-1,0), base, skin, 1, 1);
                    emitQuad(out, a001,a101,a111,a011, glm::vec3(0,0,1),  base, skin, 1, 1);
                    emitQuad(out, a100,a000,a010,a110, glm::vec3(0,0,-1), base, skin, 1, 1);
                    emitQuad(out, a101,a100,a110,a111, glm::vec3(1,0,0),  base, skin, 1, 1);
                    emitQuad(out, a000,a001,a011,a010, glm::vec3(-1,0,0), base, skin, 1, 1);
                };
                auto emitRing = [&](std::vector<float>& out, double cx,double cz,double R,double ri,
                                    double y0,double y1, glm::vec3 base, bool skin){
                    const int N  = skin ? 220 : 72;   // angular segments
                    const int NR = skin ? 8   : 1;    // radial bands across the ring width
                    for (int i=0;i<N;++i){
                        const double t0=2*kPiL*i/N, t1=2*kPiL*(i+1)/N;
                        const double c0=std::cos(t0), s0=std::sin(t0), c1=std::cos(t1), s1=std::sin(t1);
                        auto P=[&](double rad,double cs,double sn,double y){ return glm::vec3(cx+rad*cs, y, cz+rad*sn); };
                        for (int k=0;k<NR;++k){
                            const double r0=ri+(R-ri)*k/NR, r1=ri+(R-ri)*(k+1)/NR;
                            const double rm=0.5*(r0+r1), tm=0.5*(t0+t1);
                            const glm::vec3 col = skin ? skinColAt(cx+rm*std::cos(tm), cz+rm*std::sin(tm), base) : base;
                            emitTriC(out, P(r0,c0,s0,y1), P(r1,c0,s0,y1), P(r1,c1,s1,y1), glm::vec3(0,1,0), col);
                            emitTriC(out, P(r0,c0,s0,y1), P(r1,c1,s1,y1), P(r0,c1,s1,y1), glm::vec3(0,1,0), col);
                            emitTriC(out, P(r0,c0,s0,y0), P(r1,c1,s1,y0), P(r1,c0,s0,y0), glm::vec3(0,-1,0), col);
                            emitTriC(out, P(r0,c0,s0,y0), P(r0,c1,s1,y0), P(r1,c1,s1,y0), glm::vec3(0,-1,0), col);
                        }
                        const glm::vec3 nrm(float(c0),0.0f,float(s0));
                        const glm::vec3 colO = skin ? skinColAt(cx+R*c0, cz+R*s0, base) : base;
                        emitTriC(out, P(R,c0,s0,y0), P(R,c1,s1,y0), P(R,c1,s1,y1), nrm, colO);
                        emitTriC(out, P(R,c0,s0,y0), P(R,c1,s1,y1), P(R,c0,s0,y1), nrm, colO);
                        if (ri > 1e-9) {
                            const glm::vec3 colI = skin ? skinColAt(cx+ri*c0, cz+ri*s0, base) : base;
                            emitTriC(out, P(ri,c0,s0,y0), P(ri,c1,s1,y1), P(ri,c1,s1,y0), -nrm, colI);
                            emitTriC(out, P(ri,c0,s0,y0), P(ri,c0,s0,y1), P(ri,c1,s1,y1), -nrm, colI);
                        }
                    }
                };

                for (const MicrostripSim::Prim& p : gMicro->prims()) {
                    const bool copper = (p.mat==MicrostripSim::Pec && p.ymin > 1e-9);
                    const bool diel   = (p.mat==MicrostripSim::Dielectric);
                    const glm::vec3 base = copper ? glm::vec3(0.82f,0.55f,0.24f)
                                         : diel   ? glm::vec3(0.35f,0.78f,0.47f)
                                                  : glm::vec3(0.66f,0.66f,0.69f);
                    std::vector<float>& dst = diel ? tr : op;
                    skinYs = p.ymin - 1.5*(p.ymax - p.ymin);
                    if (p.kind == 0)      emitBox(dst, p.xmin,p.xmax,p.ymin,p.ymax,p.zmin,p.zmax, base, copper);
                    else if (p.kind == 2) emitRing(dst, p.cx,p.cz,p.radius,p.rinner,p.ymin,p.ymax, base, copper);
                }
                renderer.updateStructure(op, tr);
            }
            // Opaque pass first: writes depth so the cloud/spheres get occluded.
            renderer.drawStructure(view, proj, /*translucent=*/false, 1.0f);
        }

        // Baked field cloud (point sprites), animated on the GPU by the phase.
        if (showCloud && domainHasField) {
            // Point radius from the measured mean spacing of the accepted points
            // (set in bakeCloud). More samples -> tighter spacing -> smaller dots.
            const float worldR = gCloudMeanSpacing * 1.5f;
            const float pointScale = worldR * proj[1][1] * float(fbH) * 0.5f;
            renderer.drawCloud(view, proj, float(phase), gCloudInvPeak, pointScale,
                               gCloudOpaque);
        }

        // ---- Source/sense port markers: real lit 3D spheres (microstrip) ----
        // Positions are centred (physical - domain centre) to match the cloud's
        // world space; per-sphere radius (half the strip width) rides in intensity.
        if (gUseMicro && gMicro) {
            const Bounds pb = gMicro->bounds();
            const glm::vec3 pc(pb.width*0.5f, pb.height*0.5f, pb.depth*0.5f);
            auto mkPort = [&](const std::array<double,4>& mk, float cr, float cg, float cb){
                Particle p; p.x=float(mk[0])-pc.x; p.y=float(mk[1])-pc.y; p.z=float(mk[2])-pc.z;
                p.r=cr; p.g=cg; p.b=cb; p.intensity=float(mk[3]); return p;
            };
            const std::vector<Particle> ports = {
                mkPort(gMicro->sourceMarker(), 0.20f, 0.62f, 0.95f),   // source: cyan
                mkPort(gMicro->senseMarker(),  0.95f, 0.62f, 0.10f),   // sense:  amber
            };
            renderer.drawSpheres(view, proj, ports);
            // Translucent substrate last (blended, no depth write).
            renderer.drawStructure(view, proj, /*translucent=*/true, 0.30f);
        }

        // ---- Microstrip structure overlay (ground / substrate / trace) ----
        // Drawn as projected wireframe boxes so the physical schematic is visible
        // around the field: gray ground plate, green translucent substrate slab,
        // filled copper trace on top.
        if (gUseMicro && gMicro) {
            const Bounds mb = gMicro->bounds();
            const glm::vec3 ctr(mb.width*0.5f, mb.height*0.5f, mb.depth*0.5f);
            const glm::vec4 vpM(0.0f, 0.0f, float(fbW), float(fbH));
            ImDrawList* bg = ImGui::GetBackgroundDrawList();
            auto projM = [&](const glm::vec3& w, ImVec2& out) -> bool {
                glm::vec3 s = glm::project(w - ctr, view, proj, vpM);
                if (s.z < 0.0f || s.z > 1.0f) return false;
                out = ImVec2(s.x, float(fbH) - s.y); return true;
            };
            // The ground / substrate / copper are now real 3D geometry drawn in the
            // GL pass (see the solid-structure block above), so nothing is projected
            // here any more. Port markers are 3D spheres too; only their labels are
            // drawn as overlay text below.
            // Port markers are drawn as real 3D spheres in the GL pass; here we
            // just label them at the projected marker positions.
            auto portLabel = [&](const std::array<double,4>& mk, ImU32 col, const char* tag){
                ImVec2 s2;
                if (projM(glm::vec3(float(mk[0]), float(mk[1]), float(mk[2])), s2))
                    bg->AddText(ImVec2(s2.x + 8.0f, s2.y - 8.0f), col, tag);
            };
            portLabel(gMicro->sourceMarker(), IM_COL32(80,200,255,255), "source");
            portLabel(gMicro->senseMarker(),  IM_COL32(255,205,70,255), "sense");
        }

        // ---- [DEV] Rotating magnetic dipole (pulsar) ---------------------------
        // m(t) = m[ sinA cos(wt), sinA sin(wt), cosA ], w = 2*pi/P. The magnetic
        // moment precesses around the spin axis z at obliquity A; the transverse
        // (rotating) part is what radiates. Drawn as: star sphere, spin axis,
        // precessing m-vector + its cone, and the rotating dipole field lines.
        if (devDipole) {
            const double A = double(uiDipAlphaDeg) * 3.14159265358979 / 180.0;
            const double w = 2.0 * 3.14159265358979 / std::max(1e-3f, uiDipPeriod);
            const double phi = w * gDipoleT;
            const double sA = std::sin(A), cA = std::cos(A);
            const glm::vec3 mhat(float(sA*std::cos(phi)), float(sA*std::sin(phi)), float(cA));

            const glm::vec4 vpD(0.0f, 0.0f, float(fbW), float(fbH));
            ImDrawList* bg = ImGui::GetBackgroundDrawList();
            auto projD = [&](const glm::vec3& w3, ImVec2& out) -> bool {
                glm::vec3 s = glm::project(w3, view, proj, vpD);
                if (s.z < 0.0f || s.z > 1.0f) return false;
                out = ImVec2(s.x, float(fbH) - s.y); return true;
            };
            // Depth cue: nearness in [0,1] from view-space z (1 = closest to eye).
            // Far elements fade toward the background (aerial perspective) and get
            // thinner; near elements pop -> reads as 3D depth on the flat overlay.
            const float vz0 = (view * glm::vec4(0.0f,0.0f,0.0f,1.0f)).z;
            auto nearOf = [&](const glm::vec3& w)->float {
                const float vz = (view * glm::vec4(w,1.0f)).z;
                return std::clamp(0.5f + 0.5f*(vz - vz0)/1.8f, 0.0f, 1.0f);
            };
            auto shadeA = [&](ImU32 col, float k)->ImU32 {
                const int a = int(((col >> IM_COL32_A_SHIFT) & 0xFF) * std::clamp(k,0.0f,1.0f));
                return (col & ~IM_COL32_A_MASK) | (ImU32(a) << IM_COL32_A_SHIFT);
            };
            auto line3 = [&](const glm::vec3& a, const glm::vec3& b, ImU32 col, float th){
                ImVec2 pa, pb; if (projD(a,pa) && projD(b,pb)) {
                    const float n = 0.5f*(nearOf(a)+nearOf(b));
                    bg->AddLine(pa,pb, shadeA(col, 0.32f+0.68f*n), th*(0.65f+0.7f*n));
                } };

            const float R = 0.35f;      // star radius
            const float Lm = 1.15f;     // m-vector length
            // Star: a spherical wireframe grid that SPINS about z with the star
            // (the meridians rotate by the spin angle phi, revealing the rotation).
            const float cP = float(std::cos(phi)), sP = float(std::sin(phi));
            auto spun = [&](const glm::vec3& p){ return glm::vec3(cP*p.x - sP*p.y, sP*p.x + cP*p.y, p.z); };
            const ImU32 gcol = IM_COL32(70,95,150,200);
            const int NC = 44;
            // Soft shaded fill: gives the star volume, lit from the upper-left
            // (concentric disks from dark rim -> bright core, offset toward light).
            const glm::vec3 camRight(view[0][0], view[1][0], view[2][0]);
            { ImVec2 cS, rS;
              if (projD(glm::vec3(0.0f), cS) && projD(R*camRight, rS)) {
                  const float rad = std::max(4.0f, float(std::hypot(rS.x-cS.x, rS.y-cS.y)));
                  const ImVec2 lite(cS.x - rad*0.35f, cS.y - rad*0.35f);
                  const int steps = 14;
                  for (int s=steps; s>=1; --s){ const float f=float(s)/steps, t=1.0f-f;
                      auto L=[&](int lo,int hi){ return int(lo + (hi-lo)*t); };
                      const ImU32 col = IM_COL32(L(40,175), L(58,200), L(98,235), 255);
                      const ImVec2 c(cS.x+(lite.x-cS.x)*t, cS.y+(lite.y-cS.y)*t);
                      bg->AddCircleFilled(c, rad*f, col, 40);
                  }
                  bg->AddCircle(cS, rad, IM_COL32(30,42,72,255), 44, 1.5f);
              } }
            for (int lat=1; lat<6; ++lat) {                 // latitude circles
                const double th=3.14159265*lat/6.0, rr=R*std::sin(th), zz=R*std::cos(th);
                glm::vec3 prev; bool have=false;
                for (int i=0;i<=NC;++i){ const double a=2.0*3.14159265*i/NC;
                    glm::vec3 pnt=spun(glm::vec3(float(rr*std::cos(a)),float(rr*std::sin(a)),float(zz)));
                    if (have) line3(prev,pnt,gcol,1.0f); prev=pnt; have=true; }
            }
            for (int lon=0; lon<8; ++lon) {                 // meridians (show spin)
                const double ph=2.0*3.14159265*lon/8.0;
                glm::vec3 prev; bool have=false;
                for (int i=0;i<=NC;++i){ const double th=3.14159265*i/NC;
                    glm::vec3 base(float(R*std::sin(th)*std::cos(ph)),float(R*std::sin(th)*std::sin(ph)),float(R*std::cos(th)));
                    glm::vec3 pnt=spun(base);
                    if (have) line3(prev,pnt,gcol,1.0f); prev=pnt; have=true; }
            }
            // Mark the magnetic poles on the surface (where the m axis pierces it).
            { ImVec2 pN,pS;
              if (projD( R*mhat,pN)) bg->AddCircleFilled(pN,4.0f,IM_COL32(230,80,80,255));
              if (projD(-R*mhat,pS)) bg->AddCircleFilled(pS,4.0f,IM_COL32(80,120,230,255)); }

            // Spin axis z (with arrowheads) + label anchor.
            line3(glm::vec3(0,0,-1.5f), glm::vec3(0,0,1.6f), IM_COL32(90,90,100,220), 1.5f);
            line3(glm::vec3(0,0,1.6f), glm::vec3(0.05f,0,1.45f), IM_COL32(90,90,100,220), 1.5f);
            line3(glm::vec3(0,0,1.6f), glm::vec3(-0.05f,0,1.45f), IM_COL32(90,90,100,220), 1.5f);
            { ImVec2 t; if (projD(glm::vec3(0,0,1.7f),t)) bg->AddText(t, IM_COL32(60,60,70,255), "z (Omega)"); }

            // Precession cone: circle traced by the m-tip + a few cone rulings.
            { glm::vec3 prev; bool have=false;
              for (int i=0;i<=64;++i){ const double a=2.0*3.14159265*i/64;
                  glm::vec3 pnt(float(Lm*sA*std::cos(a)), float(Lm*sA*std::sin(a)), float(Lm*cA));
                  if (have) line3(prev,pnt,IM_COL32(80,160,230,150),1.2f); prev=pnt; have=true; }
              for (int i=0;i<8;++i){ const double a=2.0*3.14159265*i/8;
                  line3(glm::vec3(0,0,0), glm::vec3(float(Lm*sA*std::cos(a)),float(Lm*sA*std::sin(a)),float(Lm*cA)),
                        IM_COL32(80,160,230,60),1.0f); } }

            // Magnetic dipole field lines: they LEAVE the star near one pole,
            // loop out, and RE-ENTER near the other pole. r = L sin^2(theta)
            // (theta from the m axis), clipped to r >= R so the footpoints sit on
            // the star surface. Arrowheads show the B direction (out of N, into S).
            if (uiDipFieldLines) {
                const glm::vec3 up = (std::abs(mhat.z) < 0.9f) ? glm::vec3(0,0,1) : glm::vec3(1,0,0);
                const glm::vec3 u = glm::normalize(glm::cross(up, mhat));
                const glm::vec3 v = glm::normalize(glm::cross(mhat, u));
                const int NAZ = std::max(3, uiDipFieldAz), NTH = 48;
                const int NSH = std::max(1, uiDipShellsN);
                const ImU32 faint  = IM_COL32(210,90,90,70);
                const ImU32 bright = IM_COL32(240,110,90,235);
                const double flow  = gDipoleT * 2.2;  // dashes march N -> S
                const double Ndash = 4.0;
                for (int si=0; si<NSH; ++si) {
                    const float L = R * 1.5f * std::pow(1.42f, float(si)); // nested shells
                    if (L <= R*1.02f) continue;
                    const double th1 = std::asin(std::min(1.0, std::sqrt(double(R)/L)));
                    for (int az=0; az<NAZ; ++az) {
                        const double psi = 2.0*3.14159265*az/NAZ;
                        const glm::vec3 rad = float(std::cos(psi))*u + float(std::sin(psi))*v;
                        glm::vec3 prev; bool have=false;
                        for (int i=0;i<=NTH;++i) {
                            const double g  = double(i)/NTH;   // 0 = N footpoint, 1 = S footpoint
                            const double th = th1 + (3.14159265 - 2.0*th1)*g;
                            const double r  = L*std::sin(th)*std::sin(th);
                            const glm::vec3 pnt = float(r*std::sin(th))*rad + float(r*std::cos(th))*mhat;
                            if (have) {
                                line3(prev, pnt, faint, 1.0f);            // faint continuous line
                                double ph = g*Ndash - flow; ph -= std::floor(ph);
                                if (ph < 0.42) line3(prev, pnt, bright, 2.3f); // moving dash (flow)
                            }
                            prev=pnt; have=true;
                        }
                    }
                }
            }

            // Pulsar radiation cones: one at each magnetic pole, tip on the pole,
            // axis along +/- m, half-angle proportional to the dipole moment |m|,
            // filled with a Monte-Carlo point cloud that streams outward (the
            // radiation being expelled). Everything rotates with m(t).
            if (uiDipCones) {
                const glm::vec3 up = (std::abs(mhat.z) < 0.9f) ? glm::vec3(0,0,1) : glm::vec3(1,0,0);
                const glm::vec3 u = glm::normalize(glm::cross(up, mhat));
                const glm::vec3 v = glm::normalize(glm::cross(mhat, u));
                const double thc = std::clamp(double(uiDipMag) * 0.35, 0.05, 0.9); // half-angle ~ |m|
                const float  Lc = uiDipConeLen, tanc = float(std::tan(thc));

                // MC samples in cone-local coords: s = axial fraction, u = radial^2,
                // a = azimuth. Initialized once; s advances each frame (expulsion).
                static std::vector<float> cS, cU, cA; static int cN = 0;
                const int N = 380;
                if (cN != N) { cN = N; cS.resize(N); cU.resize(N); cA.resize(N);
                    std::mt19937 rng(1234); std::uniform_real_distribution<float> uni(0.f,1.f);
                    for (int i=0;i<N;++i){ cS[i]=uni(rng); cU[i]=uni(rng); cA[i]=uni(rng)*6.2831853f; } }
                const float adv = float(dt) * 0.35f * std::max(0.1f, uiDipSpeed);
                for (int i=0;i<N;++i){ cS[i]+=adv; if (cS[i]>1.0f) cS[i]-=1.0f; }

                for (int sgn=0; sgn<2; ++sgn) {
                    const glm::vec3 axis = (sgn ? -1.0f : 1.0f) * mhat;
                    const glm::vec3 tip  = float(R) * axis;      // magnetic pole
                    const float rimR = Lc * tanc;
                    // Cone outline: rim circle + rulings from the tip.
                    glm::vec3 prev; bool have=false;
                    for (int i=0;i<=40;++i){ const double a=2.0*3.14159265*i/40;
                        glm::vec3 pnt = tip + Lc*axis + rimR*(float(std::cos(a))*u + float(std::sin(a))*v);
                        if (have) line3(prev,pnt, IM_COL32(120,50,190,150), 1.2f); prev=pnt; have=true; }
                    for (int i=0;i<8;++i){ const double a=2.0*3.14159265*i/8;
                        glm::vec3 rim = tip + Lc*axis + rimR*(float(std::cos(a))*u + float(std::sin(a))*v);
                        line3(tip, rim, IM_COL32(120,50,190,90), 1.0f); }
                    // Monte-Carlo radiation points (saturated violet, readable on
                    // a light background; brightness fades toward the beam end).
                    for (int i=0;i<N;++i){
                        const float s = cS[i];
                        const float rr = s*rimR*std::sqrt(cU[i]);
                        const glm::vec3 p = tip + s*Lc*axis + rr*(std::cos(cA[i])*u + std::sin(cA[i])*v);
                        ImVec2 sc; if (projD(p, sc)) {
                            const float n = nearOf(p);                       // depth cue
                            const int al = std::clamp(int((150.0f*(1.0f - s) + 105.0f) * (0.55f + 0.45f*n)), 0, 255);
                            bg->AddCircleFilled(sc, 1.3f + 1.1f*n, IM_COL32(140, 40, 210, al));
                        }
                    }
                }
            }

            // The dipole moment vector m(t) (thick arrow).
            const glm::vec3 mtip = Lm * mhat;
            line3(glm::vec3(0,0,0), mtip, IM_COL32(230,60,60,255), 3.0f);
            { glm::vec3 up=(std::abs(mhat.z)<0.9f)?glm::vec3(0,0,1):glm::vec3(1,0,0);
              glm::vec3 s=glm::normalize(glm::cross(up,mhat));
              line3(mtip, mtip - 0.14f*mhat + 0.06f*s, IM_COL32(230,60,60,255), 3.0f);
              line3(mtip, mtip - 0.14f*mhat - 0.06f*s, IM_COL32(230,60,60,255), 3.0f);
              ImVec2 t; if (projD(mtip + 0.08f*mhat, t)) bg->AddText(t, IM_COL32(200,40,40,255), "m(t)"); }
        }

        // ---- Geometry builder 3D preview (CSG wireframe) ----
        // Each primitive is drawn as a wireframe, cyan for additive (union) and
        // orange for subtractive (difference), centered on the shape's AABB.
        if (gBuilderOn) {
            const Geometry g = buildGeo();
            const Aabb bb = g.bounds();
            const glm::vec3 ctr(float(0.5*(bb.xmin+bb.xmax)),
                                float(0.5*(bb.ymin+bb.ymax)),
                                float(0.5*(bb.zmin+bb.zmax)));
            const glm::vec4 vpB(0.0f, 0.0f, float(fbW), float(fbH));
            ImDrawList* bg = ImGui::GetBackgroundDrawList();
            auto projW = [&](const glm::vec3& w, ImVec2& out) -> bool {
                glm::vec3 s = glm::project(w - ctr, view, proj, vpB);
                if (s.z < 0.0f || s.z > 1.0f) return false;
                out = ImVec2(s.x, float(fbH) - s.y); return true;
            };
            auto edge = [&](const glm::vec3& A, const glm::vec3& B, ImU32 col) {
                ImVec2 a, b; if (projW(A,a) && projW(B,b)) bg->AddLine(a,b,col,1.7f);
            };
            for (const BuildStep& s : gSteps) {
                const ImU32 col = (s.op == 0) ? IM_COL32(90,200,255,235)
                                              : IM_COL32(255,130,70,235);
                if (s.type == 0) {
                    const double cx=s.p[0]/1e3, cy=s.p[1]/1e3, cz=s.p[2]/1e3;
                    const double hx=s.p[3]/2e3, hy=s.p[4]/2e3, hz=s.p[5]/2e3;
                    auto corner = [&](int i){ return glm::vec3(
                        float(cx + ((i&1)?hx:-hx)), float(cy + ((i&2)?hy:-hy)),
                        float(cz + ((i&4)?hz:-hz))); };
                    for (int i=0;i<8;++i)
                        for (int bit : {1,2,4})
                            if (!(i&bit)) edge(corner(i), corner(i|bit), col);
                } else {
                    const double px=s.p[0]/1e3, py=s.p[1]/1e3;
                    const double z0=s.p[2]/1e3, z1=s.p[3]/1e3, r=s.p[4]/1e3;
                    const int N=28;
                    glm::vec3 prev0, prev1;
                    for (int i=0;i<=N;++i){
                        const double ang=2.0*3.14159265358979*i/N;
                        const glm::vec3 c0(float(px+r*std::cos(ang)), float(py+r*std::sin(ang)), float(z0));
                        const glm::vec3 c1(float(px+r*std::cos(ang)), float(py+r*std::sin(ang)), float(z1));
                        if (i>0){ edge(prev0,c0,col); edge(prev1,c1,col); }
                        if (i % 7 == 0) edge(c0,c1,col);
                        prev0=c0; prev1=c1;
                    }
                }
            }
        }

        // ---- 3D slice planes (mirror the open cross-section slices) ----
        // Drawn as translucent projected quads/discs so you can see exactly
        // where each cross-section cuts the 3D field. In manual-slice mode the
        // sliders move these and the plots follow.
        // ---- FDTD port openings: color-coded, click an opening to cycle role ----
        if (gPortConfig && gDomain == SimDomain::Geometry) {
            if (gPorts.empty()) gPorts = detectPorts(buildGeo(), 40);
            const Aabb pbb = buildGeo().bounds();
            const float hw = float(pbb.sizeX()*0.5), hh = float(pbb.sizeY()*0.5), hd = float(pbb.sizeZ()*0.5);
            const glm::vec4 vpF(0.0f, 0.0f, float(fbW), float(fbH));
            ImDrawList* bg = ImGui::GetBackgroundDrawList();
            // Map a port's (u,v) transverse model coord to centered world.
            auto portCorner = [&](const FdtdPort& p, double u, double v) -> glm::vec3 {
                const float ax = p.side ? ( (p.axis==0)?hw:(p.axis==1)?hh:hd )
                                        : (-((p.axis==0)?hw:(p.axis==1)?hh:hd));
                if (p.axis == 0) return glm::vec3(ax, float(u - hh), float(v - hd)); // u=y,v=z
                if (p.axis == 1) return glm::vec3(float(u - hw), ax, float(v - hd)); // u=x,v=z
                return glm::vec3(float(u - hw), float(v - hh), ax);                  // u=x,v=y
            };
            const char* axName[3] = {"x","y","z"};
            for (std::size_t pi = 0; pi < gPorts.size(); ++pi) {
                const FdtdPort& p = gPorts[pi];
                const glm::vec3 c3[4] = { portCorner(p,p.uMin,p.vMin), portCorner(p,p.uMax,p.vMin),
                                         portCorner(p,p.uMax,p.vMax), portCorner(p,p.uMin,p.vMax) };
                ImVec2 s[4]; bool ok = true;
                for (int c = 0; c < 4; ++c) { glm::vec3 pr = glm::project(c3[c], view, proj, vpF);
                    if (pr.z<0.0f||pr.z>1.0f){ ok=false; break; } s[c]=ImVec2(pr.x, float(fbH)-pr.y); }
                if (!ok) continue;
                ImU32 fill, border; const char* role;
                if (p.role==1){ fill=IM_COL32(60,220,90,70); border=IM_COL32(90,240,120,235); role="IN"; }
                else if (p.role==2){ fill=IM_COL32(240,90,80,70); border=IM_COL32(255,120,100,235); role="OUT"; }
                else { fill=IM_COL32(150,150,160,40); border=IM_COL32(180,180,190,170); role="off"; }
                bg->AddConvexPolyFilled(s, 4, fill);
                bg->AddPolyline(s, 4, border, ImDrawFlags_Closed, 2.0f);
                const ImVec2 ctr((s[0].x+s[1].x+s[2].x+s[3].x)*0.25f, (s[0].y+s[1].y+s[2].y+s[3].y)*0.25f);
                char lbl[40]; std::snprintf(lbl, sizeof(lbl), "P%d %s%s: %s", int(pi)+1, axName[p.axis], p.side?"+":"-", role);
                bg->AddText(ctr, border, lbl);
            }
            if (gPortClickPending) {
                gPortClickPending = false;
                const glm::vec3 o = glm::unProject(glm::vec3(float(gPortClickX), float(fbH)-float(gPortClickY), 0.0f), view, proj, vpF);
                const glm::vec3 fp = glm::unProject(glm::vec3(float(gPortClickX), float(fbH)-float(gPortClickY), 1.0f), view, proj, vpF);
                const glm::vec3 d = fp - o;
                int best = -1; double bestT = 1e30;
                for (std::size_t pi = 0; pi < gPorts.size(); ++pi) {
                    const FdtdPort& p = gPorts[pi];
                    const int axis = p.axis;
                    const float planeC = p.side ? ((axis==0)?hw:(axis==1)?hh:hd) : (-((axis==0)?hw:(axis==1)?hh:hd));
                    if (std::abs(d[axis]) < 1e-12) continue;
                    const double t = (planeC - o[axis]) / d[axis];
                    if (t <= 0.0 || t >= bestT) continue;
                    const glm::vec3 h = o + float(t)*d;
                    double u, v;
                    if (axis==0){ u=h.y+hh; v=h.z+hd; } else if (axis==1){ u=h.x+hw; v=h.z+hd; } else { u=h.x+hw; v=h.y+hh; }
                    if (u>=p.uMin && u<=p.uMax && v>=p.vMin && v<=p.vMax) { bestT=t; best=int(pi); }
                }
                if (best >= 0) gPorts[best].role = (gPorts[best].role + 1) % 3;
            }
        }

        if (gShowSlicePlanes && !designMode && !devDipole && hasFieldForDomain()) {
            const Bounds sb = activeBounds();
            const float hw = 0.5f * sb.width;
            const float hh = 0.5f * sb.height;
            const float hd = 0.5f * sb.depth;
            const glm::vec4 vpS(0.0f, 0.0f, float(fbW), float(fbH));
            const float rad3 = (uiGeometry == 1) ? float(cylModel.radius()) : 0.0f;
            ImDrawList* bg = ImGui::GetBackgroundDrawList();

            struct PlaneStyle { ImU32 fill, border; const char* tag; };
            const PlaneStyle styles[3] = {
                { IM_COL32( 70,150,255, 45), IM_COL32( 90,170,255,225), "XY" },
                { IM_COL32(255,160, 40, 45), IM_COL32(255,180, 70,225), "ZX" },
                { IM_COL32( 80,220,120, 45), IM_COL32(110,235,150,225), "ZY" },
            };
            auto projPt = [&](const glm::vec3& w, ImVec2& out) -> bool {
                glm::vec3 s = glm::project(w, view, proj, vpS);
                if (s.z < 0.0f || s.z > 1.0f) return false;
                out = ImVec2(s.x, float(fbH) - s.y);
                return true;
            };

            for (int pl = 0; pl < 3; ++pl) {
                if (!secShown[pl]) continue;
                const float c = float(secSlicePosC[pl]); // centered perp coord
                std::vector<glm::vec3> poly;
                if (pl == 0 && uiGeometry == 1) {        // XY disc at z=c
                    const int seg = 48;
                    for (int i = 0; i < seg; ++i) {
                        const float a = 2.0f * 3.14159265358979f * i / seg;
                        poly.emplace_back(rad3*std::cos(a), rad3*std::sin(a), c);
                    }
                } else if (pl == 0) {                    // XY rect at z=c
                    poly = { {-hw,-hh,c}, { hw,-hh,c}, { hw, hh,c}, {-hw, hh,c} };
                } else if (pl == 1) {                    // ZX rect at y=c
                    poly = { {-hw,c,-hd}, { hw,c,-hd}, { hw,c, hd}, {-hw,c, hd} };
                } else {                                 // ZY rect at x=c
                    poly = { {c,-hh,-hd}, {c, hh,-hd}, {c, hh, hd}, {c,-hh, hd} };
                }
                std::vector<ImVec2> pts; pts.reserve(poly.size());
                bool ok = true;
                for (const glm::vec3& w : poly) {
                    ImVec2 s; if (!projPt(w, s)) { ok = false; break; }
                    pts.push_back(s);
                }
                if (!ok || pts.size() < 3) continue;
                bg->AddConvexPolyFilled(pts.data(), int(pts.size()), styles[pl].fill);
                bg->AddPolyline(pts.data(), int(pts.size()), styles[pl].border,
                                ImDrawFlags_Closed, 2.0f);
                ImVec2 ctr(0,0);
                for (const ImVec2& s : pts) { ctr.x += s.x; ctr.y += s.y; }
                ctr.x /= float(pts.size()); ctr.y /= float(pts.size());
                bg->AddText(ctr, styles[pl].border, styles[pl].tag);
            }
        }

        // ---- 3D field lines (streamlines through the volume) ----
        // Traced through the real 3D field at the current phase and drawn as
        // depth-projected polylines. An occupancy voxel grid keeps the lines
        // evenly spaced and bounds the cost.
        if (gView3D != 0 && !designMode && !devDipole && hasFieldForDomain()) {
            const Bounds fb3 = activeBounds();
            const double a3 = fb3.width, b3 = fb3.height, d3 = fb3.depth;
            // Ask the source where its domain actually is rather than branching
            // on uiGeometry here: that hand-rolled invariant is what bakeCloud
            // got wrong for the cylinder.
            double xLo, xHi, yLo, yHi, zLo, zHi;
            active()->domain(xLo, xHi, yLo, yHi, zLo, zHi);
            const double xc3 = 0.5*(xLo+xHi), yc3 = 0.5*(yLo+yHi), zc3 = 0.5*(zLo+zHi);
            const double peak = activePeak();
            if (peak > 0.0) {
                const glm::vec4 vpS(0.0f, 0.0f, float(fbW), float(fbH));
                ImDrawList* bg = ImGui::GetBackgroundDrawList();
                auto fireCol = [&](double t) -> ImU32 {
                    float r, g, b; fireColor(float(t), r, g, b);
                    return ImGui::ColorConvertFloat4ToU32(ImVec4(r, g, b, 0.85f));
                };
                auto fieldAt = [&](double x,double y,double z) -> std::array<double,3> {
                    return active()->fieldVector(x,y,z,phase);
                };
                auto toWorld = [&](double x,double y,double z) -> glm::vec3 {
                    return glm::vec3(float(x-xc3), float(y-yc3), float(z-zc3));
                };
                auto projW = [&](const glm::vec3& w, ImVec2& out) -> bool {
                    glm::vec3 s = glm::project(w, view, proj, vpS);
                    if (s.z < 0.0f || s.z > 1.0f) return false;
                    out = ImVec2(s.x, float(fbH) - s.y); return true;
                };

                const float lnDens = std::clamp(gFieldLineDensity, 0.3f, 2.5f);
                // Finer occupancy than before so strong-field regions (seeded more
                // heavily by the importance sampler below) can pack more lines.
                const int OX=std::max(5,int(24*lnDens)), OY=std::max(4,int(15*lnDens)), OZ=std::max(10,int(56*lnDens));
                std::vector<int> occ(size_t(OX)*OY*OZ, 0);
                auto vox = [&](double x,double y,double z) -> int {
                    const int ix=int((x-xLo)/(xHi-xLo)*OX);
                    const int iy=int((y-yLo)/(yHi-yLo)*OY);
                    const int iz=int((z-zLo)/(zHi-zLo)*OZ);
                    if(ix<0||ix>=OX||iy<0||iy>=OY||iz<0||iz>=OZ) return -1;
                    return (iz*OY+iy)*OX+ix;
                };

                const double diag = std::sqrt((xHi-xLo)*(xHi-xLo) +
                                              (yHi-yLo)*(yHi-yLo) +
                                              (zHi-zLo)*(zHi-zLo));
                const double hstep = diag/240.0;
                const double eps = 1e-3*peak;
                const int maxSteps = 600;
                // The occupancy grid is the real cost bound (each voxel is
                // marked once); the budget is only a safety ceiling for
                // pathological cases and must stay well ABOVE the work needed to
                // fill the volume, or the far end of the guide is left untraced.
                int evals = 0; const int budget = int(320000 * lnDens * lnDens);
                int sid = 0;

                // Copper occlusion: now that the conductor is opaque, a field point
                // is hidden when the straight line from it to the eye crosses the
                // copper sheet (it sits behind the trace from this viewpoint).
                const bool msOcc = (gUseMicro && gMicro);
                const glm::vec3 eyeC = gCamera.position();  // centered world
                double copperY = -1.0;
                if (msOcc) for (const MicrostripSim::Prim& p : gMicro->prims())
                    if (p.mat==MicrostripSim::Pec && ((p.kind==0 && p.ymin>1e-9) || p.kind==2))
                        { copperY = 0.5*(p.ymin+p.ymax); break; }
                auto inCopper = [&](double xd,double zd)->bool {
                    for (const MicrostripSim::Prim& p : gMicro->prims()) {
                        if (p.mat!=MicrostripSim::Pec) continue;
                        if (p.kind==0) { if (p.ymin<=1e-9) continue;   // skip the ground plane
                            if (xd>=p.xmin&&xd<=p.xmax&&zd>=p.zmin&&zd<=p.zmax) return true; }
                        else if (p.kind==2) { const double dx=xd-p.cx,dz=zd-p.cz,d2=dx*dx+dz*dz;
                            if (d2<=p.radius*p.radius && d2>=p.rinner*p.rinner) return true; }
                    }
                    return false;
                };
                auto occluded = [&](double x,double y,double z)->bool {
                    if (!msOcc || copperY<0.0) return false;
                    const glm::vec3 Pc = toWorld(x,y,z);       // centered
                    const double denom = double(eyeC.y) - Pc.y;
                    if (std::fabs(denom) < 1e-9) return false;
                    const double t = (copperY-0.5*b3 - Pc.y) / denom;   // hit copper plane
                    if (t <= 0.0 || t >= 1.0) return false;    // copper not between point and eye
                    const double xi = Pc.x + t*(double(eyeC.x)-Pc.x) + 0.5*a3;
                    const double zi = Pc.z + t*(double(eyeC.z)-Pc.z) + 0.5*d3;
                    return inCopper(xi, zi);
                };

                // Color/seed normalizer. On a GRIDDED source (microstrip / FDTD)
                // peak_ is the global instantaneous max of the grid, dominated by
                // the driven source cell (injectSource) and the conductor-edge
                // singularity -- values 1-2 orders above the field in the
                // substrate/air the streamlines actually traverse. Using it pins
                // every bulk segment to m/peak ~ 0.05, i.e. the blue floor of the
                // fire palette, so all lines look blue. There we take a high
                // percentile of |F| sampled across the interior instead: the
                // palette spans what the lines see and the lone hot source cell
                // cannot crush the scale.
                //
                // On the ANALYTIC sources (waveguide / cavity / cylindrical) there
                // is no singular cell, and the percentile is actively harmful: it
                // is a THIRD scale, differing from the 3D cloud (peakField, see
                // uInvPeak in cloud.vert), from the cross sections and from the
                // "Color scale" legend, so the same point came out one colour as a
                // cloud dot and another as a line through it. Keep peakField there.
                const bool gridSource = (gUseMicro && gMicro) || gUseFdtd;
                double colorNorm = peak;
                if (gridSource) {
                    std::mt19937 crng(0x1234ABCDu);
                    std::uniform_real_distribution<double> Cx(xLo,xHi),Cy(yLo,yHi),Cz(zLo,zHi);
                    std::vector<double> mags; mags.reserve(4096);
                    for (int s=0; s<6000 && int(mags.size())<4000; ++s) {
                        const double x=Cx(crng), y=Cy(crng), z=Cz(crng);
                        if (!active()->inside(x,y,z)) continue;
                        auto F = fieldAt(x,y,z);
                        const double m = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
                        if (m > eps) mags.push_back(m);
                    }
                    if (mags.size() >= 20) {
                        const size_t q = size_t(0.95 * double(mags.size()-1));
                        std::nth_element(mags.begin(), mags.begin()+q, mags.end());
                        colorNorm = std::max(mags[q], 1e-30);
                    }
                }
                const double invNorm = 1.0 / colorNorm;

                auto trace3 = [&](double sx,double sy,double sz) {
                    ++sid;
                    for (int dir=-1; dir<=1; dir+=2) {
                        double x=sx,y=sy,z=sz;
                        glm::vec3 wprev = toWorld(x,y,z);
                        int since=0, ac=0;
                        for (int s=0; s<maxSteps && evals<budget; ++s) {
                            auto F = fieldAt(x,y,z); ++evals;
                            double m = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
                            if (m < eps) break;
                            const double ux=F[0]/m, uy=F[1]/m, uz=F[2]/m;
                            auto F2 = fieldAt(x+dir*ux*hstep*0.5, y+dir*uy*hstep*0.5,
                                              z+dir*uz*hstep*0.5); ++evals;
                            double m2 = std::sqrt(F2[0]*F2[0]+F2[1]*F2[1]+F2[2]*F2[2]);
                            if (m2 < eps) break;
                            const double nx=x+dir*F2[0]/m2*hstep;
                            const double ny=y+dir*F2[1]/m2*hstep;
                            const double nz=z+dir*F2[2]/m2*hstep;
                            if (nx<xLo||nx>xHi||ny<yLo||ny>yHi||nz<zLo||nz>zHi) break;
                            if (!active()->inside(nx,ny,nz)) break;
                            const int cv = vox(nx,ny,nz);
                            if (cv < 0) break;
                            if (occ[cv]!=0 && occ[cv]!=sid) break;
                            const bool wasEmpty = (occ[cv]==0);
                            occ[cv] = sid;
                            const glm::vec3 wcur = toWorld(nx,ny,nz);
                            ImVec2 pa, pb;
                            if (!occluded(nx,ny,nz) && projW(wprev,pa) && projW(wcur,pb)) {
                                // sqrt spreads the mid-tones so the palette reads
                                // as a gradient rather than saturating; normalize
                                // by the interior percentile, not peak_.
                                const ImU32 col = fireCol(std::sqrt(std::min(1.0, m*invNorm)));
                                bg->AddLine(pa, pb, col, 1.5f);
                                if (++ac >= 26) {           // arrowhead in FIELD direction
                                    ac = 0;
                                    // The backward half (dir<0) integrates against
                                    // the field, so multiply by dir to recover the
                                    // true field direction; tip goes on the
                                    // downstream end.
                                    const float sgn = float(dir);
                                    const float dx=(pb.x-pa.x)*sgn, dy=(pb.y-pa.y)*sgn;
                                    const float ln=std::sqrt(dx*dx+dy*dy);
                                    if (ln>1e-3f){ const float hx=dx/ln,hy=dy/ln,px=-hy,py=hx,hl=4.5f;
                                        const ImVec2 tip = (dir>0)? pb : pa;
                                        bg->AddLine(tip,ImVec2(tip.x-hx*hl+px*hl*0.5f,tip.y-hy*hl+py*hl*0.5f),col,1.5f);
                                        bg->AddLine(tip,ImVec2(tip.x-hx*hl-px*hl*0.5f,tip.y-hy*hl-py*hl*0.5f),col,1.5f); }
                                }
                            }
                            wprev = wcur; x=nx; y=ny; z=nz;
                            if (wasEmpty) since=0; else if (++since>24) break;
                        }
                    }
                };

                // Monte-Carlo seeding, importance-sampled by field intensity:
                // candidate points are drawn at random across the whole volume
                // (so transverse TE/TM fields still get lines along all z) and
                // accepted as streamline seeds with probability (|F|/peak)^0.7 --
                // so seeds, and thus lines, concentrate where the field is strong,
                // exactly like the point cloud. A FIXED seed keeps the line set
                // stable as the field animates. The occupancy grid still stops any
                // two lines from overdrawing the same voxel.
                std::mt19937 rng(0x5EED3Du);
                std::uniform_real_distribution<double> U01(0.0, 1.0);
                std::uniform_real_distribution<double> Ux(xLo, xHi), Uy(yLo, yHi), Uz(zLo, zHi);
                const int seedAttempts = std::max(2000, int(9000 * lnDens * lnDens));
                for (int s=0; s<seedAttempts && evals<budget; ++s) {
                    const double x=Ux(rng), y=Uy(rng), z=Uz(rng);
                    if (!active()->inside(x,y,z)) continue;
                    auto F = fieldAt(x,y,z); ++evals;
                    const double m = std::sqrt(F[0]*F[0]+F[1]*F[1]+F[2]*F[2]);
                    if (m < eps) continue;
                    const double prob = std::pow(std::min(1.0, m*invNorm), 0.7);
                    if (U01(rng) > prob) continue;          // reject -> denser where |F| is large
                    const int cv = vox(x,y,z);
                    if (cv < 0 || occ[cv]!=0) continue;
                    trace3(x,y,z);
                }
            }
        }

        // ---- Dimension labels on waveguide edges ----
        // Same guard as the 3D mesh (analyticShape): the a/b/L cotas belong to the
        // analytic rect/cyl waveguide only. Without !gUseMicro && !devDipole they
        // leaked onto the microstrip/dipole scenes when switching domains.
        if (analyticShape) {
            const Bounds eb = activeBounds();
            const float hw = 0.5f * eb.width;
            const float hh = 0.5f * eb.height;
            const float hd = 0.5f * eb.depth;
            const float off = std::max({eb.width, eb.height, eb.depth}) * 0.05f;
            const glm::vec4 vp(0.0f, 0.0f, float(fbW), float(fbH));

            auto worldToScreen = [&](const glm::vec3& wp) -> glm::vec3 {
                return glm::project(wp, view, proj, vp);
            };

            auto drawDimLine = [&](const glm::vec3& w0, const glm::vec3& w1,
                                    const glm::vec3& ofs, const char* label)
            {
                glm::vec3 s0 = worldToScreen(w0 + ofs);
                glm::vec3 s1 = worldToScreen(w1 + ofs);
                if (s0.z < 0.0f || s0.z > 1.0f ||
                    s1.z < 0.0f || s1.z > 1.0f) return;
                ImVec2 p0(s0.x, float(fbH) - s0.y);
                ImVec2 p1(s1.x, float(fbH) - s1.y);
                ImVec2 mid((p0.x + p1.x) * 0.5f, (p0.y + p1.y) * 0.5f);

                ImDrawList* fg = ImGui::GetForegroundDrawList();
                const ImU32 lineCol = IM_COL32(220, 220, 180, 200);
                const ImU32 textCol = IM_COL32(255, 255, 220, 255);

                // Dimension line
                fg->AddLine(p0, p1, lineCol, 1.0f);

                // End ticks (perpendicular)
                float dx = p1.x - p0.x, dy = p1.y - p0.y;
                float len = std::sqrt(dx * dx + dy * dy);
                if (len > 1e-3f) {
                    float nx = -dy / len * 5.0f, ny = dx / len * 5.0f;
                    fg->AddLine(ImVec2(p0.x + nx, p0.y + ny),
                                ImVec2(p0.x - nx, p0.y - ny), lineCol, 1.0f);
                    fg->AddLine(ImVec2(p1.x + nx, p1.y + ny),
                                ImVec2(p1.x - nx, p1.y - ny), lineCol, 1.0f);
                }

                // Label with background
                ImVec2 ts = ImGui::CalcTextSize(label);
                ImVec2 tp(mid.x - ts.x * 0.5f, mid.y - ts.y - 4.0f);
                fg->AddRectFilled(ImVec2(tp.x - 2, tp.y - 1),
                                  ImVec2(tp.x + ts.x + 2, tp.y + ts.y + 1),
                                  IM_COL32(0, 0, 0, 180));
                fg->AddText(tp, textCol, label);
            };

            char wBuf[32], hBuf[32], dBuf[32];
            std::snprintf(dBuf, sizeof(dBuf), "L = %.1f mm", eb.depth * 1000.0f);

            if (uiGeometry == 0) {
                std::snprintf(wBuf, sizeof(wBuf), "a = %.2f mm",
                              eb.width * 1000.0f);
                std::snprintf(hBuf, sizeof(hBuf), "b = %.2f mm",
                              eb.height * 1000.0f);
                // Width: bottom-front edge, offset down
                drawDimLine({-hw, -hh, -hd}, {+hw, -hh, -hd},
                            {0, -off, 0}, wBuf);
                // Height: front-left edge, offset left
                drawDimLine({-hw, -hh, -hd}, {-hw, +hh, -hd},
                            {-off, 0, 0}, hBuf);
                // Depth: bottom-right edge, offset down
                drawDimLine({+hw, -hh, -hd}, {+hw, -hh, +hd},
                            {0, -off, 0}, dBuf);
            } else {
                float rad = hw;
                std::snprintf(wBuf, sizeof(wBuf), "D = %.2f mm",
                              eb.width * 1000.0f);
                // Diameter: across front face, offset forward
                drawDimLine({-rad, 0.0f, -hd}, {+rad, 0.0f, -hd},
                            {0, 0, -off}, wBuf);
                // Depth: along side, offset outward
                drawDimLine({rad, 0.0f, -hd}, {rad, 0.0f, +hd},
                            {off, 0, 0}, dBuf);
            }
        }

        // ---- Export button (bottom-right) ----
        static ExportOptions expOpts;
        static bool openExportModal = false;
        static std::string lastExportMsg;
        {
            int fbW5, fbH5;
            glfwGetFramebufferSize(window, &fbW5, &fbH5);
            const float bw = 170.0f, bh = 34.0f;
            ImGui::SetNextWindowPos(ImVec2(float(fbW5) - bw - 10.0f,
                                           float(fbH5) - bh - 10.0f));
            ImGui::SetNextWindowSize(ImVec2(bw, bh));
            ImGui::Begin("##exportBtn", nullptr,
                         ImGuiWindowFlags_NoDecoration |
                         ImGuiWindowFlags_NoMove |
                         ImGuiWindowFlags_NoSavedSettings |
                         ImGuiWindowFlags_NoFocusOnAppearing);
            if (ImGui::Button("Exportar simulacao", ImVec2(bw - 16.0f, bh - 10.0f)))
                openExportModal = true;
            ImGui::End();
        }
        if (openExportModal) {
            ImGui::OpenPopup("Exportar simulacao");
            openExportModal = false;
        }
        if (ImGui::BeginPopupModal("Exportar simulacao", nullptr,
                                    ImGuiWindowFlags_AlwaysAutoResize))
        {
            // 0 = single stacked GIF, 1 = single stacked PNG, 2 = LaTeX project
            static int fmt = 2;
            ImGui::RadioButton("Projeto LaTeX", &fmt, 2);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                "Cria uma pasta com o nome do arquivo contendo:\n"
                "  main.tex      documento com tabelas e figuras\n"
                "  pictures/     uma imagem PNG por item selecionado\n"
                "  <nome>.gif    a animacao, se GIF estiver marcado\n"
                "Compile com: pdflatex main.tex");
            ImGui::RadioButton("GIF (imagem unica)", &fmt, 0); ImGui::SameLine();
            ImGui::RadioButton("PNG (imagem unica)", &fmt, 1);
            ImGui::Separator();
            ImGui::Text("Selecione o que exportar:");
            ImGui::Checkbox("Simulacao 3D",       &expOpts.incScene3D);
            ImGui::Checkbox("Corte XY",           &expOpts.incXY);
            ImGui::Checkbox("Corte ZX",           &expOpts.incZX);
            ImGui::Checkbox("Corte ZY",           &expOpts.incZY);
            ImGui::Checkbox("Espectro / cutoff",  &expOpts.incSpectrum);
            ImGui::Checkbox("Escala de cores",    &expOpts.incColorBar);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                "Anexa a barra de escala ao lado de CADA figura de campo, em vez\n"
                "de gerar uma figura separada. Assim uma figura levada sozinha\n"
                "para um artigo continua dizendo o que suas cores significam.");
            ImGui::Separator();
            static bool latexGif = false;   // also write the animation next to main.tex
            static int  nPhases  = 1;       // phase series for the animated items
            static char titleBuf[256]  = "";
            static char authorBuf[256] = "";
            if (fmt == 2) {
                ImGui::SliderInt("Fases por figura", &nPhases, 1, 8);
                if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                    "Itens animados (3D e cortes) sao exportados nesta quantidade\n"
                    "de fases igualmente espacadas em um ciclo, cada uma como uma\n"
                    "figura propria. E assim que se mostra um modo viajante no papel.\n"
                    "Itens estaticos (espectro, escala) saem uma vez so.");
                ImGui::Checkbox("Gerar tambem o GIF", &latexGif);
                ImGui::SliderInt("Largura (px)", &expOpts.imageW, 480, 2400);
                ImGui::SliderInt("Altura (px)",  &expOpts.imageH, 360, 1800);
                if (ImGui::IsItemHovered()) ImGui::SetTooltip(
                    "Caixa alvo de cada figura. Cada item usa o maximo dela que\n"
                    "sua proporcao real permitir: um corte ZX longo ocupa toda a\n"
                    "largura, um corte XY quadrado e limitado pela altura.");
                ImGui::InputText("Titulo",  titleBuf,  IM_ARRAYSIZE(titleBuf));
                if (ImGui::IsItemHovered()) ImGui::SetTooltip("Vazio = gerado do modo.");
                ImGui::InputText("Autor",   authorBuf, IM_ARRAYSIZE(authorBuf));
            }
            if (fmt == 0 || (fmt == 2 && latexGif)) {
                ImGui::SliderInt("Frames",      &expOpts.frames,  8, 120);
                ImGui::SliderInt("Delay (cs)",  &expOpts.delayCs, 2, 20);
            }
            static char fnameBuf[256] = "waveguide_export";
            ImGui::InputText(fmt == 2 ? "Pasta" : "Arquivo", fnameBuf, IM_ARRAYSIZE(fnameBuf));
            if (!lastExportMsg.empty())
                ImGui::TextColored(ImVec4(0.6f,1.0f,0.6f,1.0f), "%s", lastExportMsg.c_str());

            if (ImGui::Button("Exportar", ImVec2(120, 0))) {
                expOpts.latexProject = (fmt == 2);
                expOpts.isGif = (fmt == 0) || (fmt == 2 && latexGif);
                expOpts.filename = fnameBuf;
                expOpts.title  = titleBuf;
                expOpts.author = authorBuf;
                expOpts.phasesDeg.clear();
                for (int i = 0; i < std::max(1, nPhases); ++i)
                    expOpts.phasesDeg.push_back(360.0 * i / std::max(1, nPhases));
                ExportContext exc;
                exc.rect     = &rectModel;
                exc.cyl      = &cylModel;
                exc.renderer = &renderer;
                exc.geometry = uiGeometry;
                exc.fieldKind= uiField;
                exc.modeType = uiModeType;
                exc.modeM    = uiM;
                exc.modeN    = uiN;
                exc.modeL    = uiL;
                exc.structure= uiStructure;
                exc.freqHz   = double(uiFreqGHz) * 1e9;
                exc.mediumName = gMediumName;
                exc.powerW     = double(uiPowerW);
                exc.view3D     = gView3D;
                exc.lineDensity= gFieldLineDensity;
                exc.cloudInvPeak     = gCloudInvPeak;
                exc.cloudMeanSpacing = gCloudMeanSpacing;
                exc.cloudOpaque      = gCloudOpaque;
                exc.sectionStreamlines = gFieldLines;
                exc.sectionLineDensity = gSecLineDensity;
                const bool ok = runExport(exc, expOpts);
                lastExportMsg = ok
                    ? (expOpts.latexProject
                        ? std::string("OK: ") + expOpts.filename +
                          "/main.tex (+ pictures/). Compile: pdflatex main.tex"
                        : std::string("OK: ") + expOpts.filename +
                          (expOpts.isGif ? ".gif" : ".png"))
                    : std::string("Falha ao exportar.");
                // Nothing to restore: the exporter renders the same baked cloud
                // the live view uses and never touches the sphere instance
                // buffer. Re-uploading particles here (a leftover from when the
                // exporter DID clobber that buffer) was what made the old
                // solid-sphere plot appear over the live 3D view after an
                // export, and stay there.
            }
            ImGui::SameLine();
            if (ImGui::Button("Fechar", ImVec2(120, 0))) {
                lastExportMsg.clear();
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }

        // ---- Axis gizmo overlay (X=red, Y=green, Z=blue) ----
        if (showAxis) {
            const float gSize = 120.0f;
            int fbW4, fbH4;
            glfwGetFramebufferSize(window, &fbW4, &fbH4);
            ImGui::SetNextWindowPos(ImVec2(10.0f, float(fbH4) - gSize - 10.0f));
            ImGui::SetNextWindowSize(ImVec2(gSize, gSize));
            ImGui::Begin("##axisGizmo", nullptr,
                         ImGuiWindowFlags_NoDecoration |
                         ImGuiWindowFlags_NoBackground |
                         ImGuiWindowFlags_NoMove |
                         ImGuiWindowFlags_NoInputs |
                         ImGuiWindowFlags_NoSavedSettings |
                         ImGuiWindowFlags_NoFocusOnAppearing);

            const ImVec2 wp = ImGui::GetWindowPos();
            const ImVec2 ctr(wp.x + gSize * 0.5f, wp.y + gSize * 0.5f);
            ImDrawList* dl3 = ImGui::GetWindowDrawList();
            const float L = gSize * 0.35f;

            // Use view matrix rotation to project world axes onto the screen.
            const glm::mat4 v3 = view;
            auto project = [&](const glm::vec3& w) -> ImVec2 {
                const glm::vec4 e = v3 * glm::vec4(w, 0.0f);
                return ImVec2(ctr.x + e.x * L, ctr.y - e.y * L);
            };

            struct Ax { glm::vec3 dir; ImU32 col; const char* label; };
            const Ax axes[3] = {
                {{1,0,0}, IM_COL32(230, 60, 60,255), "X"},
                {{0,1,0}, IM_COL32(60, 220, 60,255), "Y"},
                {{0,0,1}, IM_COL32(60, 120,255,255), "Z"},
            };
            // Depth-sort so the front axis draws on top
            int order[3] = {0,1,2};
            float depth[3];
            for (int k = 0; k < 3; ++k)
                depth[k] = (v3 * glm::vec4(axes[k].dir, 0.0f)).z;
            std::sort(order, order + 3, [&](int a, int b){ return depth[a] < depth[b]; });

            for (int k = 0; k < 3; ++k) {
                const Ax& a = axes[order[k]];
                const ImVec2 tip = project(a.dir);
                dl3->AddLine(ctr, tip, a.col, 2.0f);
                // arrowhead
                const float dx = tip.x - ctr.x, dy = tip.y - ctr.y;
                const float len = std::sqrt(dx*dx + dy*dy);
                if (len > 1e-3f) {
                    const float ux = dx/len, uy = dy/len;
                    const float px_ = -uy, py_ = ux;
                    const float h = 6.0f;
                    dl3->AddTriangleFilled(
                        tip,
                        ImVec2(tip.x - ux*h + px_*h*0.5f, tip.y - uy*h + py_*h*0.5f),
                        ImVec2(tip.x - ux*h - px_*h*0.5f, tip.y - uy*h - py_*h*0.5f),
                        a.col);
                }
                dl3->AddText(ImVec2(tip.x + 4.0f, tip.y - 6.0f), a.col, a.label);
            }
            ImGui::End();
        }

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return EXIT_SUCCESS;
}
