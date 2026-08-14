#pragma once

#include "TEmnModel.hpp"
#include "gl_loader.h"

#include <glm/glm.hpp>
#include <string>
#include <vector>

namespace waveguide {

// Renders a particle cloud in the visual style of kavan010/Atoms:
//   - low-poly sphere mesh, instanced once per particle
//   - flat (unlit) color from the fire heatmap
//   - white grid floor
//   - (optional) wireframe box marking the waveguide walls
class Renderer {
public:
    Renderer() = default;
    ~Renderer();

    Renderer(const Renderer&)            = delete;
    Renderer& operator=(const Renderer&) = delete;

    // Compiles shaders, builds the sphere mesh, uploads particles + grid.
    bool init(const std::string& shaderDir,
              const std::vector<Particle>& particles,
              const Bounds& bounds);

    // Replace the particle buffer (e.g. when sampling count changes).
    void updateParticles(const std::vector<Particle>& particles);
    void updateBounds(const Bounds& bounds);
    void updateCylinder(float radius, float length);

    // Baked point-sprite cloud: uploaded once (on change), animated on the GPU.
    void updateCloud(const std::vector<CloudPoint>& points);
    // `opaque` makes near dots OCCLUDE far ones (depth write on, no blending),
    // so the cloud reads as a solid body instead of a translucent haze you can
    // see straight through. Translucent is still available: it reveals interior
    // structure at the cost of every dot summing with everything behind it.
    void drawCloud(const glm::mat4& view, const glm::mat4& proj,
                   float phase, float invPeak, float pointScale,
                   bool opaque = true);

    // Draw a small set of solid, diffuse-lit 3D spheres. Per-sphere world radius
    // is taken from Particle.intensity (color from r,g,b). Used for the
    // source/sense port markers; positions must be in the same centred world
    // space as the cloud.
    void drawSpheres(const glm::mat4& view, const glm::mat4& proj,
                     const std::vector<Particle>& spheres);

    // Solid structure mesh (microstrip ground / substrate / copper). Vertices are
    // interleaved [px,py,pz, nx,ny,nz, r,g,b] in the same centred world space as
    // the cloud. The opaque set writes depth (so copper occludes the cloud and the
    // port spheres); the translucent set (substrate) is blended without depth write.
    void updateStructure(const std::vector<float>& opaqueVerts,
                         const std::vector<float>& translucentVerts);
    void drawStructure(const glm::mat4& view, const glm::mat4& proj,
                       bool translucentPass, float alpha);

    void draw(const glm::mat4& view,
              const glm::mat4& proj,
              float sphereRadius,
              bool showGrid = true,
              bool showBox  = true,
              bool showCylinder = false);

    // What an offscreen pass should contain. Defaults reproduce the live view:
    // the baked point-sprite cloud, which is what drawCloud() puts on screen.
    //
    // The legacy instanced-sphere path (useCloud = false) draws the particle
    // buffer with each model's own colormap, which is NOT the palette the UI or
    // the colour bar uses. Keep it only for callers that specifically want it.
    struct OffscreenOpts {
        bool  showGrid        = true;
        bool  showBox         = true;
        bool  showCylinder    = false;
        bool  useCloud        = true;
        float cloudPhase      = 0.0f;
        float cloudInvPeak    = 1.0f;
        float cloudPointScale = 1.0f;
        bool  cloudOpaque     = true;
        float sphereRadius    = 0.0f;   // sphere path only
    };

    // Render the scene into an offscreen FBO and read back RGBA8 pixels.
    // Returns true on success. `outPixels` is resized to w*h*4.
    bool renderOffscreen(int w, int h,
                         const glm::mat4& view,
                         const glm::mat4& proj,
                         const OffscreenOpts& opts,
                         std::vector<unsigned char>& outPixels);

private:
    GLuint loadProgram(const std::string& vertPath, const std::string& fragPath);
    static std::string readFile(const std::string& path);
    static GLuint      compileShader(GLenum type, const std::string& src, const std::string& tag);

    void buildSphereMesh(int stacks, int sectors);
    void buildGrid(float size, int divisions);
    void buildBoxWireframe(const Bounds& b);
    void buildCylinderMesh(float radius, float length);

    // Sphere instancing
    GLuint sphereProgram_ = 0;
    GLuint sphereVao_     = 0;
    GLuint sphereVbo_     = 0;   // mesh vertices (shared)
    GLuint instanceVbo_   = 0;   // per-particle Particle structs
    GLsizei sphereVertexCount_ = 0;
    GLsizei instanceCount_     = 0;

    // Dedicated instanced buffer for the port marker spheres (shares the mesh).
    GLuint portVao_ = 0;
    GLuint portVbo_ = 0;

    // Solid structure mesh (opaque + translucent), pos/normal/colour interleaved.
    GLuint solidProgram_ = 0;
    GLuint structVao_ = 0,  structVbo_ = 0;
    GLuint structTVao_ = 0, structTVbo_ = 0;
    GLsizei structCount_ = 0, structTCount_ = 0;

    // Baked point-sprite cloud
    GLuint cloudProgram_ = 0;
    GLuint cloudVao_     = 0;
    GLuint cloudVbo_     = 0;
    GLsizei cloudCount_  = 0;

    // Grid floor
    GLuint lineProgram_ = 0;
    GLuint gridVao_     = 0;
    GLuint gridVbo_     = 0;
    GLsizei gridVertexCount_ = 0;

    // Wireframe box (waveguide walls)
    GLuint boxVao_ = 0;
    GLuint boxVbo_ = 0;
    GLsizei boxVertexCount_ = 0;

    // Semi-transparent cylinder (cylindrical waveguide)
    GLuint cylVao_ = 0;
    GLuint cylVbo_ = 0;
    GLsizei cylVertexCount_ = 0;
};

} // namespace waveguide
