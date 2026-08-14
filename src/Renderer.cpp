#include "Renderer.hpp"

#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <vector>

namespace waveguide {

namespace {
constexpr float kPi = 3.14159265358979323846f;
}

Renderer::~Renderer() {
    if (sphereVbo_)     glDeleteBuffers(1, &sphereVbo_);
    if (instanceVbo_)   glDeleteBuffers(1, &instanceVbo_);
    if (sphereVao_)     glDeleteVertexArrays(1, &sphereVao_);
    if (sphereProgram_) glDeleteProgram(sphereProgram_);

    if (gridVbo_)       glDeleteBuffers(1, &gridVbo_);
    if (gridVao_)       glDeleteVertexArrays(1, &gridVao_);

    if (boxVbo_)        glDeleteBuffers(1, &boxVbo_);
    if (boxVao_)        glDeleteVertexArrays(1, &boxVao_);

    if (cylVbo_)        glDeleteBuffers(1, &cylVbo_);
    if (cylVao_)        glDeleteVertexArrays(1, &cylVao_);

    if (cloudVbo_)      glDeleteBuffers(1, &cloudVbo_);
    if (cloudVao_)      glDeleteVertexArrays(1, &cloudVao_);
    if (cloudProgram_)  glDeleteProgram(cloudProgram_);

    if (lineProgram_)   glDeleteProgram(lineProgram_);
}

std::string Renderer::readFile(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "Failed to open shader file: %s\n", path.c_str());
        return {};
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

GLuint Renderer::compileShader(GLenum type, const std::string& src, const std::string& tag) {
    GLuint s = glCreateShader(type);
    const char* cstr = src.c_str();
    glShaderSource(s, 1, &cstr, nullptr);
    glCompileShader(s);

    GLint ok = GL_FALSE;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetShaderInfoLog(s, sizeof(log), nullptr, log);
        std::fprintf(stderr, "[%s] shader compile error:\n%s\n", tag.c_str(), log);
        glDeleteShader(s);
        return 0;
    }
    return s;
}

GLuint Renderer::loadProgram(const std::string& vertPath, const std::string& fragPath) {
    const std::string vs = readFile(vertPath);
    const std::string fs = readFile(fragPath);
    if (vs.empty() || fs.empty()) return 0;

    GLuint v = compileShader(GL_VERTEX_SHADER,   vs, vertPath);
    GLuint f = compileShader(GL_FRAGMENT_SHADER, fs, fragPath);
    if (!v || !f) {
        if (v) glDeleteShader(v);
        if (f) glDeleteShader(f);
        return 0;
    }

    GLuint prog = glCreateProgram();
    glAttachShader(prog, v);
    glAttachShader(prog, f);
    glLinkProgram(prog);

    GLint ok = GL_FALSE;
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetProgramInfoLog(prog, sizeof(log), nullptr, log);
        std::fprintf(stderr, "Program link error (%s, %s):\n%s\n",
                     vertPath.c_str(), fragPath.c_str(), log);
        glDeleteProgram(prog);
        prog = 0;
    }

    glDeleteShader(v);
    glDeleteShader(f);
    return prog;
}

// ---------------------- mesh construction ----------------------

// Low-poly UV sphere (stacks/sectors). Unit radius.
void Renderer::buildSphereMesh(int stacks, int sectors) {
    std::vector<float> verts;
    verts.reserve(size_t(stacks * sectors * 6 * 3));

    auto pos = [](float t, float p) {
        return glm::vec3(std::sin(t) * std::cos(p),
                         std::cos(t),
                         std::sin(t) * std::sin(p));
    };

    for (int i = 0; i < stacks; ++i) {
        const float t1 = float(i)     / float(stacks) * kPi;
        const float t2 = float(i + 1) / float(stacks) * kPi;
        for (int j = 0; j < sectors; ++j) {
            const float p1 = float(j)     / float(sectors) * 2.0f * kPi;
            const float p2 = float(j + 1) / float(sectors) * 2.0f * kPi;

            const glm::vec3 v1 = pos(t1, p1);
            const glm::vec3 v2 = pos(t1, p2);
            const glm::vec3 v3 = pos(t2, p1);
            const glm::vec3 v4 = pos(t2, p2);

            // Two triangles per quad.
            verts.insert(verts.end(), {v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z});
            verts.insert(verts.end(), {v2.x, v2.y, v2.z, v4.x, v4.y, v4.z, v3.x, v3.y, v3.z});
        }
    }

    sphereVertexCount_ = GLsizei(verts.size() / 3);

    glGenVertexArrays(1, &sphereVao_);
    glGenBuffers(1, &sphereVbo_);

    glBindVertexArray(sphereVao_);

    // --- mesh positions (location 0) ---
    glBindBuffer(GL_ARRAY_BUFFER, sphereVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(verts.size() * sizeof(float)),
                 verts.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);

    // --- per-instance attributes (locations 1 = pos, 2 = color) ---
    glGenBuffers(1, &instanceVbo_);
    glBindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
    // data uploaded later in updateParticles(); reserve nothing for now.

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Particle),
                          (void*)offsetof(Particle, x));
    glVertexAttribDivisor(1, 1);

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Particle),
                          (void*)offsetof(Particle, r));
    glVertexAttribDivisor(2, 1);

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(Particle),
                          (void*)offsetof(Particle, intensity));
    glVertexAttribDivisor(3, 1);

    glBindVertexArray(0);

    // --- separate instanced VAO for the port marker spheres (reuses the mesh) ---
    glGenVertexArrays(1, &portVao_);
    glGenBuffers(1, &portVbo_);
    glBindVertexArray(portVao_);
    glBindBuffer(GL_ARRAY_BUFFER, sphereVbo_);            // shared mesh (location 0)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindBuffer(GL_ARRAY_BUFFER, portVbo_);              // per-instance Particle
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)offsetof(Particle, x));
    glVertexAttribDivisor(1, 1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)offsetof(Particle, r));
    glVertexAttribDivisor(2, 1);
    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, sizeof(Particle), (void*)offsetof(Particle, intensity));
    glVertexAttribDivisor(3, 1);
    glBindVertexArray(0);
}

void Renderer::updateStructure(const std::vector<float>& opaqueVerts,
                               const std::vector<float>& translucentVerts)
{
    structCount_  = GLsizei(opaqueVerts.size() / 9);
    structTCount_ = GLsizei(translucentVerts.size() / 9);
    glBindBuffer(GL_ARRAY_BUFFER, structVbo_);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(opaqueVerts.size() * sizeof(float)),
                 opaqueVerts.empty() ? nullptr : opaqueVerts.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, structTVbo_);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(translucentVerts.size() * sizeof(float)),
                 translucentVerts.empty() ? nullptr : translucentVerts.data(), GL_DYNAMIC_DRAW);
}

void Renderer::drawStructure(const glm::mat4& view, const glm::mat4& proj,
                             bool translucentPass, float alpha)
{
    const GLsizei n = translucentPass ? structTCount_ : structCount_;
    if (n <= 0) return;
    glUseProgram(solidProgram_);
    glUniformMatrix4fv(glGetUniformLocation(solidProgram_, "uView"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(solidProgram_, "uProj"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniform1f(glGetUniformLocation(solidProgram_, "uAlpha"), alpha);
    glEnable(GL_DEPTH_TEST);
    if (translucentPass) {                 // blend, keep depth read but no write
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
    } else {
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    }
    glBindVertexArray(translucentPass ? structTVao_ : structVao_);
    glDrawArrays(GL_TRIANGLES, 0, n);
    glBindVertexArray(0);
    if (translucentPass) { glDepthMask(GL_TRUE); glDisable(GL_BLEND); }
}

void Renderer::drawSpheres(const glm::mat4& view, const glm::mat4& proj,
                           const std::vector<Particle>& spheres)
{
    if (spheres.empty()) return;
    glBindBuffer(GL_ARRAY_BUFFER, portVbo_);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(spheres.size() * sizeof(Particle)),
                 spheres.data(), GL_DYNAMIC_DRAW);
    glUseProgram(sphereProgram_);
    glUniformMatrix4fv(glGetUniformLocation(sphereProgram_, "uView"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(sphereProgram_, "uProj"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniform1f(glGetUniformLocation(sphereProgram_, "uSphereRadius"), 1.0f); // radius from intensity
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glBindVertexArray(portVao_);
    glDrawArraysInstanced(GL_TRIANGLES, 0, sphereVertexCount_, GLsizei(spheres.size()));
    glBindVertexArray(0);
}

// A simple XZ-plane grid of white lines, centered on origin.
void Renderer::buildGrid(float size, int divisions) {
    std::vector<float> verts;
    const float half = size * 0.5f;
    const float step = size / float(divisions);

    // Lines parallel to Z
    for (int i = 0; i <= divisions; ++i) {
        const float x = -half + float(i) * step;
        verts.push_back(x); verts.push_back(0.0f); verts.push_back(-half);
        verts.push_back(x); verts.push_back(0.0f); verts.push_back( half);
    }
    // Lines parallel to X
    for (int i = 0; i <= divisions; ++i) {
        const float z = -half + float(i) * step;
        verts.push_back(-half); verts.push_back(0.0f); verts.push_back(z);
        verts.push_back( half); verts.push_back(0.0f); verts.push_back(z);
    }

    gridVertexCount_ = GLsizei(verts.size() / 3);

    glGenVertexArrays(1, &gridVao_);
    glGenBuffers(1, &gridVbo_);
    glBindVertexArray(gridVao_);
    glBindBuffer(GL_ARRAY_BUFFER, gridVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(verts.size() * sizeof(float)),
                 verts.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);
}

// Wireframe box matching the waveguide, centered on origin.
void Renderer::buildBoxWireframe(const Bounds& b) {
    const float x0 = -0.5f * b.width,  x1 =  0.5f * b.width;
    const float y0 = -0.5f * b.height, y1 =  0.5f * b.height;
    const float z0 = -0.5f * b.depth,  z1 =  0.5f * b.depth;

    const float c[8][3] = {
        {x0,y0,z0},{x1,y0,z0},{x1,y1,z0},{x0,y1,z0},
        {x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1},
    };
    const int edges[12][2] = {
        {0,1},{1,2},{2,3},{3,0},
        {4,5},{5,6},{6,7},{7,4},
        {0,4},{1,5},{2,6},{3,7},
    };

    std::vector<float> verts;
    verts.reserve(12 * 2 * 3);
    for (auto& e : edges) {
        for (int idx : e) {
            verts.push_back(c[idx][0]);
            verts.push_back(c[idx][1]);
            verts.push_back(c[idx][2]);
        }
    }
    boxVertexCount_ = GLsizei(verts.size() / 3);

    glGenVertexArrays(1, &boxVao_);
    glGenBuffers(1, &boxVbo_);
    glBindVertexArray(boxVao_);
    glBindBuffer(GL_ARRAY_BUFFER, boxVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(verts.size() * sizeof(float)),
                 verts.data(),
                 GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);
}

// ---------------------- cylinder mesh ----------------------

void Renderer::buildCylinderMesh(float radius, float length) {
    const int segs = 64;
    std::vector<float> verts;
    const float halfL = length * 0.5f;

    // Side surface
    for (int i = 0; i < segs; ++i) {
        const float a0 = 2.0f * kPi * float(i)     / float(segs);
        const float a1 = 2.0f * kPi * float(i + 1) / float(segs);
        const float c0 = std::cos(a0), s0 = std::sin(a0);
        const float c1 = std::cos(a1), s1 = std::sin(a1);

        verts.insert(verts.end(), {radius*c0, radius*s0, -halfL,
                                   radius*c1, radius*s1, -halfL,
                                   radius*c0, radius*s0,  halfL});
        verts.insert(verts.end(), {radius*c1, radius*s1, -halfL,
                                   radius*c1, radius*s1,  halfL,
                                   radius*c0, radius*s0,  halfL});
    }
    // Top cap
    for (int i = 0; i < segs; ++i) {
        const float a0 = 2.0f * kPi * float(i)     / float(segs);
        const float a1 = 2.0f * kPi * float(i + 1) / float(segs);
        verts.insert(verts.end(), {0.0f, 0.0f, halfL,
                                   radius*std::cos(a0), radius*std::sin(a0), halfL,
                                   radius*std::cos(a1), radius*std::sin(a1), halfL});
    }
    // Bottom cap
    for (int i = 0; i < segs; ++i) {
        const float a0 = 2.0f * kPi * float(i)     / float(segs);
        const float a1 = 2.0f * kPi * float(i + 1) / float(segs);
        verts.insert(verts.end(), {0.0f, 0.0f, -halfL,
                                   radius*std::cos(a1), radius*std::sin(a1), -halfL,
                                   radius*std::cos(a0), radius*std::sin(a0), -halfL});
    }

    cylVertexCount_ = GLsizei(verts.size() / 3);

    if (!cylVao_) {
        glGenVertexArrays(1, &cylVao_);
        glGenBuffers(1, &cylVbo_);
    }
    glBindVertexArray(cylVao_);
    glBindBuffer(GL_ARRAY_BUFFER, cylVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(verts.size() * sizeof(float)),
                 verts.data(), GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glBindVertexArray(0);
}

void Renderer::updateCylinder(float radius, float length) {
    buildCylinderMesh(radius, length);
}

// ---------------------- public API ----------------------

bool Renderer::init(const std::string& shaderDir,
                    const std::vector<Particle>& particles,
                    const Bounds& bounds)
{
    sphereProgram_ = loadProgram(shaderDir + "/point.vert", shaderDir + "/point.frag");
    lineProgram_   = loadProgram(shaderDir + "/line.vert",  shaderDir + "/line.frag");
    cloudProgram_  = loadProgram(shaderDir + "/cloud.vert", shaderDir + "/cloud.frag");
    solidProgram_  = loadProgram(shaderDir + "/solid.vert", shaderDir + "/solid.frag");
    if (!sphereProgram_ || !lineProgram_ || !cloudProgram_ || !solidProgram_) return false;

    // Solid structure VAOs (opaque + translucent): pos(3) normal(3) colour(3).
    {
        const GLsizei stride = 9 * sizeof(float);
        GLuint* vaos[2] = { &structVao_,  &structTVao_ };
        GLuint* vbos[2] = { &structVbo_,  &structTVbo_ };
        for (int i = 0; i < 2; ++i) {
            glGenVertexArrays(1, vaos[i]);
            glGenBuffers(1, vbos[i]);
            glBindVertexArray(*vaos[i]);
            glBindBuffer(GL_ARRAY_BUFFER, *vbos[i]);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
            glEnableVertexAttribArray(2);
            glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, (void*)(6 * sizeof(float)));
            glBindVertexArray(0);
        }
    }

    buildSphereMesh(/*stacks=*/10, /*sectors=*/10);

    // Baked point-sprite cloud VAO (pos + Fourier terms).
    glGenVertexArrays(1, &cloudVao_);
    glGenBuffers(1, &cloudVbo_);
    glBindVertexArray(cloudVao_);
    glBindBuffer(GL_ARRAY_BUFFER, cloudVbo_);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(CloudPoint), (void*)offsetof(CloudPoint, x));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(CloudPoint), (void*)offsetof(CloudPoint, dc));
    glBindVertexArray(0);

    // Grid scale: 4x the waveguide's largest dimension.
    const float gridSize = std::max({bounds.width, bounds.height, bounds.depth}) * 4.0f;
    buildGrid(gridSize, /*divisions=*/20);

    buildBoxWireframe(bounds);

    updateParticles(particles);
    return true;
}

void Renderer::updateBounds(const Bounds& b) {
    // Re-upload the 12 edges of the box wireframe to the existing VBO.
    const float x0 = -0.5f * b.width,  x1 =  0.5f * b.width;
    const float y0 = -0.5f * b.height, y1 =  0.5f * b.height;
    const float z0 = -0.5f * b.depth,  z1 =  0.5f * b.depth;
    const float c[8][3] = {
        {x0,y0,z0},{x1,y0,z0},{x1,y1,z0},{x0,y1,z0},
        {x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1},
    };
    const int edges[12][2] = {
        {0,1},{1,2},{2,3},{3,0},
        {4,5},{5,6},{6,7},{7,4},
        {0,4},{1,5},{2,6},{3,7},
    };
    std::vector<float> verts;
    verts.reserve(12 * 2 * 3);
    for (auto& e : edges)
        for (int idx : e) {
            verts.push_back(c[idx][0]);
            verts.push_back(c[idx][1]);
            verts.push_back(c[idx][2]);
        }
    boxVertexCount_ = GLsizei(verts.size() / 3);
    glBindBuffer(GL_ARRAY_BUFFER, boxVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(verts.size() * sizeof(float)),
                 verts.data(),
                 GL_DYNAMIC_DRAW);
}

void Renderer::updateParticles(const std::vector<Particle>& particles) {
    glBindBuffer(GL_ARRAY_BUFFER, instanceVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(particles.size() * sizeof(Particle)),
                 particles.data(),
                 GL_DYNAMIC_DRAW);
    instanceCount_ = GLsizei(particles.size());
}

void Renderer::updateCloud(const std::vector<CloudPoint>& points) {
    glBindBuffer(GL_ARRAY_BUFFER, cloudVbo_);
    glBufferData(GL_ARRAY_BUFFER,
                 GLsizeiptr(points.size() * sizeof(CloudPoint)),
                 points.data(),
                 GL_DYNAMIC_DRAW);
    cloudCount_ = GLsizei(points.size());
}

void Renderer::drawCloud(const glm::mat4& view, const glm::mat4& proj,
                         float phase, float invPeak, float pointScale,
                         bool opaque) {
    if (cloudCount_ == 0) return;
    glEnable(GL_PROGRAM_POINT_SIZE);
    if (opaque) {
        // Depth write ON and no blending: the nearest dot wins each pixel, so
        // the cloud has a visible surface. The shader still discards outside
        // the sprite disc, and a discarded fragment writes no depth, so the
        // dots stay round rather than becoming squares.
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
    } else {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
    }
    glUseProgram(cloudProgram_);
    // float, not int: gl_loader only resolves glUniform1f/3fv/Matrix4fv.
    glUniform1f(glGetUniformLocation(cloudProgram_, "uOpaque"), opaque ? 1.0f : 0.0f);
    glUniformMatrix4fv(glGetUniformLocation(cloudProgram_, "uView"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(cloudProgram_, "uProj"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniform1f(glGetUniformLocation(cloudProgram_, "uPhase"), phase);
    glUniform1f(glGetUniformLocation(cloudProgram_, "uInvPeak"), invPeak);
    glUniform1f(glGetUniformLocation(cloudProgram_, "uPointScale"), pointScale);
    glBindVertexArray(cloudVao_);
    glDrawArrays(GL_POINTS, 0, cloudCount_);
    glBindVertexArray(0);
    glUseProgram(0);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
}

void Renderer::draw(const glm::mat4& view,
                    const glm::mat4& proj,
                    float sphereRadius,
                    bool showGrid,
                    bool showBox,
                    bool showCylinder)
{
    glUseProgram(lineProgram_);
    glUniformMatrix4fv(glGetUniformLocation(lineProgram_, "uView"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(lineProgram_, "uProj"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniform1f(glGetUniformLocation(lineProgram_, "uAlpha"), 1.0f);

    // --- grid (under everything) ---
    if (showGrid) {
        const float gridColor[3] = { 0.62f, 0.63f, 0.66f };
        glUniform3fv(glGetUniformLocation(lineProgram_, "uColor"), 1, gridColor);
        glLineWidth(1.0f);
        glBindVertexArray(gridVao_);
        glDrawArrays(GL_LINES, 0, gridVertexCount_);
    }

    // --- wireframe box ---
    if (showBox) {
        const float boxColor[3] = { 0.9f, 0.25f, 0.25f };
        glUniform3fv(glGetUniformLocation(lineProgram_, "uColor"), 1, boxColor);
        glLineWidth(2.0f);
        glBindVertexArray(boxVao_);
        glDrawArrays(GL_LINES, 0, boxVertexCount_);
    }

    // --- instanced spheres ---
    glUseProgram(sphereProgram_);
    glUniformMatrix4fv(glGetUniformLocation(sphereProgram_, "uView"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(sphereProgram_, "uProj"), 1, GL_FALSE, glm::value_ptr(proj));
    glUniform1f(glGetUniformLocation(sphereProgram_, "uSphereRadius"), sphereRadius);

    glBindVertexArray(sphereVao_);
    glDrawArraysInstanced(GL_TRIANGLES, 0, sphereVertexCount_, instanceCount_);

    // --- semi-transparent cylinder (after opaque objects) ---
    if (showCylinder && cylVertexCount_ > 0) {
        glUseProgram(lineProgram_);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        const float cylColor[3] = { 0.5f, 0.7f, 1.0f };
        glUniform3fv(glGetUniformLocation(lineProgram_, "uColor"), 1, cylColor);
        glUniform1f(glGetUniformLocation(lineProgram_, "uAlpha"), 0.15f);
        glBindVertexArray(cylVao_);
        glDrawArrays(GL_TRIANGLES, 0, cylVertexCount_);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
    }

    glBindVertexArray(0);
    glUseProgram(0);
}

bool Renderer::renderOffscreen(int w, int h,
                                const glm::mat4& view,
                                const glm::mat4& proj,
                                const OffscreenOpts& opts,
                                std::vector<unsigned char>& outPixels)
{
    if (w <= 0 || h <= 0) return false;

    GLuint fbo = 0, colorTex = 0, depthRbo = 0;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);

    glGenTextures(1, &colorTex);
    glBindTexture(GL_TEXTURE_2D, colorTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,
                           GL_TEXTURE_2D, colorTex, 0);

    glGenRenderbuffers(1, &depthRbo);
    glBindRenderbuffer(GL_RENDERBUFFER, depthRbo);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                              GL_RENDERBUFFER, depthRbo);

    bool ok = (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
    if (ok) {
        GLint prevVp[4];
        glGetIntegerv(GL_VIEWPORT, prevVp);
        glViewport(0, 0, w, h);
        glClearColor(0.90f, 0.91f, 0.93f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        // Scene furniture (floor grid, guide outline). Passing a zero sphere
        // radius keeps the instanced-sphere pass from contributing when the
        // cloud is what we actually want.
        draw(view, proj, opts.useCloud ? 0.0f : opts.sphereRadius,
             opts.showGrid, opts.showBox, opts.showCylinder);
        // Exactly the call the live view makes, so an exported figure and the
        // screen show the same geometry AND the same palette.
        if (opts.useCloud)
            drawCloud(view, proj, opts.cloudPhase, opts.cloudInvPeak,
                      opts.cloudPointScale, opts.cloudOpaque);

        outPixels.assign(size_t(w) * size_t(h) * 4, 0);
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, outPixels.data());
        // Flip vertically (OpenGL origin is bottom-left)
        const int stride = w * 4;
        std::vector<unsigned char> row(stride);
        for (int y = 0; y < h / 2; ++y) {
            unsigned char* a = outPixels.data() + y * stride;
            unsigned char* b = outPixels.data() + (h - 1 - y) * stride;
            std::copy(a, a + stride, row.begin());
            std::copy(b, b + stride, a);
            std::copy(row.begin(), row.end(), b);
        }
        glViewport(prevVp[0], prevVp[1], prevVp[2], prevVp[3]);
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDeleteRenderbuffers(1, &depthRbo);
    glDeleteTextures(1, &colorTex);
    glDeleteFramebuffers(1, &fbo);
    return ok;
}

} // namespace waveguide
