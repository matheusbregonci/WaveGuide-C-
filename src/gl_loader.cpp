#include "gl_loader.h"

#include <cstdio>

// Definitions of the function pointers declared in gl_loader.h.
PFNGLCREATESHADERPROC      glCreateShader      = nullptr;
PFNGLDELETESHADERPROC      glDeleteShader      = nullptr;
PFNGLSHADERSOURCEPROC      glShaderSource      = nullptr;
PFNGLCOMPILESHADERPROC     glCompileShader     = nullptr;
PFNGLGETSHADERIVPROC       glGetShaderiv       = nullptr;
PFNGLGETSHADERINFOLOGPROC  glGetShaderInfoLog  = nullptr;

PFNGLCREATEPROGRAMPROC     glCreateProgram     = nullptr;
PFNGLDELETEPROGRAMPROC     glDeleteProgram     = nullptr;
PFNGLATTACHSHADERPROC      glAttachShader      = nullptr;
PFNGLLINKPROGRAMPROC       glLinkProgram       = nullptr;
PFNGLGETPROGRAMIVPROC      glGetProgramiv      = nullptr;
PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog = nullptr;
PFNGLUSEPROGRAMPROC        glUseProgram        = nullptr;

PFNGLGENBUFFERSPROC               glGenBuffers              = nullptr;
PFNGLDELETEBUFFERSPROC            glDeleteBuffers           = nullptr;
PFNGLBINDBUFFERPROC               glBindBuffer              = nullptr;
PFNGLBUFFERDATAPROC               glBufferData              = nullptr;

PFNGLGENVERTEXARRAYSPROC          glGenVertexArrays         = nullptr;
PFNGLDELETEVERTEXARRAYSPROC       glDeleteVertexArrays      = nullptr;
PFNGLBINDVERTEXARRAYPROC          glBindVertexArray         = nullptr;

PFNGLENABLEVERTEXATTRIBARRAYPROC  glEnableVertexAttribArray = nullptr;
PFNGLVERTEXATTRIBPOINTERPROC      glVertexAttribPointer     = nullptr;
PFNGLVERTEXATTRIBDIVISORPROC      glVertexAttribDivisor     = nullptr;
PFNGLDRAWARRAYSINSTANCEDPROC      glDrawArraysInstanced     = nullptr;

PFNGLGENFRAMEBUFFERSPROC         glGenFramebuffers         = nullptr;
PFNGLDELETEFRAMEBUFFERSPROC      glDeleteFramebuffers      = nullptr;
PFNGLBINDFRAMEBUFFERPROC         glBindFramebuffer         = nullptr;
PFNGLFRAMEBUFFERTEXTURE2DPROC    glFramebufferTexture2D    = nullptr;
PFNGLFRAMEBUFFERRENDERBUFFERPROC glFramebufferRenderbuffer = nullptr;
PFNGLCHECKFRAMEBUFFERSTATUSPROC  glCheckFramebufferStatus  = nullptr;
PFNGLGENRENDERBUFFERSPROC        glGenRenderbuffers        = nullptr;
PFNGLDELETERENDERBUFFERSPROC     glDeleteRenderbuffers     = nullptr;
PFNGLBINDRENDERBUFFERPROC        glBindRenderbuffer        = nullptr;
PFNGLRENDERBUFFERSTORAGEPROC     glRenderbufferStorage     = nullptr;

PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation = nullptr;
PFNGLUNIFORMMATRIX4FVPROC   glUniformMatrix4fv   = nullptr;
PFNGLUNIFORM1FPROC          glUniform1f          = nullptr;
PFNGLUNIFORM3FVPROC         glUniform3fv         = nullptr;

namespace {
// Resolve one pointer and record failures.
template <typename Fn>
bool grab(GLLoadProc proc, const char* name, Fn& out) {
    out = reinterpret_cast<Fn>(proc(name));
    if (!out) {
        std::fprintf(stderr, "OpenGL loader: missing entry point '%s'\n", name);
        return false;
    }
    return true;
}
} // namespace

int loadOpenGLFunctions(GLLoadProc proc) {
    bool ok = true;

    ok &= grab(proc, "glCreateShader",      glCreateShader);
    ok &= grab(proc, "glDeleteShader",      glDeleteShader);
    ok &= grab(proc, "glShaderSource",      glShaderSource);
    ok &= grab(proc, "glCompileShader",     glCompileShader);
    ok &= grab(proc, "glGetShaderiv",       glGetShaderiv);
    ok &= grab(proc, "glGetShaderInfoLog",  glGetShaderInfoLog);

    ok &= grab(proc, "glCreateProgram",     glCreateProgram);
    ok &= grab(proc, "glDeleteProgram",     glDeleteProgram);
    ok &= grab(proc, "glAttachShader",      glAttachShader);
    ok &= grab(proc, "glLinkProgram",       glLinkProgram);
    ok &= grab(proc, "glGetProgramiv",      glGetProgramiv);
    ok &= grab(proc, "glGetProgramInfoLog", glGetProgramInfoLog);
    ok &= grab(proc, "glUseProgram",        glUseProgram);

    ok &= grab(proc, "glGenBuffers",        glGenBuffers);
    ok &= grab(proc, "glDeleteBuffers",     glDeleteBuffers);
    ok &= grab(proc, "glBindBuffer",        glBindBuffer);
    ok &= grab(proc, "glBufferData",        glBufferData);

    ok &= grab(proc, "glGenVertexArrays",   glGenVertexArrays);
    ok &= grab(proc, "glDeleteVertexArrays", glDeleteVertexArrays);
    ok &= grab(proc, "glBindVertexArray",   glBindVertexArray);

    ok &= grab(proc, "glEnableVertexAttribArray", glEnableVertexAttribArray);
    ok &= grab(proc, "glVertexAttribPointer",     glVertexAttribPointer);
    ok &= grab(proc, "glVertexAttribDivisor",     glVertexAttribDivisor);
    ok &= grab(proc, "glDrawArraysInstanced",     glDrawArraysInstanced);

    ok &= grab(proc, "glGenFramebuffers",         glGenFramebuffers);
    ok &= grab(proc, "glDeleteFramebuffers",      glDeleteFramebuffers);
    ok &= grab(proc, "glBindFramebuffer",         glBindFramebuffer);
    ok &= grab(proc, "glFramebufferTexture2D",    glFramebufferTexture2D);
    ok &= grab(proc, "glFramebufferRenderbuffer", glFramebufferRenderbuffer);
    ok &= grab(proc, "glCheckFramebufferStatus",  glCheckFramebufferStatus);
    ok &= grab(proc, "glGenRenderbuffers",        glGenRenderbuffers);
    ok &= grab(proc, "glDeleteRenderbuffers",     glDeleteRenderbuffers);
    ok &= grab(proc, "glBindRenderbuffer",        glBindRenderbuffer);
    ok &= grab(proc, "glRenderbufferStorage",     glRenderbufferStorage);

    ok &= grab(proc, "glGetUniformLocation", glGetUniformLocation);
    ok &= grab(proc, "glUniformMatrix4fv",   glUniformMatrix4fv);
    ok &= grab(proc, "glUniform1f",          glUniform1f);
    ok &= grab(proc, "glUniform3fv",         glUniform3fv);

    return ok ? 1 : 0;
}
