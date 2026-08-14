#pragma once

// Minimal OpenGL 3.3 function loader - pure C++, no Python, no glad.
//
// Windows' opengl32.dll only statically exports OpenGL 1.1 symbols, so any
// 2.0+ function (shaders, VBOs, VAOs, uniforms) must be loaded at runtime
// via wglGetProcAddress. GLFW gives us `glfwGetProcAddress` which wraps
// that detail, so we only need to:
//   1. declare function-pointer variables for every modern GL entry point
//      we actually call
//   2. populate them once after context creation
//
// The standard `PFN...PROC` typedefs come from <GL/glext.h>, which ships
// with every MinGW-w64 / MSYS2 distribution.

#if defined(_WIN32)
  #ifndef WIN32_LEAN_AND_MEAN
  #define WIN32_LEAN_AND_MEAN
  #endif
  #ifndef NOMINMAX
  #define NOMINMAX
  #endif
  #include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glext.h>

#ifdef __cplusplus
extern "C" {
#endif

// ---- Shader / program ----
extern PFNGLCREATESHADERPROC      glCreateShader;
extern PFNGLDELETESHADERPROC      glDeleteShader;
extern PFNGLSHADERSOURCEPROC      glShaderSource;
extern PFNGLCOMPILESHADERPROC     glCompileShader;
extern PFNGLGETSHADERIVPROC       glGetShaderiv;
extern PFNGLGETSHADERINFOLOGPROC  glGetShaderInfoLog;

extern PFNGLCREATEPROGRAMPROC     glCreateProgram;
extern PFNGLDELETEPROGRAMPROC     glDeleteProgram;
extern PFNGLATTACHSHADERPROC      glAttachShader;
extern PFNGLLINKPROGRAMPROC       glLinkProgram;
extern PFNGLGETPROGRAMIVPROC      glGetProgramiv;
extern PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog;
extern PFNGLUSEPROGRAMPROC        glUseProgram;

// ---- Buffers / VAO ----
extern PFNGLGENBUFFERSPROC               glGenBuffers;
extern PFNGLDELETEBUFFERSPROC            glDeleteBuffers;
extern PFNGLBINDBUFFERPROC               glBindBuffer;
extern PFNGLBUFFERDATAPROC               glBufferData;

extern PFNGLGENVERTEXARRAYSPROC          glGenVertexArrays;
extern PFNGLDELETEVERTEXARRAYSPROC       glDeleteVertexArrays;
extern PFNGLBINDVERTEXARRAYPROC          glBindVertexArray;

extern PFNGLENABLEVERTEXATTRIBARRAYPROC  glEnableVertexAttribArray;
extern PFNGLVERTEXATTRIBPOINTERPROC      glVertexAttribPointer;
extern PFNGLVERTEXATTRIBDIVISORPROC      glVertexAttribDivisor;
extern PFNGLDRAWARRAYSINSTANCEDPROC      glDrawArraysInstanced;

// ---- Framebuffers / renderbuffers / textures (for offscreen capture) ----
extern PFNGLGENFRAMEBUFFERSPROC         glGenFramebuffers;
extern PFNGLDELETEFRAMEBUFFERSPROC      glDeleteFramebuffers;
extern PFNGLBINDFRAMEBUFFERPROC         glBindFramebuffer;
extern PFNGLFRAMEBUFFERTEXTURE2DPROC    glFramebufferTexture2D;
extern PFNGLFRAMEBUFFERRENDERBUFFERPROC glFramebufferRenderbuffer;
extern PFNGLCHECKFRAMEBUFFERSTATUSPROC  glCheckFramebufferStatus;
extern PFNGLGENRENDERBUFFERSPROC        glGenRenderbuffers;
extern PFNGLDELETERENDERBUFFERSPROC     glDeleteRenderbuffers;
extern PFNGLBINDRENDERBUFFERPROC        glBindRenderbuffer;
extern PFNGLRENDERBUFFERSTORAGEPROC     glRenderbufferStorage;

// ---- Uniforms ----
extern PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
extern PFNGLUNIFORMMATRIX4FVPROC   glUniformMatrix4fv;
extern PFNGLUNIFORM1FPROC          glUniform1f;
extern PFNGLUNIFORM3FVPROC         glUniform3fv;

// ---- Loader entry point ----
typedef void* (*GLLoadProc)(const char* name);

// Populates all pointers above. Returns non-zero on success, 0 if any
// required function could not be resolved.
int loadOpenGLFunctions(GLLoadProc proc);

#ifdef __cplusplus
} // extern "C"
#endif
