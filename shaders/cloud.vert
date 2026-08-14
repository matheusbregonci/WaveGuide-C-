#version 330 core

// Baked, GPU-animated field cloud. Each point carries the Fourier terms of
// |field|^2 in the animation phase; the instantaneous intensity is reconstructed
// here, so the CPU never resamples the field per frame.

layout(location = 0) in vec3 aPos; // centered world position
layout(location = 1) in vec3 aF;   // dc, a2, b2 of |field|^2(phase)

uniform mat4  uView;
uniform mat4  uProj;
uniform float uPhase;      // global animation phase
uniform float uInvPeak;    // 1 / peak|field|
uniform float uPointScale; // pixels ~ uPointScale * (0.4+t) / dist

out float vT;

void main() {
    float c = cos(2.0 * uPhase);
    float s = sin(2.0 * uPhase);
    float m2   = aF.x + aF.y * c + aF.z * s;         // instantaneous |field|^2
    float ac   = length(vec2(aF.y, aF.z));
    float env  = sqrt(max(aF.x + ac, 0.0));          // envelope (peak over time)
    float inst = sqrt(max(m2, 0.0));                 // instantaneous |field|
    // Blend a little envelope in so nodes never fully vanish (no hard blink).
    float t = clamp((0.2 * env + 0.8 * inst) * uInvPeak, 0.0, 1.0);
    vT = t;

    vec4 vp = uView * vec4(aPos, 1.0);
    gl_Position = uProj * vp;
    float dist = max(0.02, -vp.z);
    gl_PointSize = clamp(uPointScale * (0.4 + t) / dist, 1.5, 48.0);
}
