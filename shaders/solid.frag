#version 330 core

// Diffuse-lit solid surface. uAlpha lets the substrate slab render translucent
// while the ground plane and copper stay opaque.

in vec3 vNormal;
in vec3 vColor;

uniform float uAlpha;

out vec4 FragColor;

void main() {
    vec3 N = normalize(vNormal);
    vec3 L = normalize(vec3(-0.35, 0.85, 0.40));
    // Two-sided: light the back faces too (slabs are viewed from both sides).
    float diff  = max(abs(dot(N, L)), 0.0);
    float shade = 0.38 + 0.62 * diff;
    FragColor = vec4(clamp(vColor * shade, 0.0, 1.0), uAlpha);
}
