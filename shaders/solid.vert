#version 330 core

// Solid triangle mesh with per-vertex normal + colour, used for the microstrip
// structure (ground plane, substrate slab, copper trace). Writes depth so the
// field cloud and the port spheres are properly occluded by the copper.

layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in vec3 aColor;

uniform mat4 uView;
uniform mat4 uProj;

out vec3 vNormal;
out vec3 vColor;

void main() {
    gl_Position = uProj * uView * vec4(aPos, 1.0);
    vNormal = aNormal;
    vColor  = aColor;
}
