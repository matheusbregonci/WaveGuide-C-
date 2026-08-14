#version 330 core

// Instanced low-poly sphere. Per-instance radius = uSphereRadius * aIntensity
// (encode a world radius in aIntensity with uSphereRadius = 1). Outputs the
// world-space normal so the fragment shader can apply simple diffuse lighting.

layout(location = 0) in vec3  aVertex;      // sphere mesh vertex (unit sphere, local)
layout(location = 1) in vec3  aInstancePos; // per-instance world position
layout(location = 2) in vec3  aInstanceCol; // per-instance color
layout(location = 3) in float aIntensity;   // per-instance radius scale

uniform mat4  uView;
uniform mat4  uProj;
uniform float uSphereRadius;

out vec3 vColor;
out vec3 vNormal;

void main() {
    float radius = uSphereRadius * aIntensity;
    vec3 worldPos = aInstancePos + aVertex * radius;
    gl_Position = uProj * uView * vec4(worldPos, 1.0);
    vColor  = aInstanceCol;
    vNormal = normalize(aVertex);   // unit-sphere vertex = world normal (no rotation)
}
