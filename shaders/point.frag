#version 330 core

// Diffuse-lit sphere: ambient + Lambert diffuse from a fixed light, plus a small
// specular highlight, so the instanced spheres read as solid 3D balls.

in  vec3 vColor;
in  vec3 vNormal;
out vec4 FragColor;

void main() {
    vec3 N = normalize(vNormal);
    vec3 L = normalize(vec3(-0.35, 0.85, 0.40));   // light from up / up-left
    float diff  = max(dot(N, L), 0.0);
    float shade = 0.32 + 0.68 * diff;              // ambient 0.32
    float spec  = pow(diff, 24.0) * 0.5;           // subtle highlight
    vec3  col   = vColor * shade + vec3(spec);
    FragColor = vec4(clamp(col, 0.0, 1.0), 1.0);
}
