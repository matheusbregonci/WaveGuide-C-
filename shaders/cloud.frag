#version 330 core

// Soft round point sprite. Light-background heatmap: blue (low) -> dark red
// (high), so every level stays visible on the pale viewport (a white top end
// would disappear).

in  float vT;
out vec4  FragColor;

// >0.5 = the cloud is drawn as a solid body (depth write on, blending off), so
// a soft alpha rim would have nothing to blend against. Harden the disc instead.
// Carried as a float because gl_loader only resolves glUniform1f.
uniform float uOpaque;

vec3 fire(float t) {
    t = clamp(t, 0.0, 1.0) * 5.0;
    int i = int(t);
    float l = t - float(i);
    if (i < 1) return mix(vec3(0.20,0.45,0.95), vec3(0.10,0.72,0.85), l);
    if (i < 2) return mix(vec3(0.10,0.72,0.85), vec3(0.20,0.72,0.25), l);
    if (i < 3) return mix(vec3(0.20,0.72,0.25), vec3(0.98,0.70,0.10), l);
    if (i < 4) return mix(vec3(0.98,0.70,0.10), vec3(0.92,0.25,0.10), l);
    return          mix(vec3(0.92,0.25,0.10),   vec3(0.55,0.00,0.08), l);
}

void main() {
    vec2 d = gl_PointCoord - vec2(0.5);
    float r2 = dot(d, d);
    if (r2 > 0.25 || vT < 0.02) discard;   // round sprite, drop faint nodes
    if (uOpaque > 0.5) {
        // Trim the ragged outermost ring so the discs still read as circles
        // without an alpha ramp, then shade fully opaque.
        if (r2 > 0.22) discard;
        FragColor = vec4(fire(vT), 1.0);
        return;
    }
    float a = smoothstep(0.25, 0.03, r2);  // soft edge
    FragColor = vec4(fire(vT), a);
}
