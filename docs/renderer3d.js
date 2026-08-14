// WebGL2 point-sprite cloud, the same one the desktop draws.
//
// The GLSL is ported from shaders/cloud.{vert,frag}: GLES 3.0 rather than
// GL 3.3 core, which means a #version 300 es header and explicit precision
// qualifiers. Two desktop-only details had to go:
//   - glEnable(GL_PROGRAM_POINT_SIZE) does not exist here; writing gl_PointSize
//     is always allowed in GLES and the enable would be an INVALID_ENUM.
//   - lineWidth() above 1.0 is ignored by every browser, so the guide outline is
//     drawn as plain 1px lines and anything needing thickness lives on the 2D
//     canvas instead.
//
// The colour ramp is the same blue -> dark red used by the desktop cloud, the
// cross sections and the colour bar. It is duplicated here for the same reason
// it is duplicated in cloud.frag: shaders cannot include Colormap.hpp.

const VERT = `#version 300 es
precision highp float;
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aF;    // dc, a2, b2 of |field|^2(phase)

uniform mat4  uView;
uniform mat4  uProj;
uniform float uPhase;
uniform float uInvPeak;
uniform float uPointScale;

out float vT;

void main() {
    float c = cos(2.0 * uPhase);
    float s = sin(2.0 * uPhase);
    float m2   = aF.x + aF.y * c + aF.z * s;
    float ac   = length(vec2(aF.y, aF.z));
    float env  = sqrt(max(aF.x + ac, 0.0));
    float inst = sqrt(max(m2, 0.0));
    // A little envelope mixed in so nodes never fully vanish (no hard blink).
    float t = clamp((0.2 * env + 0.8 * inst) * uInvPeak, 0.0, 1.0);
    vT = t;

    vec4 vp = uView * vec4(aPos, 1.0);
    gl_Position = uProj * vp;
    float dist = max(0.02, -vp.z);
    gl_PointSize = clamp(uPointScale * (0.4 + t) / dist, 1.5, 48.0);
}`;

const FRAG = `#version 300 es
precision highp float;
in  float vT;
out vec4  FragColor;

vec3 fire(float t) {
    t = clamp(t, 0.0, 1.0) * 5.0;
    int i = int(t);
    float l = t - float(i);
    if (i < 1) return mix(vec3(0.20,0.45,0.95), vec3(0.10,0.72,0.85), l);
    if (i < 2) return mix(vec3(0.10,0.72,0.85), vec3(0.20,0.72,0.25), l);
    if (i < 3) return mix(vec3(0.20,0.72,0.25), vec3(0.98,0.70,0.10), l);
    if (i < 4) return mix(vec3(0.98,0.70,0.10), vec3(0.92,0.25,0.10), l);
    return          mix(vec3(0.92,0.25,0.10),  vec3(0.55,0.00,0.08), l);
}

void main() {
    vec2 d = gl_PointCoord - vec2(0.5);
    float r2 = dot(d, d);
    // Opaque discs: depth writing is on, so near dots occlude far ones and the
    // cloud reads as a body instead of a haze summed along the line of sight.
    if (r2 > 0.22 || vT < 0.02) discard;
    FragColor = vec4(fire(vT), 1.0);
}`;

const LINE_VERT = `#version 300 es
precision highp float;
layout(location = 0) in vec3 aPos;
uniform mat4 uView; uniform mat4 uProj;
void main() { gl_Position = uProj * uView * vec4(aPos, 1.0); }`;

const LINE_FRAG = `#version 300 es
precision highp float;
uniform vec3 uColor;
out vec4 FragColor;
void main() { FragColor = vec4(uColor, 1.0); }`;

// ---- minimal mat4 (column-major, same convention as glm) ----
export function perspective(fovy, aspect, near, far) {
  const f = 1 / Math.tan(fovy / 2), nf = 1 / (near - far);
  return new Float32Array([
    f / aspect, 0, 0, 0,
    0, f, 0, 0,
    0, 0, (far + near) * nf, -1,
    0, 0, 2 * far * near * nf, 0]);
}
export function lookAt(eye, center, up) {
  const z = norm(sub(eye, center));
  const x = norm(cross(up, z));
  const y = cross(z, x);
  return new Float32Array([
    x[0], y[0], z[0], 0,
    x[1], y[1], z[1], 0,
    x[2], y[2], z[2], 0,
    -dot(x, eye), -dot(y, eye), -dot(z, eye), 1]);
}
const sub = (a, b) => [a[0]-b[0], a[1]-b[1], a[2]-b[2]];
const dot = (a, b) => a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
const cross = (a, b) => [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
const norm = a => { const l = Math.hypot(...a) || 1; return [a[0]/l, a[1]/l, a[2]/l]; };

function compile(gl, vsrc, fsrc) {
  const mk = (type, src) => {
    const s = gl.createShader(type);
    gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
      throw new Error(gl.getShaderInfoLog(s) + '\n' + src);
    return s;
  };
  const p = gl.createProgram();
  gl.attachShader(p, mk(gl.VERTEX_SHADER, vsrc));
  gl.attachShader(p, mk(gl.FRAGMENT_SHADER, fsrc));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS))
    throw new Error(gl.getProgramInfoLog(p));
  return p;
}

export class Renderer3D {
  constructor(canvas) {
    const gl = canvas.getContext('webgl2', { antialias: true, depth: true });
    if (!gl) throw new Error('WebGL2 indisponivel neste navegador');
    this.gl = gl;
    this.canvas = canvas;

    this.progCloud = compile(gl, VERT, FRAG);
    this.progLine = compile(gl, LINE_VERT, LINE_FRAG);

    this.vaoCloud = gl.createVertexArray();
    this.vboCloud = gl.createBuffer();
    gl.bindVertexArray(this.vaoCloud);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vboCloud);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 24, 0);
    gl.enableVertexAttribArray(1);
    gl.vertexAttribPointer(1, 3, gl.FLOAT, false, 24, 12);
    gl.bindVertexArray(null);

    this.vaoLine = gl.createVertexArray();
    this.vboLine = gl.createBuffer();
    gl.bindVertexArray(this.vaoLine);
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vboLine);
    gl.enableVertexAttribArray(0);
    gl.vertexAttribPointer(0, 3, gl.FLOAT, false, 12, 0);
    gl.bindVertexArray(null);

    this.count = 0;
    this.lineCount = 0;
    this.invPeak = 1;
    this.pointScale = 40;

    // Orbit camera
    this.yaw = 0.6; this.pitch = 0.42; this.dist = 1;
    this._drag = false;
    canvas.addEventListener('pointerdown', e => {
      this._drag = true; this._px = e.clientX; this._py = e.clientY;
      canvas.setPointerCapture(e.pointerId);
    });
    canvas.addEventListener('pointerup', e => {
      this._drag = false;
      try { canvas.releasePointerCapture(e.pointerId); } catch {}
    });
    canvas.addEventListener('pointermove', e => {
      if (!this._drag) return;
      this.yaw -= (e.clientX - this._px) * 0.008;
      this.pitch = Math.max(-1.5, Math.min(1.5,
                    this.pitch + (e.clientY - this._py) * 0.008));
      this._px = e.clientX; this._py = e.clientY;
    });
    canvas.addEventListener('wheel', e => {
      e.preventDefault();
      this.dist *= Math.exp(e.deltaY * 0.0012);
    }, { passive: false });
  }

  // points: Float32Array of [x,y,z,dc,a2,b2] per vertex, already centred.
  setCloud(points, invPeak, meanSpacing) {
    const gl = this.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vboCloud);
    gl.bufferData(gl.ARRAY_BUFFER, points, gl.STATIC_DRAW);
    this.count = points.length / 6;
    this.invPeak = invPeak;
    this.meanSpacing = meanSpacing;
  }

  // Wireframe of the guide, in the same centred frame as the cloud.
  setOutline(dom, geometry) {
    const cx = 0.5 * (dom.x0 + dom.x1), cy = 0.5 * (dom.y0 + dom.y1),
          cz = 0.5 * (dom.z0 + dom.z1);
    const hx = 0.5 * (dom.x1 - dom.x0), hy = 0.5 * (dom.y1 - dom.y0),
          hz = 0.5 * (dom.z1 - dom.z0);
    this.radius = Math.hypot(hx, hy, hz);
    const v = [];
    if (geometry === 0) {
      const C = [];
      for (let i = 0; i < 8; i++)
        C.push([(i & 1 ? hx : -hx), (i & 2 ? hy : -hy), (i & 4 ? hz : -hz)]);
      const E = [[0,1],[2,3],[4,5],[6,7],[0,2],[1,3],[4,6],[5,7],
                 [0,4],[1,5],[2,6],[3,7]];
      for (const [a, b] of E) v.push(...C[a], ...C[b]);
    } else {
      const R = hx, seg = 64;
      for (const zz of [-hz, hz])
        for (let i = 0; i < seg; i++) {
          const a0 = 2*Math.PI*i/seg, a1 = 2*Math.PI*(i+1)/seg;
          v.push(R*Math.cos(a0), R*Math.sin(a0), zz,
                 R*Math.cos(a1), R*Math.sin(a1), zz);
        }
      for (let i = 0; i < 4; i++) {
        const a = 2*Math.PI*i/4;
        v.push(R*Math.cos(a), R*Math.sin(a), -hz,
               R*Math.cos(a), R*Math.sin(a),  hz);
      }
    }
    const gl = this.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.vboLine);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(v), gl.STATIC_DRAW);
    this.lineCount = v.length / 3;
    this._offset = [cx, cy, cz];
  }

  frame(phase) {
    const gl = this.gl, c = this.canvas;
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const w = Math.max(1, Math.round(c.clientWidth * dpr));
    const h = Math.max(1, Math.round(c.clientHeight * dpr));
    if (c.width !== w || c.height !== h) { c.width = w; c.height = h; }
    gl.viewport(0, 0, w, h);
    gl.clearColor(0.925, 0.933, 0.945, 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.enable(gl.DEPTH_TEST);
    gl.disable(gl.BLEND);

    const r = (this.radius || 1) * 2.2 * this.dist;
    const eye = [r * Math.cos(this.pitch) * Math.sin(this.yaw),
                 r * Math.sin(this.pitch),
                 r * Math.cos(this.pitch) * Math.cos(this.yaw)];
    const view = lookAt(eye, [0, 0, 0], [0, 1, 0]);
    const proj = perspective(Math.PI / 4, w / h, (this.radius || 1) * 0.01,
                             (this.radius || 1) * 40);

    if (this.lineCount) {
      gl.useProgram(this.progLine);
      gl.uniformMatrix4fv(gl.getUniformLocation(this.progLine, 'uView'), false, view);
      gl.uniformMatrix4fv(gl.getUniformLocation(this.progLine, 'uProj'), false, proj);
      gl.uniform3f(gl.getUniformLocation(this.progLine, 'uColor'), 0.55, 0.18, 0.18);
      gl.bindVertexArray(this.vaoLine);
      gl.drawArrays(gl.LINES, 0, this.lineCount);
    }

    if (this.count) {
      gl.useProgram(this.progCloud);
      gl.uniformMatrix4fv(gl.getUniformLocation(this.progCloud, 'uView'), false, view);
      gl.uniformMatrix4fv(gl.getUniformLocation(this.progCloud, 'uProj'), false, proj);
      gl.uniform1f(gl.getUniformLocation(this.progCloud, 'uPhase'), phase);
      gl.uniform1f(gl.getUniformLocation(this.progCloud, 'uInvPeak'), this.invPeak);
      // Same expression the desktop uses: dot size from the mean sample spacing,
      // carried into pixels through the projection.
      const ps = (this.meanSpacing || 1e-3) * 1.5 * proj[5] * h * 0.5;
      gl.uniform1f(gl.getUniformLocation(this.progCloud, 'uPointScale'), ps);
      gl.bindVertexArray(this.vaoCloud);
      gl.drawArrays(gl.POINTS, 0, this.count);
    }
    gl.bindVertexArray(null);
  }
}
