// Cross-section plot on a 2D canvas rather than in WebGL.
//
// Deliberate: browsers ignore gl.lineWidth() above 1.0, and this plot is almost
// entirely lines, arrowheads and text. Canvas2D gives real stroke widths, cheap
// text for the axes, and none of it is hot enough to need the GPU.
//
// Geometry arrives from WASM in PHYSICAL coordinates with two pre-normalized
// scalars per point, exactly as the desktop consumes it:
//   inPlane -> arrow / head LENGTH  (readability)
//   total   -> COLOUR on the global peak scale (matches the 3D view and the bar)

const STOPS = [
  [0.20, 0.45, 0.95], [0.10, 0.72, 0.85], [0.20, 0.72, 0.25],
  [0.98, 0.70, 0.10], [0.92, 0.25, 0.10], [0.55, 0.00, 0.08],
];

export function fireColor(t) {
  t = Math.max(0, Math.min(1, t)) * 5;
  const i = Math.min(4, Math.floor(t)), l = t - i;
  const a = STOPS[i], b = STOPS[i + 1];
  const r = Math.round(255 * (a[0] + l * (b[0] - a[0])));
  const g = Math.round(255 * (a[1] + l * (b[1] - a[1])));
  const c = Math.round(255 * (a[2] + l * (b[2] - a[2])));
  return `rgb(${r},${g},${c})`;
}

const M = { l: 58, r: 96, t: 14, b: 40 };   // right margin holds the scale bar

// Plane coordinates -> model coordinates, mirroring FieldViz.cpp's toModel().
// The highlight functions in theory.js are written against x/y/z, so this is
// where the cut's own (u, v) get translated back.
function toModel(plane, u, v, slice) {
  if (plane === 0) return { x: u, y: v, z: slice };
  if (plane === 1) return { x: v, y: slice, z: u };
  return { x: slice, y: v, z: u };
}

function niceStep(span, want) {
  if (span <= 0) return 1;
  const raw = span / Math.max(1, want);
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const n = raw / mag;
  return (n < 1.5 ? 1 : n < 3.5 ? 2 : n < 7.5 ? 5 : 10) * mag;
}

export function drawSection(canvas, info, data, opts) {
  const dpr = Math.min(2, window.devicePixelRatio || 1);
  const W = Math.max(1, Math.round(canvas.clientWidth * dpr));
  const H = Math.max(1, Math.round(canvas.clientHeight * dpr));
  if (canvas.width !== W || canvas.height !== H) { canvas.width = W; canvas.height = H; }
  const g = canvas.getContext('2d');
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  const w = W / dpr, h = H / dpr;

  g.fillStyle = '#12131a';
  g.fillRect(0, 0, w, h);

  const x0 = M.l, y0 = M.t, x1 = w - M.r, y1 = h - M.b;
  const pw = Math.max(1, x1 - x0), ph = Math.max(1, y1 - y0);
  const uSpan = info.uMax - info.uMin, vSpan = info.vMax - info.vMin;
  const sx = u => x0 + (u - info.uMin) / uSpan * pw;
  const sy = v => y0 + (1 - (v - info.vMin) / vSpan) * ph;

  // Outline: a circle when looking down a cylinder, a rectangle otherwise.
  g.strokeStyle = '#8b90a0'; g.lineWidth = 1;
  if (opts.geometry === 1 && opts.plane === 0) {
    g.beginPath();
    g.arc((x0 + x1) / 2, (y0 + y1) / 2, Math.min(pw, ph) / 2, 0, 2 * Math.PI);
    g.stroke();
  } else {
    g.strokeRect(x0, y0, pw, ph);
  }

  // ---- axes, in millimetres, with the coordinate origin marked ----
  g.fillStyle = '#c9cdd8'; g.font = '11px system-ui, sans-serif';
  g.textAlign = 'center'; g.textBaseline = 'top';
  const du = niceStep((info.uMax - info.uMin) * 1e3, 6);
  for (let t = Math.ceil(info.uMin * 1e3 / du) * du; t <= info.uMax * 1e3 + 1e-9; t += du) {
    const px = sx(t / 1e3), origin = Math.abs(t) < 1e-9;
    g.strokeStyle = origin ? '#eef' : '#8b90a0';
    g.beginPath(); g.moveTo(px, y1); g.lineTo(px, y1 + (origin ? 7 : 4)); g.stroke();
    g.fillStyle = origin ? '#fff' : '#c9cdd8';
    g.fillText(String(+t.toFixed(3)), px, y1 + 9);
  }
  g.textAlign = 'right'; g.textBaseline = 'middle';
  const dv = niceStep((info.vMax - info.vMin) * 1e3, 4);
  for (let t = Math.ceil(info.vMin * 1e3 / dv) * dv; t <= info.vMax * 1e3 + 1e-9; t += dv) {
    const py = sy(t / 1e3), origin = Math.abs(t) < 1e-9;
    g.strokeStyle = origin ? '#eef' : '#8b90a0';
    g.beginPath(); g.moveTo(x0 - (origin ? 7 : 4), py); g.lineTo(x0, py); g.stroke();
    g.fillStyle = origin ? '#fff' : '#c9cdd8';
    g.fillText(String(+t.toFixed(3)), x0 - 10, py);
  }
  g.fillStyle = '#c9cdd8'; g.textAlign = 'center'; g.textBaseline = 'top';
  g.fillText(opts.uLabel, (x0 + x1) / 2, y1 + 22);
  g.textAlign = 'left'; g.fillText(opts.vLabel, 4, 2);

  // ---- field ----
  // With a term hovered in the theory column, everything the factor does not
  // control is dimmed. Showing the factor as a separate heat map would compete
  // with the field's own colours; dimming keeps one image and one scale.
  const hl = opts.highlight;
  const weight = (u, v) => {
    if (!hl) return 1;
    const w = hl(toModel(opts.plane, u, v, info.slice));
    return 0.10 + 0.90 * Math.max(0, Math.min(1, w));
  };

  if (data.kind === 'lines') {
    const { verts, starts } = data;
    for (let s = 0; s + 1 < starts.length; ++s) {
      const a = starts[s], b = starts[s + 1];
      if (b - a < 2) continue;
      let since = 21;
      for (let i = a; i + 1 < b; ++i) {
        const ux = sx(verts[i*4]),     uy = sy(verts[i*4+1]);
        const vx = sx(verts[(i+1)*4]), vy = sy(verts[(i+1)*4+1]);
        g.globalAlpha = weight(verts[i*4], verts[i*4+1]);
        g.strokeStyle = fireColor(verts[i*4+3]);
        g.lineWidth = 1.5;
        g.beginPath(); g.moveTo(ux, uy); g.lineTo(vx, vy); g.stroke();

        const dx = vx - ux, dy = vy - uy, seg = Math.hypot(dx, dy);
        since += seg;
        // Heads carry the SIGN, which dashes alone cannot: the field reverses
        // every half cycle and a plain line looks identical either way.
        if (since >= 42 && seg > 1e-3) {
          since = 0;
          const nx = dx / seg, ny = dy / seg;
          // headLen, not hl: `hl` is the highlight function in the enclosing
          // scope and shadowing it here worked only by accident of block scope.
          const headLen = 8 * (0.30 + 0.70 * Math.sqrt(Math.max(0, Math.min(1, verts[i*4+2]))));
          g.beginPath();
          g.moveTo(vx, vy); g.lineTo(vx - nx*headLen - ny*headLen*0.55, vy - ny*headLen + nx*headLen*0.55);
          g.moveTo(vx, vy); g.lineTo(vx - nx*headLen + ny*headLen*0.55, vy - ny*headLen - nx*headLen*0.55);
          g.stroke();
        }
      }
    }
  } else {
    const a = data.arrows;
    const cell = Math.min(pw / opts.nu, ph / opts.nv);
    const amax = 0.9 * cell;
    for (let i = 0; i < a.length; i += 6) {
      const cx = sx(a[i]), cy = sy(a[i + 1]);
      const len = Math.sqrt(a[i + 4]) * amax;
      const dx = a[i + 2] * len, dy = -a[i + 3] * len;   // canvas y grows down
      g.globalAlpha = weight(a[i], a[i + 1]);
      g.strokeStyle = fireColor(a[i + 5]);
      g.lineWidth = 1.4;
      const ax = cx - dx / 2, ay = cy - dy / 2, bx = cx + dx / 2, by = cy + dy / 2;
      g.beginPath(); g.moveTo(ax, ay); g.lineTo(bx, by);
      const l = Math.hypot(dx, dy);
      if (l > 1e-3) {
        const nx = dx / l, ny = dy / l, hh = 3.2;
        g.moveTo(bx, by); g.lineTo(bx - nx*hh - ny*hh*0.6, by - ny*hh + nx*hh*0.6);
        g.moveTo(bx, by); g.lineTo(bx - nx*hh + ny*hh*0.6, by - ny*hh - nx*hh*0.6);
      }
      g.stroke();
    }
  }

  g.globalAlpha = 1;

  // ---- colour scale, beside the plot ----
  const bx = x1 + 18, bw = 16, by0 = y0 + 4, by1 = y1 - 4;
  for (let y = by0; y < by1; ++y) {
    g.fillStyle = fireColor(1 - (y - by0) / (by1 - by0));
    g.fillRect(bx, y, bw, 1);
  }
  g.strokeStyle = '#8b90a0'; g.lineWidth = 1;
  g.strokeRect(bx, by0, bw, by1 - by0);
  g.fillStyle = '#c9cdd8'; g.textAlign = 'left'; g.textBaseline = 'middle';
  g.font = '10px system-ui, sans-serif';
  for (let k = 0; k <= 4; ++k) {
    const f = k / 4, py = by1 - f * (by1 - by0);
    g.fillText(fmt(opts.peak * f), bx + bw + 4, py);
  }
  g.textBaseline = 'bottom';
  g.fillText(opts.unit, bx, by0 - 3);
}

function fmt(v) {
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e4 || a < 1e-2) return v.toExponential(1);
  return String(+v.toPrecision(3));
}
