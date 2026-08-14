// Wiring: read the controls, hand them to the WASM module, hand its buffers to
// the two renderers. No physics here — that all lives in the same C++ the
// desktop build uses, so the two cannot disagree.

import createWaveguide from './dist/waveguide.js';
import { Renderer3D } from './renderer3d.js';
import { drawSection } from './section2d.js';
import { buildTheory, factorFor, termNote } from './theory.js';

const $ = id => document.getElementById(id);
const statusEl = $('status');

const NUM = ['widthMM','heightMM','radiusMM','depthMM','freqGHz','epsR','muR','powerW'];
const INT = ['geometry','structure','modeType','field','modeM','modeN','modeL'];

let Module, guide, r3d;
let phase = 0, lastT = 0;
let highlight = null;      // active factor from the theory column, or null
let dirty = true;          // geometry needs rebuilding from WASM
let cache = null;

function readConfig() {
  const c = {};
  for (const k of NUM) c[k] = parseFloat($(k).value);
  for (const k of INT) c[k] = parseInt($(k).value, 10);
  return c;
}

// The same validity rules the desktop enforces, surfaced BEFORE the module
// silently substitutes a mode: TE allows one zero index but not both, TM needs
// both >= 1, and the circular radial index counts Bessel roots from 1.
function modeNote(c) {
  if (c.geometry === 1) {
    if (c.modeM < 1) return 'Guia circular: o indice radial m conta raizes de Bessel e comeca em 1 — usando m = 1.';
    return '';
  }
  if (c.modeType === 1 && (c.modeM < 1 || c.modeN < 1))
    return 'TM exige M>=1 e N>=1 (Ez ~ sin·sin se anula) — usando TM'
           + Math.max(1, c.modeM) + Math.max(1, c.modeN) + '.';
  if (c.modeType === 0 && c.modeM === 0 && c.modeN === 0)
    return 'TE00 nao existe (kc = 0, todo o campo se anula) — usando TE10.';
  return '';
}

function rebuild() {
  const c = readConfig();
  guide.configure(c);
  $('modeNote').textContent = modeNote(c);

  const dom = guide.domain();
  r3d.setOutline(dom, c.geometry);

  const pts = guide.cloud(parseInt($('points').value, 10), 8, $('cutaway').checked);
  r3d.setCloud(pts, guide.cloudInvPeak(), guide.cloudMeanSpacing());

  const plane = parseInt($('plane').value, 10);
  const info = guide.sectionInfo(plane, -1);

  const phys = guide.physicalUnits();
  const unit = phys ? (c.field === 0 ? 'V/m' : 'A/m') : 'u.a.';
  cache = { c, dom, plane, info, phys, unit, peak: guide.peakField(), nPts: pts.length / 6 };

  renderTheory(c);
  readout(cache);
  dirty = false;
}

// Rebuild the equation column. Re-rendered on every parameter change because
// the formulas themselves change with TE/TM, E/H and the geometry -- this is
// the panel, not just the numbers in it.
function renderTheory(cfg) {
  const kc = 2 * Math.PI * guide.cutoffHz() /
             (299792458 / Math.sqrt(cfg.epsR * cfg.muR));
  const v = 299792458 / Math.sqrt(cfg.epsR * cfg.muR);
  const k = 2 * Math.PI * cfg.freqGHz * 1e9 / v;
  const disc = k * k - kc * kc;
  const { ctx, eqs } = buildTheory(cfg, {
    kc, alpha: disc < 0 ? Math.sqrt(-disc) : 0,
  });

  const kind = cfg.modeType === 0 ? 'TE' : 'TM';
  const idx = cfg.geometry === 0
    ? `${ctx.m}${ctx.n}${cfg.structure === 1 ? cfg.modeL : ''}`
    : `${ctx.n}${Math.max(1, cfg.modeM)}`;
  $('theoryTitle').innerHTML =
    `${kind}<sub>${idx}</sub> &middot; campo ${cfg.field === 0 ? 'elétrico' : 'magnético'}`;

  $('equations').innerHTML = eqs.map(e =>
    `<div class="eq"><span class="lhs">${e.lhs} =</span>`
    + `<math display="inline">${e.mathml}</math></div>`).join('');

  // Hover a factor -> dim everything it does not control.
  for (const el of $('equations').querySelectorAll('.term')) {
    const id = el.dataset.term;
    const enter = () => {
      highlight = factorFor(ctx, id);
      $('termNote').textContent = termNote(id);
      $('termNote').classList.add('active');
      el.classList.add('on');
    };
    const leave = () => {
      highlight = null;
      $('termNote').textContent = 'Passe o mouse sobre um fator para ver que regiao do corte ele controla.';
      $('termNote').classList.remove('active');
      el.classList.remove('on');
    };
    el.addEventListener('pointerenter', enter);
    el.addEventListener('pointerleave', leave);
    // Touch has no hover: tapping toggles, which is how a phone gets the
    // same affordance instead of silently losing the feature.
    el.addEventListener('click', () => (highlight ? leave() : enter()));
  }
}

function readout(s) {
  const fc = guide.cutoffHz() / 1e9;
  const beta = guide.beta();
  const prop = guide.propagating();
  const lg = prop ? 2 * Math.PI / beta * 1e3 : 0;
  const rows = [
    ['Pico', `${fmt(s.peak)} ${s.unit}`],
    [s.c.structure === 1 ? 'f_res' : 'f_c', `${fc.toFixed(4)} GHz`],
    ['Regime', prop ? 'propagante' : 'evanescente'],
  ];
  if (prop) {
    rows.push(['β', `${beta.toFixed(1)} rad/m`]);
    rows.push(['λ_g', `${lg.toFixed(2)} mm`]);
  }
  if (s.phys) rows.push(['Amplitude A', fmt(guide.amplitude())]);
  rows.push(['Pontos', String(s.nPts)]);

  const t = $('readout');
  t.innerHTML = '';
  for (const [k, v] of rows) {
    const tr = document.createElement('tr');
    if (k === 'Regime' && !prop) tr.className = 'warn';
    tr.innerHTML = `<td>${k}</td><td>${v}</td>`;
    t.appendChild(tr);
  }
  if (!s.phys) {
    const tr = document.createElement('tr');
    tr.className = 'warn';
    tr.innerHTML = `<td colspan="2">Amplitude arbitraria: cavidade e modo evanescente nao transportam potencia.</td>`;
    t.appendChild(tr);
  }
}

const fmt = v => (v === 0 ? '0'
  : Math.abs(v) >= 1e4 || Math.abs(v) < 1e-2 ? v.toExponential(3)
  : String(+v.toPrecision(5)));

const PLANE_LABELS = [
  { u: 'x (mm)', v: 'y (mm)', name: 'XY (transversal)' },
  { u: 'z (mm)', v: 'x (mm)', name: 'ZX' },
  { u: 'z (mm)', v: 'y (mm)', name: 'ZY' },
];

function frame(t) {
  requestAnimationFrame(frame);
  const dt = lastT ? (t - lastT) / 1000 : 0;
  lastT = t;
  if ($('animate').checked) phase = (phase + dt * 2.5) % (2 * Math.PI);
  if (dirty) { try { rebuild(); } catch (e) { statusEl.textContent = 'erro: ' + e.message; dirty = false; } }
  if (!cache) return;

  r3d.frame(phase);

  const lab = PLANE_LABELS[cache.plane];
  const stream = $('streamlines').checked;
  const data = stream
    ? Object.assign({ kind: 'lines' }, guide.sectionLines(cache.plane, phase, 14, 9))
    : { kind: 'arrows', arrows: guide.sectionArrows(cache.plane, phase) };

  drawSection($('sec'), cache.info, data, {
    geometry: cache.c.geometry, plane: cache.plane,
    uLabel: lab.u, vLabel: lab.v,
    nu: cache.plane === 0 ? 18 : 44, nv: cache.plane === 0 ? 18 : 14,
    peak: cache.peak, unit: cache.unit,
    highlight,
  });

  const fr = (cache.info.wMax > cache.info.wMin)
    ? (cache.info.slice - cache.info.wMin) / (cache.info.wMax - cache.info.wMin) : 0.5;
  $('secCaption').textContent =
    `${lab.name} — corte em ${(fr * 100).toFixed(0)} % do eixo perpendicular · `
    + `cor = |campo| total, seta = componente no plano`;
}

function markDirty() { dirty = true; }

async function main() {
  try {
    Module = await createWaveguide();
  } catch (e) {
    statusEl.textContent = 'falha ao carregar o WASM: ' + e.message
      + ' — sirva a pasta por HTTP (file:// nao funciona).';
    return;
  }
  guide = new Module.Guide();
  try {
    r3d = new Renderer3D($('gl'));
  } catch (e) {
    statusEl.textContent = e.message;
    return;
  }

  for (const k of [...NUM, ...INT, 'points', 'cutaway', 'plane'])
    $(k).addEventListener('input', markDirty);

  // Geometry/structure switches change which controls are meaningful.
  const sync = () => {
    const g = parseInt($('geometry').value, 10);
    document.querySelectorAll('[data-geom]').forEach(
      el => { el.hidden = parseInt(el.dataset.geom, 10) !== g; });
    document.querySelectorAll('[data-cavity]').forEach(
      el => { el.hidden = $('structure').value !== '1'; });
  };
  $('geometry').addEventListener('input', sync);
  $('structure').addEventListener('input', sync);
  sync();

  statusEl.textContent = 'modulo carregado — a fisica e o mesmo C++ do aplicativo desktop';
  requestAnimationFrame(frame);
}

main();
