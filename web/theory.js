// The equations of the current mode, and the link between a term and the
// region of the plot it controls.
//
// These are transcribed from TEmnModel.cpp / CylindricalModel.cpp, so the page
// shows the same expressions the solver actually evaluates. If the two ever
// disagree, the page is wrong -- treat it as documentation OF the code, not a
// second independent derivation. (One of these very factors, Ex carrying
// sin(n pi y/b) rather than cos, was wrong in the C++ for a long time: the
// tangential E has to vanish on the y walls, and cos puts a maximum there.)
//
// Rendering uses native MathML. No KaTeX, no MathJax, no build step -- and no
// 300 kB of font payload on a student's phone.

const mi = s => `<mi>${s}</mi>`;
const mn = s => `<mn>${s}</mn>`;
const mo = s => `<mo>${s}</mo>`;
const row = (...x) => `<mrow>${x.join('')}</mrow>`;
const frac = (a, b) => `<mfrac>${a}${b}</mfrac>`;
const sub = (a, b) => `<msub>${a}${b}</msub>`;
const sup = (a, b) => `<msup>${a}${b}</msup>`;
// A hoverable factor: the id ties it to a highlight function below.
const term = (id, inner) => `<mrow class="term" data-term="${id}">${inner}</mrow>`;

const PI = mi('&#960;');
const kc2 = sup(sub(mi('k'), mi('c')), mn('2'));

// sin(m pi x / a) and friends, as hoverable groups.
function trig(fn, idxSym, varSym, denSym, id) {
  return term(id, row(
    mi(fn), mo('&#8289;'), mo('('),
    frac(row(mi(idxSym), PI, mi(varSym)), mi(denSym)),
    mo(')')));
}

// ---------------------------------------------------------------------------
// Spatial factors. Each returns |factor| in [0,1] at a MODEL point, so the plot
// can dim everything the hovered term does not control.
// ---------------------------------------------------------------------------
const F = {
  sin_mx: (c, x) => Math.abs(Math.sin(c.m * Math.PI * x.x / c.a)),
  cos_mx: (c, x) => Math.abs(Math.cos(c.m * Math.PI * x.x / c.a)),
  sin_ny: (c, x) => Math.abs(Math.sin(c.n * Math.PI * x.y / c.b)),
  cos_ny: (c, x) => Math.abs(Math.cos(c.n * Math.PI * x.y / c.b)),
  // e^{-j beta z}: the travelling factor. Its MAGNITUDE is flat above cutoff,
  // so highlighting it would light the whole plot; below cutoff it decays and
  // that is exactly what is worth seeing.
  expz: (c, x) => (c.alpha > 0 ? Math.exp(-c.alpha * x.z) : 1),
  Jn:   (c, x) => {
    const rho = Math.hypot(x.x, x.y);
    return Math.min(1, Math.abs(besselJ(c.n, c.kc * rho)) / (c.jnMax || 1));
  },
};

// Same series as Bessel.hpp -- see there for why std::cyl_bessel_j is unusable.
export function besselJ(n, x) {
  if (x < 0) return (n & 1 ? -1 : 1) * besselJ(n, -x);
  const h = 0.5 * x;
  let t = 1;
  for (let i = 1; i <= n; ++i) t *= h / i;
  let s = t;
  const h2 = h * h;
  for (let k = 1; k < 80; ++k) {
    t *= -h2 / (k * (n + k));
    s += t;
    if (Math.abs(t) < 1e-18 * Math.abs(s)) break;
  }
  return s;
}

const NOTE = {
  sin_mx: 'Zero nas paredes x = 0 e x = a. E tangencial tem de se anular no condutor.',
  cos_mx: 'Maximo nas paredes x = 0 e x = a, zero no meio.',
  sin_ny: 'Zero nas paredes y = 0 e y = b.',
  cos_ny: 'Maximo nas paredes y = 0 e y = b, zero no meio.',
  expz:   'Propagacao ao longo de z. Acima do corte a amplitude nao cai; abaixo, decai como e^(-alpha z).',
  Jn:     'Perfil radial. O corte sai da raiz de Bessel: e ela que fixa k_c = p/R.',
};

// ---------------------------------------------------------------------------
// Formula sets. Component -> { mathml, factors[] }
// ---------------------------------------------------------------------------
function rectFormulas(c) {
  const A = mi('A');
  const j = mi('j');
  const beta = mi('&#946;');
  const omega = mi('&#969;');
  const mu = mi('&#956;');
  const eps = mi('&#949;');
  const mpi = frac(row(mi('m'), PI), mi('a'));
  const npi = frac(row(mi('n'), PI), mi('b'));
  const ez = term('expz', sup(mi('e'), row(mo('&#8722;'), j, beta, mi('z'))));

  const S_mx = trig('sin', 'm', 'x', 'a', 'sin_mx');
  const C_mx = trig('cos', 'm', 'x', 'a', 'cos_mx');
  const S_ny = trig('sin', 'n', 'y', 'b', 'sin_ny');
  const C_ny = trig('cos', 'n', 'y', 'b', 'cos_ny');

  const eq = (lhs, body) => ({ lhs, mathml: row(...body) });

  if (c.modeType === 0 && c.field === 0)       // TE, electric
    return [
      eq('E_x', [j, omega, mu, frac(npi, kc2), A, C_mx, S_ny, ez]),
      eq('E_y', [mo('&#8722;'), j, omega, mu, frac(mpi, kc2), A, S_mx, C_ny, ez]),
      eq('E_z', [mn('0')]),
    ];
  if (c.modeType === 0 && c.field === 1)       // TE, magnetic
    return [
      eq('H_x', [j, beta, frac(mpi, kc2), A, S_mx, C_ny, ez]),
      eq('H_y', [j, beta, frac(npi, kc2), A, C_mx, S_ny, ez]),
      eq('H_z', [A, C_mx, C_ny, ez]),
    ];
  if (c.modeType === 1 && c.field === 0)       // TM, electric
    return [
      eq('E_x', [mo('&#8722;'), j, beta, frac(mpi, kc2), A, C_mx, S_ny, ez]),
      eq('E_y', [mo('&#8722;'), j, beta, frac(npi, kc2), A, S_mx, C_ny, ez]),
      eq('E_z', [A, S_mx, S_ny, ez]),
    ];
  return [                                      // TM, magnetic
    eq('H_x', [j, omega, eps, frac(npi, kc2), A, S_mx, C_ny, ez]),
    eq('H_y', [mo('&#8722;'), j, omega, eps, frac(mpi, kc2), A, C_mx, S_ny, ez]),
    eq('H_z', [mn('0')]),
  ];
}

function cylFormulas(c) {
  const A = mi('A'), j = mi('j'), beta = mi('&#946;');
  const rho = mi('&#961;'), phi = mi('&#966;');
  const Jn = term('Jn', row(sub(mi('J'), mi('n')), mo('('),
                            sub(mi('k'), mi('c')), rho, mo(')')));
  const Jnp = term('Jn', row(sup(sub(mi('J'), mi('n')), mo('&#8242;')), mo('('),
                             sub(mi('k'), mi('c')), rho, mo(')')));
  const ang = row(mi('cos'), mo('&#8289;'), mo('('), mi('n'), phi, mo(')'));
  const ez = term('expz', sup(mi('e'), row(mo('&#8722;'), j, beta, mi('z'))));
  const eq = (lhs, body) => ({ lhs, mathml: row(...body) });

  if (c.field === 1 && c.modeType === 0)
    return [
      eq('H_z', [A, Jn, ang, ez]),
      eq('H_&#961;', [mo('&#8722;'), j, frac(beta, sub(mi('k'), mi('c'))), A, Jnp, ang, ez]),
      eq('H_&#966;', [mo('&#8722;'), j, frac(row(beta, mi('n')),
        row(kc2, rho)), A, Jn, mi('sin'), mo('('), mi('n'), phi, mo(')'), ez]),
    ];
  return [
    eq('E_&#961;', [mo('&#8722;'), j, frac(mi('&#969;&#956;'), sub(mi('k'), mi('c'))),
                    A, Jnp, ang, ez]),
    eq('E_&#966;', [j, frac(mi('&#969;&#956;'), kc2), A, Jn, ang, ez]),
  ];
}

// ---------------------------------------------------------------------------
export function buildTheory(cfg, derived) {
  const c = {
    m: Math.max(cfg.geometry === 0 && cfg.modeType === 0 && cfg.modeM === 0 && cfg.modeN === 0
        ? 1 : cfg.modeM, cfg.modeType === 1 ? 1 : 0),
    n: Math.max(cfg.modeN, cfg.modeType === 1 ? 1 : 0),
    a: cfg.widthMM / 1000, b: cfg.heightMM / 1000,
    kc: derived.kc, alpha: derived.alpha,
    modeType: cfg.modeType, field: cfg.field, geometry: cfg.geometry,
  };
  // Normalizer so the Bessel highlight spans [0,1] over the disc.
  if (cfg.geometry === 1) {
    let mx = 1e-9;
    for (let i = 0; i <= 64; ++i)
      mx = Math.max(mx, Math.abs(besselJ(c.n, c.kc * (cfg.radiusMM / 1000) * i / 64)));
    c.jnMax = mx;
  }
  const eqs = cfg.geometry === 0 ? rectFormulas(c) : cylFormulas(c);
  return { ctx: c, eqs };
}

export function factorFor(ctx, id) {
  const f = F[id];
  if (!f) return null;
  return p => f(ctx, p);
}

export const termNote = id => NOTE[id] || '';
