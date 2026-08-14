// Exercises the JS<->WASM boundary in the exact order app.js uses it, and
// checks the SHAPES the renderers rely on. Syntax checks and the numeric
// acceptance test both pass without touching this: a stride mistake or a
// renamed field would sail through them and only show up as a blank canvas.
//
//   node web/integration.mjs
import createWaveguide from './dist/waveguide.js';

let fail = 0;
const ok = (cond, msg) => {
  console.log(`  ${cond ? 'ok  ' : 'FALHA'}  ${msg}`);
  if (!cond) fail++;
};

const Module = await createWaveguide();
const guide = new Module.Guide();

for (const [name, cfg] of [
  ['retangular TE10', { geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83,
      depthMM:300, freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1,
      modeType:0, field:0, structure:0, powerW:1 }],
  ['cilindrica TE01', { geometry:1, widthMM:22.86, heightMM:10.16, radiusMM:23.83,
      depthMM:300, freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1,
      modeType:0, field:1, structure:0, powerW:1 }],
]) {
  console.log(`\n=== ${name} ===`);
  guide.configure(cfg);

  const dom = guide.domain();
  ok(['x0','x1','y0','y1','z0','z1'].every(k => typeof dom[k] === 'number'),
     'domain() traz as seis bordas');
  ok(dom.x1 > dom.x0 && dom.z1 > dom.z0, 'dominio nao degenerado');
  // The cylinder is centred on its axis; the rectangle starts at the corner.
  ok(cfg.geometry === 1 ? dom.x0 < 0 : dom.x0 === 0,
     `referencial ${cfg.geometry === 1 ? 'centrado' : 'no canto'}`);

  const cloud = guide.cloud(20000, 8, false);
  ok(cloud instanceof Float32Array, 'cloud() devolve Float32Array');
  ok(cloud.length % 6 === 0 && cloud.length > 0,
     `stride 6 respeitado (${cloud.length / 6} pontos)`);
  ok(Number.isFinite(guide.cloudInvPeak()) && guide.cloudInvPeak() > 0,
     `invPeak = ${guide.cloudInvPeak().toExponential(3)}`);
  ok(guide.cloudMeanSpacing() > 0,
     `espacamento = ${guide.cloudMeanSpacing().toExponential(3)} m`);
  // Every point must sit inside the centred bounding box the outline draws.
  let outside = 0;
  const hx = (dom.x1 - dom.x0) / 2, hy = (dom.y1 - dom.y0) / 2, hz = (dom.z1 - dom.z0) / 2;
  for (let i = 0; i < cloud.length; i += 6)
    if (Math.abs(cloud[i]) > hx + 1e-9 || Math.abs(cloud[i+1]) > hy + 1e-9 ||
        Math.abs(cloud[i+2]) > hz + 1e-9) outside++;
  ok(outside === 0, `nuvem dentro do contorno (${outside} fora)`);

  for (let pl = 0; pl < 3; ++pl) {
    const info = guide.sectionInfo(pl, -1);
    ok(info.uMax > info.uMin && info.vMax > info.vMin, `plano ${pl}: span valido`);
    ok(info.slice >= info.wMin - 1e-12 && info.slice <= info.wMax + 1e-12,
       `plano ${pl}: slice dentro do eixo perpendicular`);

    const a = guide.sectionArrows(pl, 0.7);
    ok(a instanceof Float32Array && a.length % 6 === 0,
       `plano ${pl}: setas stride 6 (${a.length / 6})`);
    let bad = 0;
    for (let i = 0; i < a.length; i += 6) {
      const dirLen = Math.hypot(a[i+2], a[i+3]);
      if (Math.abs(dirLen - 1) > 1e-3) bad++;                 // direcao unitaria
      if (a[i+4] < 0 || a[i+4] > 1.001) bad++;                // inPlane em [0,1]
      if (a[i+5] < 0 || a[i+5] > 1.001) bad++;                // total em [0,1]
    }
    ok(bad === 0, `plano ${pl}: direcoes unitarias e escalares normalizados`);

    const ln = guide.sectionLines(pl, 0.7, 14, 9);
    ok(ln.verts instanceof Float32Array && ln.starts instanceof Float32Array,
       `plano ${pl}: linhas trazem verts+starts`);
    ok(ln.verts.length % 4 === 0, `plano ${pl}: linhas stride 4`);
    const n = ln.starts.length - 1;
    ok(n >= 0 && ln.starts[ln.starts.length - 1] === ln.verts.length / 4,
       `plano ${pl}: sentinela final fecha (${n} linhas)`);
    let mono = true;
    for (let i = 1; i < ln.starts.length; ++i)
      if (ln.starts[i] < ln.starts[i-1]) mono = false;
    ok(mono, `plano ${pl}: offsets crescentes`);
  }
}

// Reconfiguring must not leak the previous mode's cached span/reference.
console.log('\n=== troca de modo ===');
guide.configure({ geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83,
  depthMM:300, freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1,
  modeType:0, field:0, structure:0, powerW:1 });
const before = guide.sectionInfo(0, -1).reference;
guide.configure({ geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83,
  depthMM:300, freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1,
  modeType:0, field:0, structure:0, powerW:100 });
// Fields scale as sqrt(power): 100x the power is 10x the field.
const after = guide.sectionInfo(0, -1).reference;
ok(Math.abs(after / before - 10) < 1e-3,
   `referencia acompanha sqrt(P): ${before.toPrecision(5)} -> ${after.toPrecision(5)} (x${(after/before).toFixed(3)})`);

// Asking for section geometry WITHOUT sectionInfo first must still work.
const fresh = new Module.Guide();
fresh.configure({ geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83,
  depthMM:300, freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1,
  modeType:0, field:0, structure:0, powerW:1 });
ok(fresh.sectionArrows(0, 0.7).length > 0, 'sectionArrows sem sectionInfo antes');
fresh.delete();

guide.delete();
console.log(fail === 0 ? '\nOK: fronteira JS<->WASM consistente' : `\n${fail} FALHA(S)`);
process.exit(fail === 0 ? 0 : 1);
