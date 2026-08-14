// Acceptance test for the WebAssembly build: the numbers it produces must match
// the desktop binary. Anything else means the browser would be showing different
// physics from the app the results were validated with.
//
// Expected values come from the native build (the FieldViz probe) and, for the
// TE11 case, from an independent calculation done while comparing against Ansys
// HFSS -- so this pins the whole chain, not just C++-to-WASM.
//
//   node web/check.mjs
import createWaveguide from './dist/waveguide.js';

const CASES = [
  { name: 'TE10 WR-90 12 GHz E (1 W)',
    cfg: { geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83, depthMM:300,
           freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1, modeType:0,
           field:0, structure:0, powerW:1 },
    peak: 2783.34, amp: 4.037 },

  { name: 'TE11 WR-90 17 GHz H (10 W)',
    cfg: { geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83, depthMM:50,
           freqGHz:17, epsR:1, muR:1, modeM:1, modeN:1, modeL:1, modeType:0,
           field:1, structure:0, powerW:10 },
    peak: 51.3186, amp: 51.32 },

  { name: 'Cavidade TE101 E',
    cfg: { geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83, depthMM:300,
           freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1, modeType:0,
           field:0, structure:1, powerW:1 },
    peak: 689.441 },

  { name: 'Cilindrica TE01 (n=0) H',
    cfg: { geometry:1, widthMM:22.86, heightMM:10.16, radiusMM:23.83, depthMM:300,
           freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1, modeType:0,
           field:1, structure:0, powerW:1 },
    peak: 1.0 },
];

const rel = (a, b) => Math.abs(a - b) / Math.max(1e-30, Math.abs(b));

createWaveguide().then(Module => {
  const g = new Module.Guide();
  let fail = 0;

  for (const c of CASES) {
    g.configure(c.cfg);
    console.log(`\n=== ${c.name} ===`);

    const pk = g.peakField();
    const ePk = rel(pk, c.peak);
    console.log(`  peakField = ${pk.toPrecision(9)}  (desktop ${c.peak})`
                + `  err ${(ePk * 100).toExponential(2)} %${ePk > 1e-4 ? '  <== DIVERGE' : ''}`);
    if (ePk > 1e-4) fail++;

    if (c.amp !== undefined) {
      const am = g.amplitude(), eAm = rel(am, c.amp);
      console.log(`  amplitude = ${am.toPrecision(9)}  (desktop ${c.amp})`
                  + `  err ${(eAm * 100).toExponential(2)} %${eAm > 1e-3 ? '  <== DIVERGE' : ''}`);
      if (eAm > 1e-3) fail++;
    }
    console.log(`  fisico=${g.physicalUnits()}  propaga=${g.propagating()}`
                + `  beta=${g.beta().toFixed(2)}  fc=${(g.cutoffHz() / 1e9).toFixed(4)} GHz`);

    // No buffer may come back empty: that is how every silent failure in this
    // project has presented itself so far.
    const cloud = g.cloud(20000, 8, false);
    const nPts = cloud.length / 6;
    console.log(`  nuvem: ${nPts} pontos  invPeak ${g.cloudInvPeak().toExponential(4)}`
                + (nPts === 0 ? '  <== VAZIA' : ''));
    if (nPts === 0) fail++;

    for (let pl = 0; pl < 3; ++pl) {
      const info = g.sectionInfo(pl, -1);
      const ar = g.sectionArrows(pl, 0.7);
      const ln = g.sectionLines(pl, 0.7, 12, 8);
      console.log(`  plano ${pl}: slice ${info.slice.toExponential(3)}`
                  + `  ref ${info.reference.toPrecision(4)}`
                  + `  | setas ${ar.length / 6}`
                  + `  linhas ${ln.starts.length - 1} (${ln.verts.length / 4} vert)`);
    }
  }

  g.delete();
  console.log(fail === 0 ? '\nOK: WASM bate com o desktop'
                         : `\n${fail} DIVERGENCIA(S)`);
  process.exit(fail === 0 ? 0 : 1);
}).catch(e => { console.error('falha ao carregar o modulo:', e); process.exit(2); });
