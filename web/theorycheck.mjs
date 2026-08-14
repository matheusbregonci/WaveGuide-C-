// Checks the theory column without a browser: that every mode produces
// well-formed MathML with hoverable terms, and that a highlight factor really
// peaks where the physics says it should.
//
//   node web/theorycheck.mjs
import { buildTheory, factorFor, besselJ, } from './theory.js';

let fail = 0;
const ok = (c, m) => { console.log(`  ${c ? 'ok  ' : 'FALHA'}  ${m}`); if (!c) fail++; };

const BASE = { geometry:0, widthMM:22.86, heightMM:10.16, radiusMM:23.83,
  depthMM:300, freqGHz:12, epsR:1, muR:1, modeM:1, modeN:0, modeL:1,
  modeType:0, field:0, structure:0, powerW:1 };

const countTag = (s, tag) => (s.split('<' + tag).length - 1);

for (const [mt, fl, name] of [[0,0,'TE-E'], [0,1,'TE-H'], [1,0,'TM-E'], [1,1,'TM-H']]) {
  const cfg = { ...BASE, modeType: mt, field: fl, modeN: mt === 1 ? 1 : 0 };
  const { ctx, eqs } = buildTheory(cfg, { kc: 137.44, alpha: 0 });
  console.log(`\n=== ${name} ===`);
  ok(eqs.length === 3, `${eqs.length} componentes`);
  let terms = 0, balanced = true;
  for (const e of eqs) {
    terms += (e.mathml.split('data-term=').length - 1);
    if (countTag(e.mathml, 'mrow') !== (e.mathml.split('</mrow>').length - 1)) balanced = false;
    if (countTag(e.mathml, 'mfrac') !== (e.mathml.split('</mfrac>').length - 1)) balanced = false;
  }
  ok(balanced, 'tags MathML balanceadas');
  ok(terms > 0, `${terms} fatores hoveraveis`);
}

console.log('\n=== o realce aponta para a regiao certa? ===');
{
  const { ctx } = buildTheory({ ...BASE, modeM: 1, modeN: 0 }, { kc: 137.44, alpha: 0 });
  const a = BASE.widthMM / 1000;
  const s = factorFor(ctx, 'sin_mx');
  const c = factorFor(ctx, 'cos_mx');
  // TE10: sin(pi x/a) zera nas paredes e e maximo no meio; cos e o oposto.
  ok(s({ x: 0, y: 0, z: 0 }) < 1e-9, 'sin(pi x/a) = 0 na parede x=0');
  ok(s({ x: a, y: 0, z: 0 }) < 1e-9, 'sin(pi x/a) = 0 na parede x=a');
  ok(Math.abs(s({ x: a / 2, y: 0, z: 0 }) - 1) < 1e-9, 'sin(pi x/a) = 1 no meio');
  ok(Math.abs(c({ x: 0, y: 0, z: 0 }) - 1) < 1e-9, 'cos(pi x/a) = 1 na parede');
  ok(c({ x: a / 2, y: 0, z: 0 }) < 1e-9, 'cos(pi x/a) = 0 no meio');
}
{
  // Abaixo do corte o fator de propagacao tem de DECAIR ao longo de z.
  const { ctx } = buildTheory(BASE, { kc: 137.44, alpha: 226.4 });
  const e = factorFor(ctx, 'expz');
  ok(e({ x:0, y:0, z:0 }) === 1 && e({ x:0, y:0, z:0.02 }) < 0.02,
     'e^(-alpha z) decai abaixo do corte');
}
{
  const { ctx } = buildTheory({ ...BASE, geometry: 1, modeN: 0, modeM: 1, field: 1 },
                              { kc: 160.79, alpha: 0 });
  const J = factorFor(ctx, 'Jn');
  const R = BASE.radiusMM / 1000;
  ok(Math.abs(J({ x: 0, y: 0, z: 0 }) - 1) < 1e-6, 'J0 maximo no eixo (rho=0)');
  ok(J({ x: R, y: 0, z: 0 }) < J({ x: 0, y: 0, z: 0 }), 'J0 menor na parede');
}
ok(Math.abs(besselJ(0, 2.4048255577)) < 1e-9, 'besselJ do JS bate com a raiz tabelada');

// The cylindrical set was wrong on the page for a while: E_rho and E_phi had
// J_n swapped with J_n-prime, the prefactors swapped with them, and TM had no
// branch at all, so it silently rendered the TE electric equations.
console.log('\n=== guia circular: as quatro combinacoes existem e diferem ===');
{
  const seen = new Map();
  for (const [mt, fl, name] of [[0,0,'TE-E'], [0,1,'TE-H'], [1,0,'TM-E'], [1,1,'TM-H']]) {
    const cfg = { ...BASE, geometry: 1, modeType: mt, field: fl, modeN: 1, modeM: 1 };
    const { eqs } = buildTheory(cfg, { kc: 160.79, alpha: 0 });
    const key = eqs.map(e => e.lhs + e.mathml).join('|');
    ok(eqs.length === 3, name + ': 3 componentes');
    ok(!seen.has(key), name + ': conjunto proprio (nao repete outro modo)');
    seen.set(key, name);
  }
}
{
  // Table 3.5: the axial component and whatever comes from d/drho carry the
  // derivative J_n-prime; anything born from d/dphi carries plain J_n AND the
  // n/(kc^2 rho) factor.
  const { eqs } = buildTheory({ ...BASE, geometry: 1, modeType: 0, field: 0, modeN: 1 },
                              { kc: 160.79, alpha: 0 });
  const Erho = eqs.find(e => e.lhs.indexOf('961') >= 0);
  const Ephi = eqs.find(e => e.lhs.indexOf('966') >= 0);
  const prime = t => t.indexOf('8242') >= 0;
  ok(!prime(Erho.mathml) && Erho.mathml.indexOf('mfrac') >= 0,
     'E_rho usa J_n (sem linha) com o fator n/(kc^2 rho)');
  ok(prime(Ephi.mathml), 'E_phi usa J_n-prime');
  ok(Erho.mathml.indexOf('data-term="gphi"') >= 0, 'E_rho carrega a combinacao g(phi)');
  ok(Ephi.mathml.indexOf('data-term="fphi"') >= 0, 'E_phi carrega a combinacao f(phi)');
}
{
  const { ctx } = buildTheory({ ...BASE, geometry: 1, modeN: 0, modeM: 1, field: 0 },
                              { kc: 160.79, alpha: 0 });
  ok(typeof factorFor(ctx, 'gphi') === 'function', 'g(phi) tem funcao de realce');
}

console.log(fail === 0 ? '\nOK: coluna de teoria consistente' : `\n${fail} FALHA(S)`);
process.exit(fail === 0 ? 0 : 1);
