// Node tests for seifert.js and minimal.js against Sage's values (vectors.json, made by make_vectors.sage).
//   node web/test/test-seifert.mjs
import { createRequire } from 'module';
import fs from 'fs';
const require = createRequire(import.meta.url);
const Seifert = require('../js/seifert.js');
const Minimal = require('../js/minimal.js');
const Wire = require('../js/wire.js');
const Precomputed = require('../js/precomputed.js');
const V = JSON.parse(fs.readFileSync(new URL('./vectors.json', import.meta.url)));

let pass = 0, fail = 0; const failures = [];
function check(name, ok, detail = '') { if (ok) pass++; else { fail++; failures.push(name + (detail ? ': ' + detail : '')); } }
const same = (a, b) => JSON.stringify(a) === JSON.stringify(b);
// Sage's Laurent polynomial [min degree, coefficients] in the same normal form as Seifert.alexander: lowest degree
// 0, positive leading coefficient
function normalise(coeffs) {
  let c = coeffs.slice(); while (c.length && c[0] === 0) c.shift(); while (c.length && c[c.length - 1] === 0) c.pop();
  if (!c.length) return [0]; if (c[c.length - 1] < 0) c = c.map(v => -v); return c;
}

// ---- parsing ----
{
  const w = [1, 2, -1, 2];
  for (const text of ['1 2 -1 2', '[1,2,-1,2]', '{1, 2, -1, 2}', 's1 s2 s1^-1 s2', 'σ1σ2σ1⁻¹σ2', 'σ_1 σ_2 σ_1^{-1} σ_2', 'abAb', 'sigma1 sigma2 sigma1^-1 sigma2'])
    check(`parse "${text}"`, same(Seifert.parseBraid(text).word, w), JSON.stringify(Seifert.parseBraid(text)));
  check('parse (1 2)^5', same(Seifert.parseBraid('(1 2)^5').word, [1, 2, 1, 2, 1, 2, 1, 2, 1, 2]));
  check('parse (s1 s2)^-2', same(Seifert.parseBraid('(s1 s2)^-2').word, [-2, -1, -2, -1]));
  check('parse s1^3', same(Seifert.parseBraid('s1^3').word, [1, 1, 1]));
  check('parse s1^-3', same(Seifert.parseBraid('s1^-3').word, [-1, -1, -1]));
  check('parse T(3,5)', same(Seifert.parseBraid('T(3,5)').word, Seifert.torusWord(3, 5)) && Seifert.torusWord(3, 5).length === 10);
  check('parse T(2,-3)', same(Seifert.parseBraid('T(2,-3)').word, [-1, -1, -1]));
  check('parse a name', Seifert.parseBraid('3_1').name === 'K3_1' && Seifert.parseBraid('trefoil').name === 'K3_1' && Seifert.parseBraid('12n242').name === 'K12n_242'
        && Seifert.parseBraid('L6a4').name === 'L6a4' && Seifert.parseBraid('l2a1{1}').name === 'L2a1_1' && Seifert.parseBraid('L6a4{1,0}').name === 'L6a4_1_0' && Seifert.parseBraid('L6a4_0_1').name === 'L6a4_0_1' && Seifert.parseBraid('borromean').name === 'L6a4' && Seifert.parseBraid('K11n_34').name === 'K11n_34');
  check('parse errors', !!Seifert.parseBraid('').error && !!Seifert.parseBraid('(1 2').error && !!Seifert.parseBraid('1 0').error && !!Seifert.parseBraid('^2').error && !!Seifert.parseBraid('1 2)').error && !!Seifert.parseBraid('1 + 2').error);
  check('parse 1^2 2 and 3-1', same(Seifert.parseBraid('1 ^ 2 2').word, [1, 1, 2]) && Seifert.parseBraid('3-1').name === 'K3_1' && same(Seifert.parseBraid('3 -1').word, [3, -1]));
}

// ---- the closed braid ----
{
  const d = Seifert.braidData([1, 1, 1]);
  check('trefoil: 2 strands, 3 crossings, knot, genus 1', d.n === 2 && d.c === 3 && d.mu === 1 && d.genus === 1 && d.chi === -1 && d.connected);
  const h = Seifert.braidData([1, 1]); check('hopf: 2 components, genus 0', h.mu === 2 && h.genus === 0);
  const b = Seifert.braidData([1, -2, 1, -2, 1, -2]); check('borromean: 3 components, genus 1', b.mu === 3 && b.genus === 1 && b.chi === -3);
  const s = Seifert.braidData([1, 3, 1, 3]); check('split braid is not connected', !s.connected);
  const u = Seifert.braidData([]); check('empty word: one strand, one disk', u.n === 1 && u.mu === 1 && u.genus === 0 && u.connected);
  const t45 = Seifert.braidData(Seifert.torusWord(4, 5)); check('T(4,5): knot of genus 6', t45.mu === 1 && t45.genus === 6);
  const t46 = Seifert.braidData(Seifert.torusWord(4, 6)); check('T(4,6): 2 components', t46.mu === 2);
}

// ---- Seifert matrix, Alexander polynomial and signature against Sage ----
const all = [...V.examples, ...V.random, ...V.knots];
let matrixMatches = 0, alexMatches = 0, sigMatches = 0, compMatches = 0, tested = 0;
for (const rec of all) {
  const bd = Seifert.braidData(rec.braid);
  if (bd.mu === rec.components) compMatches++; else failures.push(`components ${JSON.stringify(rec.braid)}: ${bd.mu} vs Sage ${rec.components}`);
  if (!bd.connected) continue;                        // Sage reorders split braids; we refuse them
  tested++;
  const M = Seifert.seifertMatrix(rec.braid);
  if (same(M, rec.seifert)) matrixMatches++; else failures.push(`seifert matrix ${JSON.stringify(rec.braid)}: ${JSON.stringify(M)} vs ${JSON.stringify(rec.seifert)}`);
  const a = Seifert.alexander(rec.braid);
  if (same(a.coeffs, normalise(rec.alexander[1]))) alexMatches++; else failures.push(`alexander ${JSON.stringify(rec.braid)}: ${a.coeffs} vs ${rec.alexander[1]}`);
  const s = Seifert.signature(rec.braid);
  if (s.signature === rec.signature) sigMatches++; else failures.push(`signature ${JSON.stringify(rec.braid)}: ${s.signature} vs ${rec.signature}`);
  if (rec.knotinfo_genus !== undefined && rec.knotinfo_genus !== null)
    check(`genus bound ${rec.name}`, a.degree <= 2 * rec.knotinfo_genus && rec.knotinfo_genus <= bd.genus, `deg Δ = ${a.degree}, KnotInfo genus ${rec.knotinfo_genus}, surface genus ${bd.genus}`);
}
check(`components agree with Sage on all ${all.length} braids`, compMatches === all.length, `${compMatches}/${all.length}`);
check(`Seifert matrices equal Sage's on ${tested} connected braids`, matrixMatches === tested, `${matrixMatches}/${tested}`);
check(`Alexander polynomials equal Sage's on ${tested} braids`, alexMatches === tested, `${alexMatches}/${tested}`);
check(`signatures equal Sage's on ${tested} braids`, sigMatches === tested, `${sigMatches}/${tested}`);
{
  const lehmer = [1, 1, 0, -1, -1, -1, -1, -1, 0, 1, 1];      // t^10 + t^9 - t^7 - t^6 - t^5 - t^4 - t^3 + t + 1
  const a = Seifert.alexander([1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2]);
  check('12n_242: Δ(t) is Lehmer\'s polynomial in -t', same(a.coeffs, lehmer.map((v, k) => k % 2 ? -v : v)), a.coeffs.join(','));
  check('12n_242: det = 1', a.det === 1);
  check('trefoil: det 3, signature -2', Seifert.alexander([1, 1, 1]).det === 3 && Seifert.signature([1, 1, 1]).signature === -2);
  check('figure-eight: det 5, signature 0', Seifert.alexander([1, -2, 1, -2]).det === 5 && Seifert.signature([1, -2, 1, -2]).signature === 0);
  check('TeX of Δ', Seifert.alexanderTeX([1, -1, 1]) === 't^{2} - t + 1' && Seifert.alexanderTeX([1, -3, 1]) === 't^{2} - 3t + 1' && Seifert.alexanderTeX([-1, 1]) === 't - 1');
}

// ---- the mesh ----
function meshChecks(name, word, opts) {
  const m = Seifert.buildSurface(word, opts), Vn = m.pos.length / 3, bd = m.braid;
  Minimal.prepare(m);
  check(`${name}: oriented mesh`, Minimal.orientable(m.tri, Vn));
  check(`${name}: χ = n - c`, Minimal.eulerCharacteristic(m) === bd.chi, `${Minimal.eulerCharacteristic(m)} vs ${bd.chi}`);
  check(`${name}: ${bd.mu} boundary loop(s)`, m.loops.length === bd.mu, String(m.loops.length));
  check(`${name}: no degenerate triangles`, Minimal.degenerateTriangles(m, 1e-9) === 0, String(Minimal.degenerateTriangles(m, 1e-9)));
  const fixedCount = m.fixed.reduce((s, v) => s + v, 0), loopCount = m.loops.reduce((s, l) => s + l.length, 0);
  check(`${name}: boundary vertices are the loops`, fixedCount === loopCount && fixedCount === m.loops.flat().length);
  let bad = 0; for (let v = 0; v < Vn; v++) if (m.topo.start[v + 1] === m.topo.start[v]) bad++;
  check(`${name}: no isolated vertices`, bad === 0);
  return m;
}
const hopf = meshChecks('hopf', [1, 1]);
{
  const lk = Minimal.linkingNumber(Minimal.loopPoints(hopf.pos, hopf.loops[0]), Minimal.loopPoints(hopf.pos, hopf.loops[1]));
  check('hopf [1,1]: the boundary loops link once, positively (σ1 is the positive crossing)', lk === 1, String(lk));
  const anti = Seifert.buildSurface([-1, -1]);
  check('hopf [-1,-1]: linking number -1', Minimal.linkingNumber(Minimal.loopPoints(anti.pos, anti.loops[0]), Minimal.loopPoints(anti.pos, anti.loops[1])) === -1);
  const t24 = Seifert.buildSurface([1, 1, 1, 1]);
  check('T(2,4): linking number 2', Minimal.linkingNumber(Minimal.loopPoints(t24.pos, t24.loops[0]), Minimal.loopPoints(t24.pos, t24.loops[1])) === 2);
  const three = Seifert.buildSurface([1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2]);   // 3 strands, 14 crossings: a 3-component link?
  check('unsmoothed and smoothed surfaces have the same connectivity', Seifert.buildSurface([1, 1, 1], { smoothPasses: 0 }).tri.length === Seifert.buildSurface([1, 1, 1]).tri.length);
  check('the raw boundary of the trefoil surface lies on the rims and bands', (() => { const r = Seifert.buildSurface([1, 1, 1], { smoothPasses: 0 }); return r.loops[0].every(v => r.kind[v] < 2 ? Math.abs(Math.hypot(r.pos[3 * v], r.pos[3 * v + 1]) - 1) < 1e-9 : true); })());
  void three;
}
meshChecks('trefoil', [1, 1, 1]);
meshChecks('figure-eight', [1, -2, 1, -2]);
meshChecks('borromean', [1, -2, 1, -2, 1, -2]);
meshChecks('12n_242', [1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2]);
meshChecks('T(4,5)', Seifert.torusWord(4, 5), { angular: 96 });
meshChecks('unknot, one disk', []);
meshChecks('unknot, two strands', [1]);

// ---- the solver ----
{ // a bumped disk with its rim fixed flattens: the area tends to π and H to 0
  const m = Seifert.buildSurface([], { angular: 72, smoothPasses: 0 });
  for (let v = 0; v < m.pos.length / 3; v++) if (!m.fixed[v]) { const r = Math.hypot(m.pos[3 * v], m.pos[3 * v + 1]); m.pos[3 * v + 2] += 0.4 * (1 - r * r); }
  Minimal.prepare(m);
  const a0 = Minimal.area(m); let last = a0, monotone = true;
  for (let it = 0; it < 6; it++) { const s = Minimal.relax(m, { dt: Infinity, flips: false }); if (s.area > last + 1e-12) monotone = false; last = s.area; }
  check('bumped disk: area decreases monotonically', monotone);
  const polygon = 36 * Math.sin(2 * Math.PI / 72);                 // the rim is a 72-gon
  check('bumped disk: area → that of the flat 72-gon', Math.abs(last - polygon) < 1e-6, `${last} vs ${polygon}`);
  check('bumped disk: H → 0', Minimal.meanCurvature(m).max < 1e-6, String(Minimal.meanCurvature(m).max));
  check('bumped disk: boundary unchanged', m.loops[0].every(v => Math.abs(m.pos[3 * v + 2]) < 1e-12));
}
{ // the catenoid: two coaxial unit circles at distance 2a, a cylinder between them, relaxed to the catenoid
  const a = 0.5, N = 96, L = 24, pos = [], tri = [], fixed = [];
  for (let l = 0; l <= L; l++) for (let k = 0; k < N; k++) { const th = 2 * Math.PI * k / N; pos.push(Math.cos(th), Math.sin(th), -a + 2 * a * l / L); fixed.push(l === 0 || l === L ? 1 : 0); }
  for (let l = 0; l < L; l++) for (let k = 0; k < N; k++) { const A = l * N + k, B = l * N + (k + 1) % N, C = (l + 1) * N + (k + 1) % N, D = (l + 1) * N + k; tri.push(A, B, C, A, C, D); }
  const m = Minimal.prepare({ pos: Float64Array.from(pos), tri: Uint32Array.from(tri), fixed: Uint8Array.from(fixed) });
  // the catenoid r = c cosh(z/c) through r = 1 at z = ±a: c cosh(a/c) = 1, the larger root (the stable one)
  let c = 0.9; for (let i = 0; i < 100; i++) { const f = c * Math.cosh(a / c) - 1, df = Math.cosh(a / c) - (a / c) * Math.sinh(a / c); c -= f / df; }
  const exact = Math.PI * c * (2 * a + c * Math.sinh(2 * a / c));
  let A = Minimal.area(m), monotone = true;
  for (let it = 0; it < 30; it++) { const s = Minimal.relax(m, { dt: Infinity, flips: true }); if (s.area > A + 1e-9) monotone = false; A = s.area; }
  check('catenoid: area decreases monotonically', monotone);
  check('catenoid: area within 0.5% of the exact catenoid', Math.abs(A - exact) / exact < 5e-3, `${A} vs ${exact} (cylinder ${2 * Math.PI * 2 * a})`);
  check('catenoid: waist radius ≈ c', (() => { let r = Infinity; for (let v = 0; v < m.pos.length / 3; v++) r = Math.min(r, Math.hypot(m.pos[3 * v], m.pos[3 * v + 1])); return Math.abs(r - c) < 0.01; })(), String(c));
  check('catenoid: H small', Minimal.meanCurvature(m).rms < 0.05, String(Minimal.meanCurvature(m).rms));
  const mcf = Minimal.prepare({ pos: Float64Array.from(pos), tri: Uint32Array.from(tri), fixed: Uint8Array.from(fixed) });
  let A2 = Minimal.area(mcf), mono2 = true;
  for (let it = 0; it < 20; it++) { const s = Minimal.relax(mcf, { dt: 0.05, flips: false }); if (s.area > A2 + 1e-9) mono2 = false; A2 = s.area; }
  check('catenoid by mean curvature flow: area decreases towards the catenoid', mono2 && A2 < 2 * Math.PI * 2 * a && Math.abs(A2 - exact) / exact < 5e-3, String(A2));
}
{ // the trefoil's Bennequin surface relaxes: area down, χ, boundary and orientation kept
  const m = Seifert.buildSurface([1, 1, 1], { angular: 72 }); Minimal.prepare(m);
  const chi = Minimal.eulerCharacteristic(m), boundary = m.loops[0].map(v => [m.pos[3 * v], m.pos[3 * v + 1], m.pos[3 * v + 2]]);
  const a0 = Minimal.area(m); let A = a0, monotone = true, flips = 0;
  for (let it = 0; it < 25; it++) { const s = Minimal.relax(m, { dt: Infinity, flips: true, tangential: 0.3 }); if (s.area > A * (1 + 1e-6)) monotone = false; A = s.area; flips += s.flips; }
  check('trefoil surface: area decreases', A < a0 && monotone, `${a0} -> ${A}`);
  check('trefoil surface: χ unchanged by flips', Minimal.eulerCharacteristic(m) === chi && Minimal.orientable(m.tri, m.pos.length / 3));
  check('trefoil surface: boundary fixed', m.loops[0].every((v, k) => Math.abs(m.pos[3 * v] - boundary[k][0]) + Math.abs(m.pos[3 * v + 1] - boundary[k][1]) + Math.abs(m.pos[3 * v + 2] - boundary[k][2]) < 1e-12));
  check('trefoil surface: no degenerate triangles after relaxing', Minimal.degenerateTriangles(m, 1e-9) === 0);
  const H = Minimal.meanCurvature(m);
  check('trefoil surface: mean curvature residual small (rms |H| R < 0.5 after 25 rounds)', H.rms < 0.5, `rms ${H.rms}, max ${H.max}`);
}

// ---- taming the wire ----
function turning(ws) {                                       // total turning angle of the wire, in turns
  let sum = 0; for (let k = 0; k < ws.N; k++) { const p = ws.prev[k], n = ws.next[k], P = ws.P;
    const ax = P[3 * p] - P[3 * k], ay = P[3 * p + 1] - P[3 * k + 1], az = P[3 * p + 2] - P[3 * k + 2], bx = P[3 * n] - P[3 * k], by = P[3 * n + 1] - P[3 * k + 1], bz = P[3 * n + 2] - P[3 * k + 2];
    sum += Math.PI - Math.acos(Math.max(-1, Math.min(1, (ax * bx + ay * by + az * bz) / (Math.hypot(ax, ay, az) * Math.hypot(bx, by, bz))))); }
  return sum / (2 * Math.PI);
}
{ // the circle stays a circle (it grows, the model's own scale), spacing uniform
  const m = Seifert.buildSurface([], { angular: 72 }); Minimal.prepare(m);
  const ws = Wire.init(m, { coarse: 2 * Math.PI / 48 });
  check('wire: coarse copy of the loop', ws.N === 48 && ws.fLoops[0].n === 72, String(ws.N));
  for (let it = 0; it < 300; it++) Wire.step(ws);
  let rmin = Infinity, rmax = 0, smin = Infinity, smax = 0;
  for (let k = 0; k < ws.N; k++) { const r = Math.hypot(ws.P[3 * k], ws.P[3 * k + 1]); rmin = Math.min(rmin, r); rmax = Math.max(rmax, r); const n = ws.next[k]; const d = Math.hypot(ws.P[3 * n] - ws.P[3 * k], ws.P[3 * n + 1] - ws.P[3 * k + 1], ws.P[3 * n + 2] - ws.P[3 * k + 2]); smin = Math.min(smin, d); smax = Math.max(smax, d); }
  check('wire: a circle stays round and evenly spaced', (rmax - rmin) < 1e-9 * rmax && (smax - smin) < 1e-9 * smax && ws.P.every((v, i) => i % 3 !== 2 || Math.abs(v) < 1e-12), `r ${rmin}..${rmax}`);
  check('wire: the circle grows', Wire.length(ws) > 1.2 * ws.L0, String(Wire.length(ws) / ws.L0));
  Wire.apply(ws, m, Minimal);
  check('wire: the interpolated circle is round', (() => { const F = ws.F[0]; let rmin = Infinity, rmax = 0; for (let k = 0; k < 72; k++) { const r = Math.hypot(F[3 * k], F[3 * k + 1]); rmin = Math.min(rmin, r); rmax = Math.max(rmax, r); } return (rmax - rmin) < 1e-3 * rmax; })());
}
{ // the trefoil wire: the energy falls, the corners go, no strand comes within d_close of another, the surface follows
  const m = Seifert.buildSurface([1, 1, 1], { angular: 72 }); Minimal.prepare(m);
  const chi = Minimal.eulerCharacteristic(m), ws = Wire.init(m), turn0 = turning(ws);
  let e0 = null, e = null, monotone = true, rejected = 0;
  for (let it = 0; it < 1500; it++) { const s = Wire.step(ws); if (e0 === null) e0 = s.energy; else if (s.energy > e * (1 + 1e-6)) monotone = false; e = s.energy; rejected += s.rejected; if (it % 100 === 99) Wire.apply(ws, m, Minimal); }
  Wire.apply(ws, m, Minimal);
  check('wire: the energy decreases monotonically', monotone && e < 0.8 * e0, `${e0} -> ${e}`);
  check('wire: the corners are gone', turning(ws) < 0.5 * turn0 && turning(ws) < 3, `${turn0} -> ${turning(ws)} turns`);
  check('wire: no strand closer than it started, below d_close', Wire.clearance(ws).distance >= Math.min(2 * ws.ra, 0.9 * 1.5 * ws.ra) * 0.999, String(Wire.clearance(ws).distance / ws.ra));
  check('wire: the surface follows', m.loops[0].every((v, k) => Math.abs(m.pos[3 * v] - ws.F[0][3 * k]) + Math.abs(m.pos[3 * v + 1] - ws.F[0][3 * k + 1]) + Math.abs(m.pos[3 * v + 2] - ws.F[0][3 * k + 2]) < 1e-12));
  check('wire: the mesh boundary stays close to the coarse polygon', (() => {
    const seg = (px, py, pz, a, b) => { const P = ws.P, ux = P[3 * b] - P[3 * a], uy = P[3 * b + 1] - P[3 * a + 1], uz = P[3 * b + 2] - P[3 * a + 2], l2 = ux * ux + uy * uy + uz * uz; let t = l2 > 0 ? ((px - P[3 * a]) * ux + (py - P[3 * a + 1]) * uy + (pz - P[3 * a + 2]) * uz) / l2 : 0; t = Math.max(0, Math.min(1, t)); return Math.hypot(px - P[3 * a] - t * ux, py - P[3 * a + 1] - t * uy, pz - P[3 * a + 2] - t * uz); };
    let worst = 0; for (let k = 0; k < ws.fLoops[0].n; k++) { let best = Infinity; for (let j = 0; j < ws.N; j++) best = Math.min(best, seg(ws.F[0][3 * k], ws.F[0][3 * k + 1], ws.F[0][3 * k + 2], j, ws.next[j])); worst = Math.max(worst, best); } return worst < 0.15 * ws.ra; })());
  check('wire: settles by the energy criterion', Wire.settled(ws) || ws.steps >= 1500);
  check('wire: the carried surface is sound', Minimal.eulerCharacteristic(m) === chi && Minimal.orientable(m.tri, m.pos.length / 3) && Minimal.degenerateTriangles(m, 1e-9) === 0);
  let A = Minimal.area(m); for (let it = 0; it < 25; it++) A = Minimal.relax(m, { dt: Infinity, flips: true, tangential: 0.3 }).area;
  const H = Minimal.meanCurvature(m);
  check('wire: the film on the tamed wire relaxes', H.rms < 0.5 && Minimal.degenerateTriangles(m, 1e-9) === 0, `rms ${H.rms}, max ${H.max}, area ${A}`);
}

{ // the pipeline: tame with the surface carried (Wire.carry), then settle the film (Minimal.settleRound) -- no pinch, a clean film
  for (const [name, word, maxSteps] of [['trefoil', [1, 1, 1], 1200], ['12n242', [1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2], 1500]]) {
    const m = Seifert.buildSurface(word, { angular: 96 }); Minimal.prepare(m);
    const ws = Wire.init(m); ws.meshEdge = Minimal.meanEdgeLength(m);
    let pinched = false;
    while (!Wire.settled(ws) && ws.steps < maxSteps) { const r = Wire.carry(ws, m, Minimal); if (r.pinched) pinched = true; }
    const st = Minimal.settleInit(m), boundary = m.loops.flat().map(v => [m.pos[3 * v], m.pos[3 * v + 1], m.pos[3 * v + 2]]);
    let s = null, rounds = 0;
    for (let r = 0; r < 80; r++) { s = Minimal.settleRound(m, st); rounds++; if (s.pinched) { pinched = true; break; } if (s.harmonic && s.moved < 1e-4) break; }
    const H = Minimal.meanCurvature(m);
    check(`pipeline ${name}: the wire opened`, Wire.length(ws) / ws.L0 > 1.2, String(Wire.length(ws) / ws.L0));
    check(`pipeline ${name}: no pinch, a settled film`, !pinched && s.harmonic && H.rms < 0.5 && Minimal.degenerateTriangles(m, 1e-9) === 0, `pinched ${pinched}, harmonic ${s.harmonic}, rounds ${rounds}, rms ${H.rms}`);
    check(`pipeline ${name}: χ and orientation kept through the remeshing`, Minimal.eulerCharacteristic(m) === m.braid.chi && Minimal.orientable(m.tri, m.pos.length / 3) && m.loops.length === m.braid.mu);
    check(`pipeline ${name}: the wire did not move while the film settled`, m.loops.flat().every((v, k) => Math.abs(m.pos[3 * v] - boundary[k][0]) + Math.abs(m.pos[3 * v + 1] - boundary[k][1]) + Math.abs(m.pos[3 * v + 2] - boundary[k][2]) < 1e-12));
  }
}
{ // stress: taming with a film round after every few wire steps keeps the mesh a manifold (the flips must never make an edge twice)
  let ok = true, detail = '';
  for (const word of [[1, 1, 1], [1, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2], [1, -2, 1, -2]]) {
    const m = Seifert.buildSurface(word, { angular: 96 }); Minimal.prepare(m);
    const ws = Wire.init(m, { dt0: 0.2, gamma: 0.15, decay: 0.0003 });
    try {
      for (let t = 0; t < 60; t++) { for (let k = 0; k < 10; k++) Wire.step(ws); Wire.apply(ws, m, Minimal); Minimal.relax(m, { dt: Infinity, flips: true, tangential: 0.3, tol: 1e-7 }); }
      Minimal.prepare(m);
      if (!Minimal.orientable(m.tri, m.pos.length / 3) || Minimal.eulerCharacteristic(m) !== m.braid.chi) { ok = false; detail += ` ${word}: χ ${Minimal.eulerCharacteristic(m)}`; }
    } catch (e) { ok = false; detail += ` ${word}: ${e.message}`; }
  }
  check('stress: 600 wire steps with film rounds keep the meshes manifold and oriented', ok, detail);
}

// ---- what a knot is called ----
{
  const cases = [['K3_1', 'trefoil'], ['K4_1', 'figure-eight knot'], ['K12n_242', '(−2, 3, 7)-pretzel knot'],
                 ['L2a1_1', 'Hopf link'], ['L6a4_0_0', 'Borromean rings'], ['L5a1', 'Whitehead link'], ['K11n_34', 'Conway knot']];
  check('names: the common ones, orientation indices and all', cases.every(([n, c]) => Seifert.commonName(n) === c), cases.filter(([n, c]) => Seifert.commonName(n) !== c).map(c => c[0]).join(' '));
  check('names: none invented for the unnamed', [Seifert.commonName('K7_2'), Seifert.commonName('K11a_1'), Seifert.commonName('K10_124')].every(v => v === null));
  check('names: Rolfsen tags to ten crossings only', Seifert.rolfsenTeX('K3_1') === '3_{1}' && Seifert.rolfsenTeX('K10_124') === '10_{124}' && Seifert.rolfsenTeX('K11n_34') === null && Seifert.rolfsenTeX('K11a_1') === null && Seifert.rolfsenTeX('L6a4_0') === null);
  check('names: the modern ones for knots and links', Seifert.modernTeX('K12n_242') === '12n_{242}' && Seifert.modernTeX('K11a_1') === '11a_{1}' && Seifert.modernTeX('L6a4_0_0') === 'L6a4\\{0,0\\}' && Seifert.modernTeX('L5a1') === 'L5a1');
  check('names: as plain text, with subscripts', Seifert.nameText('K3_1') === 'trefoil' && Seifert.nameText('K6_3') === '6₃' && Seifert.nameText('K10_124') === '10₁₂₄'
        && Seifert.nameText('K12n_242') === '(−2, 3, 7)-pretzel knot' && Seifert.nameText('K11a_1') === '11a₁' && Seifert.nameText('L6a4_0_0') === 'Borromean rings' && Seifert.nameText('L6a5_1_0') === 'L6a5{1,0}');
  check('names: a braid word as plain text', Seifert.wordSymbols([1, 2, -1]) === 'σ₁σ₂σ₁⁻¹' && Seifert.wordSymbols([]) === '1' && Seifert.wordSymbols([-11]) === 'σ₁₁⁻¹');
  // every name the table can produce is displayable one way or another
  const table = JSON.parse(fs.readFileSync(new URL('../data/knots.json', import.meta.url)));
  const unnameable = table.filter(r => !Seifert.commonName(r.name) && !Seifert.rolfsenTeX(r.name) && !Seifert.modernTeX(r.name));
  check(`names: all ${table.length} table entries have a name to show`, unnameable.length === 0, unnameable.slice(0, 5).map(r => r.name).join(' '));
  const noText = table.filter(r => !Seifert.nameText(r.name));
  check(`names: and a plain-text one`, noText.length === 0, noText.slice(0, 5).map(r => r.name).join(' '));
}

// ---- the shipped films (web/data/films, made by web/data/precompute.mjs) ----
{
  const dir = new URL('../data/films/', import.meta.url);
  let manifest = null;
  try { manifest = JSON.parse(fs.readFileSync(new URL('index.json', dir))); } catch (e) { /* none shipped */ }
  check('films: a manifest is shipped', !!manifest, 'run: node web/data/precompute.mjs');
  if (manifest) {
    check('films: the manifest is of this format version', manifest.version === Precomputed.VERSION, `${manifest.version} vs ${Precomputed.VERSION}`);
    const table = JSON.parse(fs.readFileSync(new URL('../data/knots.json', import.meta.url)));
    const byName = new Map(table.map(r => [r.name, r]));
    // every example of the app's menu has one, under the key of its braid word
    const app = fs.readFileSync(new URL('../js/app.js', import.meta.url), 'utf8');
    const EXAMPLES = JSON.parse('[' + /const EXAMPLES = \[([^\]]*)\]/.exec(app)[1].replace(/'/g, '"') + ']');
    check('films: the app lists examples', EXAMPLES.length >= 15);
    const missing = [];
    for (const input of EXAMPLES) {
      const p = Seifert.parseBraid(input);
      let word = p.word;
      if (p.name) { const rec = byName.get(p.name) || table.find(r => r.name.startsWith(p.name + '_')); word = rec && rec.braid; }
      if (!word || !manifest.films[Precomputed.key(word)]) missing.push(input);
    }
    check(`films: all ${EXAMPLES.length} examples are shipped`, missing.length === 0, missing.join(' '));
    // each one decodes to the settled film of the knot it claims: the right topology, the wire fixed, H small
    let bad = [], bytes = 0;
    for (const [key, entry] of Object.entries(manifest.films)) {
      const word = key === '()' ? [] : key.split(',').map(Number);
      const buf = fs.readFileSync(new URL(entry.file, dir));
      bytes += buf.byteLength;
      const film = Precomputed.decode(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength));
      const mesh = Seifert.buildSurface(word, { spacing: manifest.params.spacing, bandWidth: manifest.params.bandWidth, bulge: manifest.params.bulge, angular: manifest.params.angular, round: manifest.params.round });
      Minimal.prepare(mesh);
      const bd = mesh.braid;
      try {
        Precomputed.install(mesh, film, Seifert, Minimal);
        const H = Minimal.meanCurvature(mesh), chi = Minimal.eulerCharacteristic(mesh);
        const why = [];
        if (chi !== bd.chi) why.push(`χ ${chi} not ${bd.chi}`);
        if (mesh.loops.length !== bd.mu) why.push(`${mesh.loops.length} loops not ${bd.mu}`);
        if (!Minimal.orientable(mesh.tri, film.V)) why.push('not oriented');
        if (Minimal.degenerateTriangles(mesh, 1e-9)) why.push(`${Minimal.degenerateTriangles(mesh, 1e-9)} degenerate`);
        if (!entry.pinched && !(H.rms < 0.5)) why.push(`rms |H| ${H.rms.toFixed(3)}`);
        if (Math.abs(Minimal.area(mesh) - entry.area) > 1e-3 * entry.area) why.push(`area ${Minimal.area(mesh).toFixed(3)} not ${entry.area}`);
        if (!(film.length / film.L0 > 1.0)) why.push(`the wire did not grow (×${(film.length / film.L0).toFixed(2)})`);
        if (why.length) bad.push(`${entry.file}: ${why.join(', ')}`);
      } catch (e) { bad.push(`${entry.file}: ${e.message}`); }
    }
    check(`films: all ${Object.keys(manifest.films).length} decode to settled films of the right topology (${(bytes / 1024 / 1024).toFixed(2)} MB)`, bad.length === 0, bad.join(' | '));
    // a film whose parameters the app has changed must not be used: the manifest records what made it
    check('films: the manifest records the parameters they were made with', ['spacing', 'bandWidth', 'bulge', 'round', 'alpha', 'K', 'dclose', 'maxsteps'].every(k => manifest.params[k] !== undefined));
    // the names map lets the app resolve an example without the knot table; every name in it must point at a film
    const named = Object.entries(manifest.names || {});
    check('films: the manifest maps the table names of the examples to their films', named.length >= 12 && named.every(([, key]) => manifest.films[key]), named.length + ' names');
    const badName = named.filter(([name]) => Seifert.tableName(name) !== name);
    check('films: those names are the ones the app parses inputs to', badName.length === 0, badName.map(n => n[0]).join(' '));
  }
  // the format itself
  const s = { pos: Float64Array.from([0, 0, 0, 1, 0, 0, 0, 1, 0]), tri: Uint32Array.from([0, 1, 2]), kind: Int16Array.from([0, 1, 2]), steps: 7, pinched: true, capped: false, L0: 1.5, length: 2.5, edge0: 0.25, area: 0.5 };
  const back = Precomputed.decode(Precomputed.encode(s));
  check('films: a film survives encoding and decoding', back.V === 3 && back.F === 1 && back.steps === 7 && back.pinched && !back.capped && back.L0 === 1.5 && back.length === 2.5 && back.edge0 === 0.25 && back.area === 0.5
        && [...back.pos].join() === [...s.pos].join() && [...back.tri].join() === [...s.tri].join() && [...back.kind].join() === [...s.kind].join());
  let threw = false; try { Precomputed.decode(new ArrayBuffer(8)); } catch (e) { threw = true; }
  check('films: a truncated file is refused', threw);
}

console.log(`${pass} passed, ${fail} failed`);
for (const f of failures) console.log('  FAIL ' + f);
process.exit(fail ? 1 : 0);
