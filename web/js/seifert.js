/* seifert.js -- braids, the Bennequin surface of a closed braid, and the invariants read off it.
   Node-compatible (module.exports) and browser (window.Seifert).  No dependencies.

   A braid on n strands with c crossings closes to a link.  Seifert's algorithm on the closed braid gives the
   n Seifert circles (one per strand: the oriented smoothing of every crossing joins strand i to strand i), so the
   surface is n disks, one per strand, stacked along the braid axis, joined by one half-twisted band per crossing.
   This is the Bennequin surface of the closed braid; its Euler characteristic is n - c, and for a link with mu
   components its genus is (c - n + 2 - mu)/2.  The Seifert matrix is computed by the algorithm of Collins
   (2013, section 3.3) exactly as Sage's Link.seifert_matrix does, so the tests can compare against Sage. */
(function () {
'use strict';
const Seifert = {};

// ------------------------------------------------------------------ braid words
// Accepted: integer lists "1 2 -1 2", "[1,2,-1,2]"; sigma words "s1 s2 s1^-1", "σ1σ2σ1⁻¹", "σ_1^3"; letters "abA"
// (a = σ1, A = σ1^-1: the reverse of SeifertView's case); groups with powers "(1 2)^5"; torus links "T(3,5)"; a table name ("3_1", "12n242", "L6a4")
// or an alias ("trefoil") is returned as { name } for the app to look up.
const SUP = { '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9', '⁻': '-' };
const ALIASES = {
  unknot: 'K0_1', trefoil: 'K3_1', 'figure-eight': 'K4_1', 'figure eight': 'K4_1', figureeight: 'K4_1', cinquefoil: 'K5_1',
  'three-twist': 'K5_2', 'three twist': 'K5_2', stevedore: 'K6_1', 'miller institute': 'K6_2', conway: 'K11n_34',
  'kinoshita-terasaka': 'K11n_42', 'kinoshita terasaka': 'K11n_42', 'pretzel(-2,3,7)': 'K12n_242',
  'p(-2,3,7)': 'K12n_242', '(-2,3,7)': 'K12n_242', hopf: 'L2a1_1', whitehead: 'L5a1', borromean: 'L6a4', solomon: 'L4a1_1',
};
Seifert.ALIASES = ALIASES;
// What to call a knot or link, in the order one would say it: its common name, else its Rolfsen tag (the knot
// tables' own numbering, which Rolfsen's tables carry to ten crossings), else the modern name of the
// Hoste-Thistlethwaite tables.  A link's orientation indices are not part of its common name.
const COMMON = {
  K0_1: 'unknot', K3_1: 'trefoil', K4_1: 'figure-eight knot', K5_1: 'cinquefoil', K5_2: 'three-twist knot',
  K6_1: 'stevedore knot', K6_2: 'Miller Institute knot', K7_4: 'endless knot', K8_18: 'carrick bend',
  K11n_34: 'Conway knot', K11n_42: 'Kinoshita–Terasaka knot', K12n_242: '(−2, 3, 7)-pretzel knot',
  L2a1: 'Hopf link', L4a1: "Solomon's link", L5a1: 'Whitehead link', L6a4: 'Borromean rings',
};
Seifert.commonName = name => COMMON[name] || COMMON[String(name).replace(/(_\d+)+$/, '')] || null;
// the Rolfsen tag as TeX, for the knots that have one (to ten crossings, where the tables need no a/n letter)
Seifert.rolfsenTeX = name => { const m = /^K(\d{1,2})_(\d+)$/.exec(name || ''); return m && +m[1] <= 10 ? `${m[1]}_{${m[2]}}` : null; };
// the modern name as TeX: 11n_{34}, L6a4\{0,0\}
Seifert.modernTeX = name => {
  let m = /^K(\d+)([an]?)_(\d+)$/.exec(name || '');
  if (m) return `${m[1]}${m[2]}_{${m[3]}}`;
  m = /^L(\d+)([an])(\d+)((?:_\d+)*)$/.exec(name || '');
  if (m) return `L${m[1]}${m[2]}${m[3]}` + (m[4] ? `\\{${m[4].slice(1).split('_').join(',')}\\}` : '');
  return null;
};
// KnotInfo names: knots K3_1, K11n_34, K12a_5; LinkInfo names L2a1_0, L6a4_1_0 (one orientation index per component
// after the first).  Accepts 3_1, 3-1, 11n34, 12n_242, L6a4, L6a4{1,0}, L6a4_1_0; a link without indices is returned
// bare (L6a4) and the app takes the table's first entry for it.
function tableName(text) {
  let s = text.trim().toLowerCase().replace(/\s+/g, ' ');
  if (ALIASES[s]) return ALIASES[s];
  if (/\s/.test(s)) return null;                           // "1 2 -1 2" is a braid word, not 12_12
  let m = s.match(/^l(\d+)([an])(\d+)(?:[_{]([\d,_]+)\}?)?$/);
  if (m) return `L${m[1]}${m[2]}${m[3]}` + (m[4] ? '_' + m[4].split(/[,_]/).join('_') : '');
  m = s.match(/^k?(\d+)([an]?)[_-]?(\d+)$/);
  if (m) return `K${m[1]}${m[2]}_${m[3]}`;
  return null;
}
Seifert.tableName = tableName;

function torusWord(p, q) {                          // T(p, q) as the closure of (σ1 σ2 ... σ_{p-1})^q on p strands
  const word = [], s = q < 0 ? -1 : 1;
  for (let r = 0; r < Math.abs(q); r++) for (let i = 1; i < p; i++) word.push(s * i);
  return word;
}
Seifert.torusWord = torusWord;

function parseBraid(text) {
  let s = String(text || '').trim();
  if (!s) return { error: 'empty braid word' };
  const name = tableName(s);
  if (name) return { name, text: s };
  let m = s.match(/^t\s*\(\s*(-?\d+)\s*,\s*(-?\d+)\s*\)$/i);
  if (m) {
    const p = +m[1], q = +m[2];
    if (p < 2) return { error: 'T(p, q) needs p ≥ 2' };
    if (p > 12 || Math.abs(q) * (p - 1) > 400) return { error: 'that torus link is too big to draw' };
    return { word: torusWord(p, q), text: s, label: `T(${p}, ${q})` };
  }
  s = s.replace(/[⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+/g, sup => '^' + [...sup].map(ch => SUP[ch]).join(''));
  s = s.replace(/^(braid|b)\s*[:=]\s*/i, '').replace(/[\[\]{}]/g, ' ').replace(/,/g, ' ');
  const tokens = [];
  const re = /\s+|\(|\)|\^\s*\{?\s*(-?\d+)\s*\}?|(?:σ|s|sigma)_?\s*(\d+)|(-?\d+)|([a-zA-Z])/gy;
  let pos = 0, tk;
  while (pos < s.length) {
    re.lastIndex = pos; tk = re.exec(s);
    if (!tk || tk.index !== pos) return { error: `cannot read the braid word at "${s.slice(pos, pos + 12)}"` };
    pos = re.lastIndex;
    const t = tk[0];
    if (/^\s+$/.test(t)) continue;
    if (t === '(' || t === ')') tokens.push({ t });
    else if (tk[1] !== undefined) tokens.push({ t: '^', n: +tk[1] });
    else if (tk[2] !== undefined) tokens.push({ t: 'g', g: +tk[2] });
    else if (tk[3] !== undefined) tokens.push({ t: 'g', g: +tk[3] });
    else { const ch = tk[4], lower = ch.toLowerCase(); const g = lower.charCodeAt(0) - 96; tokens.push({ t: 'g', g: ch === lower ? g : -g }); }
  }
  let k = 0;
  function expr() { const out = []; while (k < tokens.length && tokens[k].t !== ')') out.push(...term()); return out; }
  function term() {
    let atom;
    if (tokens[k].t === '(') { k++; atom = expr(); if (k >= tokens.length || tokens[k].t !== ')') throw new Error('unbalanced parentheses'); k++; }
    else if (tokens[k].t === 'g') { atom = [tokens[k].g]; k++; }
    else throw new Error('unexpected "' + (tokens[k].t === '^' ? '^' : tokens[k].t) + '"');
    if (k < tokens.length && tokens[k].t === '^') {
      const n = tokens[k].n; k++;
      const inv = n < 0 ? atom.slice().reverse().map(g => -g) : atom, out = [];
      for (let r = 0; r < Math.abs(n); r++) out.push(...inv);
      atom = out;
    }
    return atom;
  }
  let word;
  try { word = expr(); if (k < tokens.length) throw new Error('unbalanced parentheses'); } catch (e) { return { error: e.message }; }
  if (word.some(g => g === 0)) return { error: 'σ0 is not a generator' };
  if (word.length > 400) return { error: 'braid words up to 400 crossings' };
  return { word, text: s };
}
Seifert.parseBraid = parseBraid;
const wordText = word => word.join(' ');
Seifert.wordText = wordText;
Seifert.wordTeX = word => word.length ? word.map(g => `\\sigma_{${Math.abs(g)}}${g < 0 ? '^{-1}' : ''}`).join('') : '1';

// ------------------------------------------------------------------ the closed braid
// n strands (max |σ_i| + 1), c crossings, the permutation, the components (cycles), Euler characteristic n - c,
// genus (c - n + 2 - mu)/2 of the Bennequin surface, and whether that surface is connected (every σ_i occurs).
function braidData(word) {
  const c = word.length, n = c ? Math.max(...word.map(Math.abs)) + 1 : 1;
  const perm = [...Array(n).keys()];                   // perm[p] = the position where the strand starting at p ends
  const pos = [...Array(n).keys()];                    // pos[strand] = current position
  const at = [...Array(n).keys()];                     // at[position] = strand
  for (const g of word) { const i = Math.abs(g) - 1, a = at[i], b = at[i + 1]; at[i] = b; at[i + 1] = a; pos[a] = i + 1; pos[b] = i; }
  for (let s = 0; s < n; s++) perm[s] = pos[s];
  const comp = new Array(n).fill(-1); let mu = 0;
  for (let s = 0; s < n; s++) if (comp[s] < 0) { let t = s; while (comp[t] < 0) { comp[t] = mu; t = perm[t]; } mu++; }
  const used = new Set(word.map(Math.abs)); let connected = true;
  for (let i = 1; i < n; i++) if (!used.has(i)) connected = false;
  const chi = n - c, genus = (c - n + 2 - mu) / 2;
  return { n, c, perm, comp, mu, chi, genus, connected, positive: word.every(g => g > 0), negative: word.every(g => g < 0) };
}
Seifert.braidData = braidData;

// ------------------------------------------------------------------ Seifert matrix (Collins 2013 §3.3, as in Sage)
// The generators of H_1 of the Bennequin surface: for each crossing j whose generator σ_i occurs again at a later
// crossing h(j), the loop up band j, along disk i+1, down band h(j), back along disk i.  Sage's Link.seifert_matrix
// verbatim; for a split braid (some σ_i missing) Sage reorders the letters first, we simply refuse.
function seifertMatrix(word) {
  const x = word, len = x.length;
  const h = [];
  for (let j = 0; j < len - 1; j++) { const a = Math.abs(x[j]); let hj = 0; for (let i = j + 1; i < len; i++) if (Math.abs(x[i]) === a) { hj = i; break; } h.push(hj); }
  const indices = []; for (let i = 0; i < h.length; i++) if (h[i]) indices.push(i);
  const N = indices.length, A = Array.from({ length: N }, () => new Array(N).fill(0));
  const sign = v => v > 0 ? 1 : v < 0 ? -1 : 0;
  for (let ni = 0; ni < N; ni++) {
    const i = indices[ni], hi = h[i];
    A[ni][ni] = -sign(x[i] + x[hi]);
    for (let nj = ni + 1; nj < N; nj++) {
      const j = indices[nj];
      if (hi > h[j] || hi < j) continue;
      if (hi === j) { if (x[j] > 0) A[nj][ni] = 1; else A[ni][nj] = -1; }
      else if (Math.abs(x[i]) - Math.abs(x[j]) === 1) A[nj][ni] = -1;
      else if (Math.abs(x[j]) - Math.abs(x[i]) === 1) A[ni][nj] = 1;
    }
  }
  return A;
}
Seifert.seifertMatrix = seifertMatrix;

// ---- exact integer linear algebra (BigInt): Bareiss determinant, and det(V - t V^T) by evaluation and interpolation
function detBigInt(M) {                                // fraction-free Gaussian elimination; M is an array of BigInt rows (copied)
  const n = M.length; if (!n) return 1n;
  const a = M.map(r => r.slice()); let sgn = 1n, prev = 1n;
  for (let k = 0; k < n - 1; k++) {
    if (a[k][k] === 0n) {
      let p = -1; for (let i = k + 1; i < n; i++) if (a[i][k] !== 0n) { p = i; break; }
      if (p < 0) return 0n;
      [a[k], a[p]] = [a[p], a[k]]; sgn = -sgn;
    }
    for (let i = k + 1; i < n; i++) {
      for (let j = k + 1; j < n; j++) a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) / prev;
      a[i][k] = 0n;
    }
    prev = a[k][k];
  }
  return sgn * a[n - 1][n - 1];
}
Seifert.detBigInt = detBigInt;
const bigAbs = v => v < 0n ? -v : v;
function bigGcd(a, b) { a = bigAbs(a); b = bigAbs(b); while (b) { [a, b] = [b, a % b]; } return a; }
// det(V - t V^T) as integer coefficients [a_0, ..., a_N] (index = degree), by evaluating at t = 0..N and interpolating
function alexanderCoefficients(V) {
  const N = V.length; if (!N) return [1n];
  const ys = [];
  for (let t = 0; t <= N; t++) {
    const M = V.map((row, i) => row.map((v, j) => BigInt(v) - BigInt(t) * BigInt(V[j][i])));
    ys.push(detBigInt(M));
  }
  // Lagrange: p(t) = Σ_i y_i Π_{j≠i} (t - j) / D_i, D_i = Π_{j≠i} (i - j); accumulate over the common denominator L
  const D = []; for (let i = 0; i <= N; i++) { let d = 1n; for (let j = 0; j <= N; j++) if (j !== i) d *= BigInt(i - j); D.push(d); }
  let L = 1n; for (const d of D) L = L / bigGcd(L, d) * bigAbs(d);
  const acc = new Array(N + 1).fill(0n);
  for (let i = 0; i <= N; i++) {
    if (ys[i] === 0n) continue;
    let poly = [1n];                                   // Π_{j≠i} (t - j)
    for (let j = 0; j <= N; j++) if (j !== i) { const next = new Array(poly.length + 1).fill(0n); for (let k = 0; k < poly.length; k++) { next[k + 1] += poly[k]; next[k] -= BigInt(j) * poly[k]; } poly = next; }
    const f = ys[i] * (L / D[i]);
    for (let k = 0; k <= N; k++) acc[k] += f * poly[k];
  }
  return acc.map(v => v / L);
}
// The Alexander polynomial of the closed braid from its Seifert matrix: det(V - t V^T), then made a polynomial of
// lowest degree 0 with a positive leading coefficient (as KnotInfo prints it).  Returns { coeffs (Number[], index =
// degree), degree, det (|Δ(-1)|), raw (BigInt[]) }.  Zero for a split link.
function alexander(word) {
  const V = seifertMatrix(word), raw = alexanderCoefficients(V);
  let lo = 0; while (lo < raw.length && raw[lo] === 0n) lo++;
  if (lo === raw.length) return { coeffs: [0], degree: 0, det: 0, raw, zero: true };
  let hi = raw.length - 1; while (raw[hi] === 0n) hi--;
  let c = raw.slice(lo, hi + 1); if (c[c.length - 1] < 0n) c = c.map(v => -v);
  let det = 0n; for (let k = 0; k < c.length; k++) det += (k % 2 ? -c[k] : c[k]);
  return { coeffs: c.map(Number), degree: c.length - 1, det: Number(bigAbs(det)), raw };
}
Seifert.alexander = alexander;
Seifert.alexanderTeX = function (coeffs, variable = 't') {
  const parts = [];
  for (let k = coeffs.length - 1; k >= 0; k--) {
    const a = coeffs[k]; if (!a) continue;
    const mag = Math.abs(a), mono = k === 0 ? '' : k === 1 ? variable : `${variable}^{${k}}`;
    const coef = (mag === 1 && k > 0) ? '' : String(mag);
    parts.push(parts.length ? (a < 0 ? ' - ' : ' + ') + coef + mono : (a < 0 ? '-' : '') + coef + mono);
  }
  return parts.length ? parts.join('') : '0';
};

// ---- signature of the symmetrised Seifert form V + V^T (Jacobi eigenvalues; the matrix is small)
function signature(word) {
  const V = seifertMatrix(word), N = V.length;
  if (!N) return { signature: 0, nullity: 0, rank: 0 };
  const S = V.map((row, i) => row.map((v, j) => v + V[j][i]));
  let norm = 0; for (const row of S) for (const v of row) norm = Math.max(norm, Math.abs(v));
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0; for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) off += S[i][j] * S[i][j];
    if (off < 1e-24 * Math.max(1, norm * norm)) break;
    for (let p = 0; p < N; p++) for (let q = p + 1; q < N; q++) {
      if (Math.abs(S[p][q]) < 1e-300) continue;
      const theta = (S[q][q] - S[p][p]) / (2 * S[p][q]), t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const cs = 1 / Math.sqrt(t * t + 1), sn = t * cs;
      for (let k = 0; k < N; k++) { const a = S[k][p], b = S[k][q]; S[k][p] = cs * a - sn * b; S[k][q] = sn * a + cs * b; }
      for (let k = 0; k < N; k++) { const a = S[p][k], b = S[q][k]; S[p][k] = cs * a - sn * b; S[q][k] = sn * a + cs * b; }
    }
  }
  let pos = 0, neg = 0; const tol = 1e-9 * Math.max(1, norm);
  for (let i = 0; i < N; i++) { if (S[i][i] > tol) pos++; else if (S[i][i] < -tol) neg++; }
  return { signature: pos - neg, nullity: N - pos - neg, rank: pos + neg };
}
Seifert.signature = signature;

// ------------------------------------------------------------------ the Bennequin surface as a triangle mesh
// Disk p (p = 0..n-1) is horizontal at height z_p = (p - (n-1)/2) h, radius R, normal +z, triangulated by rings
// whose point counts grow with the radius; the rim has M = c q points, so every band sector has q of them.
// Crossing k (between disks i and i+1, sign ε) is a band attached to q_b + 1 rim points of each disk, centred on
// θ_k, of width W = q_b × (rim segment) ≈ bandWidth × h, parametrised by u ∈ [0,1] (height) and v ∈ [-1,1] (across):
//     angle  = θ_k + v (w/2) cos φ(u),   radius = R + b sin(πu) + v (W/2) sin φ(u),   z = z_i + u h,
// with φ(u) = TWIST ε π (3u² - 2u³) a smooth half turn, w = W/R, and b the band's outward bulge.  At u = 0 this is
// the rim arc of disk i, at u = 1 the rim arc of disk i+1 traversed backwards, at u = ½ a radial segment, so the
// two edges of the band swap heights over the sector and pass each other radially, one in front of the other: a
// braid crossing seen from outside the cylinder.  TWIST = +1 makes ε = +1 the positive crossing (checked by the
// linking number of the Hopf link in the tests).  The mesh is consistently oriented, every rim edge not under a
// band and every band edge is boundary, and the boundary loops are the closed braid.  Finally the boundary
// polygon is rounded (binomial smoothing with radius `round` × R) and resampled at equal arclength: this is the
// wire, and it is fixed from then on.
Seifert.TWIST = 1;
const DEFAULTS = { R: 1, spacing: 0.45, bandWidth: 1.2, bulge: 0.7, angular: 144, round: 0.1 };
Seifert.DEFAULTS = DEFAULTS;
function buildSurface(word, opts = {}) {
  const o = Object.assign({}, DEFAULTS, opts);
  const bd = braidData(word), n = bd.n, c = bd.c, R = o.R, h = o.spacing;
  const q = c ? Math.max(6, 2 * Math.round(o.angular / (2 * c))) : 0, M = c ? c * q : Math.max(24, o.angular);
  const seg = 2 * Math.PI * R / M;                     // rim segment length: the target edge length everywhere
  const qb = c ? Math.max(2, Math.min(q - 2, Math.round(Math.min(o.bandWidth * h, 0.8 * 2 * Math.PI * R / c) / seg))) : 0;
  const rings = Math.max(3, Math.round(R / seg));
  const W = qb * seg, w = W / R, b = o.bulge * W;
  const bandLen = Math.hypot(h, 2 * b) + 0.5 * W;      // rough length of a band edge, for the level count
  const levels = Math.max(2, Math.round(bandLen / seg));
  const pos = [], kind = [];
  const P = (x, y, z, k) => { pos.push(x, y, z); kind.push(k); return pos.length / 3 - 1; };
  const zOf = p => (p - (n - 1) / 2) * h;
  // disks
  const rim = [];                                      // rim[p][m] = vertex index of disk p at angle 2π m / M
  const tri = [];
  for (let p = 0; p < n; p++) {
    const z = zOf(p), center = P(0, 0, z, p), ringIdx = [];
    for (let j = 1; j <= rings; j++) {
      const cnt = j === rings ? M : Math.max(6, Math.round(M * j / rings)), r = R * j / rings, start = pos.length / 3;
      for (let m = 0; m < cnt; m++) { const a = 2 * Math.PI * m / cnt; P(r * Math.cos(a), r * Math.sin(a), z, p); }
      ringIdx.push({ start, cnt });
    }
    // fan around the centre, then the strips between rings (merge by angle), all counterclockwise from +z
    const r1 = ringIdx[0];
    for (let m = 0; m < r1.cnt; m++) tri.push(center, r1.start + m, r1.start + (m + 1) % r1.cnt);
    for (let j = 1; j < rings; j++) {
      const A = ringIdx[j - 1], B = ringIdx[j];
      let i = 0, k = 0;
      while (i < A.cnt || k < B.cnt) {
        const ai = (i + 1) / A.cnt, bk = (k + 1) / B.cnt;
        if (i < A.cnt && (k >= B.cnt || ai <= bk)) { tri.push(A.start + i, B.start + (k % B.cnt), A.start + (i + 1) % A.cnt); i++; }
        else { tri.push(A.start + (i % A.cnt), B.start + k, B.start + (k + 1) % B.cnt); k++; }
      }
    }
    const last = ringIdx[rings - 1]; rim.push([...Array(M).keys()].map(m => last.start + m));
  }
  // bands
  const bands = [];
  for (let k = 0; k < c; k++) {
    const g = word[k], i = Math.abs(g) - 1, eps = g > 0 ? 1 : -1, zi = zOf(i);
    const m0 = Math.round(k * q + (q - qb) / 2), theta = 2 * Math.PI * (m0 + qb / 2) / M;
    const V = (l, s) => {
      if (l === 0) return rim[i][(m0 + s) % M];
      if (l === levels) return rim[i + 1][(m0 + qb - s) % M];
      return bands[k].start + (l - 1) * (qb + 1) + s;
    };
    const start = pos.length / 3;
    for (let l = 1; l < levels; l++) {
      const u = l / levels, phi = Seifert.TWIST * eps * Math.PI * (3 * u * u - 2 * u * u * u), cp = Math.cos(phi), sp = Math.sin(phi);
      for (let s = 0; s <= qb; s++) {
        const v = -1 + 2 * s / qb, ang = theta + v * (w / 2) * cp, rad = R + b * Math.sin(Math.PI * u) + v * (W / 2) * sp;
        P(rad * Math.cos(ang), rad * Math.sin(ang), zi + u * h, n + k);
      }
    }
    bands.push({ start, i, eps, theta });
    for (let l = 0; l < levels; l++) for (let s = 0; s < qb; s++) {
      const A = V(l, s), B = V(l + 1, s), C = V(l + 1, s + 1), D = V(l, s + 1);
      tri.push(A, B, C, A, C, D);
    }
  }
  const mesh = { pos: Float64Array.from(pos), tri: Uint32Array.from(tri), kind: Int16Array.from(kind), n, c, braid: bd,
                 params: { R, h, q, M, qb, rings, levels, w, W, b, seg }, word: word.slice() };
  mesh.loops = boundaryLoops(mesh.tri, mesh.pos.length / 3);
  mesh.fixed = new Uint8Array(mesh.pos.length / 3);
  for (const loop of mesh.loops) for (const v of loop) mesh.fixed[v] = 1;
  const passes = o.smoothPasses !== undefined ? o.smoothPasses : Math.round(2 * Math.pow(o.round * R / seg, 2));
  mesh.params.passes = passes;
  if (passes > 0) for (const loop of mesh.loops) { smoothLoop(mesh.pos, loop, passes); resampleLoop(mesh.pos, loop); }
  mesh.initial = mesh.pos.slice();
  return mesh;
}
Seifert.buildSurface = buildSurface;

// The boundary of an oriented triangle mesh as loops of vertex indices, each traversed as its triangles induce
// (edge a->b of a triangle with no partner b->a is boundary, and the loop follows a->b).
function boundaryLoops(tri, V) {
  const seen = new Map();                              // directed edge a->b => true
  for (let f = 0; f < tri.length; f += 3) for (let e = 0; e < 3; e++) seen.set(tri[f + e] * V + tri[f + (e + 1) % 3], 1);
  const next = new Int32Array(V).fill(-1); let count = 0;
  for (let f = 0; f < tri.length; f += 3) for (let e = 0; e < 3; e++) {
    const a = tri[f + e], b = tri[f + (e + 1) % 3];
    if (!seen.has(b * V + a)) { if (next[a] !== -1) throw new Error('non-manifold boundary at vertex ' + a); next[a] = b; count++; }
  }
  const loops = [], done = new Uint8Array(V);
  for (let s = 0; s < V; s++) {
    if (next[s] < 0 || done[s]) continue;
    const loop = []; let v = s;
    while (!done[v]) { done[v] = 1; loop.push(v); v = next[v]; if (v < 0) throw new Error('open boundary chain'); }
    loops.push(loop);
  }
  return loops;
}
Seifert.boundaryLoops = boundaryLoops;
// Binomial smoothing of a closed polygon in place, p_i <- (p_{i-1} + 2 p_i + p_{i+1}) / 4, `passes` times: rounds the
// corners where a rim arc turns into a band edge, at the scale of a few segments.
function smoothLoop(pos, loop, passes) {
  const L = loop.length; if (L < 4) return;
  const tmp = new Float64Array(3 * L);
  for (let it = 0; it < passes; it++) {
    for (let k = 0; k < L; k++) {
      const a = loop[(k + L - 1) % L], b = loop[k], c = loop[(k + 1) % L];
      for (let d = 0; d < 3; d++) tmp[3 * k + d] = 0.25 * (pos[3 * a + d] + 2 * pos[3 * b + d] + pos[3 * c + d]);
    }
    for (let k = 0; k < L; k++) for (let d = 0; d < 3; d++) pos[3 * loop[k] + d] = tmp[3 * k + d];
  }
}
Seifert.smoothLoop = smoothLoop;
// Move the vertices of a closed polygon along it to equal arclength spacing (the same number of vertices), so
// the smoothing leaves no bunched-up corners.
function resampleLoop(pos, loop) {
  const L = loop.length; if (L < 4) return;
  const P = loop.map(v => [pos[3 * v], pos[3 * v + 1], pos[3 * v + 2]]), cum = [0];
  for (let k = 0; k < L; k++) { const a = P[k], b = P[(k + 1) % L]; cum.push(cum[k] + Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2])); }
  const total = cum[L]; let j = 0;
  for (let k = 0; k < L; k++) {
    const target = total * k / L;
    while (j < L - 1 && cum[j + 1] < target) j++;
    const a = P[j], b = P[(j + 1) % L], t = cum[j + 1] > cum[j] ? (target - cum[j]) / (cum[j + 1] - cum[j]) : 0;
    for (let d = 0; d < 3; d++) pos[3 * loop[k] + d] = a[d] + t * (b[d] - a[d]);
  }
}
Seifert.resampleLoop = resampleLoop;

if (typeof module !== 'undefined' && module.exports) module.exports = Seifert; else window.Seifert = Seifert;
})();
