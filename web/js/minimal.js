/* minimal.js -- discrete minimal surfaces with a fixed boundary, on a triangle mesh.  Node-compatible and browser.

   A surface X is minimal iff its coordinate functions are harmonic for its own Laplace-Beltrami operator,
   Δ_X X = 2H n = 0: the elliptic Dirichlet problem with the wire as boundary data, quasilinear because the metric
   depends on X.  Pinkall and Polthier (Experiment. Math. 2, 1993) linearise it on a mesh: freeze the metric of the
   current mesh, solve the cotangent Laplace equation L x = 0 for the interior vertices with the boundary fixed, move
   the vertices there, repeat.  Every step decreases the area, and the fixed points are the discrete minimal
   surfaces (each vertex at the cotangent-weighted mean of its neighbours, the discrete H = 0).  Between the
   harmonic step and no step lies implicit mean curvature flow, (M + dt L) x_new = M x_old (Desbrun, Meyer,
   Schröder, Barr 1999), which the same solver does with a finite dt.  The sparse SPD systems are solved by
   conjugate gradients with a Jacobi preconditioner, warm-started.  Edge flips to the intrinsic Delaunay
   triangulation and tangential smoothing keep the mesh usable while the vertices drift, as in Brakke's Surface
   Evolver (equiangulation and vertex averaging). */
(function () {
'use strict';
const Minimal = {};

// ------------------------------------------------------------------ topology
// edges: for each undirected edge (a < b), the two opposite vertices (or -1); a CSR neighbour list per vertex with
// the edge index of each neighbour; the edge ids of each triangle's three edges (ab, bc, ca).
function topology(tri, V) {
  const F = tri.length / 3, key = (a, b) => a < b ? a * V + b : b * V + a;
  const map = new Map(), ea = [], eb = [], oppA = [], oppB = [], triEdge = new Int32Array(3 * F);
  for (let f = 0; f < F; f++) for (let e = 0; e < 3; e++) {
    const a = tri[3 * f + e], b = tri[3 * f + (e + 1) % 3], c = tri[3 * f + (e + 2) % 3], k = key(a, b);
    let id = map.get(k);
    if (id === undefined) { id = ea.length; map.set(k, id); ea.push(Math.min(a, b)); eb.push(Math.max(a, b)); oppA.push(-1); oppB.push(-1); }
    if (oppA[id] < 0) oppA[id] = c; else if (oppB[id] < 0) oppB[id] = c; else throw new Error('edge with three triangles');
    triEdge[3 * f + e] = id;
  }
  const E = ea.length, deg = new Int32Array(V);
  for (let e = 0; e < E; e++) { deg[ea[e]]++; deg[eb[e]]++; }
  const start = new Int32Array(V + 1); for (let v = 0; v < V; v++) start[v + 1] = start[v] + deg[v];
  const nbr = new Int32Array(2 * E), nbrEdge = new Int32Array(2 * E), fill = start.slice(0, V);
  for (let e = 0; e < E; e++) { const a = ea[e], b = eb[e]; nbr[fill[a]] = b; nbrEdge[fill[a]++] = e; nbr[fill[b]] = a; nbrEdge[fill[b]++] = e; }
  const boundaryEdge = new Uint8Array(E); for (let e = 0; e < E; e++) boundaryEdge[e] = oppB[e] < 0 ? 1 : 0;
  return { V, F, E, ea: Int32Array.from(ea), eb: Int32Array.from(eb), oppA: Int32Array.from(oppA), oppB: Int32Array.from(oppB), triEdge, start, nbr, nbrEdge, boundaryEdge, edgeMap: map };
}
Minimal.topology = topology;
// The mesh object: { pos, tri, fixed } plus the topology, rebuilt after every change of connectivity.
function prepare(mesh) { mesh.topo = topology(mesh.tri, mesh.pos.length / 3); return mesh; }
Minimal.prepare = prepare;

// ------------------------------------------------------------------ geometry
function triArea(p, a, b, c) {
  const ux = p[3 * b] - p[3 * a], uy = p[3 * b + 1] - p[3 * a + 1], uz = p[3 * b + 2] - p[3 * a + 2];
  const vx = p[3 * c] - p[3 * a], vy = p[3 * c + 1] - p[3 * a + 1], vz = p[3 * c + 2] - p[3 * a + 2];
  const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
  return 0.5 * Math.sqrt(nx * nx + ny * ny + nz * nz);
}
function area(mesh) { let A = 0; const t = mesh.tri; for (let f = 0; f < t.length; f += 3) A += triArea(mesh.pos, t[f], t[f + 1], t[f + 2]); return A; }
Minimal.area = area;
function meanEdgeLength(mesh) {
  const T = mesh.topo || topology(mesh.tri, mesh.pos.length / 3), p = mesh.pos; let s = 0;
  for (let e = 0; e < T.E; e++) { const a = T.ea[e], b = T.eb[e]; s += Math.hypot(p[3 * a] - p[3 * b], p[3 * a + 1] - p[3 * b + 1], p[3 * a + 2] - p[3 * b + 2]); }
  return s / T.E;
}
Minimal.meanEdgeLength = meanEdgeLength;
// cotangent weights w_e = (cot α + cot β)/2 per edge, and the lumped (barycentric) vertex areas
function weights(mesh) {
  const T = mesh.topo, p = mesh.pos, w = new Float64Array(T.E), mass = new Float64Array(T.V);
  for (let f = 0; f < T.F; f++) {
    const a = mesh.tri[3 * f], b = mesh.tri[3 * f + 1], c = mesh.tri[3 * f + 2];
    const cot = (i, j, k) => {                       // cot of the angle at i in triangle (i, j, k)
      const ux = p[3 * j] - p[3 * i], uy = p[3 * j + 1] - p[3 * i + 1], uz = p[3 * j + 2] - p[3 * i + 2];
      const vx = p[3 * k] - p[3 * i], vy = p[3 * k + 1] - p[3 * i + 1], vz = p[3 * k + 2] - p[3 * i + 2];
      const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx, cr = Math.sqrt(nx * nx + ny * ny + nz * nz);
      return cr > 1e-300 ? (ux * vx + uy * vy + uz * vz) / cr : 0;
    };
    w[T.triEdge[3 * f]] += 0.5 * cot(c, a, b);       // edge ab, opposite c
    w[T.triEdge[3 * f + 1]] += 0.5 * cot(a, b, c);   // edge bc, opposite a
    w[T.triEdge[3 * f + 2]] += 0.5 * cot(b, c, a);   // edge ca, opposite b
    const A = triArea(p, a, b, c) / 3; mass[a] += A; mass[b] += A; mass[c] += A;
  }
  return { w, mass };
}
Minimal.weights = weights;
// the mean curvature vector per vertex, (1/A_i) Σ_j w_ij (x_i - x_j) = 2 H n, and |H| for the interior vertices
function meanCurvature(mesh) {
  const T = mesh.topo, p = mesh.pos, { w, mass } = weights(mesh), H = new Float64Array(T.V);
  let max = 0, sum2 = 0, cnt = 0;
  for (let i = 0; i < T.V; i++) {
    if (mesh.fixed[i]) continue;
    let x = 0, y = 0, z = 0;
    for (let k = T.start[i]; k < T.start[i + 1]; k++) { const j = T.nbr[k], wij = w[T.nbrEdge[k]]; x += wij * (p[3 * i] - p[3 * j]); y += wij * (p[3 * i + 1] - p[3 * j + 1]); z += wij * (p[3 * i + 2] - p[3 * j + 2]); }
    const h = Math.sqrt(x * x + y * y + z * z) / (2 * Math.max(mass[i], 1e-300));
    H[i] = h; max = Math.max(max, h); sum2 += h * h; cnt++;
  }
  return { H, max, rms: cnt ? Math.sqrt(sum2 / cnt) : 0 };
}
Minimal.meanCurvature = meanCurvature;

// ------------------------------------------------------------------ the elliptic solve
// Solve (M/dt + L) x = M x0/dt + (boundary terms) for the interior vertices, x prescribed on the fixed vertices.
// `target` (3V) holds the prescribed values on the fixed vertices and, on the interior, x0 and the warm start.
// dt = Infinity is the harmonic problem L x = 0.  Returns { x (3V: the solution on the interior, target on the
// boundary), iterations, residual }.
function solve(mesh, target, opts = {}) {
  const dt = opts.dt === undefined ? Infinity : opts.dt, invdt = dt === Infinity ? 0 : 1 / dt;
  const T = mesh.topo, V = T.V, { w, mass } = weights(mesh), out = Float64Array.from(target);
  const idx = new Int32Array(V).fill(-1); let n = 0;
  for (let v = 0; v < V; v++) if (!mesh.fixed[v]) idx[v] = n++;
  if (!n) return { x: out, iterations: 0, residual: 0 };
  const diag = new Float64Array(n), verts = new Int32Array(n);
  for (let v = 0; v < V; v++) if (idx[v] >= 0) {
    verts[idx[v]] = v; let s = mass[v] * invdt;
    for (let k = T.start[v]; k < T.start[v + 1]; k++) s += w[T.nbrEdge[k]];
    diag[idx[v]] = s;
  }
  const matvec = (x, res) => {
    for (let r = 0; r < n; r++) {
      const v = verts[r]; let s = diag[r] * x[r];
      for (let k = T.start[v]; k < T.start[v + 1]; k++) { const j = T.nbr[k]; if (idx[j] >= 0) s -= w[T.nbrEdge[k]] * x[idx[j]]; }
      res[r] = s;
    }
  };
  const tol = opts.tol || 1e-9, maxIter = opts.maxIter || 2000;
  let iterations = 0, residual = 0;
  const x = new Float64Array(n), b = new Float64Array(n), r = new Float64Array(n), z = new Float64Array(n), q = new Float64Array(n), d = new Float64Array(n);
  const pre = new Float64Array(n); for (let i = 0; i < n; i++) pre[i] = 1 / Math.max(diag[i], 1e-12);
  for (let c = 0; c < 3; c++) {
    for (let i = 0; i < n; i++) {
      const v = verts[i]; x[i] = target[3 * v + c]; let s = mass[v] * invdt * target[3 * v + c];
      for (let k = T.start[v]; k < T.start[v + 1]; k++) { const j = T.nbr[k]; if (idx[j] < 0) s += w[T.nbrEdge[k]] * target[3 * j + c]; }
      b[i] = s;
    }
    // preconditioned conjugate gradients from the warm start
    matvec(x, q); let rz = 0, bnorm = 0;
    for (let i = 0; i < n; i++) { r[i] = b[i] - q[i]; z[i] = pre[i] * r[i]; d[i] = z[i]; rz += r[i] * z[i]; bnorm += b[i] * b[i]; }
    bnorm = Math.sqrt(bnorm) || 1; let it = 0, rn = 0;
    for (; it < maxIter; it++) {
      rn = 0; for (let i = 0; i < n; i++) rn += r[i] * r[i]; rn = Math.sqrt(rn);
      if (rn <= tol * bnorm) break;
      matvec(d, q); let dq = 0; for (let i = 0; i < n; i++) dq += d[i] * q[i];
      if (!(dq > 0)) break;                          // not SPD numerically: stop with what we have
      const alpha = rz / dq; let rz2 = 0;
      for (let i = 0; i < n; i++) { x[i] += alpha * d[i]; r[i] -= alpha * q[i]; z[i] = pre[i] * r[i]; rz2 += r[i] * z[i]; }
      const beta = rz2 / rz; rz = rz2;
      for (let i = 0; i < n; i++) d[i] = z[i] + beta * d[i];
    }
    iterations = Math.max(iterations, it); residual = Math.max(residual, rn / bnorm);
    for (let i = 0; i < n; i++) out[3 * verts[i] + c] = x[i];
  }
  return { x: out, iterations, residual };
}
Minimal.solve = solve;
// One step of the surface with the wire fixed: the harmonic step (dt = Infinity) or implicit mean curvature flow.
// Returns { iterations, residual, moved } (the largest vertex displacement).
function step(mesh, opts = {}) {
  const r = solve(mesh, mesh.pos, opts), p = mesh.pos; let moved = 0;
  for (let i = 0; i < p.length; i++) { moved = Math.max(moved, Math.abs(r.x[i] - p[i])); p[i] = r.x[i]; }
  return { iterations: r.iterations, residual: r.residual, moved };
}
Minimal.step = step;
// Carry the surface along with a displacement of the wire, as a rubber sheet: the harmonic extension of the
// boundary displacement `disp` (3V, read on the fixed vertices) is added to every vertex.
function extend(mesh, disp, opts = {}) {
  const target = new Float64Array(disp.length);
  for (let v = 0; v < mesh.fixed.length; v++) if (mesh.fixed[v]) for (let c = 0; c < 3; c++) target[3 * v + c] = disp[3 * v + c];
  const r = solve(mesh, target, Object.assign({ dt: Infinity }, opts)), p = mesh.pos; let moved = 0;
  for (let i = 0; i < p.length; i++) { moved = Math.max(moved, Math.abs(r.x[i])); p[i] += r.x[i]; }
  return { iterations: r.iterations, residual: r.residual, moved };
}
Minimal.extend = extend;

// ------------------------------------------------------------------ mesh maintenance
// Edge flips towards the intrinsic Delaunay triangulation: an interior edge whose two opposite angles sum to more
// than π is flipped, unless the new edge exists already or the flip would fold the quad.  Passes until no flips.
function flipDelaunay(mesh, maxPasses = 10) {
  let total = 0;
  for (let pass = 0; pass < maxPasses; pass++) {
    const T = mesh.topo, p = mesh.pos, tri = mesh.tri, V = T.V;
    const angleAt = (i, j, k) => {                   // angle at i in triangle (i, j, k)
      const ux = p[3 * j] - p[3 * i], uy = p[3 * j + 1] - p[3 * i + 1], uz = p[3 * j + 2] - p[3 * i + 2];
      const vx = p[3 * k] - p[3 * i], vy = p[3 * k + 1] - p[3 * i + 1], vz = p[3 * k + 2] - p[3 * i + 2];
      const dot = ux * vx + uy * vy + uz * vz, nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
      return Math.atan2(Math.sqrt(nx * nx + ny * ny + nz * nz), dot);
    };
    // which triangle holds each directed edge a->b (there is one, the mesh being oriented)
    const faceOf = new Map();
    for (let f = 0; f < T.F; f++) for (let e = 0; e < 3; e++) faceOf.set(tri[3 * f + e] * V + tri[3 * f + (e + 1) % 3], f);
    const touched = new Uint8Array(T.F); let flips = 0;
    for (let e = 0; e < T.E; e++) {
      if (T.boundaryEdge[e]) continue;
      const a = T.ea[e], b = T.eb[e], c = T.oppA[e], d = T.oppB[e];
      if (angleAt(c, a, b) + angleAt(d, a, b) <= Math.PI + 1e-6) continue;
      if (T.edgeMap.has(c < d ? c * V + d : d * V + c)) continue;
      const f1 = faceOf.get(a * V + b), f2 = faceOf.get(b * V + a);
      if (f1 === undefined || f2 === undefined || touched[f1] || touched[f2]) continue;
      // f1 holds a->b with third vertex c1, f2 holds b->a with third vertex d1 (c1, d1 are c, d in some order)
      const c1 = tri[3 * f1] + tri[3 * f1 + 1] + tri[3 * f1 + 2] - a - b, d1 = tri[3 * f2] + tri[3 * f2 + 1] + tri[3 * f2 + 2] - a - b;
      // the new triangles (a, d1, c1) and (b, c1, d1) must not fold: their normals should agree with the old ones
      const nrm = (i, j, k) => { const ux = p[3 * j] - p[3 * i], uy = p[3 * j + 1] - p[3 * i + 1], uz = p[3 * j + 2] - p[3 * i + 2], vx = p[3 * k] - p[3 * i], vy = p[3 * k + 1] - p[3 * i + 1], vz = p[3 * k + 2] - p[3 * i + 2]; return [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]; };
      const o1 = nrm(a, b, c1), o2 = nrm(b, a, d1), n1 = nrm(a, d1, c1), n2 = nrm(b, c1, d1), old = [o1[0] + o2[0], o1[1] + o2[1], o1[2] + o2[2]];
      if (n1[0] * old[0] + n1[1] * old[1] + n1[2] * old[2] <= 0 || n2[0] * old[0] + n2[1] * old[1] + n2[2] * old[2] <= 0) continue;
      tri[3 * f1] = a; tri[3 * f1 + 1] = d1; tri[3 * f1 + 2] = c1;
      tri[3 * f2] = b; tri[3 * f2 + 1] = c1; tri[3 * f2 + 2] = d1;
      touched[f1] = touched[f2] = 1; flips++;
    }
    total += flips;
    if (!flips) break;
    prepare(mesh);
  }
  return total;
}
Minimal.flipDelaunay = flipDelaunay;
// Tangential smoothing: each interior vertex moves a fraction λ of the way to the centroid of its neighbours, with
// the normal component of that displacement removed, so the shape stays and the vertices spread out.
function tangentialSmooth(mesh, lambda = 0.5) {
  const T = mesh.topo, p = mesh.pos, V = T.V, nrm = new Float64Array(3 * V), disp = new Float64Array(3 * V);
  for (let f = 0; f < T.F; f++) {
    const a = mesh.tri[3 * f], b = mesh.tri[3 * f + 1], c = mesh.tri[3 * f + 2];
    const ux = p[3 * b] - p[3 * a], uy = p[3 * b + 1] - p[3 * a + 1], uz = p[3 * b + 2] - p[3 * a + 2];
    const vx = p[3 * c] - p[3 * a], vy = p[3 * c + 1] - p[3 * a + 1], vz = p[3 * c + 2] - p[3 * a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;   // area-weighted
    for (const v of [a, b, c]) { nrm[3 * v] += nx; nrm[3 * v + 1] += ny; nrm[3 * v + 2] += nz; }
  }
  for (let i = 0; i < V; i++) {
    if (mesh.fixed[i]) continue;
    let cx = 0, cy = 0, cz = 0, k0 = T.start[i], k1 = T.start[i + 1];
    if (k1 === k0) continue;
    for (let k = k0; k < k1; k++) { const j = T.nbr[k]; cx += p[3 * j]; cy += p[3 * j + 1]; cz += p[3 * j + 2]; }
    cx = cx / (k1 - k0) - p[3 * i]; cy = cy / (k1 - k0) - p[3 * i + 1]; cz = cz / (k1 - k0) - p[3 * i + 2];
    let nx = nrm[3 * i], ny = nrm[3 * i + 1], nz = nrm[3 * i + 2]; const nl = Math.hypot(nx, ny, nz);
    if (nl > 0) { nx /= nl; ny /= nl; nz /= nl; const dn = cx * nx + cy * ny + cz * nz; cx -= dn * nx; cy -= dn * ny; cz -= dn * nz; }
    disp[3 * i] = lambda * cx; disp[3 * i + 1] = lambda * cy; disp[3 * i + 2] = lambda * cz;
  }
  for (let i = 0; i < 3 * V; i++) p[i] += disp[i];
}
Minimal.tangentialSmooth = tangentialSmooth;
// one relaxation round: the solve, then the optional flips and tangential smoothing
function relax(mesh, opts = {}) {
  if (!mesh.topo) prepare(mesh);
  const s = step(mesh, opts);
  let flips = 0;
  if (opts.flips !== false) flips = flipDelaunay(mesh, opts.flipPasses || 3);
  if (opts.tangential) tangentialSmooth(mesh, typeof opts.tangential === 'number' ? opts.tangential : 0.5);
  return Object.assign(s, { flips, area: area(mesh) });
}
Minimal.relax = relax;

// ------------------------------------------------------------------ checks and diagnostics
function eulerCharacteristic(mesh) { const T = mesh.topo || topology(mesh.tri, mesh.pos.length / 3); return T.V - T.E + T.F; }
Minimal.eulerCharacteristic = eulerCharacteristic;
// every interior edge is traversed once in each direction, every boundary edge once: a consistently oriented surface
function orientable(tri, V) {
  const seen = new Map();
  for (let f = 0; f < tri.length; f += 3) for (let e = 0; e < 3; e++) { const k = tri[f + e] * V + tri[f + (e + 1) % 3]; if (seen.has(k)) return false; seen.set(k, 1); }
  return true;
}
Minimal.orientable = orientable;
function degenerateTriangles(mesh, tol = 1e-12) { let n = 0; const t = mesh.tri; for (let f = 0; f < t.length; f += 3) if (triArea(mesh.pos, t[f], t[f + 1], t[f + 2]) < tol) n++; return n; }
Minimal.degenerateTriangles = degenerateTriangles;
// Gauss linking number of two closed polygons (arrays of [x, y, z]) by counting signed crossings of their
// projections along a generic direction; the two directions used must agree.
function linkingNumber(A, B) {
  const dirs = [[0.21, 0.17, 0.96], [0.53, -0.31, 0.79]];
  const results = dirs.map(dz => {
    const dx = [dz[1], -dz[0], 0], l = Math.hypot(dx[0], dx[1]); dx[0] /= l; dx[1] /= l;
    const dy = [dz[1] * dx[2] - dz[2] * dx[1], dz[2] * dx[0] - dz[0] * dx[2], dz[0] * dx[1] - dz[1] * dx[0]];
    const proj = P => P.map(p => [p[0] * dx[0] + p[1] * dx[1] + p[2] * dx[2], p[0] * dy[0] + p[1] * dy[1] + p[2] * dy[2], p[0] * dz[0] + p[1] * dz[1] + p[2] * dz[2]]);
    const a = proj(A), b = proj(B); let s = 0;
    for (let i = 0; i < a.length; i++) {
      const p0 = a[i], p1 = a[(i + 1) % a.length], ux = p1[0] - p0[0], uy = p1[1] - p0[1];
      for (let j = 0; j < b.length; j++) {
        const q0 = b[j], q1 = b[(j + 1) % b.length], vx = q1[0] - q0[0], vy = q1[1] - q0[1];
        const den = ux * vy - uy * vx; if (Math.abs(den) < 1e-15) continue;
        const wx = q0[0] - p0[0], wy = q0[1] - p0[1];
        const t = (wx * vy - wy * vx) / den, u = (wx * uy - wy * ux) / den;
        if (t < 0 || t >= 1 || u < 0 || u >= 1) continue;
        const za = p0[2] + t * (p1[2] - p0[2]), zb = q0[2] + u * (q1[2] - q0[2]);
        s += (za > zb ? 1 : -1) * Math.sign(den);   // sign((over × under)·z): positive crossing = +1
      }
    }
    return s / 2;
  });
  if (Math.abs(results[0] - results[1]) > 1e-9) throw new Error('linking number: the two projections disagree');
  return results[0];
}
Minimal.linkingNumber = linkingNumber;
Minimal.loopPoints = (pos, loop) => loop.map(v => [pos[3 * v], pos[3 * v + 1], pos[3 * v + 2]]);

if (typeof module !== 'undefined' && module.exports) module.exports = Minimal; else window.Minimal = Minimal;
})();
