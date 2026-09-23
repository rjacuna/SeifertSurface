/* wire.js -- taming the wire: Scharein's KnotPlot relaxation of the knot, as van Wijk and Cohen adopt it
   (IEEE TVCG 2006, §5.3), applied to the boundary loops of the surface mesh.  Node-compatible and browser.

   Every point of the wire is a unit mass, attracted by its two neighbours along the loop with the generalised Hooke
   force F_a(r) = H r^(1+β) and repelled by every other point of the wire with the generalised electrostatic force
   F_r(r) = K r^(-(2+α)), distances measured in units of r_a, the initial spacing of the points (β = 1, α = 0:
   inverse-square repulsion; α = 4 is KnotPlot's "strong repulsion").  Damped explicit Euler, v <- (1 - γ) v + F dt,
   x <- x + v dt with the displacement clamped to d_max, and a move is refused when it would bring the point within
   d_close of a segment of the wire not adjacent to it: Scharein's argument that d_close > d_max then keeps the knot
   from passing through itself.  KnotPlot's d_close is a fraction of the spacing; here it is 2 spacings, since a
   film between two strands closer than a few mesh edges cannot be resolved and pinches.  The step dt decays by (1 - μ) each step, so a cycle settles; a new cycle restarts
   from dt0.  The wire grows while it relaxes (the repulsion wins over the springs until the spacing is about
   1.5 r_a), the model's own scale; nothing here depends on it.

   Two things are ours.  The dynamics runs on a coarse copy of each boundary loop, resampled at equal arclength with
   spacing `coarse` (0.1 disk radii, under the wire's smallest features), since the forces cost N² and the mesh
   boundary has several hundred points per loop; the mesh boundary is then interpolated from the coarse polygon by a
   Catmull-Rom spline at equal arclength.  And the
   surface is not given springs of its own: the app carries it along as a rubber sheet, by the harmonic extension of
   the wire's displacement (Minimal.extend), and settles the film on it. */
(function () {
'use strict';
const Wire = {};
Wire.DEFAULTS = { alpha: 0, beta: 1, H: 1, K: 1, gamma: 0.15, dt0: 0.2, decay: 0.0003, dmax: 0.25, dclose: 2, coarse: 0.1 };   // coarse: the spacing of the coarse polygon, in the mesh's units (R = 1)

// ---- closed polygons
function polyLength(P, n, off = 0) { let L = 0; for (let k = 0; k < n; k++) { const a = off + 3 * k, b = off + 3 * ((k + 1) % n); L += Math.hypot(P[b] - P[a], P[b + 1] - P[a + 1], P[b + 2] - P[a + 2]); } return L; }
// m points at equal arclength along the closed polygon (P, n), starting at its first vertex
function resample(P, n, m, off = 0) {
  const cum = [0]; for (let k = 0; k < n; k++) { const a = off + 3 * k, b = off + 3 * ((k + 1) % n); cum.push(cum[k] + Math.hypot(P[b] - P[a], P[b + 1] - P[a + 1], P[b + 2] - P[a + 2])); }
  const total = cum[n], out = new Float64Array(3 * m); let j = 0;
  for (let k = 0; k < m; k++) {
    const target = total * k / m; while (j < n - 1 && cum[j + 1] < target) j++;
    const a = off + 3 * j, b = off + 3 * ((j + 1) % n), t = cum[j + 1] > cum[j] ? (target - cum[j]) / (cum[j + 1] - cum[j]) : 0;
    for (let d = 0; d < 3; d++) out[3 * k + d] = P[a + d] + t * (P[b + d] - P[a + d]);
  }
  return out;
}
// the closed Catmull-Rom spline through (P, n), sampled `per` times per segment
function spline(P, n, per, off = 0) {
  const out = new Float64Array(3 * n * per);
  for (let i = 0; i < n; i++) {
    const p0 = off + 3 * ((i + n - 1) % n), p1 = off + 3 * i, p2 = off + 3 * ((i + 1) % n), p3 = off + 3 * ((i + 2) % n);
    for (let s = 0; s < per; s++) {
      const t = s / per, t2 = t * t, t3 = t2 * t;
      for (let d = 0; d < 3; d++) out[3 * (i * per + s) + d] = 0.5 * (2 * P[p1 + d] + (-P[p0 + d] + P[p2 + d]) * t + (2 * P[p0 + d] - 5 * P[p1 + d] + 4 * P[p2 + d] - P[p3 + d]) * t2 + (-P[p0 + d] + 3 * P[p1 + d] - 3 * P[p2 + d] + P[p3 + d]) * t3);
    }
  }
  return out;
}
Wire.resample = resample; Wire.spline = spline;

// The state of a relaxation.  Coarse: loops of points `coarse` apart (P, V, prev, next, loopOf).  Fine: for each
// loop the mesh vertices in order (vert) and the interpolated positions F.
function init(mesh, opts = {}) {
  const o = Object.assign({}, Wire.DEFAULTS, opts), loops = mesh.loops, cLoops = [], fLoops = [];
  const P = [], prev = [], next = [], loopOf = [];
  loops.forEach((L, l) => {
    const n = L.length, fine = new Float64Array(3 * n);
    for (let k = 0; k < n; k++) for (let d = 0; d < 3; d++) fine[3 * k + d] = mesh.pos[3 * L[k] + d];
    const m = Math.max(12, Math.min(n, Math.round(polyLength(fine, n) / o.coarse))), base = P.length / 3, Q = resample(fine, n, m);
    for (let k = 0; k < m; k++) { P.push(Q[3 * k], Q[3 * k + 1], Q[3 * k + 2]); prev.push(base + (k + m - 1) % m); next.push(base + (k + 1) % m); loopOf.push(l); }
    cLoops.push({ base, m }); fLoops.push({ vert: Int32Array.from(L), n });
  });
  const N = P.length / 3;
  const ws = { o, N, P: Float64Array.from(P), V: new Float64Array(3 * N), prev: Int32Array.from(prev), next: Int32Array.from(next), loopOf: Int32Array.from(loopOf),
               cLoops, fLoops, dt: o.dt0, steps: 0, energies: [] };
  ws.ra = length(ws) / N; ws.L0 = length(ws);
  ws.F = interpolate(ws);
  return ws;
}
Wire.init = init;
function length(ws) { return polyLength(ws.P, ws.N); }
Wire.length = length;
// The fine loops from the coarse ones: the spline, sampled at equal arclength with the loop's point count.  The
// coarse points slide along the loop as they relax (nothing in the forces pins the parametrisation), and the mesh
// boundary must not slide with them, or the surface next to it is sheared without end; so the first fine point is
// anchored at the point of the new curve nearest to where it was, and the others follow at equal arclength.
function interpolate(ws) {
  const out = [];
  ws.cLoops.forEach((c, l) => {
    const f = ws.fLoops[l], per = Math.max(2, Math.ceil(4 * f.n / c.m)), dense = spline(ws.P, c.m, per, 3 * c.base), nd = c.m * per;
    let start = 0;
    if (ws.F && ws.F[l]) {
      const x = ws.F[l][0], y = ws.F[l][1], z = ws.F[l][2]; let best = Infinity;
      for (let k = 0; k < nd; k++) { const d = (dense[3 * k] - x) ** 2 + (dense[3 * k + 1] - y) ** 2 + (dense[3 * k + 2] - z) ** 2; if (d < best) { best = d; start = k; } }
    }
    const rolled = start ? Float64Array.from([...dense.subarray(3 * start), ...dense.subarray(0, 3 * start)]) : dense;
    out.push(resample(rolled, nd, f.n));
  });
  return out;
}
Wire.interpolate = interpolate;
// distance from point p to segment ab
function pointSegment(px, py, pz, ax, ay, az, bx, by, bz) {
  const ux = bx - ax, uy = by - ay, uz = bz - az, l2 = ux * ux + uy * uy + uz * uz;
  let t = l2 > 0 ? ((px - ax) * ux + (py - ay) * uy + (pz - az) * uz) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
  return Math.hypot(px - ax - t * ux, py - ay - t * uy, pz - az - t * uz);
}
// the smallest distance from a point of the coarse wire to a segment not adjacent to it (in R), and the pair
function clearance(ws) {
  const { N, P, prev, next } = ws; let best = Infinity, pair = [-1, -1];
  for (let k = 0; k < N; k++) {
    const skip = new Set([k, prev[k], prev[prev[k]], next[k]]);
    for (let j = 0; j < N; j++) {
      if (skip.has(j)) continue;
      const n = next[j], d = pointSegment(P[3 * k], P[3 * k + 1], P[3 * k + 2], P[3 * j], P[3 * j + 1], P[3 * j + 2], P[3 * n], P[3 * n + 1], P[3 * n + 2]);
      if (d < best) { best = d; pair = [k, j]; }
    }
  }
  return { distance: best, pair };
}
Wire.clearance = clearance;

// One step of the coarse wire.  Returns { moved (largest displacement, in R), rejected, energy, length }.
function step(ws) {
  const { o, N, P, V, prev, next } = ws, ra = ws.ra, dt = ws.dt;
  const F = new Float64Array(3 * N); let energy = 0;
  const expR = 3 + o.alpha, plain = o.alpha === 0 && o.beta === 1;   // F_a = H r^(1+β) d/r = H r^β d;  F_r = K r^(-(2+α)) d/r = K d / r^(3+α)
  for (let k = 0; k < N; k++) {
    const xk = P[3 * k] / ra, yk = P[3 * k + 1] / ra, zk = P[3 * k + 2] / ra;
    for (const j of [prev[k], next[k]]) {
      const dx = P[3 * j] / ra - xk, dy = P[3 * j + 1] / ra - yk, dz = P[3 * j + 2] / ra - zk, r = Math.hypot(dx, dy, dz);
      if (r < 1e-12) continue;
      const f = o.H * (plain ? r : Math.pow(r, o.beta));
      F[3 * k] += f * dx; F[3 * k + 1] += f * dy; F[3 * k + 2] += f * dz;
      if (j > k) energy += o.H * Math.pow(r, o.beta + 2) / (o.beta + 2);
    }
    for (let j = k + 1; j < N; j++) {
      if (j === prev[k] || j === next[k]) continue;
      const dx = xk - P[3 * j] / ra, dy = yk - P[3 * j + 1] / ra, dz = zk - P[3 * j + 2] / ra, r2 = dx * dx + dy * dy + dz * dz;
      if (r2 < 1e-24) continue;
      const r = Math.sqrt(r2), f = plain ? o.K / (r2 * r) : o.K / Math.pow(r, expR);
      F[3 * k] += f * dx; F[3 * k + 1] += f * dy; F[3 * k + 2] += f * dz;
      F[3 * j] -= f * dx; F[3 * j + 1] -= f * dy; F[3 * j + 2] -= f * dz;
      energy += plain ? o.K / r : o.K * Math.pow(r, -(1 + o.alpha)) / (1 + o.alpha);
    }
  }
  // damped Euler in units of ra, the displacement clamped to dmax
  const Q = new Float64Array(3 * N); let rejected = 0, moved = 0;
  for (let k = 0; k < N; k++) {
    for (let c = 0; c < 3; c++) V[3 * k + c] = (1 - o.gamma) * V[3 * k + c] + F[3 * k + c] * dt;
    let dx = V[3 * k] * dt, dy = V[3 * k + 1] * dt, dz = V[3 * k + 2] * dt; const l = Math.hypot(dx, dy, dz);
    if (l > o.dmax) { const s = o.dmax / l; dx *= s; dy *= s; dz *= s; V[3 * k] *= s; V[3 * k + 1] *= s; V[3 * k + 2] *= s; }
    Q[3 * k] = P[3 * k] + dx * ra; Q[3 * k + 1] = P[3 * k + 1] + dy * ra; Q[3 * k + 2] = P[3 * k + 2] + dz * ra;
    moved = Math.max(moved, Math.hypot(dx, dy, dz) * ra);
  }
  // Scharein's rejection: a point may not come within dclose of a segment that is not adjacent to it (unless the
  // move takes it further from that segment than it was, so a point that starts too close can still get away)
  const dclose = o.dclose * ra;
  for (let k = 0; k < N; k++) {
    const skip = new Set([k, prev[k], prev[prev[k]], next[k]]);
    let ok = true;
    for (let j = 0; j < N && ok; j++) {
      if (skip.has(j)) continue;
      const n = next[j];
      const d = pointSegment(Q[3 * k], Q[3 * k + 1], Q[3 * k + 2], P[3 * j], P[3 * j + 1], P[3 * j + 2], P[3 * n], P[3 * n + 1], P[3 * n + 2]);
      if (d < dclose && d < pointSegment(P[3 * k], P[3 * k + 1], P[3 * k + 2], P[3 * j], P[3 * j + 1], P[3 * j + 2], P[3 * n], P[3 * n + 1], P[3 * n + 2])) ok = false;
    }
    if (!ok) { rejected++; for (let c = 0; c < 3; c++) { Q[3 * k + c] = P[3 * k + c]; V[3 * k + c] = 0; } }
  }
  for (let i = 0; i < 3 * N; i++) P[i] = Q[i];
  ws.dt *= 1 - o.decay; ws.steps++; ws.energies.push(energy);
  return { moved, rejected, energy, length: length(ws) };
}
Wire.step = step;
// settled: the energy fell by less than `tol` of itself over the last `window` steps
Wire.settled = (ws, window = 100, tol = 1e-3) => { const e = ws.energies, n = e.length; return n > window && (e[n - 1 - window] - e[n - 1]) / Math.abs(e[n - 1]) < tol; };
// Write the wire into the mesh: the fine loops interpolated from the coarse ones, the surface carried along as a
// rubber sheet (Minimal.extend).  Returns the largest displacement of any vertex.
function apply(ws, mesh, Minimal) {
  ws.F = interpolate(ws);
  const disp = new Float64Array(mesh.pos.length); let maxWire = 0;
  ws.fLoops.forEach((f, l) => { const Fl = ws.F[l]; for (let k = 0; k < f.n; k++) { const v = f.vert[k]; for (let c = 0; c < 3; c++) { disp[3 * v + c] = Fl[3 * k + c] - mesh.pos[3 * v + c]; maxWire = Math.max(maxWire, Math.abs(disp[3 * v + c])); } } });
  if (maxWire === 0) return { moved: 0, wire: 0 };
  const r = Minimal.extend(mesh, disp, { tol: 1e-8 });
  ws.fLoops.forEach((f, l) => { const Fl = ws.F[l]; for (let k = 0; k < f.n; k++) { const v = f.vert[k]; for (let c = 0; c < 3; c++) mesh.pos[3 * v + c] = Fl[3 * k + c]; } });   // exactly, not up to the solver
  return { moved: r.moved, wire: maxWire };
}
Wire.apply = apply;
Wire.restart = ws => { ws.dt = ws.o.dt0; ws.V.fill(0); ws.energies = []; };
// A copy of everything taming changes, so that a tick which pinches the film can be undone: the wire is then as
// far open as it can be while the surface still follows it, which is where taming should stop.
Wire.snapshot = (ws, mesh) => ({ P: ws.P.slice(), V: ws.V.slice(), dt: ws.dt, steps: ws.steps, energies: ws.energies.slice(), F: ws.F.map(f => f.slice()),
                                 fLoops: ws.fLoops.map(f => ({ vert: f.vert.slice(), n: f.n })),
                                 pos: mesh.pos.slice(), tri: mesh.tri.slice(), kind: mesh.kind.slice(), fixed: mesh.fixed.slice(), loops: mesh.loops.map(l => l.slice()) });
Wire.rollback = (ws, mesh, snap, Minimal) => {
  ws.P.set(snap.P); ws.V.set(snap.V); ws.dt = snap.dt; ws.steps = snap.steps; ws.energies = snap.energies.slice();
  ws.F = snap.F.map(f => f.slice()); ws.fLoops = snap.fLoops.map(f => ({ vert: f.vert.slice(), n: f.n }));
  mesh.pos = snap.pos.slice(); mesh.tri = snap.tri.slice(); mesh.kind = snap.kind.slice(); mesh.fixed = snap.fixed.slice(); mesh.loops = snap.loops.map(l => l.slice());
  Minimal.prepare(mesh);
};
// One tick of taming with the surface carried along: wire steps until some coarse point has moved `threshold` r_a
// (or `maxSteps` steps), then the surface follows (apply), is remeshed to the edge length `target` scaled by the
// wire's growth, and takes a gentle film round.  Returns { steps, moved, remesh, film, pinched }.
Wire.carry = (ws, mesh, Minimal, opts = {}) => {
  const threshold = (opts.threshold || 0.5) * ws.ra, maxSteps = opts.maxSteps || 25, P0 = ws.P.slice();
  let steps = 0, moved = 0;
  while (steps < maxSteps && moved < threshold) {
    step(ws); steps++;
    for (let k = 0; k < ws.N; k++) moved = Math.max(moved, Math.abs(ws.P[3 * k] - P0[3 * k]), Math.abs(ws.P[3 * k + 1] - P0[3 * k + 1]), Math.abs(ws.P[3 * k + 2] - P0[3 * k + 2]));
  }
  apply(ws, mesh, Minimal);
  const target = (opts.target || ws.meshEdge || (ws.meshEdge = Minimal.meanEdgeLength(mesh))) * length(ws) / ws.L0;
  const r = Minimal.remesh(mesh, target, { tangential: opts.tangential === undefined ? 0.3 : opts.tangential });
  Wire.renumber(ws, r.map);
  // then a gentle film round (a few edge lengths squared of mean curvature flow with positive weights), so the
  // surface tracks the soap film as the wire moves instead of drifting away from it as a rubber sheet
  let film = null;
  if (opts.film !== false) { const f = Minimal.relax(mesh, { dt: (opts.filmStep || 4) * target * target, flips: true, tangential: 0.3, positive: true, tol: 1e-7 }); if (!f.failed) film = f; }
  return { steps, moved, remesh: r, film, pinched: Minimal.bunched(mesh) };
};
// after a remeshing that renumbered the mesh's vertices (Minimal.remesh returns the map)
Wire.renumber = (ws, map) => { for (const f of ws.fLoops) for (let k = 0; k < f.n; k++) f.vert[k] = map[f.vert[k]]; };

if (typeof module !== 'undefined' && module.exports) module.exports = Wire; else window.Wire = Wire;
})();
