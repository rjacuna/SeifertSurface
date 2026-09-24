/* wire.js -- taming the wire: Scharein's KnotPlot relaxation of the knot, as van Wijk and Cohen adopt it
   (IEEE TVCG 2006, §5.3), applied to the boundary loops of the surface mesh.  Node-compatible and browser.

   Every point of the wire is a unit mass, attracted by its two neighbours along the loop with the generalised Hooke
   force F_a(r) = H r^(1+β) and repelled by every other point of the wire with the generalised electrostatic force
   F_r(r) = K r^(-(2+α)), distances measured in units of r_a, the initial spacing of the points (β = 1, α = 0:
   inverse-square repulsion; α = 4 is KnotPlot's "strong repulsion").  Damped explicit Euler, v <- (1 - γ) v + F dt,
   x <- x + v dt with the displacement clamped to d_max, and a move is refused when it would bring the point within
   d_close (KnotPlot's 0.5 r_a, `guard`) of a segment of the wire not adjacent to it: Scharein's argument that
   d_close > d_max then keeps the knot from passing through itself.  The step dt decays by (1 - μ) each step, so a
   cycle settles; a new cycle restarts from dt0.  The wire grows while it relaxes (the repulsion wins over the springs
   until the spacing is about 1.5 r_a), the model's own scale; nothing here depends on it.

   Three things are ours.  First, there are two knots.  The one tamed is a fat knot, a solid torus of radius
   `radius` whose core is the coarse polygon; it is never drawn and has no surface.  KnotPlot's forces act on the
   solid torus: two of its points repel across the gap between their cross-sections, the core distance less the
   tube's diameter (less (2/π)·s for two points s apart along one strand, what the tube bent at its own radius leaves),
   so the strands keep a diameter apart, and every crossing is prescribed that separation, however the knot relaxes.
   The other is the thin knot, the wire: the boundary of the mesh, drawn and spanned by the film.  It lies inside the
   fat knot and follows its core, record by record, as fast as the film can follow it; the film has no say in the
   fat knot.
   Second, the dynamics runs on a coarse copy of each boundary loop, resampled at equal arclength with spacing
   `coarse` (0.1 disk radii, under the wire's smallest features), since the forces cost N² and the mesh boundary has
   several hundred points per loop; the thin knot is interpolated from the fat knot's core by a Catmull-Rom spline at
   equal arclength.  And third, the surface is not given springs of its own: it is carried along as a rubber sheet,
   by the harmonic extension of the thin knot's displacement (Minimal.extend), and the film settled on it. */
(function () {
'use strict';
const Wire = {};
// coarse: the spacing of the coarse polygon; radius: the tube's; both in the mesh's units (R = 1).  guard and dmax
// are in units of r_a.
Wire.DEFAULTS = { alpha: 0, beta: 1, H: 1, K: 1, gamma: 0.15, dt0: 0.2, decay: 0.0003, dmax: 0.25, guard: 0.5, radius: 0.25, gmin: 0.1, threshold: 0.5, ramp: 300, coarse: 0.1 };

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
// opts.L0: the wire's length when its taming began, for a wire already tamed (a shipped film): the coarse copy then
// has as many points as it was tamed with, and r_a is what it was, so the forces are where they were.
function init(mesh, opts = {}) {
  const o = Object.assign({}, Wire.DEFAULTS, opts), loops = mesh.loops, cLoops = [], fLoops = [];
  const P = [], prev = [], next = [], loopOf = [];
  const fines = loops.map(L => { const n = L.length, fine = new Float64Array(3 * n); for (let k = 0; k < n; k++) for (let d = 0; d < 3; d++) fine[3 * k + d] = mesh.pos[3 * L[k] + d]; return fine; });
  const total = fines.reduce((t, f, l) => t + polyLength(f, loops[l].length), 0), coarse = o.coarse * (opts.L0 ? total / opts.L0 : 1);
  loops.forEach((L, l) => {
    const n = L.length, fine = fines[l];
    const m = Math.max(12, Math.min(n, Math.round(polyLength(fine, n) / coarse))), base = P.length / 3, Q = resample(fine, n, m);
    for (let k = 0; k < m; k++) { P.push(Q[3 * k], Q[3 * k + 1], Q[3 * k + 2]); prev.push(base + (k + m - 1) % m); next.push(base + (k + 1) % m); loopOf.push(l); }
    cLoops.push({ base, m }); fLoops.push({ vert: Int32Array.from(L), n });
  });
  const N = P.length / 3;
  const ws = { o, N, P: Float64Array.from(P), V: new Float64Array(3 * N), prev: Int32Array.from(prev), next: Int32Array.from(next), loopOf: Int32Array.from(loopOf),
               cLoops, fLoops, dt: o.dt0, steps: 0, energies: [] };
  ws.L0 = opts.L0 || length(ws); ws.ra = ws.L0 / N;
  ws.F = interpolate(ws); ws.thin = ws.P.slice();
  return ws;
}
Wire.init = init;
function length(ws) { return polyLength(ws.P, ws.N); }
Wire.length = length;
// The fine loops from the coarse ones: the spline, sampled at equal arclength with the loop's point count.  The
// coarse points slide along the loop as they relax (nothing in the forces pins the parametrisation), and the mesh
// boundary must not slide with them, or the surface next to it is sheared without end; so the first fine point is
// anchored at the point of the new curve nearest to where it was, and the others follow at equal arclength.
function interpolate(ws, P = ws.P) {
  const out = [];
  ws.cLoops.forEach((c, l) => {
    const f = ws.fLoops[l], per = Math.max(2, Math.ceil(4 * f.n / c.m)), dense = spline(P, c.m, per, 3 * c.base), nd = c.m * per;
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

// The position of every coarse point along its loop, in units of ra, and the distance of two points of a loop along
// it, the shorter way round.
function arclengths(ws) {
  const { P, cLoops, loopOf } = ws, ra = ws.ra, at = new Float64Array(ws.N), total = cLoops.map(({ base, m }) => {
    let t = 0;
    for (let k = 0; k < m; k++) { const a = base + k, b = base + (k + 1) % m; at[a] = t; t += Math.hypot(P[3 * b] - P[3 * a], P[3 * b + 1] - P[3 * a + 1], P[3 * b + 2] - P[3 * a + 2]) / ra; }
    return t;
  });
  return { at, along: (k, j) => { const d = Math.abs(at[k] - at[j]), L = total[loopOf[k]]; return Math.min(d, L - d); } };
}
Wire.arclengths = arclengths;
// The tube's clearance: the smallest distance between two points of the wire that are not neighbours on the tube
// (further apart along it than half its circumference, or on different loops), in R, and the pair.
function thickness(ws) {
  const { N, P, loopOf } = ws, s = arclengths(ws), window = Math.PI * ws.o.radius / ws.ra; let best = Infinity, pair = [-1, -1];
  for (let k = 0; k < N; k++) for (let j = k + 1; j < N; j++) {
    if (loopOf[k] === loopOf[j] && s.along(k, j) <= window) continue;
    const d = Math.hypot(P[3 * k] - P[3 * j], P[3 * k + 1] - P[3 * j + 1], P[3 * k + 2] - P[3 * j + 2]);
    if (d < best) { best = d; pair = [k, j]; }
  }
  return { distance: best, pair };
}
Wire.thickness = thickness;

// One step of the coarse wire.  Returns { moved (largest displacement, in R), rejected, energy, length }.
function step(ws) {
  const { o, N, P, V, prev, next, loopOf } = ws, ra = ws.ra, dt = ws.dt;
  const F = new Float64Array(3 * N); let energy = 0;
  const expR = 3 + o.alpha, plain = o.alpha === 0 && o.beta === 1;   // F_a = H r^(1+β) d/r = H r^β d;  F_r = K r^(-(2+α)) d/r = K d / r^(3+α)
  // the fat knot's diameter, in units of ra, and where along it each point is
  const diam = 2 * o.radius / ra, twoOverPi = 2 / Math.PI, s = arclengths(ws);
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
      // The fat knot repels itself across the gap between its surfaces: the distance of the two points' cross-
      // sections of the tube, the core distance less what the tube takes up between them.  That is its diameter for
      // points on different strands, and (2/π) s for two points s apart along one strand, up to the diameter at
      // s = π·radius: a tube bent as tightly as it can be, at its own radius, leaves exactly that (Jordan's
      // inequality, sin x ≥ 2x/π), so the gap stays positive on any curve the tube can follow.
      const r = Math.sqrt(r2), taken = loopOf[k] !== loopOf[j] ? diam : Math.min(diam, twoOverPi * s.along(k, j));
      const g = Math.max(r - taken, o.gmin);
      const f = (plain ? o.K / (g * g) : o.K / Math.pow(g, expR - 1)) / r;
      energy += plain ? o.K / g : o.K * Math.pow(g, -(1 + o.alpha)) / (1 + o.alpha);
      F[3 * k] += f * dx; F[3 * k + 1] += f * dy; F[3 * k + 2] += f * dz;
      F[3 * j] -= f * dx; F[3 * j + 1] -= f * dy; F[3 * j + 2] -= f * dz;
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
  // move takes it further from that segment than it was, so a point that starts too close can still get away).  It
  // only keeps the knot from passing through itself; the fat knot's repulsion keeps the strands apart.
  const dclose = o.guard * ra;
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
// ---- the two knots
// The fat knot is the state above: the coarse polygon (P, V), the core of a solid torus of radius o.radius, tamed
// by step().  It is never drawn and carries no surface, and nothing the film does moves it back.  The thin knot is
// the wire: the mesh's boundary, drawn and spanned by the film; it lies inside the fat knot and follows its core.
// The fat knot runs ahead: each time its core has moved `threshold` r_a since the last record it is recorded
// (ws.keys), and the thin knot goes through the records one at a time, as fast as the film can follow.  Two
// records are that close and the tube keeps the strands a diameter apart, so going from one to the next the thin
// knot cannot pass through itself.  ws.thin is the record the thin knot is at.

// The thin knot onto the core of the fat knot at P (its current state by default): the fine loops interpolated
// from P, the surface carried along as a rubber sheet (Minimal.extend).  Returns the largest displacement.
function apply(ws, mesh, Minimal, P = ws.P) {
  // the mesh's own loops, which every remeshing renumbers, not a copy of them
  ws.fLoops = mesh.loops.map(L => ({ vert: Int32Array.from(L), n: L.length }));
  ws.F = interpolate(ws, P); ws.thin = Float64Array.from(P);
  const disp = new Float64Array(mesh.pos.length); let maxWire = 0;
  ws.fLoops.forEach((f, l) => { const Fl = ws.F[l]; for (let k = 0; k < f.n; k++) { const v = f.vert[k]; for (let c = 0; c < 3; c++) { disp[3 * v + c] = Fl[3 * k + c] - mesh.pos[3 * v + c]; maxWire = Math.max(maxWire, Math.abs(disp[3 * v + c])); } } });
  if (maxWire === 0) return { moved: 0, wire: 0 };
  const r = Minimal.extend(mesh, disp, { tol: 1e-8 });
  ws.fLoops.forEach((f, l) => { const Fl = ws.F[l]; for (let k = 0; k < f.n; k++) { const v = f.vert[k]; for (let c = 0; c < 3; c++) mesh.pos[3 * v + c] = Fl[3 * k + c]; } });   // exactly, not up to the solver
  return { moved: r.moved, wire: maxWire };
}
Wire.apply = apply;
// (Re)start taming: the fat knot from where the thin knot is, at rest, with the options in force.  It inflates: it
// starts as thick as the wire already allows (half the closest approach of two strands) and grows to the radius set
// over the first `ramp` steps, so that no strand is flung off another and the thin knot, and the film, can follow.
Wire.begin = (ws, opts = {}) => {
  Object.assign(ws.o, opts);
  ws.P.set(ws.thin); ws.V.fill(0); ws.dt = ws.o.dt0; ws.energies = []; ws.steps = 0;
  ws.keys = []; ws.last = ws.P.slice(); ws.done = false;
  ws.radius = ws.o.radius; ws.radius0 = Math.min(ws.radius, 0.5 * thickness(ws).distance);
};
// The fat knot's own taming, at most `budget` steps of it: records its core as it goes, and is done when its
// energy has settled or after `maxSteps` steps.  Returns the number of steps taken.
Wire.tame = (ws, budget = 50, maxSteps = 2500) => {
  const threshold = ws.o.threshold * ws.ra; let n = 0;
  if (!ws.keys) Wire.begin(ws);
  while (n < budget && !ws.done) {
    const inflating = ws.steps < ws.o.ramp;
    ws.o.radius = ws.radius0 + (ws.radius - ws.radius0) * Math.min(1, ws.steps / ws.o.ramp);
    step(ws); n++;
    if (inflating) ws.energies = [];                          // the energy of a tube still inflating says nothing
    let moved = 0; for (let i = 0; i < 3 * ws.N; i++) moved = Math.max(moved, Math.abs(ws.P[i] - ws.last[i]));
    ws.done = (!inflating && Wire.settled(ws)) || ws.steps >= maxSteps;
    if (moved >= threshold || (ws.done && moved > 0)) { ws.keys.push(ws.P.slice()); ws.last = ws.P.slice(); }
  }
  return n;
};
// Taming is over when the fat knot is done and the thin knot has caught up with it.
Wire.caughtUp = ws => !!ws.done && !ws.keys.length;
// The thin knot one record further, the surface with it: carried (apply), remeshed to the edge length `target`
// scaled by the knot's growth, and a gentle film round.  If that pinches the film the step is undone and the thin
// knot stays where it was; the fat knot is not touched either way.  Returns null when there is no record to go to,
// else { remesh, film, pinched }.
Wire.follow = (ws, mesh, Minimal, opts = {}) => {
  if (!ws.keys || !ws.keys.length) return null;
  const key = ws.keys[0], undo = Wire.snapshot(ws, mesh);
  apply(ws, mesh, Minimal, key);
  const target = (opts.target || ws.meshEdge || (ws.meshEdge = Minimal.meanEdgeLength(mesh))) * polyLength(key, ws.N) / ws.L0;
  const r = Minimal.remesh(mesh, target, { tangential: opts.tangential === undefined ? 0.3 : opts.tangential });
  ws.fLoops = mesh.loops.map(L => ({ vert: Int32Array.from(L), n: L.length }));
  // then a gentle film round (a few edge lengths squared of mean curvature flow with positive weights), so the
  // surface tracks the soap film as the wire moves instead of drifting away from it as a rubber sheet
  let film = null;
  if (opts.film !== false) { const f = Minimal.relax(mesh, { dt: (opts.filmStep || 4) * target * target, flips: true, tangential: 0.3, positive: true, tol: 1e-7 }); if (!f.failed) film = f; }
  if (Minimal.bunched(mesh)) { Wire.rollback(ws, mesh, undo, Minimal); return { remesh: r, film, pinched: true }; }
  ws.keys.shift();
  return { remesh: r, film, pinched: false };
};
// The thin knot and its surface, so that a step the film cannot take can be undone, or the thin knot walked back
// to where the film could still settle.  The fat knot is not part of it.
Wire.snapshot = (ws, mesh) => ({ thin: ws.thin.slice(), F: ws.F.map(f => f.slice()), fLoops: ws.fLoops.map(f => ({ vert: f.vert.slice(), n: f.n })),
                                 pos: mesh.pos.slice(), tri: mesh.tri.slice(), kind: mesh.kind.slice(), fixed: mesh.fixed.slice(), loops: mesh.loops.map(l => l.slice()) });
Wire.rollback = (ws, mesh, snap, Minimal) => {
  ws.thin = snap.thin.slice(); ws.F = snap.F.map(f => f.slice()); ws.fLoops = snap.fLoops.map(f => ({ vert: f.vert.slice(), n: f.n }));
  mesh.pos = snap.pos.slice(); mesh.tri = snap.tri.slice(); mesh.kind = snap.kind.slice(); mesh.fixed = snap.fixed.slice(); mesh.loops = snap.loops.map(l => l.slice());
  Minimal.prepare(mesh);
};
// the thin knot's length, which is the wire's
Wire.thinLength = ws => polyLength(ws.thin, ws.N);

if (typeof module !== 'undefined' && module.exports) module.exports = Wire; else window.Wire = Wire;
})();
