/* wire.js -- taming the wire: Scharein's KnotPlot relaxation of the knot, as van Wijk and Cohen adopt it
   (IEEE TVCG 2006, §5.3), applied to the boundary loops of the surface mesh.  Node-compatible and browser.

   Every wire vertex is a point mass, attracted by its two neighbours along the loop with the generalised Hooke force
   F_a(r) = H r^(1+β) and repelled by every other wire vertex with the generalised electrostatic force
   F_r(r) = K r^(-(2+α)), distances measured in units of r_a, the initial spacing of the vertices (β = 0 linear
   springs, α = 0 inverse-square repulsion; α = 4 is KnotPlot's "strong repulsion").  Damped explicit Euler,
   v <- (1 - γ) v + F dt, x <- x + v dt with the displacement clamped to d_max, and a move is rejected when it would
   bring the vertex within d_close of a segment of the wire that is not adjacent to it: Scharein's argument that
   d_close > d_max then keeps the knot from passing through itself.  The step dt decays by (1 - μ) each step, so a
   cycle settles; a new cycle restarts from dt0.  The wire grows while it relaxes (the repulsion wins over the springs
   until the spacing is about 1.5 r_a), which is the model's own scale; nothing here depends on it.  One departure
   from the paper: the surface is not given springs of its own, the app carries it along as a rubber sheet, by the
   harmonic extension of the wire's displacement (Minimal.extend). */
(function () {
'use strict';
const Wire = {};
Wire.DEFAULTS = { alpha: 0, beta: 1, H: 1, K: 1, gamma: 0.15, dt0: 0.2, decay: 0.0003, dmax: 0.25, dclose: 0.5 };

// The state of a relaxation: the loops as wire indices, positions P (3N, in units of R), velocities, the step.
function init(mesh, opts = {}) {
  const o = Object.assign({}, Wire.DEFAULTS, opts), loops = mesh.loops, vert = [], loopOf = [], prev = [], next = [];
  for (let l = 0; l < loops.length; l++) {
    const L = loops[l], base = vert.length;
    for (let k = 0; k < L.length; k++) { vert.push(L[k]); loopOf.push(l); prev.push(base + (k + L.length - 1) % L.length); next.push(base + (k + 1) % L.length); }
  }
  const N = vert.length, P = new Float64Array(3 * N);
  for (let k = 0; k < N; k++) for (let c = 0; c < 3; c++) P[3 * k + c] = mesh.pos[3 * vert[k] + c];
  const ws = { o, N, vert: Int32Array.from(vert), loopOf: Int32Array.from(loopOf), prev: Int32Array.from(prev), next: Int32Array.from(next),
               P, V: new Float64Array(3 * N), dt: o.dt0, steps: 0 };
  ws.ra = length(ws) / N; ws.L0 = ws.ra * N;
  return ws;
}
Wire.init = init;
function length(ws) { let L = 0; for (let k = 0; k < ws.N; k++) { const j = ws.next[k]; L += Math.hypot(ws.P[3 * j] - ws.P[3 * k], ws.P[3 * j + 1] - ws.P[3 * k + 1], ws.P[3 * j + 2] - ws.P[3 * k + 2]); } return L; }
Wire.length = length;
// distance from point p to segment ab
function pointSegment(px, py, pz, ax, ay, az, bx, by, bz) {
  const ux = bx - ax, uy = by - ay, uz = bz - az, l2 = ux * ux + uy * uy + uz * uz;
  let t = l2 > 0 ? ((px - ax) * ux + (py - ay) * uy + (pz - az) * uz) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
  return Math.hypot(px - ax - t * ux, py - ay - t * uy, pz - az - t * uz);
}
// the smallest distance from a vertex to a segment of the wire not adjacent to it (in R), and the pair
function clearance(ws) {
  const { N, P, prev, next } = ws; let best = Infinity, pair = [-1, -1];
  for (let k = 0; k < N; k++) {
    const skip = new Set([k, prev[k], prev[prev[k]], next[k]]);          // segments (j, next j) touching k or its neighbours
    for (let j = 0; j < N; j++) {
      if (skip.has(j)) continue;
      const n = next[j], d = pointSegment(P[3 * k], P[3 * k + 1], P[3 * k + 2], P[3 * j], P[3 * j + 1], P[3 * j + 2], P[3 * n], P[3 * n + 1], P[3 * n + 2]);
      if (d < best) { best = d; pair = [k, j]; }
    }
  }
  return { distance: best, pair };
}
Wire.clearance = clearance;

// One step.  Returns { moved (largest displacement, in R), rejected, energy, length }.
function step(ws) {
  const { o, N, P, V, prev, next } = ws, ra = ws.ra, dt = ws.dt;
  const F = new Float64Array(3 * N); let energy = 0;
  const expA = 1 + o.beta, expR = 3 + o.alpha;               // F_a = H r^(1+β) along d/r = H r^β d;  F_r = K r^(-(2+α)) d/r = K d / r^(3+α)
  for (let k = 0; k < N; k++) {
    const xk = P[3 * k] / ra, yk = P[3 * k + 1] / ra, zk = P[3 * k + 2] / ra;
    for (const j of [prev[k], next[k]]) {
      const dx = P[3 * j] / ra - xk, dy = P[3 * j + 1] / ra - yk, dz = P[3 * j + 2] / ra - zk, r = Math.hypot(dx, dy, dz);
      if (r < 1e-12) continue;
      const f = o.H * Math.pow(r, o.beta);
      F[3 * k] += f * dx; F[3 * k + 1] += f * dy; F[3 * k + 2] += f * dz;
      if (j > k) energy += o.H * Math.pow(r, expA + 1) / (expA + 1);
    }
    for (let j = k + 1; j < N; j++) {
      if (j === prev[k] || j === next[k]) continue;
      const dx = xk - P[3 * j] / ra, dy = yk - P[3 * j + 1] / ra, dz = zk - P[3 * j + 2] / ra, r = Math.hypot(dx, dy, dz);
      if (r < 1e-12) continue;
      const f = o.K / Math.pow(r, expR);
      F[3 * k] += f * dx; F[3 * k + 1] += f * dy; F[3 * k + 2] += f * dz;
      F[3 * j] -= f * dx; F[3 * j + 1] -= f * dy; F[3 * j + 2] -= f * dz;
      energy += o.K * Math.pow(r, -(1 + o.alpha)) / (1 + o.alpha);
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
  // Scharein's rejection: a vertex may not come within dclose of a segment that is not adjacent to it
  const dclose = o.dclose * ra;
  for (let k = 0; k < N; k++) {
    const skip = new Set([k, prev[k], prev[prev[k]], next[k]]);
    let ok = true;
    for (let j = 0; j < N && ok; j++) {
      if (skip.has(j)) continue;
      const n = next[j];
      if (pointSegment(Q[3 * k], Q[3 * k + 1], Q[3 * k + 2], P[3 * j], P[3 * j + 1], P[3 * j + 2], P[3 * n], P[3 * n + 1], P[3 * n + 2]) < dclose) ok = false;
    }
    if (!ok) { rejected++; for (let c = 0; c < 3; c++) { Q[3 * k + c] = P[3 * k + c]; V[3 * k + c] = 0; } }
  }
  for (let i = 0; i < 3 * N; i++) P[i] = Q[i];
  ws.dt *= 1 - o.decay; ws.steps++;
  return { moved, rejected, energy, length: length(ws) };
}
Wire.step = step;
// Write the wire's positions into the mesh, carrying the surface along as a rubber sheet (Minimal.extend); returns
// the largest displacement of any vertex.
function apply(ws, mesh, Minimal) {
  const disp = new Float64Array(mesh.pos.length); let maxWire = 0;
  for (let k = 0; k < ws.N; k++) { const v = ws.vert[k]; for (let c = 0; c < 3; c++) { disp[3 * v + c] = ws.P[3 * k + c] - mesh.pos[3 * v + c]; maxWire = Math.max(maxWire, Math.abs(disp[3 * v + c])); } }
  if (maxWire === 0) return { moved: 0, wire: 0 };
  const r = Minimal.extend(mesh, disp, { tol: 1e-8 });
  for (let k = 0; k < ws.N; k++) { const v = ws.vert[k]; for (let c = 0; c < 3; c++) mesh.pos[3 * v + c] = ws.P[3 * k + c]; }   // exactly, not up to the solver
  return { moved: r.moved, wire: maxWire };
}
Wire.apply = apply;
Wire.restart = ws => { ws.dt = ws.o.dt0; ws.V.fill(0); };

if (typeof module !== 'undefined' && module.exports) module.exports = Wire; else window.Wire = Wire;
})();
