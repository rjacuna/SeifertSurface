/* precomputed.js -- the shipped films: a tamed wire and the minimal surface settled on it, computed once by
   web/data/precompute.mjs and loaded instead of running the pipeline in the browser.  Node-compatible and browser.

   A file holds one surface: the vertex positions (single precision, enough for a mesh whose edges are 0.02 or
   longer), the triangles, and the vertex kinds of the scaffold they came from (the disks-and-bands coloring), after
   a header with the counts and the few numbers the app shows -- the taming steps, the wire's length then and now,
   the scaffold's edge length, the area.  Everything else is recomputed on load: the boundary loops from the mesh
   (Seifert.boundaryLoops), the coarse wire from those loops (Wire.init), the scaffold from the braid word, which is
   what Reset goes back to.

   layout (little endian)
     0   char[4]   "SFRT"
     4   uint32    format version
     8   uint32    V, vertices
     12  uint32    F, triangles
     16  uint32    taming steps
     20  uint32    flags: bit 0 the film pinched, bit 1 the taming hit its step limit
     24  float64   L0, the wire's length before taming
     32  float64   the wire's length after
     40  float64   edge0, the scaffold's mean edge length
     48  float64   area
     56  float32[3V]  positions
     ... uint32[3F]   triangles
     ... int16[V]     kinds
*/
(function () {
'use strict';
const Precomputed = {};
const MAGIC = 0x54524653;                                   // "SFRT" little endian
const VERSION = 1, HEADER = 56;
Precomputed.VERSION = VERSION;

Precomputed.encode = function (s) {
  const V = s.pos.length / 3, F = s.tri.length / 3;
  const buf = new ArrayBuffer(HEADER + 12 * V + 12 * F + 2 * V), view = new DataView(buf);
  view.setUint32(0, MAGIC, true); view.setUint32(4, VERSION, true);
  view.setUint32(8, V, true); view.setUint32(12, F, true);
  view.setUint32(16, s.steps | 0, true);
  view.setUint32(20, (s.pinched ? 1 : 0) | (s.capped ? 2 : 0), true);
  view.setFloat64(24, s.L0, true); view.setFloat64(32, s.length, true); view.setFloat64(40, s.edge0, true); view.setFloat64(48, s.area, true);
  new Float32Array(buf, HEADER, 3 * V).set(s.pos);
  new Uint32Array(buf, HEADER + 12 * V, 3 * F).set(s.tri);
  new Int16Array(buf, HEADER + 12 * V + 12 * F, V).set(s.kind);
  return buf;
};

Precomputed.decode = function (buf) {
  if (buf.byteLength < HEADER) throw new Error('a precomputed film is truncated');
  const view = new DataView(buf);
  if (view.getUint32(0, true) !== MAGIC) throw new Error('not a precomputed film');
  const version = view.getUint32(4, true);
  if (version !== VERSION) throw new Error(`a precomputed film of version ${version}, this is version ${VERSION}`);
  const V = view.getUint32(8, true), F = view.getUint32(12, true), flags = view.getUint32(20, true);
  const want = HEADER + 12 * V + 12 * F + 2 * V;
  if (buf.byteLength !== want) throw new Error(`a precomputed film of ${buf.byteLength} bytes, expected ${want}`);
  return { V, F, steps: view.getUint32(16, true), pinched: !!(flags & 1), capped: !!(flags & 2),
           L0: view.getFloat64(24, true), length: view.getFloat64(32, true), edge0: view.getFloat64(40, true), area: view.getFloat64(48, true),
           pos: new Float32Array(buf, HEADER, 3 * V), tri: new Uint32Array(buf, HEADER + 12 * V, 3 * F), kind: new Int16Array(buf, HEADER + 12 * V + 12 * F, V) };
};

// the manifest's key for a braid word: the word itself, so that a name and the word it resolves to share an entry
Precomputed.key = word => word.join(',') || '()';

// Put a decoded film into a mesh built from the same braid word (Seifert.buildSurface), keeping the mesh's params,
// braid and word; the boundary loops and the fixed vertices are recomputed from the triangles.  Returns the mesh.
Precomputed.install = function (mesh, film, Seifert, Minimal) {
  mesh.pos = Float64Array.from(film.pos);
  mesh.tri = Uint32Array.from(film.tri);
  mesh.kind = Int16Array.from(film.kind);
  mesh.loops = Seifert.boundaryLoops(mesh.tri, film.V);
  mesh.fixed = new Uint8Array(film.V);
  for (const loop of mesh.loops) for (const v of loop) mesh.fixed[v] = 1;
  Minimal.prepare(mesh);
  return mesh;
};

if (typeof module !== 'undefined' && module.exports) module.exports = Precomputed; else window.Precomputed = Precomputed;
})();
