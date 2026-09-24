// Precompute the films the app ships: for every example in its menu, build the scaffold, tame the wire with the
// surface carried along, settle the film on it, and write the result to web/data/films/.  The app loads these
// instead of running the pipeline, so an example appears at once; a knot that is not among them, or a change to
// any of the parameters below, still runs live in the browser.
//
//   node web/data/precompute.mjs            all of them, about four minutes
//   node web/data/precompute.mjs 12n242 T(4,5)    only these
//
// The parameters used are written into the manifest, and the app compares them against its own state: they are the
// app's defaults, and `angular` is deliberately not compared, since it sets the scaffold's resolution only and the
// shipped mesh carries its own.
import { createRequire } from 'module';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
const require = createRequire(import.meta.url);
const Seifert = require('../js/seifert.js'), Minimal = require('../js/minimal.js'), Wire = require('../js/wire.js'), Precomputed = require('../js/precomputed.js');
const here = path.dirname(fileURLToPath(import.meta.url)), out = path.join(here, 'films');

// the app's defaults (web/js/app.js `state`), and its stopping rules
const PARAMS = { spacing: 0.45, bandWidth: 1.2, bulge: 0.7, round: 0.1, alpha: 0, K: 1, radius: 0.25, maxsteps: 2500, angular: 144 };
const SETTLE = { every: 5, tangential: 0.3, stop: -7, maxRounds: 400 };
// the inputs of the app's Examples menu
const EXAMPLES = ['3_1', '4_1', '5_1', '5_2', '6_1', '6_2', '6_3', '7_1', '8_19', '10_124', '11n34', '12n242',
                  'hopf', 'whitehead', 'borromean', 'T(3,6)', 'T(4,5)', '1 2 3 -1 2 -3 1 2', '()'];

const table = JSON.parse(fs.readFileSync(path.join(here, 'knots.json')));
const byName = new Map(table.map(r => [r.name, r]));
function wordOf(text) {
  const parsed = Seifert.parseBraid(text);
  if (parsed.error) throw new Error(`${text}: ${parsed.error}`);
  if (!parsed.name) return parsed.word;
  const rec = byName.get(parsed.name) || table.find(r => r.name.startsWith(parsed.name + '_'));
  if (!rec) throw new Error(`${text}: no ${parsed.name} in the table`);
  return rec.braid;
}

// the pipeline, exactly as the app runs it
function compute(word) {
  const mesh = Seifert.buildSurface(word, { spacing: PARAMS.spacing, bandWidth: PARAMS.bandWidth, bulge: PARAMS.bulge, angular: PARAMS.angular, round: PARAMS.round });
  Minimal.prepare(mesh); mesh.attrs = ['kind']; mesh.edge0 = Minimal.meanEdgeLength(mesh);
  const ws = Wire.init(mesh, { alpha: PARAMS.alpha, K: PARAMS.K, H: 1, radius: PARAMS.radius });
  ws.meshEdge = mesh.edge0;
  // Taming, a tick at a time, as the app does: the fat knot tames on (never drawn, no surface, never moved back by
  // the film) and the thin knot, the wire, follows its core one record further with the surface.  A step the film
  // cannot take is undone and the thin knot stops there; a history of the thin knot, one every HISTORY ticks, lets
  // it be walked back if the film then cannot settle.
  const HISTORY = 4, history = [];
  let rolled = false, tick = 0;
  Wire.begin(ws);
  while (!Wire.caughtUp(ws)) {
    Wire.tame(ws, 50, PARAMS.maxsteps);
    const snap = Wire.snapshot(ws, mesh), r = Wire.follow(ws, mesh, Minimal, { tangential: SETTLE.tangential });
    if (r && r.pinched) { rolled = true; break; }
    if (r && tick++ % HISTORY === 0) { history.push(snap); if (history.length > 12) history.shift(); }
  }
  const capped = !rolled && ws.steps >= PARAMS.maxsteps;
  function settleHere() {
    const st = Minimal.settleInit(mesh, { every: SETTLE.every, tangential: SETTLE.tangential });
    const areas = []; let s = null;
    for (let r = 0; r < SETTLE.maxRounds; r++) {
      s = Minimal.settleRound(mesh, st); areas.push(s.area);
      if (s.pinched) return { rounds: st.round, area: s.area, bad: true };
      if (s.harmonic && areas.length > 5 && (areas[areas.length - 6] - s.area) / s.area < Math.pow(10, SETTLE.stop)) break;
    }
    return { rounds: st.round, area: s.area, bad: false };
  }
  let pinched = false, rounds = 0, area = Minimal.area(mesh);
  for (;;) {
    const r = settleHere(); rounds += r.rounds; area = r.area;
    if (!r.bad) { pinched = false; break; }
    pinched = true;
    const back = history.pop();
    if (!back) break;
    Wire.rollback(ws, mesh, back, Minimal); rolled = true;
  }
  return { mesh, ws, pinched, capped, rolled, rounds, area, H: Minimal.meanCurvature(mesh) };
}

const want = process.argv.slice(2), inputs = want.length ? want : EXAMPLES;
fs.mkdirSync(out, { recursive: true });
const manifestPath = path.join(out, 'index.json');
const manifest = (want.length && fs.existsSync(manifestPath)) ? JSON.parse(fs.readFileSync(manifestPath)) : { version: Precomputed.VERSION, params: PARAMS, settle: SETTLE, names: {}, films: {} };
manifest.version = Precomputed.VERSION; manifest.params = PARAMS; manifest.settle = SETTLE; manifest.names = manifest.names || {};
let total = 0;
for (const input of inputs) {
  const t0 = Date.now(), word = wordOf(input), key = Precomputed.key(word), tname = Seifert.tableName(input);
  if (tname) manifest.names[tname] = key;             // so that a name finds its film without the knot table
  const name = (input.replace(/[^A-Za-z0-9]+/g, '_').replace(/^_|_$/g, '') || 'unknot') + '.bin';
  const { mesh, ws, pinched, capped, rolled, rounds, area, H } = compute(word);
  const buf = Precomputed.encode({ pos: mesh.pos, tri: mesh.tri, kind: mesh.kind, steps: ws.steps, pinched, capped,
                                   L0: ws.L0, length: Wire.thinLength(ws), edge0: mesh.edge0, area });
  fs.writeFileSync(path.join(out, name), Buffer.from(buf));
  manifest.films[key] = { file: name, bytes: buf.byteLength, V: mesh.pos.length / 3, F: mesh.tri.length / 3, steps: ws.steps,
                          growth: +(Wire.thinLength(ws) / ws.L0).toFixed(4), area: +area.toFixed(4), rms: +H.rms.toFixed(4),
                          chi: Minimal.eulerCharacteristic(mesh), loops: mesh.loops.length, pinched, capped, inputs: [input] };
  total += buf.byteLength;
  console.log(`${input.padEnd(18)} ${String(mesh.pos.length / 3).padStart(6)} vertices  ${(buf.byteLength / 1024).toFixed(0).padStart(4)} kB  ` +
              `${ws.steps} taming steps, ${rounds} film rounds, wire ×${(Wire.thinLength(ws) / ws.L0).toFixed(2)}, area ${area.toFixed(3)}, rms |H| ${H.rms.toFixed(3)}` +
              `${pinched ? ', PINCHED' : ''}${rolled ? ', taming rolled back at a pinch' : ''}${capped ? ', taming capped' : ''}  (${((Date.now() - t0) / 1000).toFixed(1)} s)`);
}
// a full run owns the directory: films no longer listed are stale and go
if (!want.length) {
  const kept = new Set(Object.values(manifest.films).map(f => f.file));
  for (const f of fs.readdirSync(out)) if (f.endsWith('.bin') && !kept.has(f)) { fs.unlinkSync(path.join(out, f)); console.log(`removed the stale ${f}`); }
  for (const [name, key] of Object.entries(manifest.names)) if (!manifest.films[key]) delete manifest.names[name];
}
fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 1) + '\n');
console.log(`\n${inputs.length} films, ${(total / 1024 / 1024).toFixed(2)} MB written, manifest with ${Object.keys(manifest.films).length} entries`);
