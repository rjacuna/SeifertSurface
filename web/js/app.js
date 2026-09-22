/* app.js -- viewer and UI for the Seifert Surface web app.  Globals: THREE, katex, Seifert, Minimal. */
(function () {
'use strict';
const $ = id => document.getElementById(id);
const isMobile = matchMedia('(max-width: 640px)').matches || navigator.maxTouchPoints > 1;

const state = { spacing: 0.45, bandwidth: 1.2, bulge: 0.7, angular: isMobile ? 96 : 144, round: 0.1,
                dt: 3.1, flips: true, tangential: 0.3, perframe: 2, stop: -7,          // dt is log10 of the step; 3.1 is ∞
                coloring: 'sides', opacity: 1, wireframe: false, wire: true, thick: 0.02, ghost: false, axes: false };
let mesh = null, record = null, lastText = '', iteration = 0, running = false, lastStats = null, sizeRadius = 2, areas = [];

// ------------------------------------------------------------------ scene
const view = $('view');
const renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
renderer.setSize(view.clientWidth, view.clientHeight);
renderer.setClearColor(0xffffff, 1);
view.appendChild(renderer.domElement);
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(40, view.clientWidth / view.clientHeight, 0.01, 200);
camera.up.set(0, 0, 1);
const controls = new THREE.OrbitControls(camera, renderer.domElement);
controls.enableDamping = true; controls.dampingFactor = 0.12;
controls.touches = { ONE: THREE.TOUCH.ROTATE, TWO: THREE.TOUCH.DOLLY_PAN };
controls.addEventListener('change', requestRender);
const keyLight = new THREE.DirectionalLight(0xffffff, 0.8); keyLight.position.set(-6, 4, 8); camera.add(keyLight); scene.add(camera);
scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const fill = new THREE.DirectionalLight(0xffffff, 0.35); fill.position.set(4, -5, -7); scene.add(fill);
let dirty = true;
function requestRender() { dirty = true; }
(function loop() {
  requestAnimationFrame(loop);
  const moved = controls.update();
  if (moved || dirty) { renderer.render(scene, camera); dirty = false; }
})();
// The relaxation runs on its own timer, not on the animation frames (which a hidden or embedded page throttles):
// `perframe` rounds per tick, the picture updated as frames come.
let ticking = false;
function tick() {
  ticking = false;
  if (!running || !mesh) return;
  for (let k = 0; k < state.perframe && running; k++) relaxOnce();
  if (running) { ticking = true; setTimeout(tick, 0); }
}
window.addEventListener('resize', () => { renderer.setSize(view.clientWidth, view.clientHeight); camera.aspect = view.clientWidth / view.clientHeight; camera.updateProjectionMatrix(); requestRender(); });
function resetView() {
  const R = sizeRadius, dist = R / Math.sin(Math.PI / 9) * 1.05;
  const el = 25 * Math.PI / 180, s3 = Math.cos(el), c3 = Math.sin(el), s5 = Math.sin(Math.PI / 5), c5 = Math.cos(Math.PI / 5);
  camera.position.set(dist * s3 * c5, -dist * s3 * s5, dist * c3);
  controls.target.set(0, 0, 0); camera.near = 0.01 * R; camera.far = 100 * R; camera.updateProjectionMatrix(); controls.update(); requestRender();
}

// ------------------------------------------------------------------ materials and scene objects
// The surface is oriented, so its two sides get two colors: the front (the side the normal points to, +z on the
// disks) in green, the back in red -- the classical way to show that a Seifert surface is two-sided.
const FRONT = 0x3a9d5d, BACK = 0xc8463a;
const frontMaterial = new THREE.MeshPhongMaterial({ color: FRONT, side: THREE.FrontSide, shininess: 30, specular: new THREE.Color(0x333333), transparent: true, opacity: 1 });
const backMaterial = new THREE.MeshPhongMaterial({ color: BACK, side: THREE.BackSide, shininess: 30, specular: new THREE.Color(0x333333), transparent: true, opacity: 1 });
const ghostMaterial = new THREE.MeshBasicMaterial({ color: 0x888888, wireframe: true, transparent: true, opacity: 0.25, depthWrite: false });
const WIRE_COLORS = [0x2d4f9e, 0xb3261e, 0xd08a00, 0x5b2a86, 0x0b7a75, 0x7a4a00];
const group = new THREE.Group(); scene.add(group);
let surfaceGeom = null, frontMesh = null, backMesh = null, ghostMesh = null, wireGroup = null, axesGroup = null;
const SERIF = '"STIX Two Text", "STIX Two Math", "Times New Roman", Times, serif';
function makeLabel(text, x, y, z) {
  const canvas = document.createElement('canvas'), ctx = canvas.getContext('2d'), pr = 2, fs = 24;
  ctx.font = `italic ${fs}px ${SERIF}`; const w = Math.ceil(ctx.measureText(text).width) + 8, h = fs + 10;
  canvas.width = w * pr; canvas.height = h * pr; ctx.scale(pr, pr);
  ctx.fillStyle = '#222'; ctx.textBaseline = 'middle'; ctx.font = `italic ${fs}px ${SERIF}`; ctx.fillText(text, 4, h / 2);
  const tex = new THREE.Texture(canvas); tex.needsUpdate = true;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, sizeAttenuation: false, depthWrite: false, transparent: true }));
  sp.position.set(x, y, z); sp.scale.set(w / 700, h / 700, 1);
  return sp;
}
function buildAxes(R) {
  const grp = new THREE.Group();
  grp.add(new THREE.AxesHelper(1.15 * R));
  grp.add(makeLabel('x', 1.3 * R, 0, 0)); grp.add(makeLabel('y', 0, 1.3 * R, 0)); grp.add(makeLabel('z', 0, 0, 1.3 * R));
  return grp;
}
function disposeObject(obj) { if (!obj) return; group.remove(obj); obj.traverse(o => { if (o.geometry) o.geometry.dispose(); }); }
function buildWire() {
  disposeObject(wireGroup); wireGroup = new THREE.Group();
  if (!mesh) return;
  mesh.loops.forEach((loop, k) => {
    const pts = loop.map(v => new THREE.Vector3(mesh.pos[3 * v], mesh.pos[3 * v + 1], mesh.pos[3 * v + 2]));
    if (pts.length < 3) return;
    const curve = new THREE.CatmullRomCurve3(pts, true, 'centripetal');
    const tube = new THREE.TubeGeometry(curve, Math.min(1200, 2 * pts.length), state.thick, 10, true);
    const mat = new THREE.MeshPhongMaterial({ color: WIRE_COLORS[k % WIRE_COLORS.length], shininess: 50, specular: new THREE.Color(0x333344) });
    wireGroup.add(new THREE.Mesh(tube, mat));
  });
  wireGroup.visible = state.wire; group.add(wireGroup);
}
function buildSurfaceObjects() {
  disposeObject(frontMesh); disposeObject(backMesh); disposeObject(ghostMesh); frontMesh = backMesh = ghostMesh = null;
  if (!mesh) return;
  surfaceGeom = new THREE.BufferGeometry();
  surfaceGeom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(mesh.pos), 3));
  surfaceGeom.setAttribute('color', new THREE.BufferAttribute(new Float32Array(mesh.pos.length), 3));
  surfaceGeom.setIndex(new THREE.BufferAttribute(new Uint32Array(mesh.tri), 1));
  surfaceGeom.computeVertexNormals();
  frontMesh = new THREE.Mesh(surfaceGeom, frontMaterial); backMesh = new THREE.Mesh(surfaceGeom, backMaterial);
  group.add(frontMesh); group.add(backMesh);
  const g0 = new THREE.BufferGeometry();
  g0.setAttribute('position', new THREE.BufferAttribute(new Float32Array(mesh.initial), 3));
  g0.setIndex(new THREE.BufferAttribute(new Uint32Array(mesh.tri), 1));
  ghostMesh = new THREE.Mesh(g0, ghostMaterial); ghostMesh.visible = state.ghost; group.add(ghostMesh);
  applyColoring(); applyOpacity();
}
function updateSurfacePositions() {                                  // after a relaxation round (connectivity may have changed by flips)
  if (!surfaceGeom) return;
  surfaceGeom.attributes.position.array.set(mesh.pos); surfaceGeom.attributes.position.needsUpdate = true;
  if (surfaceGeom.index.array.length !== mesh.tri.length || surfaceGeom.index.array[7] !== mesh.tri[7]) surfaceGeom.setIndex(new THREE.BufferAttribute(new Uint32Array(mesh.tri), 1));
  else { surfaceGeom.index.array.set(mesh.tri); surfaceGeom.index.needsUpdate = true; }
  surfaceGeom.computeVertexNormals();
  if (state.coloring === 'curvature') applyColoring();
  requestRender();
}
// coloring: two sides (plain materials), the disks and bands of the starting surface, or the discrete mean curvature
const PART_COLORS = [[0.22, 0.60, 0.36], [0.16, 0.44, 0.70], [0.80, 0.55, 0.10], [0.55, 0.25, 0.60], [0.05, 0.50, 0.48], [0.70, 0.30, 0.20], [0.45, 0.55, 0.15], [0.35, 0.35, 0.65]];
function applyColoring() {
  if (!surfaceGeom) return;
  const vertexColors = state.coloring !== 'sides';
  for (const m of [frontMaterial, backMaterial]) { m.vertexColors = vertexColors; m.needsUpdate = true; }
  frontMaterial.color.set(vertexColors ? 0xffffff : FRONT); backMaterial.color.set(vertexColors ? 0xffffff : BACK);
  if (vertexColors) {
    const col = surfaceGeom.attributes.color.array, V = mesh.pos.length / 3;
    if (state.coloring === 'parts') {
      for (let v = 0; v < V; v++) { const k = mesh.kind[v], c = k < mesh.n ? PART_COLORS[k % PART_COLORS.length] : [0.78, 0.78, 0.78]; col[3 * v] = c[0]; col[3 * v + 1] = c[1]; col[3 * v + 2] = c[2]; }
    } else {
      const H = Minimal.meanCurvature(mesh).H; let scale = 0; const sorted = Float64Array.from(H).sort();
      scale = sorted[Math.floor(0.95 * (sorted.length - 1))] || 1;   // the 95th percentile saturates
      for (let v = 0; v < V; v++) { const t = Math.min(1, H[v] / scale); col[3 * v] = 0.95 * t + 0.92 * (1 - t); col[3 * v + 1] = 0.35 * t + 0.92 * (1 - t); col[3 * v + 2] = 0.2 * t + 0.92 * (1 - t); }
    }
    surfaceGeom.attributes.color.needsUpdate = true;
  }
  requestRender();
}
function applyOpacity() {
  for (const m of [frontMaterial, backMaterial]) { m.opacity = state.opacity; m.depthWrite = state.opacity >= 1; m.wireframe = state.wireframe; }
  requestRender();
}
function rebuildDecorations() {
  disposeObject(axesGroup); axesGroup = buildAxes(sizeRadius * 0.9); axesGroup.visible = state.axes; group.add(axesGroup);
  buildWire(); requestRender();
}

// ------------------------------------------------------------------ the relaxation
const dtValue = () => state.dt >= 3.05 ? Infinity : Math.pow(10, state.dt) * (mesh ? mesh.meanEdge * mesh.meanEdge : 1);
function relaxOnce() {
  if (!mesh) return;
  const stats = Minimal.relax(mesh, { dt: dtValue(), flips: state.flips, tangential: state.tangential > 0 ? state.tangential : false, tol: 1e-8 });
  iteration++; lastStats = stats; areas.push(stats.area);
  updateSurfacePositions(); showStats();
  // stop when the area has stopped decreasing: less than 10^stop of itself over the last five rounds
  if (areas.length > 5 && (areas[areas.length - 6] - stats.area) / stats.area < Math.pow(10, state.stop)) setRunning(false);
}
function showStats() {
  $('v-iter').textContent = String(iteration);
  if (!mesh) { for (const id of ['v-area', 'v-h', 'v-moved']) $(id).textContent = '–'; return; }
  const H = Minimal.meanCurvature(mesh);
  $('v-area').textContent = Minimal.area(mesh).toFixed(4);
  $('v-h').textContent = (H.rms * mesh.params.R).toPrecision(3);
  $('v-moved').textContent = lastStats ? lastStats.moved.toExponential(1) : '–';
}
function setRunning(on) { running = !!on && !!mesh; $('run').textContent = running ? '❚❚ Pause' : '▶ Relax'; if (running && !ticking) { ticking = true; setTimeout(tick, 0); } }
function resetSurface() {
  if (!mesh) return;
  setRunning(false); mesh.pos.set(mesh.initial); mesh.tri.set(mesh.initialTri); Minimal.prepare(mesh); iteration = 0; lastStats = null; areas = [];
  updateSurfacePositions(); showStats();
}

// ------------------------------------------------------------------ text
const T = tex => katex.renderToString(tex, { throwOnError: false, output: 'html' });
const esc = t => String(t).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
const mixed = str => String(str).split('$').map((part, i) => i % 2 ? T(part) : esc(part)).join('');
for (const el of document.querySelectorAll('.tex')) katex.render(el.textContent, el, { throwOnError: false, output: 'html' });
function setInfo(html, cls) { const el = $('info'); el.innerHTML = html; el.className = cls || ''; }
function setRemarks(notes) { const el = $('remarks'); el.innerHTML = notes.map(mixed).join('<br>'); el.hidden = !notes.length; }
const prettyName = name => name.replace(/^K(\d+)([an]?)_(\d+)$/, (m, a, b, c) => `${a}${b ? b : ''}_{${c}}`).replace(/^L(\d+)([an])(\d+)((?:_\d+)+)$/, (m, a, b, c, d) => `L${a}${b}${c}\\{${d.slice(1).split('_').join(',')}\\}`);

// ------------------------------------------------------------------ the knot table (KnotInfo / LinkInfo braid words, data/knots.json)
let table = null, tableByName = null, tableByWord = null;
async function loadTable() {
  if (table) return table;
  const r = await fetch('data/knots.json'); if (!r.ok) throw new Error('the knot table could not be loaded; serve the app over http');
  table = await r.json(); tableByName = new Map(); tableByWord = new Map();
  for (const rec of table) { tableByName.set(rec.name, rec); const k = JSON.stringify(rec.braid); if (!tableByWord.has(k)) tableByWord.set(k, rec); }
  return table;
}
function nextFrame() { return new Promise(r => { let done = false; const go = () => { if (!done) { done = true; setTimeout(r, 0); } }; requestAnimationFrame(go); setTimeout(go, 100); }); }

// ------------------------------------------------------------------ build
async function build(text, opts = {}) {
  text = (text || '').trim(); if (!text) return;
  lastText = text; $('input').value = text; updateHash(text);
  $('busy').hidden = false; setRemarks([]); setRunning(false); await nextFrame();
  try {
    const parsed = Seifert.parseBraid(text);
    if (parsed.error) throw new Error(parsed.error);
    let word = parsed.word, label = parsed.label || null; record = null;
    if (parsed.name) {
      await loadTable(); record = tableByName.get(parsed.name) || table.find(r => r.name.startsWith(parsed.name + '_')) || null;
      if (!record) throw new Error(`no ${parsed.name} in the table (KnotInfo knots to 12 crossings, LinkInfo links to 11)`);
      word = record.braid;
    } else if (word.length) {
      try { await loadTable(); record = tableByWord.get(JSON.stringify(word)) || null; } catch (e) { /* names are optional */ }
    }
    const bd = Seifert.braidData(word);
    if (!bd.connected) throw new Error(`σ${[...Array(bd.n - 1).keys()].map(i => i + 1).find(i => !word.some(g => Math.abs(g) === i))} never occurs: the closed braid is split and the surface would be disconnected`);
    if (bd.n * bd.c > 6000) throw new Error('too many strands and crossings to draw');
    mesh = Seifert.buildSurface(word, { spacing: state.spacing, bandWidth: state.bandwidth, bulge: state.bulge, angular: state.angular, round: state.round });
    Minimal.prepare(mesh); mesh.initialTri = mesh.tri.slice(); mesh.meanEdge = Minimal.meanEdgeLength(mesh);
    iteration = 0; lastStats = null; areas = [];
    sizeRadius = Math.max(mesh.params.R + mesh.params.b + mesh.params.W, 0.6 * bd.n * state.spacing + 0.5);
    buildSurfaceObjects(); rebuildDecorations(); if (!opts.keepView) resetView();
    // the info line
    const desc = [];
    if (record) desc.push(`<b>${T(prettyName(record.name))}</b>${label ? ' = ' + esc(label) : ''}`); else if (label) desc.push(`<b>${esc(label)}</b>`);
    desc.push(T(Seifert.wordTeX(word)) + (word.length ? ` on ${bd.n} strands` : ' (one strand)'));
    desc.push(`${bd.c} crossing${bd.c === 1 ? '' : 's'}, ${bd.mu} component${bd.mu === 1 ? '' : 's'}`);
    desc.push(T(`\\chi = ${bd.chi},\\ g = ${bd.genus}`));
    if (word.length) {
      const a = Seifert.alexander(word), s = Seifert.signature(word);
      desc.push(T(`\\Delta(t) = ${Seifert.alexanderTeX(a.coeffs)}`));
      desc.push(T(`\\det = ${a.det},\\ \\sigma = ${s.signature}`) + (s.nullity ? mixed(` (nullity $${s.nullity}$)`) : ''));
      const bound = Math.ceil(a.degree / 2);
      if (bd.mu === 1) {
        if (record && record.genus !== undefined) desc.push(record.genus === bd.genus ? `<span class="ok">this surface has the knot's genus ${record.genus}</span>` : `knot genus ${record.genus} (KnotInfo): this surface has ${bd.genus - record.genus} extra handle${bd.genus - record.genus === 1 ? '' : 's'}`);
        else if (bound === bd.genus) desc.push(`<span class="ok">minimal genus (${mixed('$\\deg\\Delta = 2g$')})</span>`);
        else desc.push(mixed(`genus $\\ge ${bound}$ from $\\Delta$`));
      }
      if (bd.positive || bd.negative) desc.push(`${bd.positive ? 'positive' : 'negative'} braid`);
    }
    setInfo(desc.join(' | '));
    const notes = [];
    if (state.round === 0) notes.push('the wire keeps the corners of the stacked-disk construction');
    else notes.push(`wire: the boundary of the stacked disks and bands, corners rounded at $${state.round}R$, ${mesh.loops.reduce((s, l) => s + l.length, 0)} points; mesh: ${mesh.pos.length / 3} vertices, ${mesh.tri.length / 3} triangles`);
    setRemarks(notes);
    showStats();
    if (opts.autorun) setRunning(true);
  } catch (e) {
    mesh = null; record = null; disposeObject(frontMesh); disposeObject(backMesh); disposeObject(ghostMesh); disposeObject(wireGroup); frontMesh = backMesh = ghostMesh = wireGroup = null;
    setInfo(`<span class="err">${mixed(e.message)}</span>`, ''); showStats();
  } finally { $('busy').hidden = true; requestRender(); }
}
function hashFor(text) { return '#' + encodeURIComponent(text); }
function updateHash(text) { try { history.replaceState(null, '', hashFor(text)); } catch (e) {} }
function shareLink() { return location.origin + location.pathname + location.search + (lastText ? hashFor(lastText) : location.hash); }

// ------------------------------------------------------------------ UI wiring
const EXAMPLES = [
  ['trefoil  3₁ = T(2,3)', '3_1'], ['figure-eight  4₁', '4_1'], ['cinquefoil  5₁ = T(2,5)', '5_1'], ['5₂', '5_2'], ['stevedore  6₁', '6_1'], ['6₂', '6_2'], ['6₃', '6_3'],
  ['8₁₉ = T(3,4)', '8_19'], ['10₁₂₄ = T(3,5)', '10_124'], ['12n₂₄₂ = P(−2,3,7)  (Lehmer)', '12n242'], ['Conway knot  11n₃₄', '11n34'],
  ['Hopf link', 'hopf'], ['Whitehead link', 'whitehead'], ['Borromean rings', 'borromean'], ['T(4,5)', 'T(4,5)'], ['T(2,7)', 'T(2,7)'], ['T(3,6)  (3 components)', 'T(3,6)'],
  ['a 4-strand braid', '1 2 3 -1 2 -3 1 2'], ['unknot, one disk', '()'],
];
for (const [label, value] of EXAMPLES) { const o = document.createElement('option'); o.value = value; o.textContent = label; $('examples').appendChild(o); }
$('examples').addEventListener('change', e => { if (e.target.value) build(e.target.value); e.target.value = ''; });
$('build').addEventListener('click', () => build($('input').value));
$('input').addEventListener('keydown', e => { if (e.key === 'Enter') build($('input').value); });
$('run').addEventListener('click', () => setRunning(!running));
$('step-once').addEventListener('click', () => { setRunning(false); relaxOnce(); });
$('reset-surface').addEventListener('click', resetSurface);
const drawer = $('drawer'), tabs = [...drawer.querySelectorAll('.tab')];
let openSection = null;
function showSection(name) {
  openSection = name; drawer.classList.toggle('open', !!name);
  for (const t of tabs) { const on = t.dataset.section === name; t.classList.toggle('active', on); t.setAttribute('aria-expanded', String(on)); }
  if (name) for (const sec of drawer.querySelectorAll('#panel > section')) sec.hidden = sec.id !== 'section-' + name;
}
for (const t of tabs) t.addEventListener('click', () => showSection(openSection === t.dataset.section ? null : t.dataset.section));
const topbar = $('topbar');
function measureTopbar() { document.documentElement.style.setProperty('--topbar-h', topbar.offsetHeight + 'px'); }
new ResizeObserver(measureTopbar).observe(topbar); measureTopbar();
function bindRange(id, key, show, onChange) {
  const el = $(id); el.value = state[key]; $('v-' + id).textContent = show(state[key]);
  el.addEventListener('input', () => { state[key] = Number(el.value); $('v-' + id).textContent = show(state[key]); onChange('input'); });
  el.addEventListener('change', () => onChange('change'));
}
function bindCheck(id, key, onChange) { const el = $(id); el.checked = state[key]; el.addEventListener('change', () => { state[key] = el.checked; onChange(); requestRender(); }); }
const rebuild = kind => { if (kind === 'change' && lastText !== null && mesh) build(lastText, { keepView: true }); };
bindRange('spacing', 'spacing', v => v.toFixed(2), rebuild);
bindRange('bandwidth', 'bandwidth', v => v.toFixed(1), rebuild);
bindRange('bulge', 'bulge', v => v.toFixed(1), rebuild);
bindRange('angular', 'angular', v => v, rebuild);
bindRange('round', 'round', v => v.toFixed(2), rebuild);
$('rebuild').addEventListener('click', () => { if (lastText !== null) build(lastText, { keepView: true }); });
bindRange('dt', 'dt', v => v >= 3.05 ? '∞ (harmonic)' : Math.pow(10, v).toPrecision(2), () => {});
bindCheck('flips', 'flips', () => {});
bindRange('tangential', 'tangential', v => v.toFixed(2), () => {});
bindRange('perframe', 'perframe', v => v, () => {});
bindRange('stop', 'stop', v => '10^' + v, () => {});
$('coloring').value = state.coloring; $('coloring').addEventListener('change', e => { state.coloring = e.target.value; applyColoring(); });
bindRange('opacity', 'opacity', v => v.toFixed(2), applyOpacity);
bindCheck('wireframe', 'wireframe', applyOpacity);
bindCheck('wire', 'wire', () => { if (wireGroup) wireGroup.visible = state.wire; });
bindRange('thick', 'thick', v => v.toFixed(3), kind => { if (kind === 'change') buildWire(); });
bindCheck('ghost', 'ghost', () => { if (ghostMesh) ghostMesh.visible = state.ghost; });
bindCheck('axes', 'axes', () => { if (axesGroup) axesGroup.visible = state.axes; });
$('reset-view').addEventListener('click', resetView);
$('snapshot').addEventListener('click', () => { renderer.render(scene, camera); const a = document.createElement('a'); a.href = renderer.domElement.toDataURL('image/png'); a.download = 'seifert-surface.png'; document.body.appendChild(a); a.click(); a.remove(); });
$('share').addEventListener('click', async () => { const url = shareLink(); try { await navigator.clipboard.writeText(url); $('share').textContent = 'Copied'; setTimeout(() => $('share').textContent = 'Copy link', 1200); } catch (e) { prompt('Link:', url); } });
setTimeout(() => { $('hint').hidden = true; }, 9000);
window.SEIFERT_DEBUG = { state, get mesh() { return mesh; }, build, relaxOnce, setRunning, render: () => renderer.render(scene, camera), canvas: renderer.domElement, Seifert, Minimal };
const initial = decodeURIComponent((location.hash || '').slice(1));
build(initial || '3_1');
})();
