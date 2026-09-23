/* app.js -- viewer and UI for the Seifert Surface web app.  Globals: THREE, katex, Seifert, Minimal, Wire.
   Build: the scaffold (Seifert.buildSurface), then either the film that ships with the app for this braid word
   (data/films/, see precomputed.js), or the pipeline run here: the wire tamed with the film following
   (Wire.carry), then the film settled on the still wire (Minimal.settleRound).  Tame wire and Reset run it here in
   any case.  Two pictures: the soap film, thin-film interference from a per-vertex
   thickness, and the two-sided rubber sheet; the Display tab adds the disks and bands and the mean curvature. */
(function () {
'use strict';
const $ = id => document.getElementById(id);
const isMobile = matchMedia('(max-width: 640px)').matches || navigator.maxTouchPoints > 1;

const state = { spacing: 0.45, bandwidth: 1.2, bulge: 0.7, angular: isMobile ? 96 : 144, round: 0.1, auto: true,
                alpha: 0, kh: 0, dclose: 2, maxsteps: 2500,                                     // the wire: kh is log10 of K/H
                every: 5, tangential: 0.3, perframe: 2, stop: -7,                                // the film
                coloring: 'soap', opacity: 1, wireframe: false, wire: true, thick: 0.02, ghost: false, axes: false,
                soapmin: 100, soapmax: 800, soapopacity: 0.35, envbright: 1.2, wiremetal: 'gold' };
let mesh = null, scaffold = null, wire = null, record = null, lastText = '', sizeRadius = 2;
let phase = 'idle', settle = null, areas = [], lastStats = null, remarks = [];   // phase: idle | taming | settling
let surfaceGeom = null, frontMesh = null, backMesh = null, soapMesh = null, ghostMesh = null, wireGroup = null, axesGroup = null, tubeMaterials = [];
const soap = () => state.coloring === 'soap';

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
const ambient = new THREE.AmbientLight(0xffffff, 0.55); scene.add(ambient);
const fill = new THREE.DirectionalLight(0xffffff, 0.35); fill.position.set(4, -5, -7); scene.add(fill);
let dirty = true;
function requestRender() { dirty = true; }
// The film's triangles sorted back to front for the current view (the index buffer rewritten), so that its
// overlapping sheets blend in depth order.  `sortedFor` remembers the view it was done for; a moved surface clears it.
let sortedFor = null; const sortKeys = { key: null, order: null };
function sortFilm() {
  if (!surfaceGeom || !mesh) return;
  const dir = new THREE.Vector3(); camera.getWorldDirection(dir);
  const sig = camera.matrixWorld.elements.join(',') + '|' + mesh.tri.length;
  if (sortedFor === sig) return;
  sortedFor = sig;
  const F = mesh.tri.length / 3, p = mesh.pos, t = mesh.tri;
  if (!sortKeys.key || sortKeys.key.length !== F) { sortKeys.key = new Float32Array(F); sortKeys.order = new Uint32Array(F); }
  const key = sortKeys.key, order = sortKeys.order;
  for (let f = 0; f < F; f++) {
    const a = t[3 * f], b = t[3 * f + 1], c = t[3 * f + 2];
    key[f] = (p[3 * a] + p[3 * b] + p[3 * c]) * dir.x + (p[3 * a + 1] + p[3 * b + 1] + p[3 * c + 1]) * dir.y + (p[3 * a + 2] + p[3 * b + 2] + p[3 * c + 2]) * dir.z;
    order[f] = f;
  }
  order.sort((i, j) => key[j] - key[i]);                        // farthest first
  const idx = surfaceGeom.index.array;
  for (let k = 0; k < F; k++) { const f = order[k]; idx[3 * k] = t[3 * f]; idx[3 * k + 1] = t[3 * f + 1]; idx[3 * k + 2] = t[3 * f + 2]; }
  surfaceGeom.index.needsUpdate = true;
}
(function loop() {
  requestAnimationFrame(loop);
  const moved = controls.update();
  if (moved || dirty) { if (soap()) sortFilm(); renderer.render(scene, camera); dirty = false; }
})();
window.addEventListener('resize', () => { renderer.setSize(view.clientWidth, view.clientHeight); camera.aspect = view.clientWidth / view.clientHeight; camera.updateProjectionMatrix(); requestRender(); });
function resetView() {
  const R = sizeRadius, dist = R / Math.sin(Math.PI / 9) * 1.05;
  const el = 25 * Math.PI / 180, s3 = Math.cos(el), c3 = Math.sin(el), s5 = Math.sin(Math.PI / 5), c5 = Math.cos(Math.PI / 5);
  camera.position.set(dist * s3 * c5, -dist * s3 * s5, dist * c3);
  controls.target.set(0, 0, 0); camera.near = 0.01 * R; camera.far = 100 * R; camera.updateProjectionMatrix(); controls.update(); requestRender();
}
function wireRadius() {                                          // the extent of the wire, for the view
  let r = 0; if (mesh) for (const loop of mesh.loops) for (const v of loop) r = Math.max(r, Math.hypot(mesh.pos[3 * v], mesh.pos[3 * v + 1], mesh.pos[3 * v + 2]));
  return r;
}
function keepInView() {                                          // the wire grows: keep it in view, from the same direction
  const r = wireRadius();
  if (r > sizeRadius) { camera.position.multiplyScalar(r / sizeRadius); sizeRadius = r; camera.far = 100 * r; camera.updateProjectionMatrix(); }
}

// ------------------------------------------------------------------ the environment (for the soap film and the wire's metal)
// A room: a sky gradient with a few soft lights, painted on a canvas as an equirectangular map, prefiltered by PMREM.
function makeEnvironment() {
  const W = 1024, H = 512, canvas = document.createElement('canvas'); canvas.width = W; canvas.height = H;
  const ctx = canvas.getContext('2d'), g = ctx.createLinearGradient(0, 0, 0, H);
  // a studio: bright above, a mid-grey horizon, a dark floor, so that a metal reflects both light and dark
  g.addColorStop(0, '#eef1f4'); g.addColorStop(0.42, '#aeb5bd'); g.addColorStop(0.5, '#5e656e'); g.addColorStop(0.62, '#2c3036'); g.addColorStop(1, '#14171b');
  ctx.fillStyle = g; ctx.fillRect(0, 0, W, H);
  // softboxes with soft edges: a large warm key, a cool fill, a long thin strip (the highlight that reads as metal), a warm spot
  const box = (x, y, w, h, color, soft) => { const r = ctx.createRadialGradient(x, y, 0, x, y, 1); ctx.save(); ctx.translate(x, y); ctx.scale(w, h); r.addColorStop(0, color); r.addColorStop(soft, color); r.addColorStop(1, 'rgba(255,255,255,0)'); ctx.fillStyle = r; ctx.fillRect(-1, -1, 2, 2); ctx.restore(); };
  box(190, 130, 210, 110, '#fff6e6', 0.75); box(860, 150, 120, 170, '#e9f0ff', 0.7); box(500, 215, 330, 22, '#ffffff', 0.8); box(640, 330, 90, 40, '#ffe2b8', 0.6); box(330, 400, 160, 30, '#ffffff', 0.5);
  const tex = new THREE.CanvasTexture(canvas); tex.mapping = THREE.EquirectangularReflectionMapping; tex.encoding = THREE.sRGBEncoding;
  const pmrem = new THREE.PMREMGenerator(renderer); pmrem.compileEquirectangularShader();
  const env = pmrem.fromEquirectangular(tex).texture; tex.dispose(); pmrem.dispose();
  return env;
}
const envTexture = makeEnvironment();
scene.environment = envTexture;

// ------------------------------------------------------------------ materials and scene objects
// The surface is oriented, so its two sides get two colors: the front (the side the normal points to, +z on the
// disks) in green, the back in red -- the classical way to show that a Seifert surface is two-sided.
const FRONT = 0x3a9d5d, BACK = 0xc8463a;
const frontMaterial = new THREE.MeshPhongMaterial({ color: FRONT, side: THREE.FrontSide, shininess: 30, specular: new THREE.Color(0x333333), transparent: true, opacity: 1 });
const backMaterial = new THREE.MeshPhongMaterial({ color: BACK, side: THREE.BackSide, shininess: 30, specular: new THREE.Color(0x333333), transparent: true, opacity: 1 });
const ghostMaterial = new THREE.MeshBasicMaterial({ color: 0x888888, wireframe: true, transparent: true, opacity: 0.25, depthWrite: false });
// The soap film: no diffuse colour, a water-like specular response with the environment, and thin-film
// interference computed spectrally.  A free-standing film of water (n = 1.333) in air reflects, at wavelength λ and
// optical path 2 n d cos θ_f, the Airy fraction R = 2ρ²(1 − cos δ) / (1 + ρ⁴ − 2ρ² cos δ), δ = 4π n d cos θ_f / λ,
// ρ = (n − 1)/(n + 1) the amplitude reflected at each face (the two reflections differ by π, so the thinnest film
// is black).  Integrated over the visible spectrum against the CIE colour matching functions this gives the colour
// of the film as a function of its thickness, Newton's series (black, grey, white, yellow, orange, red, violet,
// blue, green, ...), tabulated once into a lookup texture and normalised to 1 at the brightest constructive
// interference.  The shader multiplies the material's specular light (whose Fresnel term the physically based
// model already provides) by that colour, looked up at the thickness `thick` given per vertex in nanometres, times
// cos θ_f for the angle of view.  Alpha is the film's opacity to what is behind it, rising at grazing angles; the
// reflected light is added unweighted (custom blending), as a film's reflection is not dimmed by its transparency.
const FILM_N = 1.333, LUT_MAX_NM = 1500, LUT_N = 512;
function newtonColours() {
  // Wyman, Sloan, Shirley (2013): analytic fits of the CIE 1931 colour matching functions
  const g = (x, mu, s1, s2) => { const t = (x - mu) / (x < mu ? s1 : s2); return Math.exp(-0.5 * t * t); };
  const xbar = l => 1.056 * g(l, 599.8, 37.9, 31.0) + 0.362 * g(l, 442.0, 16.0, 26.7) - 0.065 * g(l, 501.1, 20.4, 26.2);
  const ybar = l => 0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1);
  const zbar = l => 1.217 * g(l, 437.0, 11.8, 36.0) + 0.681 * g(l, 459.0, 26.0, 13.8);
  const rho2 = Math.pow((FILM_N - 1) / (FILM_N + 1), 2), Rmax = 4 * rho2 / Math.pow(1 + rho2, 2);
  const data = new Uint8Array(4 * LUT_N); let X0 = 0, Y0 = 0, Z0 = 0;
  for (let l = 380; l <= 730; l += 2) { X0 += xbar(l); Y0 += ybar(l); Z0 += zbar(l); }   // a flat spectrum
  for (let i = 0; i < LUT_N; i++) {
    const d = LUT_MAX_NM * i / (LUT_N - 1); let X = 0, Y = 0, Z = 0;
    for (let l = 380; l <= 730; l += 2) {
      const c = Math.cos(4 * Math.PI * FILM_N * d / l), R = 2 * rho2 * (1 - c) / (1 + rho2 * rho2 - 2 * rho2 * c) / Rmax;
      X += R * xbar(l); Y += R * ybar(l); Z += R * zbar(l);
    }
    X /= X0; Y /= Y0; Z /= Z0;                            // relative to the flat spectrum, so R ≡ 1 is white
    let r = 3.2406 * X - 1.5372 * Y - 0.4986 * Z, gg = -0.9689 * X + 1.8758 * Y + 0.0415 * Z, b = 0.0557 * X - 0.2040 * Y + 1.0570 * Z;   // XYZ (D65 white) to linear sRGB
    const m = Math.max(r, gg, b, 1);
    data[4 * i] = Math.round(255 * Math.max(0, r / m)); data[4 * i + 1] = Math.round(255 * Math.max(0, gg / m)); data[4 * i + 2] = Math.round(255 * Math.max(0, b / m)); data[4 * i + 3] = 255;
  }
  const tex = new THREE.DataTexture(data, LUT_N, 1, THREE.RGBAFormat); tex.minFilter = tex.magFilter = THREE.LinearFilter; tex.wrapS = tex.wrapT = THREE.ClampToEdgeWrapping; tex.needsUpdate = true;
  return tex;
}
const soapUniforms = { uSoapOpacity: { value: state.soapopacity }, uNewton: { value: newtonColours() }, uLutMax: { value: LUT_MAX_NM } };
const soapShader = shader => {
  Object.assign(shader.uniforms, soapUniforms);
  const v0 = shader.vertexShader, f0 = shader.fragmentShader;
  shader.vertexShader = v0
    .replace('#include <common>', '#include <common>\nattribute float thick;\nvarying float vThick;')
    .replace('#include <begin_vertex>', '#include <begin_vertex>\nvThick = thick;');
  shader.fragmentShader = f0
    .replace('#include <common>', '#include <common>\nvarying float vThick;\nuniform float uSoapOpacity;\nuniform sampler2D uNewton;\nuniform float uLutMax;')
    .replace('#include <lights_fragment_end>', '#include <lights_fragment_end>\n{\n  float cosT = clamp(abs(dot(normalize(normal), normalize(vViewPosition))), 0.0, 1.0);\n  float sin2 = (1.0 - cosT * cosT) / (1.333 * 1.333);\n  float cosF = sqrt(max(0.0, 1.0 - sin2));\n  vec3 fringe = texture2D(uNewton, vec2(clamp(vThick * cosF / uLutMax, 0.0, 1.0), 0.5)).rgb;\n  reflectedLight.directSpecular *= fringe;\n  reflectedLight.indirectSpecular *= fringe;\n}')
    .replace('#include <output_fragment>', '#include <output_fragment>\n{\n  float cosT = clamp(abs(dot(normalize(normal), normalize(vViewPosition))), 0.0, 1.0);\n  float fres = pow(1.0 - cosT, 4.0);\n  gl_FragColor.a = clamp(uSoapOpacity + (1.0 - uSoapOpacity) * fres, 0.0, 1.0);\n}');
  for (const [what, ok] of [['vertex attribute', shader.vertexShader !== v0], ['interference', shader.fragmentShader.includes('uNewton, vec2')], ['alpha', shader.fragmentShader.includes('uSoapOpacity +')]]) if (!ok) console.warn('soap film shader: the ' + what + ' hook did not match this three.js');
};
function makeSoapMaterial(side) {
  const m = new THREE.MeshPhysicalMaterial({ color: 0x000000, metalness: 0, roughness: 0.04, side, transparent: true, depthWrite: false,
    envMap: envTexture, envMapIntensity: state.envbright, reflectivity: 0.7,
    blending: THREE.CustomBlending, blendSrc: THREE.OneFactor, blendDst: THREE.OneMinusSrcAlphaFactor, blendEquation: THREE.AddEquation });
  m.onBeforeCompile = soapShader; m.customProgramCacheKey = () => 'soap-film';
  return m;
}
// One two-sided pass whose triangles are sorted back to front whenever the view or the surface changes (below),
// so that where one sheet of the film lies behind another the blend is in the right order.
const soapMaterial = makeSoapMaterial(THREE.DoubleSide);
const WIRE_COLORS = [0x2d4f9e, 0xb3261e, 0xd08a00, 0x5b2a86, 0x0b7a75, 0x7a4a00];
// the wire in the soap-film picture: gold by default (its colour is the metal's reflectance, tinting the room it
// reflects), or steel
// measured reflectances of the metals, in linear light (gold 1.00, 0.71, 0.29; steel 0.56, 0.57, 0.58), as
// MeshStandardMaterial takes the colour of a full metal for its F0
const WIRE_METALS = { gold: new THREE.MeshStandardMaterial({ color: new THREE.Color(1.0, 0.71, 0.29), metalness: 1, roughness: 0.18 }), steel: new THREE.MeshStandardMaterial({ color: new THREE.Color(0.56, 0.57, 0.58), metalness: 1, roughness: 0.25 }), dark: new THREE.MeshStandardMaterial({ color: new THREE.Color(0.05, 0.05, 0.06), metalness: 0.9, roughness: 0.35 }) };
const wireMetal = () => WIRE_METALS[state.wiremetal] || WIRE_METALS.gold;
const group = new THREE.Group(); scene.add(group);
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
  disposeObject(wireGroup); wireGroup = new THREE.Group(); tubeMaterials = [];
  if (!mesh) return;
  mesh.loops.forEach((loop, k) => {
    const pts = loop.map(v => new THREE.Vector3(mesh.pos[3 * v], mesh.pos[3 * v + 1], mesh.pos[3 * v + 2]));
    if (pts.length < 3) return;
    const curve = new THREE.CatmullRomCurve3(pts, true, 'centripetal');
    const tube = new THREE.TubeGeometry(curve, Math.min(1200, 2 * pts.length), state.thick, 16, true);
    const mat = new THREE.MeshPhongMaterial({ color: WIRE_COLORS[k % WIRE_COLORS.length], shininess: 50, specular: new THREE.Color(0x333344) });
    tubeMaterials.push(mat);
    wireGroup.add(new THREE.Mesh(tube, soap() ? wireMetal() : mat));
  });
  wireGroup.visible = state.wire; group.add(wireGroup);
}
function buildSurfaceObjects() {                                 // (re)creates the geometry: after a build, a remeshing, a reset
  disposeObject(frontMesh); disposeObject(backMesh); disposeObject(soapMesh); frontMesh = backMesh = soapMesh = null;
  if (!mesh) return;
  surfaceGeom = new THREE.BufferGeometry();
  surfaceGeom.setAttribute('position', new THREE.BufferAttribute(new Float32Array(mesh.pos), 3));
  surfaceGeom.setAttribute('color', new THREE.BufferAttribute(new Float32Array(mesh.pos.length), 3));
  surfaceGeom.setAttribute('thick', new THREE.BufferAttribute(new Float32Array(mesh.pos.length / 3).fill(state.soapmax), 1));
  surfaceGeom.setIndex(new THREE.BufferAttribute(new Uint32Array(mesh.tri), 1));
  surfaceGeom.computeVertexNormals();
  frontMesh = new THREE.Mesh(surfaceGeom, frontMaterial); backMesh = new THREE.Mesh(surfaceGeom, backMaterial);
  soapMesh = new THREE.Mesh(surfaceGeom, soapMaterial); soapMesh.renderOrder = 2;
  group.add(frontMesh); group.add(backMesh); group.add(soapMesh);
  sortedFor = null;
  applyColoring(); applyOpacity();
}
function buildGhost() {
  disposeObject(ghostMesh); ghostMesh = null;
  if (!scaffold) return;
  const g0 = new THREE.BufferGeometry();
  g0.setAttribute('position', new THREE.BufferAttribute(new Float32Array(scaffold.pos), 3));
  g0.setIndex(new THREE.BufferAttribute(new Uint32Array(scaffold.tri), 1));
  ghostMesh = new THREE.Mesh(g0, ghostMaterial); ghostMesh.visible = state.ghost; group.add(ghostMesh);
}
function updateSurfacePositions() {                              // the same vertices, moved
  if (!surfaceGeom) return;
  if (surfaceGeom.attributes.position.array.length !== mesh.pos.length || surfaceGeom.index.array.length !== mesh.tri.length) { buildSurfaceObjects(); return; }
  surfaceGeom.attributes.position.array.set(mesh.pos); surfaceGeom.attributes.position.needsUpdate = true;
  surfaceGeom.index.array.set(mesh.tri); surfaceGeom.index.needsUpdate = true; sortedFor = null;
  surfaceGeom.computeVertexNormals();
  if (state.coloring === 'curvature') applyColoring();
  requestRender();
}
// coloring: two sides (plain materials), the disks and bands of the scaffold, the discrete mean curvature, or the film
const PART_COLORS = [[0.22, 0.60, 0.36], [0.16, 0.44, 0.70], [0.80, 0.55, 0.10], [0.55, 0.25, 0.60], [0.05, 0.50, 0.48], [0.70, 0.30, 0.20], [0.45, 0.55, 0.15], [0.35, 0.35, 0.65]];
function applyColoring() {
  if (!surfaceGeom) return;
  const s = soap(), vertexColors = !s && state.coloring !== 'sides';
  frontMesh.visible = backMesh.visible = !s; soapMesh.visible = s; sortedFor = null;
  renderer.setClearColor(0xffffff, 1); scene.background = s ? envTexture : null;
  renderer.outputEncoding = s ? THREE.sRGBEncoding : THREE.LinearEncoding;   // the film and the metal are lit in linear light and encoded for the screen
  ambient.intensity = s ? 0.25 : 0.55;
  if (wireGroup) wireGroup.children.forEach((m, k) => { m.material = s ? wireMetal() : tubeMaterials[k]; });
  for (const m of [frontMaterial, backMaterial]) { m.vertexColors = vertexColors; m.needsUpdate = true; }
  frontMaterial.color.set(vertexColors ? 0xffffff : FRONT); backMaterial.color.set(vertexColors ? 0xffffff : BACK);
  if (vertexColors) {
    const col = surfaceGeom.attributes.color.array, V = mesh.pos.length / 3;
    if (state.coloring === 'parts') {
      for (let v = 0; v < V; v++) { const k = mesh.kind[v], c = k < mesh.n ? PART_COLORS[k % PART_COLORS.length] : [0.78, 0.78, 0.78]; col[3 * v] = c[0]; col[3 * v + 1] = c[1]; col[3 * v + 2] = c[2]; }
    } else {
      const H = Minimal.meanCurvature(mesh).H, sorted = Float64Array.from(H).sort(), scale = sorted[Math.floor(0.95 * (sorted.length - 1))] || 1;   // the 95th percentile saturates
      for (let v = 0; v < V; v++) { const t = Math.min(1, H[v] / scale); col[3 * v] = 0.95 * t + 0.92 * (1 - t); col[3 * v + 1] = 0.35 * t + 0.92 * (1 - t); col[3 * v + 2] = 0.2 * t + 0.92 * (1 - t); }
    }
    surfaceGeom.attributes.color.needsUpdate = true;
  }
  if (s) computeThickness();
  requestRender();
}
function applyOpacity() {
  for (const m of [frontMaterial, backMaterial]) { m.opacity = state.opacity; m.depthWrite = state.opacity >= 1; m.wireframe = state.wireframe; }
  soapMaterial.wireframe = state.wireframe; soapMaterial.envMapIntensity = state.envbright;
  soapUniforms.uSoapOpacity.value = state.soapopacity;
  requestRender();
}
// The film's thickness, in nanometres, per vertex: a drainage profile in the axis direction (thick below, thin above,
// thick at the wire, the Plateau border), a gentle unevenness, then one implicit diffusion step on the surface's own
// Laplacian so that it is smooth over many edges.
function computeThickness() {
  if (!mesh || !surfaceGeom) return;
  const V = mesh.pos.length / 3, p = mesh.pos, dmin = state.soapmin, dmax = Math.max(state.soapmax, state.soapmin + 10);
  let zmin = Infinity, zmax = -Infinity; for (let v = 0; v < V; v++) { zmin = Math.min(zmin, p[3 * v + 2]); zmax = Math.max(zmax, p[3 * v + 2]); }
  const h0 = (zmin + zmax) / 2, sc = 0.25 * (zmax - zmin) + 1e-9, d0 = new Float64Array(V);
  for (let v = 0; v < V; v++) {
    const x = p[3 * v], y = p[3 * v + 1], z = p[3 * v + 2];
    const sig = 1 / (1 + Math.exp((z - h0) / sc)), eta = Math.sin(3.1 * x + 1.7) * Math.sin(2.3 * y - 0.4) * Math.cos(2.9 * z + 0.9) + 0.5 * Math.sin(5.3 * x - 2.1 * y + 1.2 * z);
    d0[v] = mesh.fixed[v] ? dmax : dmin + (dmax - dmin) * sig + 0.05 * (dmax - dmin) * eta;
  }
  if (!mesh.topo) Minimal.prepare(mesh);
  const l = Minimal.meanEdgeLength(mesh), r = Minimal.solve(mesh, d0, { dims: 1, dt: 8 * l * l, positive: true, tol: 1e-6 });
  const a = surfaceGeom.attributes.thick.array;
  for (let v = 0; v < V; v++) a[v] = Math.max(dmin, Math.min(dmax, r.x[v]));
  surfaceGeom.attributes.thick.needsUpdate = true; requestRender();
}
function rebuildDecorations() {
  disposeObject(axesGroup); axesGroup = buildAxes(sizeRadius * 0.9); axesGroup.visible = state.axes; group.add(axesGroup);
  buildWire(); requestRender();
}

// ------------------------------------------------------------------ the pipeline: taming, then settling the film
let ticking = false;
function schedule() { if (!ticking) { ticking = true; setTimeout(tick, 0); } }
function tick() {
  ticking = false;
  if (!mesh || phase === 'idle') return;
  try {
    if (phase === 'taming') tameTick(); else if (phase === 'settling') settleTick();
  } catch (e) { phase = 'idle'; setStatus(`<span class="err">${esc(e.message)}</span>`); console.error(e); showButtons(); return; }
  showButtons();
  if (phase !== 'idle') schedule();
}
function wireOptions() { return { alpha: state.alpha, K: Math.pow(10, state.kh), H: 1, dclose: state.dclose }; }
function initWire() { wire = mesh ? Wire.init(mesh, wireOptions()) : null; if (wire) wire.meshEdge = mesh.edge0; }
function startTaming() {
  if (!mesh || !wire) return;
  Object.assign(wire.o, wireOptions()); Wire.restart(wire); wire.steps = 0; history = []; tameTicks = 0;
  phase = 'taming'; remarks = remarks.filter(r => !/^taming|^the film/.test(r)); setRemarks(remarks); schedule(); showButtons();
}
function startSettling() {
  if (!mesh) return;
  settle = Minimal.settleInit(mesh, { every: state.every, tangential: state.tangential }); areas = [];
  phase = 'settling'; schedule(); showButtons();
}
// A history of the wires the film still followed, one every few ticks.  A tick that pinches the film is undone and
// taming stops there; if the film then cannot settle on that wire either, the history is walked back until it can,
// so the wire ends as far open as the film can follow it.
const HISTORY_EVERY = 4, HISTORY_MAX = 12;
let history = [], tameTicks = 0;
function tameTick() {
  const snap = Wire.snapshot(wire, mesh);
  const r = Wire.carry(wire, mesh, Minimal, { tangential: state.tangential });
  lastStats = r;
  if (r.pinched) {
    Wire.rollback(wire, mesh, snap, Minimal);
    remarks.push(`taming stopped after ${wire.steps} steps: the next of them began to pinch the film, a handle of the surface closing, and was undone`);
    setRemarks(remarks); buildSurfaceObjects(); buildWire(); keepInView(); startSettling(); return;
  }
  if (tameTicks++ % HISTORY_EVERY === 0) { history.push(snap); if (history.length > HISTORY_MAX) history.shift(); }
  buildSurfaceObjects(); buildWire(); keepInView(); showStatus();
  if (Wire.settled(wire) || wire.steps >= state.maxsteps) startSettling();
}
function settleTick() {
  for (let k = 0; k < state.perframe && phase === 'settling'; k++) {
    const s = Minimal.settleRound(mesh, settle); lastStats = s; areas.push(s.area);
    if (s.remesh) buildSurfaceObjects(); else updateSurfacePositions();
    if (s.pinched) {
      // the film cannot settle on this wire: go back to an earlier one it could follow, and settle there instead
      const back = history.pop();
      if (back) { Wire.rollback(wire, mesh, back, Minimal); buildSurfaceObjects(); buildWire(); startSettling(); }
      else { phase = 'idle'; remarks.push('the film pinches here: a handle of the surface is closing, so no film of this genus sits on this wire'); setRemarks(remarks); }
      break;
    }
    if (s.harmonic && areas.length > 5 && (areas[areas.length - 6] - s.area) / s.area < Math.pow(10, state.stop)) { phase = 'idle'; break; }
    if (settle.round > 400) { phase = 'idle'; break; }
  }
  showStatus();
  if (phase === 'idle' && soap()) computeThickness();
}
function setStatus(html) { $('status').innerHTML = html; }
function showStatus() {
  if (!mesh) { setStatus(''); return; }
  const A = Minimal.area(mesh).toFixed(3), V = mesh.pos.length / 3;
  if (phase === 'taming') setStatus(`taming: step ${wire.steps}, wire ×${(Wire.length(wire) / wire.L0).toFixed(2)}, area ${A}`);
  else if (phase === 'settling') setStatus(`settling the film: round ${settle.round}${lastStats && lastStats.harmonic ? ' (harmonic)' : ''}, area ${A}`);
  else { const H = Minimal.meanCurvature(mesh); setStatus(`settled: area ${A}, rms ${mixed('$|H|$')} ${(H.rms * mesh.params.R).toPrecision(2)}, ${V} vertices${wire && wire.steps ? `, wire ×${(Wire.length(wire) / wire.L0).toFixed(2)}` : ''}`); }
}
function showButtons() { $('tame').textContent = phase === 'taming' ? '❚❚ Pause' : '▶ Tame wire'; }
function resetSurface() {
  if (!scaffold) return;
  phase = 'idle';
  mesh.pos = scaffold.pos.slice(); mesh.tri = scaffold.tri.slice(); mesh.kind = scaffold.kind.slice(); mesh.fixed = scaffold.fixed.slice(); mesh.loops = scaffold.loops.map(l => l.slice());
  Minimal.prepare(mesh); initWire(); lastStats = null; areas = []; remarks = remarks.filter(r => !/^taming|^the film/.test(r)); setRemarks(remarks);
  sizeRadius = scaffold.sizeRadius; buildSurfaceObjects(); buildWire(); showStatus(); showButtons();
  if (state.auto) startTaming();
}

// ------------------------------------------------------------------ the rendering: soap film or rubber sheet
// "Soap film" and "Rubber" are the two pictures of the same surface: the film as it would look, or the two-sided
// rubber sheet (green front, red back) that shows it is orientable.  The Display tab has the other colorings.
function setColoring(name) {
  state.coloring = name; $('coloring').value = name;
  for (const [id, on] of [['soap', name === 'soap'], ['rubber', name !== 'soap']]) { $(id).classList.toggle('active', on); $(id).setAttribute('aria-pressed', String(on)); }
  applyColoring();
}
function setSoap(on) { setColoring(on ? 'soap' : 'sides'); }

// ------------------------------------------------------------------ text
const T = tex => katex.renderToString(tex, { throwOnError: false, output: 'html' });
const esc = t => String(t).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
const mixed = str => String(str).split('$').map((part, i) => i % 2 ? T(part) : esc(part)).join('');
for (const el of document.querySelectorAll('.tex')) katex.render(el.textContent, el, { throwOnError: false, output: 'html' });
function setInfo(html, cls) { const el = $('info'); el.innerHTML = html; el.className = cls || ''; }
function setRemarks(notes) { const el = $('remarks'); el.innerHTML = notes.map(mixed).join('<br>'); el.hidden = !notes.length; }
let buildToken = 0;
// The info line: the braid and everything read off it, and, once the knot table has arrived, the knot's name and
// its genus.  Called again when it does.
function renderInfoLine(word, bd, label) {
  const desc = [];
  if (record) desc.push(`<b>${T(prettyName(record.name))}</b>${label ? ' = ' + esc(label) : ''}`); else if (label) desc.push(`<b>${esc(label)}</b>`);
  desc.push(T(Seifert.wordTeX(word)) + (word.length ? ` on ${bd.n} strands` : ' (one strand)'));
  desc.push(`${bd.c} crossing${bd.c === 1 ? '' : 's'}, ${bd.mu} component${bd.mu === 1 ? '' : 's'}`);
  desc.push(T(`\\chi = ${bd.chi},\\ g = ${bd.genus}`));
  if (word.length) {
    const a = Seifert.alexander(word), sg = Seifert.signature(word);
    desc.push(T(`\\Delta(t) = ${Seifert.alexanderTeX(a.coeffs)}`));
    desc.push(T(`\\det = ${a.det},\\ \\sigma = ${sg.signature}`) + (sg.nullity ? mixed(` (nullity $${sg.nullity}$)`) : ''));
    const bound = Math.ceil(a.degree / 2);
    if (bd.mu === 1) {
      if (record && record.genus !== undefined) desc.push(record.genus === bd.genus ? `<span class="ok">this surface has the knot's genus ${record.genus}</span>` : `knot genus ${record.genus} (KnotInfo): this surface has ${bd.genus - record.genus} extra handle${bd.genus - record.genus === 1 ? '' : 's'}`);
      else if (bound === bd.genus) desc.push(`<span class="ok">minimal genus (${mixed('$\\deg\\Delta = 2g$')})</span>`);
      else desc.push(mixed(`genus $\\ge ${bound}$ from $\\Delta$`));
    }
    if (bd.positive || bd.negative) desc.push(`${bd.positive ? 'positive' : 'negative'} braid`);
  }
  setInfo(desc.join(' | '));
}
const prettyName = name => name.replace(/^K(\d+)([an]?)_(\d+)$/, (m, a, b, c) => `${a}${b ? b : ''}_{${c}}`).replace(/^L(\d+)([an])(\d+)((?:_\d+)+)$/, (m, a, b, c, d) => `L${a}${b}${c}\\{${d.slice(1).split('_').join(',')}\\}`);

// ------------------------------------------------------------------ the shipped films (data/films/, made by data/precompute.mjs)
// A film is used when the parameters it was computed with are the ones in force; `angular` is not compared, as it
// sets the scaffold's resolution only and a shipped mesh carries its own.  Anything else runs the pipeline live.
let films = null;
async function loadFilms() {
  if (films !== null) return films;
  try {
    const r = await fetch('data/films/index.json');
    films = r.ok ? await r.json() : false;
    if (films && films.version !== Precomputed.VERSION) films = false;
  } catch (e) { films = false; }
  return films;
}
function filmMatches(m) {
  if (!m || !m.params) return false;
  const p = m.params;
  return p.spacing === state.spacing && p.bandWidth === state.bandwidth && p.bulge === state.bulge && p.round === state.round
      && p.alpha === state.alpha && p.K === Math.pow(10, state.kh) && p.dclose === state.dclose && p.maxsteps === state.maxsteps;
}
// the braid word of a shipped film under this KnotInfo name, without the knot table
async function filmWord(name) {
  const m = await loadFilms();
  if (!filmMatches(m) || !m.names) return null;
  const key = m.names[name];
  return key ? (key === '()' ? [] : key.split(',').map(Number)) : null;
}
// the film for a braid word, decoded, or null
async function fetchFilm(word) {
  const m = await loadFilms();
  if (!filmMatches(m)) return null;
  const entry = m.films[Precomputed.key(word)];
  if (!entry) return null;
  try {
    const r = await fetch('data/films/' + entry.file);
    if (!r.ok) return null;
    return Precomputed.decode(await r.arrayBuffer());
  } catch (e) { console.warn('the shipped film could not be read, computing it here:', e.message); return null; }
}

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
  $('busy').hidden = false; phase = 'idle'; remarks = []; setRemarks(remarks); await nextFrame();
  try {
    const parsed = Seifert.parseBraid(text);
    if (parsed.error) throw new Error(parsed.error);
    // a name is resolved from the shipped manifest when it can be, so that an example needs no knot table
    let word = parsed.word, label = parsed.label || null; record = null;
    if (parsed.name) {
      word = await filmWord(parsed.name);
      if (!word) {
        await loadTable(); record = tableByName.get(parsed.name) || table.find(r => r.name.startsWith(parsed.name + '_')) || null;
        if (!record) throw new Error(`no ${parsed.name} in the table (KnotInfo knots to 12 crossings, LinkInfo links to 11)`);
        word = record.braid;
      }
    }
    const bd = Seifert.braidData(word);
    if (!bd.connected) throw new Error(`σ${[...Array(bd.n - 1).keys()].map(i => i + 1).find(i => !word.some(g => Math.abs(g) === i))} never occurs: the closed braid is split and the surface would be disconnected`);
    if (bd.n * bd.c > 6000) throw new Error('too many strands and crossings to draw');
    mesh = Seifert.buildSurface(word, { spacing: state.spacing, bandWidth: state.bandwidth, bulge: state.bulge, angular: state.angular, round: state.round });
    Minimal.prepare(mesh); mesh.edge0 = Minimal.meanEdgeLength(mesh); mesh.attrs = ['kind'];
    sizeRadius = Math.max(mesh.params.R + mesh.params.b + mesh.params.W, 0.6 * bd.n * state.spacing + 0.5);
    scaffold = { pos: mesh.pos.slice(), tri: mesh.tri.slice(), kind: mesh.kind.slice(), fixed: mesh.fixed.slice(), loops: mesh.loops.map(l => l.slice()), sizeRadius };
    lastStats = null; areas = []; history = []; tameTicks = 0;
    const film = opts.live ? null : await fetchFilm(word);     // a shipped film, if this is one of them
    if (film) {
      Precomputed.install(mesh, film, Seifert, Minimal);
      initWire(); wire.L0 = film.L0; wire.steps = film.steps; wire.energies = [];
      sizeRadius = Math.max(sizeRadius, wireRadius() * 1.05);
    } else initWire();
    buildSurfaceObjects(); buildGhost(); rebuildDecorations(); if (!opts.keepView) resetView();
    // the info line; the name and the knot genus come from the table, which is fetched behind the picture
    const token = ++buildToken;
    renderInfoLine(word, bd, label);
    if (!record) loadTable().then(() => {
      if (token !== buildToken || !mesh) return;
      record = tableByWord.get(JSON.stringify(word)) || null;
      renderInfoLine(word, bd, label);
    }).catch(() => {});
    remarks.push(`scaffold: ${scaffold.tri.length / 3} triangles, wire of ${scaffold.loops.reduce((s, l) => s + l.length, 0)} points${state.round ? `, corners rounded at $${state.round}R$` : ''}` +
                 (film ? `; the film shipped with the app, ${film.steps} taming steps` : ''));
    if (film && film.pinched) remarks.push('the film pinched while it was computed: a handle of the surface closed, so what is drawn is not a surface of this genus');
    setRemarks(remarks); showStatus(); showButtons();
    if (!film && state.auto) startTaming();
  } catch (e) {
    mesh = null; scaffold = null; record = null; wire = null; disposeObject(frontMesh); disposeObject(backMesh); disposeObject(soapMesh); disposeObject(ghostMesh); disposeObject(wireGroup); frontMesh = backMesh = soapMesh = ghostMesh = wireGroup = null;
    setInfo(`<span class="err">${mixed(e.message)}</span>`, ''); showStatus();
  } finally { $('busy').hidden = true; requestRender(); }
}
function hashFor(text) { return '#' + encodeURIComponent(text); }
function updateHash(text) { try { history.replaceState(null, '', hashFor(text)); } catch (e) {} }
function shareLink() { return location.origin + location.pathname + location.search + (lastText ? hashFor(lastText) : location.hash); }

// ------------------------------------------------------------------ UI wiring
const EXAMPLES = [
  ['trefoil  3₁ = T(2,3)', '3_1'], ['figure-eight  4₁', '4_1'], ['cinquefoil  5₁ = T(2,5)', '5_1'], ['5₂', '5_2'], ['stevedore  6₁', '6_1'], ['6₂', '6_2'], ['6₃', '6_3'],
  ['8₁₉ = T(3,4)', '8_19'], ['10₁₂₄ = T(3,5)', '10_124'], ['12n₂₄₂ = P(−2,3,7)', '12n242'], ['Conway knot  11n₃₄', '11n34'],
  ['Hopf link', 'hopf'], ['Whitehead link', 'whitehead'], ['Borromean rings', 'borromean'], ['T(4,5)', 'T(4,5)'], ['T(2,7)', 'T(2,7)'], ['T(3,6)  (3 components)', 'T(3,6)'],
  ['a 4-strand braid', '1 2 3 -1 2 -3 1 2'], ['unknot, one disk', '()'],
];
for (const [label, value] of EXAMPLES) { const o = document.createElement('option'); o.value = value; o.textContent = label; $('examples').appendChild(o); }
$('examples').addEventListener('change', e => { if (e.target.value) build(e.target.value); e.target.value = ''; });
$('build').addEventListener('click', () => build($('input').value));
$('input').addEventListener('keydown', e => { if (e.key === 'Enter') build($('input').value); });
$('tame').addEventListener('click', () => { if (phase === 'taming') { phase = 'idle'; showButtons(); showStatus(); } else startTaming(); });
$('soap').addEventListener('click', () => setColoring('soap'));
$('rubber').addEventListener('click', () => setColoring('sides'));
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
const rebuild = kind => { if (kind === 'change' && lastText && mesh) build(lastText, { keepView: true }); };   // changed parameters no longer match a shipped film, so build() computes live
bindRange('spacing', 'spacing', v => v.toFixed(2), rebuild);
bindRange('bandwidth', 'bandwidth', v => v.toFixed(1), rebuild);
bindRange('bulge', 'bulge', v => v.toFixed(1), rebuild);
bindRange('angular', 'angular', v => v, rebuild);
bindRange('round', 'round', v => v.toFixed(2), rebuild);
bindCheck('auto', 'auto', () => {});
$('rebuild').addEventListener('click', () => { if (lastText) build(lastText, { keepView: true, live: true }); });
$('alpha').value = state.alpha; $('alpha').addEventListener('change', e => { state.alpha = Number(e.target.value); });
bindRange('kh', 'kh', v => Math.pow(10, v).toPrecision(2), () => {});
bindRange('dclose', 'dclose', v => v.toFixed(2), () => {});
bindRange('maxsteps', 'maxsteps', v => v, () => {});
bindRange('every', 'every', v => v, () => { if (settle) settle.every = state.every; });
bindRange('tangential', 'tangential', v => v.toFixed(2), () => { if (settle) settle.tangential = state.tangential; });
bindRange('perframe', 'perframe', v => v, () => {});
bindRange('stop', 'stop', v => '10^' + v, () => {});
$('coloring').addEventListener('change', e => setColoring(e.target.value));
setColoring(state.coloring);
bindRange('opacity', 'opacity', v => v.toFixed(2), applyOpacity);
bindCheck('wireframe', 'wireframe', applyOpacity);
bindRange('soapmin', 'soapmin', v => v, kind => { if (kind === 'change' && soap()) computeThickness(); });
bindRange('soapmax', 'soapmax', v => v, kind => { if (kind === 'change' && soap()) computeThickness(); });
bindRange('soapopacity', 'soapopacity', v => v.toFixed(2), applyOpacity);
bindRange('envbright', 'envbright', v => v.toFixed(1), applyOpacity);
$('wiremetal').value = state.wiremetal; $('wiremetal').addEventListener('change', e => { state.wiremetal = e.target.value; if (wireGroup && soap()) wireGroup.children.forEach(m => { m.material = wireMetal(); }); requestRender(); });
bindCheck('wire', 'wire', () => { if (wireGroup) wireGroup.visible = state.wire; });
bindRange('thick', 'thick', v => v.toFixed(3), kind => { if (kind === 'change') buildWire(); });
bindCheck('ghost', 'ghost', () => { if (ghostMesh) ghostMesh.visible = state.ghost; });
bindCheck('axes', 'axes', () => { if (axesGroup) axesGroup.visible = state.axes; });
$('reset-view').addEventListener('click', () => { if (mesh) sizeRadius = Math.max(sizeRadius * 0.5, wireRadius() * 1.05); resetView(); });
$('snapshot').addEventListener('click', () => { renderer.render(scene, camera); const a = document.createElement('a'); a.href = renderer.domElement.toDataURL('image/png'); a.download = 'seifert-surface.png'; document.body.appendChild(a); a.click(); a.remove(); });
$('share').addEventListener('click', async () => { const url = shareLink(); try { await navigator.clipboard.writeText(url); $('share').textContent = 'Copied'; setTimeout(() => $('share').textContent = 'Copy link', 1200); } catch (e) { prompt('Link:', url); } });
setTimeout(() => { $('hint').hidden = true; }, 9000);
window.SEIFERT_DEBUG = { state, get mesh() { return mesh; }, get wire() { return wire; }, get phase() { return phase; }, build, startTaming, startSettling, setSoap, render: () => renderer.render(scene, camera), canvas: renderer.domElement, Seifert, Minimal, Wire };
const initial = decodeURIComponent((location.hash || '').slice(1));
build(initial || '12n242');
})();
