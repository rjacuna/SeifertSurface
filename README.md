# SeifertSurface

A Seifert surface for a knot or link `K ⊂ ℝ³` is a compact oriented surface with boundary `K`. This project draws
one for any closed braid as the soap film on a tame wire. It builds the classical surface first, Seifert's algorithm
applied to the closed braid: one disk per strand stacked along the braid axis, one half-twisted band per crossing
(the Bennequin surface), a scaffold whose boundary is the knot. Then, on its own, it tames that wire as a knot with
Scharein's KnotPlot forces, the ones van Wijk and Cohen use, with the film following, and settles the film on the
still wire by solving the minimal surface equation on it. What is shown is always the film; bend the wire by hand
and the film follows, tame it again when it is bent too much. Alongside, it reads the invariants off the surface:
the Seifert matrix, the Alexander polynomial, the signature, the genus. A soap-film rendering shows it as a real
film would look, thin-film interference and all.

It is a single-page web app in the shape of [EllipticCurve3D](https://github.com/rjacuna/EllipticCurve3D):
`web/index.html` with `js/`, `css/` and `data/`, three.js and KaTeX vendored, no build step; the mathematics in three
dependency-free modules that also run under Node for the tests.

It is live at **https://rjacuna.github.io/SeifertSurface/**; to run it locally, serve the folder and open it.

```bash
python3 -m http.server -d ~/Projects/SeifertSurface/web 8766      # then http://localhost:8766/
```

## Layout

| path | what it is |
|---|---|
| `web/index.html` | the page |
| `web/js/seifert.js` | braid words, the closed braid, the Bennequin surface as a triangle mesh, the Seifert matrix, `Δ(t)`, signature |
| `web/js/minimal.js` | discrete minimal surfaces: cotangent Laplacian, conjugate gradients, the harmonic step and mean curvature flow, the harmonic extension that carries the surface with the wire, Delaunay flips, tangential smoothing, isotropic remeshing, the settling schedule with pinch detection, diagnostics, linking numbers |
| `web/js/wire.js` | taming the wire: Scharein's KnotPlot relaxation of a coarse copy of the boundary loops (springs, repulsion, damped Euler, strands kept apart), the mesh boundary interpolated from it, the carry of the surface |
| `web/js/precomputed.js` | the shipped films: the binary format, and installing one into a mesh |
| `web/js/app.js` | viewer and UI |
| `web/js/vendor/` | three.js r128, OrbitControls, KaTeX 0.16.11 (licenses alongside) |
| `web/css/` | the stylesheet (EllipticCurve3D's), KaTeX's with its fonts |
| `web/data/knots.json` | KnotInfo's braid words for the 2,961 knots through 12 crossings and LinkInfo's for the 4,188 links through 11, with the knots' genera |
| `web/data/braid.sage` | `sage web/data/braid.sage K12n242` (or a PD code, or `DT:[...]`): a braid word for any knot, by SnapPy |
| `web/data/films/` | the examples' films, computed once and shipped, with `index.json` |
| `web/data/precompute.mjs` | `node web/data/precompute.mjs`: recomputes them |
| `web/test/test-seifert.mjs` | `node web/test/test-seifert.mjs`: 206 checks against Sage's values in `vectors.json` (made by `make_vectors.sage`), plus the mesh, the wire, the solver and the whole pipeline |

The web app takes a braid word (`1 2 -1 2`, `s1 s2 s1^-1 s2`, `σ₁σ₂σ₁⁻¹σ₂`, `abAb`, `(1 2)^5`), a torus link `T(p, q)`,
or a name from the tables (`3_1`, `12n242`, `11n_34`, `L6a4`, `L2a1{1}`, and the aliases `trefoil`, `figure-eight`,
`hopf`, `whitehead`, `borromean`, `p(-2,3,7)`). It opens on the (−2, 3, 7)-pretzel knot. The info line gives the
braid, the strands, crossings and components,
`χ` and the genus of the surface, `Δ(t)`, `det`, `σ`, and whether the surface has the knot's genus (from KnotInfo, or
from `deg Δ = 2g`). The examples appear at once, their films having been computed in advance and shipped (below);
anything else is tamed and settled in the browser, a few seconds of animation, and so is an example once a setting
is changed. The status shows
the phase, the area, the root mean square of the discrete mean curvature `|H|` over the interior (the residual of the
equation), the vertex count and the wire's growth. **Tame wire** runs the taming again from the wire as it is (for
after bending it), **Soap film** and **Rubber** are the two pictures of the surface, **Reset** goes back to the
scaffold and starts over. The drawer has four tabs: **Surface** (the scaffold: disk spacing, band width and bulge,
resolution, the rounding of the wire's corners, whether to tame on build), **Wire** (weak or strong repulsion, its
strength, the clearance between strands, the most steps), **Film** (how often to remesh,
tangential smoothing, rounds per tick, when the film counts as settled), **Display** (two-sided coloring, or by disk
and band, or by mean curvature, or the soap film with its thickness range, opacity, light and the wire's metal;
wireframe, the wire and its thickness, the scaffold as a ghost, axes, PNG, a shareable link `#word`).

Conventions are Sage's: `σᵢ` is the positive crossing, so the closure of `σ₁³` is the right-handed trefoil, signature
`−2`. The Seifert matrix is computed by the algorithm of Collins ([Col2013]) as Sage's `Link.seifert_matrix` does,
line for line, and the tests check that the matrices agree exactly on 261 braids, the Alexander polynomials and
signatures with them. The chirality of a *named* knot is whatever the table's diagram is: KnotInfo's `3_1` has
signature `−2`, SnapPy's `3_1` is its mirror, and `braid.sage` prints the signature so one can tell. Negating a word
mirrors the knot.

## Seifert surfaces from braids

Every link is the closure of a braid (Alexander), and Seifert's algorithm on a closed braid diagram is especially
simple: the oriented smoothing of a crossing `σᵢ` joins strand `i` to strand `i` and `i+1` to `i+1`, so the Seifert
circles are the `n` strands themselves, nested around the braid axis. Realised in space, they are `n` parallel disks
stacked along the axis, all oriented the same way, and each crossing `σᵢ^{±1}` is a band with a positive or negative
half twist joining the rims of disks `i` and `i+1`. This is the surface Bennequin used in [Ben1983]; for a positive
braid it has minimal genus (Bennequin's inequality is an equality), and in general

    χ = n − c,        g = (c − n + 2 − μ) / 2

for a braid on `n` strands with `c` crossings whose closure has `μ` components. `seifert.js` builds it as a triangle
mesh: disk `p` at height `(p − (n−1)/2) h`, radius `R = 1`, rings of points whose counts grow with the radius so the
edges have one length everywhere; the rim of each disk has `M = cq` points, `q` per crossing sector; band `k` is
attached to `q_b + 1` rim points of each of its two disks, centred on the sector's middle angle `θ_k`, and is the
ruled strip

    angle = θ_k + v (w/2) cos φ(u),    radius = R + b sin(πu) + v (W/2) sin φ(u),    z = z_i + u h,

`u ∈ [0, 1]` up the band, `v ∈ [−1, 1]` across it, `φ(u) = ±π (3u² − 2u³)` a smooth half turn. At `u = 0` it is the
rim arc of the lower disk, at `u = 1` the rim arc of the upper disk traversed backwards, at `u = ½` a radial segment:
the two edges of the band swap heights over the sector and pass each other radially, one in front of the other, a
braid crossing seen from outside the cylinder. The mesh is consistently oriented; its boundary loops, followed as the
triangles induce, are the closed braid with all strands going the same way round. The test that pins the sign is the
linking number of the two boundary loops of `σ₁²`: `+1`.

The boundary polygon then has corners where a rim arc turns up a band edge. It is rounded (binomial smoothing at a
radius `0.1 R`, then resampled at equal arclength, so the corners do not bunch the points). This is the wire; the
film can be computed on it as it is, but it is the boundary of a scaffold, not a shape anyone would bend a wire into.

## Taming the wire

What makes van Wijk and Cohen's 2006 pictures look like knots rather than scaffolds is Scharein's relaxation from
KnotPlot, which they adopt in §5.3, and `wire.js` applies it to the boundary loops of the mesh: every point of the
wire is a unit mass, attracted by its two neighbours along the loop with `F_a(r) = H r^{1+β}` and repelled by every
other point of the wire with `F_r(r) = K r^{−(2+α)}`, distances in units of the initial spacing `r_a` (`β = 1`,
`α = 0` by default, inverse-square repulsion; `α = 4` is KnotPlot's strong, short-range repulsion). Damped explicit
Euler, `v ← (1 − γ) v + F dt`, `x ← x + v dt`, the displacement clamped to `d_max = 0.25 r_a`; and a move is refused
when it would bring the point within `d_close = 0.5 r_a` of a segment of the wire not adjacent to it, which with
`d_close > d_max` keeps the knot from ever passing through itself, so the knot type is preserved (the components of a
link repel each other by the same forces). The step decays, `dt ← (1 − μ) dt`, so a cycle settles; a new cycle starts
from the full step. The wire grows while it relaxes, to about 1.5 times its length, the scale at which the springs
and the repulsion balance; nothing depends on that scale. The energy falls monotonically and the corners go: on the
trefoil the total turning of the wire drops from 5.6 turns to 2.2 in a few thousand steps of 2 ms each, and the
stacked rims open into a round three-dimensional trefoil.

Two things are done differently. The forces cost `N²`, and the mesh boundary has several hundred points per loop,
so the dynamics runs on a coarse copy of each loop, resampled at equal arclength one tenth of a radius apart, and the
mesh boundary is interpolated back from it by a Catmull–Rom spline at equal arclength, its first point anchored at
the point of the new curve nearest to where it was (the coarse points slide along the loop as they relax, and the
mesh boundary must not slide with them, or the surface next to it is sheared without end). And KnotPlot's `d_close`
is a fraction of the spacing; here it is two spacings: a film between two strands closer than a few mesh edges
cannot be resolved, and pinches. With the original value the pretzel knot's wire brought two strands within two
thirds of a spacing and the film between them closed a handle.

The surface is not given springs of its own, as theirs is. Whenever the coarse wire has moved half a spacing, the
surface is carried along as a rubber sheet: the wire's displacement is extended harmonically into the interior (the
same Laplace solve as the film, with the displacement as boundary data and positive weights, so that no interior
vertex moves further than the wire does, `Minimal.extend`); the mesh is remeshed to its scaffold edge length scaled
by the wire's growth (`Minimal.remesh`: Botsch–Kobbelt isotropic remeshing, long interior edges split, short ones
collapsed under the link condition and a check that no triangle turns over, slivers and interior vertices of degree
three collapsed away, then Delaunay flips and tangential smoothing; the boundary polygon is never touched); and one
gentle round of the film is run (a few edge lengths squared of mean curvature flow with positive weights), so that
the surface tracks the film as the wire moves rather than drifting away from it as a rubber sheet. The full harmonic
step must not be used while the wire moves: it pulls interior vertices onto the wire where the surface wraps a bend
faster than any remeshing can collapse them.

A tick that pinches the film is undone, and taming stops there. Taming also keeps a short history of the wires the
film still followed, one every few ticks: if the film then cannot settle even on the wire taming stopped at, the
history is walked back until it can, so the wire ends as far open as the film is able to follow it. The 4-strand
braid in the menu needs this, its two components separating until the film between them closes; without the
rollback its film is nonsense, with it the film settles at 1,098 taming steps instead of 1,748.

Once the wire is still the film settles (`Minimal.settleRound`): implicit mean curvature flow from a step of one
edge length squared, doubled while a round moves no vertex more than half an edge, halved and the round undone when
the solve fails, positive weights and a remeshing every round while the step is small, the exact cotangent weights
once it is large, the harmonic step at the end; a vertex's mass is floored at a tenth of the mean so that a vertex
with a tiny area moves like the others and not at infinite speed. It counts as settled when the area has stopped
falling. If vertices bunch up that the remeshing cannot spread out again, a neck of the surface is closing: the film
is leaving the surface's isotopy class, a handle would be lost, and the app stops there and says so, both during
taming (the taming stops at that step) and while settling. The surface, like theirs, is not checked for
self-intersection; the wire is.

Deforming the wire by hand with the pointer was tried and taken out again: the film followed the dragged wire
correctly but did not render well while it did, and a hand drag has nothing to stop one strand from being pushed
through another. The same carry would serve it, so it is in the plan below.

Whether one tames first and then spans the wire, or carries the surface along as it is done here, the film at the
end is the same: it depends
only on the final wire and on the isotopy class carried along. The order matters for a wire that is not born from
the scaffold, a parametrised or hand-drawn knot: spanning a *given* wire needs Seifert's algorithm on its projection,
built in the wire's own geometry, which is in the plan below.

## The shipped films

Taming and settling take a few seconds, which is a few seconds of watching a scaffold turn into a knot. Worth
watching once; not every time one opens the page. So `web/data/precompute.mjs` runs the pipeline under Node for
every example of the menu, with the app's own modules and its default parameters, and writes each result to
`web/data/films/`: a header with the counts and the numbers the app shows (the taming steps, the wire's length
before and after, the scaffold's edge length, the area), then the vertex positions in single precision, the
triangles, and the scaffold kinds each vertex came from. 150 kB apiece, 2.8 MB in all, one file fetched per knot,
and the knot is on the screen in about 50 ms instead of the ten seconds it takes to compute. Everything else is recomputed when the file is read: the boundary loops from the triangles, the coarse wire
from those loops, the scaffold from the braid word, which is what **Reset** goes back to.

The manifest records the parameters the films were made with, and a film is used only while those are the
parameters in force; change the disk spacing, the repulsion or the clearance and the app computes the knot here
instead, as it does for any knot that is not among the examples. The one parameter not compared is the scaffold's
resolution, which sets the number of rim points only: the shipped mesh carries its own, having been remeshed. The
tests decode every shipped film and check that it is the settled film of the knot it claims, with the right Euler
characteristic, the right number of boundary loops, an oriented mesh with no degenerate triangles, a small mean
curvature and the area the manifest records.

## The soap-film rendering

The mathematics gives a surface of zero thickness; what makes a picture read as a soap film is a very thin layer
of water with Fresnel reflection, strong at grazing angles, a reflected environment, and thin-film interference whose
colour depends on the local thickness. So the app gives the surface a thickness field `d: M → (0, ∞)` in nanometres
(a drainage profile along the axis, thick below and thin above, thick at the wire where the Plateau border sits, a
gentle unevenness, then one implicit diffusion step on the surface's own Laplacian so that it varies over many
edges), as a vertex attribute of the mesh, and renders the surface with three.js's physically based material, its
specular light multiplied by the colour of a free-standing water film of that thickness. That colour is computed
spectrally once, into a lookup table: the film reflects the Airy fraction `R(λ) = 2ρ²(1 − cos δ) / (1 + ρ⁴ − 2ρ² cos δ)`
with `δ = 4π n d cos θ_f / λ` and `ρ = (n − 1)/(n + 1)` (the two reflections differ by π, so the thinnest film is
black, as a real one is), integrated over the visible wavelengths against the CIE colour matching functions and
converted to sRGB: Newton's series, black, grey, white, yellow, orange, red, violet, blue, green, and the pale
higher orders. Three.js's own iridescence was tried first and set aside: it models a film on a substrate of
refractive index 1.5, not a film in air, and gives a purple cast. No diffuse colour, a water-like Fresnel response
from the physically based model, a painted room as the environment (prefiltered for reflections) and as the
background, the film's opacity to what lies behind it rising with the Fresnel term at grazing angles, the reflected
light added unweighted, and the film's triangles sorted back to front for the current view whenever the camera or
the surface moves (its index buffer rewritten), so that where one sheet of the film lies behind another the blend
is in depth order rather than mesh order. In this mode the scene is lit in linear light and encoded to sRGB for the
screen, as physically based rendering needs (the two-sided view keeps the plain pipeline it always had), and the
environment is a studio: bright above, a dark floor, a warm key softbox, a cool fill, a long thin strip light. The
wire becomes gold, a full metal whose colour is gold's measured reflectance (1.00, 0.71, 0.29 in linear light) and
which therefore looks like gold only because it reflects that studio, its highlights and its dark floor; steel and
dark metal are the alternatives in the Display tab. It is rasterised, not
ray-traced: for an interactive viewer this is the right baseline, and the mesh can be exported to a path tracer for
a still. The thickness range, the opacity and the light are in the Display tab.

The Seifert form on `H₁` of this surface has a basis of loops that go up one band, along the upper disk, down the
next band of the same column and back along the lower disk; the linking numbers between them and their push-offs are
the entries `0, ±1` that Collins' rules give from the word alone ([Col2013] §3.3). `Δ(t) = det(V − tVᵀ)` is computed
exactly (BigInt Bareiss determinants at `t = 0, …, N`, interpolated), and printed with lowest degree 0 and a positive
leading coefficient, as KnotInfo does; `det = |Δ(−1)|`; `σ` is the signature of `V + Vᵀ`.

## Minimal surfaces: what is true

Can a Seifert surface be represented by a minimal surface, a soap film on the knot as a wire? Essentially yes, with
three qualifications.

**Existence.** For a smooth Jordan curve `Γ ⊂ ℝ³` the oriented Plateau problem has a solution among integral
currents (Federer–Fleming); an area-minimising 2-current in `ℝ³` is a smooth embedded surface in its interior
(Fleming, De Giorgi, Almgren: in codimension one, singularities start in dimension 7), and Hardt and Simon proved
regularity up to a `C^{1,α}` boundary [HS1979]: the minimiser is a smooth compact embedded *orientable* surface with
boundary `Γ`, and it is connected, since a closed component would be a compact minimal surface without boundary in
`ℝ³`, which does not exist. So every smooth knot bounds an embedded orientable area-minimising surface: a minimal
Seifert surface. This is the theorem behind the whole project, and it is a statement about integral currents, that is
about *oriented* surfaces with cancellation, which is exactly the class of Seifert surfaces.

**Which one depends on the wire, not the knot.** The minimiser is determined by the geometry of `Γ`, and its genus
need not be the genus of the knot. Almgren and Thurston [AT1977] built unknotted curves whose every embedded spanning
surface inside their convex hull has genus `≥ g`, for any `g`; a minimal surface lies in the convex hull of its
boundary, so the least-area surface of such an unknot has genus `≥ g`, not 0. The other way round, any Seifert surface
of a knotted wire has genus `≥ g(K)`, and the minimiser may have more. Nothing forces the area-minimising surface to
realise the knot genus; for a well-chosen wire one expects it to, and that is what the app lets one look at. The wire
may also bound several minimal surfaces, stable ones and unstable ones between them (Morse–Tompkins, Shiffman), and
a soap film realises only the stable ones; which one a numerical method finds depends on where it starts. Starting
from the Bennequin surface and never changing the mesh's topology, the app finds a minimal surface of the Bennequin
surface's genus, isotopic to it as long as the flow does not pass the surface through itself (which it does not
check). Existence of embedded minimal surfaces of *prescribed* genus spanning a curve is a harder theorem, proved by
Jost under conditions [Jos1986]; the general picture of an isotopy class realised by a least-area representative is
Meeks–Simon–Yau's, for closed surfaces in 3-manifolds.

**A physical soap film is not constrained to be orientable, or a surface.** A film is a Plateau-type minimal set:
smooth pieces meeting in threes at 120° along curves that meet in fours at tetrahedral points (Taylor [Tay1976]).
On a trefoil wire the film that forms depends on how the frame is withdrawn: the trefoil bounds a Möbius band with
three half twists, and films with triple junctions occur as well as the two-sided surface of genus 1 (Almgren's
*Plateau's problem* and Courant's experiments). Minimising area over *all* films (mod 2 currents, or Almgren's
`(M, ε, δ)`-minimal sets) can find a one-sided surface of smaller area than any Seifert surface. The oriented
minimiser above is the one that is always a Seifert surface, and a mesh-based method with fixed connectivity computes
in that class by construction.

**The equation.** A surface `X` is minimal iff its mean curvature vanishes, `H = 0`, iff its coordinate functions
are harmonic for its own Laplace–Beltrami operator, `Δ_X X = 2H n = 0`: three scalar elliptic equations, coupled and
quasilinear because the metric depends on `X`, with Dirichlet data the wire. In isothermal coordinates it is
literally `Δx = 0` plus conformality (Douglas and Radó solved it for disks; a knotted wire bounds no embedded disk, so
their solution is an immersed disk). The graph form `div(∇u/√(1+|∇u|²)) = 0` is of no use here, a Seifert surface of a
knotted wire being no graph. `minimal.js` uses the linearisation of Pinkall and Polthier [PP1993]: freeze the metric
of the current mesh, minimise the Dirichlet energy of the map to the next one, which is the cotangent Laplace equation

    Σ_j ½(cot α_ij + cot β_ij)(x_i − x_j) = 0     for every interior vertex i, x fixed on the wire,

a sparse symmetric positive definite system (the Dirichlet energy is a sum of squares of gradients whatever the sign
of individual weights), solved by conjugate gradients with a Jacobi preconditioner, warm-started; move the vertices;
repeat. The area decreases every round (`area(M_{k+1}) ≤ E_{M_k}(M_{k+1}) ≤ E_{M_k}(M_k) = area(M_k)`) and the fixed
points are the discrete minimal surfaces, every vertex at the cotangent-weighted mean of its neighbours. With a finite
step the same solve is implicit mean curvature flow, `(M + dt L) x_{k+1} = M x_k` (Desbrun, Meyer, Schröder, Barr),
which the Flow tab offers for watching. Vertices drift tangentially and triangles degenerate under either, so after
each round interior edges are flipped to the intrinsic Delaunay triangulation and each interior vertex is moved part
of the way towards the centroid of its neighbours within its tangent plane: Brakke's equiangulation and vertex
averaging [Bra1992]. The residual shown is the root mean square of the discrete `|H| = |Σ_j w_ij (x_i − x_j)| / 2A_i`
over the interior. On the trefoil at the default resolution (about 4,400 vertices) it falls below `0.1/R` away from
the wire in twenty rounds of 20–40 ms each; the catenoid test in `test-seifert.mjs` reaches the exact area within
0.5% from a cylinder.

Ken Brakke's Surface Evolver is the reference tool for soap films on wires and does all of this and much more, in C,
with a command language; it needs the starting surface of the right topology supplied, which is what `seifert.js`
produces. Nothing here competes with it; the point is a browser, a braid word, and one click.

## van Wijk's SeifertView

Jarke van Wijk and Arjeh Cohen made the two papers and the program this project starts from: *Visualization of the
genus of knots* (IEEE Visualization 2005, 567–574) and its extended version *Visualization of Seifert surfaces*
(IEEE TVCG 12(4), 2006, 485–496), implemented in SeifertView, a Windows-only executable, free for personal use, no
source, with a braid-word input (letters: uppercase a right-hand crossing, and a table of the knots through 10
crossings from Gittings' minimum braids). Their construction, in its default "stacked" style, is the one above: one
coaxial disk per strand along the braid axis, a twisted band per crossing stepping round the axis; three other styles
(split, flat, reduced) route the closure differently so that nested Seifert circles become two stacks or one plane.
Disks are ellipsoids with elliptical holes where bands attach, bands are Bézier tubes with a rotating frame, and the
knot is read off the mesh afterwards as the seam between the two sides.

They considered the soap film. Section 2.4 of both papers weighs treating the Seifert surface as a minimal surface
bounded by the knot and sets it aside: it needs a 3D embedding of the knot, a starting mesh of the right topology,
and an iterative, compute-intensive minimisation, whereas they wanted a deterministic construction from the abstract
braid. What they do instead is Catmull–Clark subdivision and vertex averaging (2005, with the remark that the
averaging "more or less" approaches a minimal surface, and a correction for the tubes it thins), and in 2006 the
physically based relaxation of Scharein's KnotPlot: the knot's vertices as point masses with a generalised Hooke
attraction between neighbours and an inverse-power repulsion between all others, the surface's vertices with the
attraction only, so the surface follows the knot (a steel rod and a rubber sheet, in their words), damped explicit
Euler with a decaying step, self-intersection of the knot prevented by rejecting close approaches. They say plainly
that this gives a smooth surface, not a minimal one, that the surface is not checked for self-intersection, and that
different starting styles relax to different results. Seifert surfaces from an arbitrary given 3D curve, and
minimal-genus surfaces, are their listed future work; the genus question was taken up by van Garderen and van Wijk
(Bridges 2013), who search braid moves for a minimal-genus presentation.

So the two objections to the soap film are exactly what this project supplies: the 3D embedding and the starting mesh
are the Bennequin surface itself, and the minimisation is one sparse linear solve per round, fast enough for a
browser; and their relaxation of the knot is kept, as the step that makes the wire worth spanning. Their pictures
remain the model for the rendering, and their letter convention is the reverse of the one here (`a` is `σ₁`, `A` its
inverse), so a SeifertView word is typed with its case swapped.

Soap films on knotted wires have been computed before, with Brakke's Surface Evolver: Brakke's own page *Soap films
on knots* has Evolver files for the trefoil and the figure-eight, the orientable films among them being numerically
computed minimal Seifert surfaces, beside Möbius films and films with triple lines; Coletti (Bridges 2024) does torus
knots and links and notes three distinct stable films on the trefoil; Stockrahm, Lahtinen, Kangas and Kotiuga (2019)
compute cut surfaces for Almgren–Thurston unknots by finite elements. Wang and Chern (SIGGRAPH 2021) minimise area
over *all* topologies at once, the boundary as a current and the problem convex, which is the way to compute the
Hardt–Simon minimiser itself rather than the minimal surface in one isotopy class; Parks (1992) and Brakke (1995) are
the theory of soap-film-like surfaces on knots.

## The (−2, 3, 7)-pretzel knot

The knot the app opens on. It is `12n_242` in the tables, the closure of the positive 3-braid `σ₁σ₂²σ₁²σ₂⁷`
(KnotInfo's word `1 2 2 1 1 2 2 2 2 2 2 2`), genus 5, signature −8, fibered, an L-space knot, and the one Fintushel
and Stern used for their surgery examples. Its Alexander polynomial is

    Δ(t) = t¹⁰ − t⁹ + t⁷ − t⁶ + t⁵ − t⁴ + t³ − t + 1 = L(−t),

Lehmer's polynomial `L(t) = t¹⁰ + t⁹ − t⁷ − t⁶ − t⁵ − t⁴ − t³ + t + 1` with `t ↦ −t`, which does not change the
Mahler measure: Lehmer's number `1.17628…` (Hironaka [Hir2001], who credits Reidemeister's book for `L(−x)`). The app
computes this `Δ` from the Bennequin surface, and the Bennequin surface of this positive braid has the knot's genus,
so its relaxed form is a minimal surface of genus 5 with the pretzel knot as boundary: the example `12n242` in the
menu. Silver and Williams relate the Mahler measure of Alexander polynomials to the growth of homology in cyclic
covers, and Seifert's theorem says every palindromic integer polynomial with `Δ(1) = ±1` is the Alexander polynomial
of some knot; the notebooks `Friedl's Algorithm.ipynb` in `~/Projects/Lehmer` build a Seifert matrix from the
coefficients of such a polynomial (Seifert's inductive construction), which is the input the planned
surface-from-a-Seifert-matrix builder below would take.

## Deployment

`.github/workflows/pages.yml` publishes `web/` to GitHub Pages on every push to `master` (the Pages source is
"GitHub Actions"), so the site is whatever the working tree's `web/` holds; there is no build step. The whole site is
about 2 MB, the largest file `web/data/knots.json` at 0.7 MB, so nothing is near GitHub's limits. Opening
`index.html` from the file system works for typed braid words, but browsers block the fetch of the knot table from
`file://`, so names such as `12n242` need the page served over http.

## Tests

```bash
node web/test/test-seifert.mjs          # 185 checks
sage web/test/make_vectors.sage         # regenerates web/data/knots.json and web/test/vectors.json (needs the KnotInfo database in Sage)
```

The checks: the shipped films, as above; parsing of every input form; strands, crossings, components and genus of closed braids; the Seifert
matrix equal to Sage's on 261 connected braids (17 named, 160 random on 2–5 strands, the 85 knots through 9
crossings), the Alexander polynomial and the signature with them, `deg Δ ≤ 2g(K) ≤ 2g(surface)` against KnotInfo's
genera; `Δ` of `12n_242` equal to `L(−t)`; the meshes oriented, with `χ = n − c`, the right number of boundary loops,
no degenerate triangles; the linking numbers of the Hopf link and `T(2,4)` from the built geometry (`+1`, `−1` for the
mirror, `2`); the solver on a bumped disk (flattens to the polygon's area, `H → 0`, boundary fixed), on the catenoid
(area within 0.5% of the exact catenoid, waist radius, monotone area under both the harmonic step and mean curvature
flow), and on the trefoil's surface (area decreases, `χ` and orientation survive the flips, the wire does not move,
the residual falls); the wire relaxation on a circle (stays round and evenly spaced, grows) and on the trefoil (the
energy falls monotonically, the total turning halves, no strand closer than it started, the surface follows and
stays sound, the film on the tamed wire relaxes); the whole pipeline on the trefoil and the pretzel knot (the wire
opens, no pinch, a settled film, `χ` and orientation kept through the remeshing, the wire still while the film
settles); a stress test of taming with film rounds.

## Planned

* **A surface from a Seifert matrix.** Given `V` (from `Friedl's Algorithm.ipynb`, or any `V` with `det(V − Vᵀ) = ±1`),
  a disk with `2g` bands whose twists and mutual linkings realise `V`, then the same relaxation: minimal surfaces for
  Lehmer's polynomial itself, `L(t)`, and for any Alexander polynomial one likes.
* **Deforming the wire by hand**, properly: a drag of the wire with the film carried along the way taming carries
  it, a check that no strand is pushed through another, and a rendering that holds up while it moves.
* **Wires that are not the Bennequin boundary.** A torus knot on a round torus, a Fourier knot, a wire drawn by hand:
  Seifert's algorithm on a generic projection of the given wire (crossings from the polygon, Seifert circles as
  chains of its arcs, disks spanning them at staggered depths, bands at the crossings), so that the surface is built
  in the wire's own geometry and the wire comes first. Van Wijk and Cohen list the same as future work.
* **Diagrams as input.** PD and DT codes through `braid.sage` (SnapPy's Vogel algorithm), as `SkeinA/data/braids.txt`
  already does for the census.
* **A self-intersection check** during the flow; order-independent transparency for the film only if sorted triangles ever prove insufficient (they are exact unless triangles intersect).
* **The global minimiser.** Wang–Chern's minimal currents for the area-minimising surface over all topologies, to compare with the fixed-topology one: when they differ, the wire's least-area surface is not the Bennequin surface's isotopy class.
* **Export** of the relaxed mesh (OBJ/STL) for printing, and GitHub Pages as for EllipticCurve3D.

## References

* [AT1977] F. J. Almgren, W. P. Thurston, *Examples of unknotted curves which bound only surfaces of high genus within their convex hulls*, Ann. of Math. 105 (1977) 527–538.
* [Ben1983] D. Bennequin, *Entrelacements et équations de Pfaff*, Astérisque 107–108 (1983) 87–161.
* [Bra1992] K. A. Brakke, *The Surface Evolver*, Experiment. Math. 1 (1992) 141–165.
* [Col2013] J. Collins, *An algorithm for computing the Seifert matrix of a link from a braid representation*, 2007, corrected 2013 (the algorithm Sage implements), https://webhomes.maths.ed.ac.uk/~v1ranick/papers/collinsseifert.pdf
* [HS1979] R. Hardt, L. Simon, *Boundary regularity and embedded solutions for the oriented Plateau problem*, Ann. of Math. 110 (1979) 439–486.
* [Hir2001] E. Hironaka, *The Lehmer polynomial and pretzel links*, Canad. Math. Bull. 44 (2001) 440–451.
* [Jos1986] J. Jost, *Existence results for embedded minimal surfaces of controlled topological type*, I–III, Ann. Scuola Norm. Sup. Pisa 13 (1986) 15–50, 401–426; 14 (1987) 165–167.
* [PP1993] U. Pinkall, K. Polthier, *Computing discrete minimal surfaces and their conjugates*, Experiment. Math. 2 (1993) 15–36.
* [Tay1976] J. E. Taylor, *The structure of singularities in soap-bubble-like and soap-film-like minimal surfaces*, Ann. of Math. 103 (1976) 489–539.
* J. J. van Wijk, A. M. Cohen, *Visualization of the genus of knots*, IEEE Visualization 2005, 567–574 (doi 10.1109/VISUAL.2005.1532843); *Visualization of Seifert surfaces*, IEEE TVCG 12 (2006) 485–496 (doi 10.1109/TVCG.2006.83). SeifertView: https://vanwijk.win.tue.nl/seifertview/
* M. van Garderen, J. J. van Wijk, *Seifert surfaces with minimal genus*, Bridges 2013, 453–456.
* K. A. Brakke, *Soap films on knots*, https://kenbrakke.com/knots/ ; *Soap films and covering spaces*, J. Geom. Anal. 5 (1995) 445–514.
* H. R. Parks, *Soap-film-like minimal surfaces spanning knots*, J. Geom. Anal. 2 (1992) 267–290.
* C. Coletti, *Volume-enclosing minimal surfaces of torus knots and links*, Bridges 2024, 463–466.
* A. Stockrahm, V. Lahtinen, J. J. J. Kangas, P. R. Kotiuga, *Cuts for 3-D magnetic scalar potentials: visualizing unintuitive surfaces arising from trivial knots*, Comput. Math. Appl. 78 (2019) 3200–3210.
* S. Wang, A. Chern, *Computing minimal surfaces with differential forms*, ACM Trans. Graph. 40(4) (2021) 113.
* T. Ekholm, B. White, D. Wienholtz, *Embeddedness of minimal surfaces with total boundary curvature at most 4π*, Ann. of Math. 155 (2002) 209–234.
* KnotInfo (C. Livingston, A. H. Moore), https://knotinfo.math.indiana.edu; LinkInfo; SnapPy (M. Culler, N. Dunfield, M. Goerner, J. Weeks).
