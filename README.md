# SeifertSurface

A Seifert surface for a knot or link `K ⊂ ℝ³` is a compact oriented surface with boundary `K`. This project draws
one for any closed braid, in two steps. First it builds the classical one, Seifert's algorithm applied to the closed
braid: one disk per strand stacked along the braid axis, one half-twisted band per crossing (the Bennequin surface).
Then it takes the boundary of that surface as a fixed wire and relaxes the surface to the minimal surface spanning
the wire, by solving the minimal surface equation on it: the soap film on a knotted wire. Alongside, it reads the
invariants off the surface: the Seifert matrix, the Alexander polynomial, the signature, the genus.

It is a single-page web app in the shape of [EllipticCurve3D](../EllipticCurve3D): `web/index.html` with `js/`, `css/`
and `data/`, three.js and KaTeX vendored, no build step; the mathematics in two dependency-free modules that also run
under Node for the tests.

```bash
python3 -m http.server -d ~/Projects/SeifertSurface/web 8766      # then http://localhost:8766/
```

## Layout

| path | what it is |
|---|---|
| `web/index.html` | the page |
| `web/js/seifert.js` | braid words, the closed braid, the Bennequin surface as a triangle mesh, the Seifert matrix, `Δ(t)`, signature |
| `web/js/minimal.js` | discrete minimal surfaces: cotangent Laplacian, conjugate gradients, the harmonic step and mean curvature flow, Delaunay flips, tangential smoothing, diagnostics, linking numbers |
| `web/js/app.js` | viewer and UI |
| `web/js/vendor/` | three.js r128, OrbitControls, KaTeX 0.16.11 (licenses alongside) |
| `web/css/` | the stylesheet (EllipticCurve3D's), KaTeX's with its fonts |
| `web/data/knots.json` | KnotInfo's braid words for the 2,961 knots through 12 crossings and LinkInfo's for the 4,188 links through 11, with the knots' genera |
| `web/data/braid.sage` | `sage web/data/braid.sage K12n242` (or a PD code, or `DT:[...]`): a braid word for any knot, by SnapPy |
| `web/test/test-seifert.mjs` | `node web/test/test-seifert.mjs`: 185 checks against Sage's values in `vectors.json` (made by `make_vectors.sage`), plus the mesh and the solver |

The web app takes a braid word (`1 2 -1 2`, `s1 s2 s1^-1 s2`, `σ₁σ₂σ₁⁻¹σ₂`, `abAb`, `(1 2)^5`), a torus link `T(p, q)`,
or a name from the tables (`3_1`, `12n242`, `11n_34`, `L6a4`, `L2a1{1}`, and the aliases `trefoil`, `figure-eight`,
`hopf`, `whitehead`, `borromean`, `lehmer`). The info line gives the braid, the strands, crossings and components,
`χ` and the genus of the surface, `Δ(t)`, `det`, `σ`, and whether the surface has the knot's genus (from KnotInfo, or
from `deg Δ = 2g`). **Relax** runs the relaxation; the flow bar shows the round, the area, the root mean square of the
discrete mean curvature `|H|` over the interior (the residual of the equation), and the largest displacement. The
drawer has three tabs: **Surface** (the starting surface: disk spacing, band width and bulge, resolution, the rounding
of the wire's corners), **Flow** (the harmonic step or mean curvature flow with a finite step, edge flips, tangential
smoothing, the stopping rule), **Display** (two-sided coloring, or by disk and band, or by mean curvature; opacity,
wireframe, the wire and its thickness, the starting surface as a ghost, axes, PNG, a shareable link `#word`).

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
radius `0.1 R`, then resampled at equal arclength, so the corners do not bunch the points), and from then on it is
the wire: it never moves again.

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
genus of knots* (IEEE Visualization 2005) and *Visualization of Seifert surfaces* (IEEE Transactions on Visualization
and Computer Graphics 12(4), 2006), implemented in SeifertView, a Windows program with a braid-word input and a very
good renderer. Their construction is the same Seifert-circles-as-disks, crossings-as-twisted-bands one, done for
general knot diagrams (nested circles at different heights, bands routed between them) and then made smooth and
pleasant to look at, with the surface's two sides in two colors to show its orientability, and the genus counted
from the diagram. The surface is a topological object, and theirs is a good picture of it. What can be done better is
to give it a canonical geometry once the wire is chosen: the soap film. That geometry is not a matter of taste, it is
the solution of an elliptic boundary value problem, and it makes the disks-and-bands scaffolding disappear into a
smooth surface that a physical wire would actually carry. Their pictures remain the model for the rendering.

## Lehmer's polynomial

The (−2, 3, 7)-pretzel knot is `12n_242` in the tables, the closure of the positive 3-braid `σ₁σ₂²σ₁²σ₂⁷` (KnotInfo's
word `1 2 2 1 1 2 2 2 2 2 2 2`), genus 5, signature −8, fibered, an L-space knot. Its Alexander polynomial is

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

## Tests

```bash
node web/test/test-seifert.mjs          # 185 checks
sage web/test/make_vectors.sage         # regenerates web/data/knots.json and web/test/vectors.json (needs the KnotInfo database in Sage)
```

The checks: parsing of every input form; strands, crossings, components and genus of closed braids; the Seifert
matrix equal to Sage's on 261 connected braids (17 named, 160 random on 2–5 strands, the 85 knots through 9
crossings), the Alexander polynomial and the signature with them, `deg Δ ≤ 2g(K) ≤ 2g(surface)` against KnotInfo's
genera; `Δ` of `12n_242` equal to `L(−t)`; the meshes oriented, with `χ = n − c`, the right number of boundary loops,
no degenerate triangles; the linking numbers of the Hopf link and `T(2,4)` from the built geometry (`+1`, `−1` for the
mirror, `2`); the solver on a bumped disk (flattens to the polygon's area, `H → 0`, boundary fixed), on the catenoid
(area within 0.5% of the exact catenoid, waist radius, monotone area under both the harmonic step and mean curvature
flow), and on the trefoil's surface (area decreases, `χ` and orientation survive the flips, the wire does not move,
the residual falls).

## Planned

* **A surface from a Seifert matrix.** Given `V` (from `Friedl's Algorithm.ipynb`, or any `V` with `det(V − Vᵀ) = ±1`),
  a disk with `2g` bands whose twists and mutual linkings realise `V`, then the same relaxation: minimal surfaces for
  Lehmer's polynomial itself, `L(t)`, and for any Alexander polynomial one likes.
* **Wires that are not the Bennequin boundary.** A torus knot on a round torus, a Fourier knot, a wire drawn by hand,
  with the Seifert surface built inside it (an isotopy from the Bennequin wire, or Seifert's algorithm on the
  projection of the given wire).
* **Diagrams as input.** PD and DT codes through `braid.sage` (SnapPy's Vogel algorithm), as `SkeinA/data/braids.txt`
  already does for the census.
* **A self-intersection check** during the flow, and isotropic remeshing (edge split and collapse) beside the flips.
* **Export** of the relaxed mesh (OBJ/STL) for printing, and GitHub Pages as for EllipticCurve3D.

## References

* [AT1977] F. J. Almgren, W. P. Thurston, *Examples of unknotted curves which bound only surfaces of high genus within their convex hulls*, Ann. of Math. 105 (1977) 527–538.
* [Ben1983] D. Bennequin, *Entrelacements et équations de Pfaff*, Astérisque 107–108 (1983) 87–161.
* [Bra1992] K. A. Brakke, *The Surface Evolver*, Experiment. Math. 1 (1992) 141–165.
* [Col2013] J. Collins, *An algorithm for computing the Seifert matrix of a link from a braid representation* (the algorithm Sage implements; cited there as [Col2013]).
* [HS1979] R. Hardt, L. Simon, *Boundary regularity and embedded solutions for the oriented Plateau problem*, Ann. of Math. 110 (1979) 439–486.
* [Hir2001] E. Hironaka, *The Lehmer polynomial and pretzel links*, Canad. Math. Bull. 44 (2001) 440–451.
* [Jos1986] J. Jost, *Existence results for embedded minimal surfaces of controlled topological type*, I–III, Ann. Scuola Norm. Sup. Pisa 13 (1986) 15–50, 401–426; 14 (1987) 165–167.
* [PP1993] U. Pinkall, K. Polthier, *Computing discrete minimal surfaces and their conjugates*, Experiment. Math. 2 (1993) 15–36.
* [Tay1976] J. E. Taylor, *The structure of singularities in soap-bubble-like and soap-film-like minimal surfaces*, Ann. of Math. 103 (1976) 489–539.
* J. J. van Wijk, A. M. Cohen, *Visualization of the genus of knots*, IEEE Visualization 2005, 567–574; *Visualization of Seifert surfaces*, IEEE TVCG 12 (2006) 485–496. SeifertView: https://www.win.tue.nl/~vanwijk/seifertview/
* KnotInfo (C. Livingston, A. H. Moore), https://knotinfo.math.indiana.edu; LinkInfo; SnapPy (M. Culler, N. Dunfield, M. Goerner, J. Weeks).
