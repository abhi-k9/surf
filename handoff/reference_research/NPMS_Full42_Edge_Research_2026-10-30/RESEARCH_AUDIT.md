# NPMS — Full 42D Davydov–Yeo Architecture and Positive Edge Integration

**Research date:** 30 October 2026  
**Status:** New exact fixed-mesh positive geometric enrichment on two rational planar fans; authentic full 42D analytic reference; proof-level local positivity and transfer; **no claim of stable spatial refinement on arbitrary meshes**.

## Executive finding

The earlier NPMS route through the 30D PS6/PS12 intersection loses compact edge-star freedom. Rather than attempt another basis in that restricted intersection, this study returns to the authentic 42D Davydov–Yeo (D–Y) Powell–Sabin-12 macro-element and reuses two independent certified data sets:

- **Two-macro full-space edge star (v1.2):** exact dimension 8. Its 7D zero-first-normal subspace admits seven independent, Bernstein-nonnegative functions; the eighth direction must change sign. This is the precise reason to keep a signed analytic layer alongside positive geometric controls.
- **Carrier-free affine fan baseline (v1.7):** 3 positive quintic functions per vertex, exact partition and affine reproduction, center's natural one-ring support, outer control support at most four sectors. On the tested valence-5 and valence-9 fans, its previous lack of radial-edge curvature freedom was remediable.

**First new result (well-scaled restricted basis):** All **3n** previously certified positive radial-edge curvature directions from v1.3 can be inserted into the v1.7 basis **with a nonnegative, exact fixed-mesh control-transfer matrix**. The insertion uses ONLY associated outer-vertex parents (never a borrowed center-star reservoir). The new curvature children retain two-sector support; each outer parent retains at most four-sector support. The resulting families contain **33** and **57** independent positive quintic functions for valences **5** and **9**, respectively, while preserving affine precision and any rational surface geometry on insertion.

The resulting child strengths are **not** tiny: every one is at least **10123/25600** (valence 5), and every one equals **1/2** (valence 9). Measured column-normalized sampled basis condition numbers are approximately **101.34** and **74.34**. These are instance diagnostics, not level-uniform bounds.

**Second new result (complete but poorly conditioned):** An exact parent-dominance screen found a compatible positive local reservoir for **every missing global Hermite direction**: 82 modes on the valence-five fan and 142 on the valence-nine fan. A constructive finite splitting then produces complete **100- and 172-function positive normalized affine-precise quintic bases** spanning the entire level-zero D–Y spline spaces on these *fixed fans*. Hence the complete positive geometric spaces themselves reproduce all P5 polynomials, but their measured column-scaled coefficient condition numbers are approximately **16,613** and **23,846**. This is a feasibility benchmark, not an acceptable stable final representation.

**Critical boundary:** Both new constructions are **control-space enrichments on the same triangles**. Neither establishes exact positive transfer under true red spatial mesh subdivision, a level-independent coefficient dual, an arbitrary-mesh theorem, or a regular-box-spline chart connection.

## 1. Mathematical reference and notation

The exact D–Y reference macro-element consists of 12 Powell–Sabin microtriangles, each carrying a total-degree-at-most-five polynomial, with 21 Bernstein coefficients per cell. Hence 252 raw coefficients are constrained to a **42-dimensional C² local space** with the prescribed additional third-order conditions. The published nested D–Y spaces have stable local determining sets, which do **not** automatically imply a positive control basis or a positive transfer between independently designed bases.

Let `n = 5` or `9`, with the exact rational fan vertices in `base_v17/fan_assembly.py` and center `O=(0,0)`. Let `B_0,...,B_(N0-1)` be the previously certified positive v1.7 quintic functions, `N0 = 3(n+1)`. In Bernstein-coefficient form each is a column in `Q^(252n)`, and:

```
B >= 0                          coefficientwise
B * 1 = 1                      exact partition of unity
B * A = X                      exact affine-coordinate precision
rank(B) = 3(n+1)               follows from nodal jet independence
```

Here `A` records the fixed 2D parameter sites; `X` is the exact coefficient matrix for the two parameter coordinates. The expression `rank` is a function-space statement, not merely an SVD.

For radial edge `i` let `E_(i,k)`, `k=0,1,2`, denote the three previously certified positive v1.3 edge modes. Their 4D local D–Y edge-nodal coefficient matrices have **first-normal component zero** and an independent three-dimensional block in the remaining edge derivatives. Each `E_(i,k)` occupies exactly the two adjacent fan sectors. Exact Bernstein conversion gives `E_(i,k) >= 0` on all microtriangles. Physical-coordinate radial C² matching was independently verified, not inferred from reference-space normals.

## 2. Why the eight-dimensional edge star matters

On the original two-macro reference square, the compact D–Y edge-star space `W` is dimension 8. All functions in `W` vanish on the common interior edge `e`. The first-normal derivative functional `L(f) = D_n f(m)` at its midpoint has rank 1 on `W`.

For any globally nonnegative `C¹` function `f` with `f|_e = 0`, every interior edge point is a local minimum. Thus `D_n f|_e = 0`. Consequently the nonnegative cone of `W` lies in `ker L`, a seven-dimensional subspace. The archived exact Bernstein certificate supplies **seven independent nonnegative directions** spanning `ker L`. Hence:

```
W = ker(L) ⊕ span{signed jet mode},     dim ker(L) = 7.
```

A fully nonnegative eight-function compact edge-star basis of `W` is **impossible**. The eighth direction must remain signed, or its first-normal effect must be provided by positive controls whose support extends beyond this compact two-macro star. The v1.2 exact `verify_committed_certificates.py` was rerun during this study; it passes 3,528 rational Bernstein positivity tests, exact 42D compatibility, edge-star dimension and rank, and fixed-patch collocation stability.

**Scope of attribution:** The original 8D/7D calculation is a previously certified local reference, not a newly discovered fact in this study. The new contribution is its structural integration with a carrier-free affine fan geometry basis and a much stronger local positive parent split.

## 3. New exact parent assignment and positive split

For each radial edge `i`, consider only the **three positive geometric controls associated with its outer endpoint** `P_i`:

```
parent candidates = B_(3(i+1)), B_(3(i+1)+1), B_(3(i+1)+2).
```

For each child mode `E_(i,k)`, form exact coefficientwise dominance ratios

```
r_(a,k) = min{ B_parent_a[c] / E_(i,k)[c] : E_(i,k)[c] > 0 },
```

setting `r_(a,k)=0` if there is any positive child coefficient where the parent Bernstein coefficient is zero. Enumerate all six bijective assignments `sigma:{0,1,2} -> {0,1,2}` and select one maximizing `min_k r_(sigma(k),k)` among positive-feasible assignments. This is **finite exact rational combinatorics**, not a floating-point optimizer and not a bilinear unknown basis/refinement problem.

Choose the certified child scale

```
epsilon_(i,k) = min(1/2, r_(sigma(k),k)/2) > 0.
```

Each outer parent is used at most once, so its fine remainder is simply

```
B'_p = B_p - epsilon_(i,k) E_(i,k) >= 0;
F_(i,k) = epsilon_(i,k) E_(i,k) >= 0;
B_p = B'_p + F_(i,k).
```

The exact Bernstein inequalities hold on **every** one of the `252n` coefficient positions for **every** new fine function. The controls become `N1 = N0 + 3n = 6n + 3`: 33 for `n=5`, 57 for `n=9`.

### Exact geometry preservation

Let `T ∈ {0,1}^{N1 × N0}` have one nonzero 1 in every row: one for every unchanged/remainder parent and one mapping each new curvature child to its associated coarse parent. By construction:

```
B_coarse = B_fine T,       T >= 0,       T*1 = 1.
```

For **arbitrary homogeneous control data** `Q_coarse = (w_i P_i, w_i)` with `w_i>0`, initialize `Q_fine = T Q_coarse`. Then

```
B_fine Q_fine = B_coarse Q_coarse
```

as an exact polynomial identity. Positive weights are retained since each new weight is copied from a parent. Therefore the rationally projected surface is unchanged at insertion. Repeating each parent's **parameter site** for its new child ensures

```
B_fine (T A_coarse) = B_coarse A_coarse = (x,y),
```

so affine-coordinate precision and normalization are preserved.

### Compact support and rank

- Each new curvature control is supported on two incident fan sectors.
- Every outer parent or remainder is supported on at most four sectors because the child support is inside the selected parent support; subtraction cannot create support outside the original parent.
- Center controls retain their natural full center vertex star, but no new center-star reservoir is borrowed.
- The original 3(n+1) geometric functions are linearly independent by their vertex value/first-jet data and global first-normal edge data.
- The 3n curvature children are independent in the remaining **three** radial-edge nodal directions and have zero first-normal edge components. The geometric base has no corresponding edge second/third-order nodal degrees of freedom. Thus the new 3n directions are independent modulo the old geometric span, giving exact rank 33/57. Numerical coefficient and collocation SVDs independently show full rank.

**Important limitation:** The geometric `B_fine` span is not the full D–Y global space; full 42D refers to the local **underlying analytic spline space**, and the positive family is a controlled restricted subspace inside it.

## 4. Reproduced exact and numerical results

| Quantity | Valence 5 | Valence 9 |
|---|---:|---:|
| Coarse positive controls | 18 | 30 |
| Fine positive controls | 33 | 57 |
| New positive edge controls | 15 | 27 |
| Minimum edge scale | `10123/25600` | `1/2` |
| Maximum edge scale | `1/2` | `1/2` |
| Negative Bernstein coefficients | 0 | 0 |
| Exact partition and affine reproduction | Pass | Pass |
| Exact nonnegative coarse/fine identity | Pass | Pass |
| Maximum outer support | 4 sectors | 4 sectors |
| Edge child support | 2 sectors | 2 sectors |
| Numerical column-scaled Bernstein coefficient condition | 59.1164 | 44.0292 |
| Numerical column-scaled sampled function condition | 101.3420 | 74.3359 |

The condition numbers refer to a fixed coefficient scaling and four interior barycentric sample sites per microtriangle. They are **not mesh-shape-independent** and must never be called uniform multilevel bounds. The v1.7 geometric span still has zero second-order jets at original fan vertices, and therefore cannot represent general quadratics.

The existing v1.3 physical-seam verifier was independently rerun: it passed **1,200** (valence 5) and **2,160** (valence 9) exact scalar equalities for values, first-normal and second-normal derivative traces. The new functions are linear combinations of these compatible source functions and the verified v1.7 base, so they inherit the same seam compatibility.

## 5. Exact P5 reproduction in the analytic layer

The signed analytic layer is the authentic full 42-dimensional reference D–Y basis, evaluated through the exact 42×42 physical Hermite nodal matrix. `verify_p5_analytic.py` independently reconstructs **all 21 monomials** `x^a y^b` with `a+b <=5`:

1. Compute each monomial's exact physical third-order vertex jets and four prescribed edge-midpoint nodal values.
2. Apply the rational inverse 42×42 nodal matrix and the complete PS12 coefficient map.
3. Independently convert the target monomial to Bernstein coefficients on each of the twelve microtriangles.
4. Compare all 252 Bernstein coefficients exactly for every monomial.

All 21 tests pass. This demonstrates local polynomial completeness of the **underlying analytic space**. The restricted **33/57** positive geometric basis is *not* polynomial complete. The subsequently constructed **100/172** positive geometric basis **is** full-dimensional and therefore contains P5, but is poorly conditioned. A well-conditioned complete positive basis with a certified stable spatial-refinement operator remains the principal missing theorem.

## 6. A restricted but genuine scale-independent edge-star block bound

The archived v1.2 seven-function positive edge-star family comes with an exact seven-site rational collocation matrix `M7`. With each function Bernstein-scaled to maximum coefficient 1, the independently rerun certificate gives

```
||Σ_j c_j ψ_j||_(L∞(Ω_ref)) >= b ||c||∞,
b = 2154636280112944149852832120012435275 /
    57028248035493506784145691689562656768
  ≈ 0.0377819125492323.
```

For **any invertible affine copy** `F` of the reference patch, set `ψ_j^F(x)=ψ_j(F^-1(x))`. The same seven-point collocation matrix is obtained at transformed sample sites `F(x_i)`, so `b` is unchanged by patch size or uniform refinement level. For a family of **pairwise interior-disjoint** affine-image edge patches (with zero boundary traces), the global block basis inherits

```
b ||c||∞ <= ||Σ_(patch,j) c_(patch,j) ψ_(patch,j)||∞ <= 7 ||c||∞.
```

This is a mathematically valid **uniform isolated-block** coefficient estimate for disjoint supports, assuming affine mapped geometry and zero extension. **It is not a theorem** for overlapping patches, a complete global basis, a mixed-level hierarchical family, or a nonnegative spatial refinement transfer. Those require separate analysis of overlap, local duals, truncation and admissibility.

## 7. What the new experiment did NOT establish

1. **True spatial subdivision:** No `B_coarse = B_fine T` was produced between geometrically red-refined PS12 macro meshes. We certified exact fixed-mesh function enrichment only.
2. **Uniform refinement stability:** The new scales are substantially stronger than old `1/32` descendants, but a positive local split at one level is not a level-independent lower coefficient bound. Overlapping multilevel supports can still be unstable.
3. **Well-conditioned P5-complete positive geometry:** The smaller 33/57 family does not contain general quadratic fields. The new complete 100/172 positive geometry family does contain P5 but has high one-level coefficient conditioning; it is not a production-quality solution.
4. **NURBS-equivalent regular grid and extraordinary-patch junction:** No exact C² nonnegative transfer to the regular triangular-lattice quintic box-spline basis was constructed.
5. **Arbitrary manifold atlas:** Only two concrete planar fan geometries, not adjacent extraordinary vertices, curved chart transitions, corners, holes or global cycle compatibility.
6. **Curvature fairness:** The newly introduced edge controls have genuine transverse second-order freedom, but no theorem comparing impulse/checkerboard/principal-curvature response against cubic/quintic NURBS.
7. **Universal guarantees for arbitrary geometry controls:** Positivity gives rational convex-hull containment, not an immersion or self-intersection-free surface for every choice of 3D controls.

## 8. Research decision / falsifiable next experiment

**Keep the full D–Y 42D analytic space.** The combined positive vertex/edge family is now a concrete local geometry basis candidate with exact affine precision and useful edge curvature controls. The natural next proof target is **spatial refinement**, not another adjustment of fixed-mesh parents:

1. Select **one common two-macro edge patch**, its uniformly red-refined descendants, and a complete rational coefficient conversion on the overlay.
2. Construct the new fine positive extraction **before** solving transfer (avoid treating the extraction and transfer simultaneously as an LP). Include geometric vertex/tangent, second-order vertex and positive edge modes.
3. Solve exact inclusion `D B_c = B_f T`, impose `T >=0` and a strict support/rank/admissibility mask. If infeasible, produce a rational Farkas dual witness for that fixed choice of fine basis.
4. Verify physical C² seams and D–Y extra per-side third-order requirements on **every refined macro interface**. Avoid equating third derivatives across macro boundaries unless the published space requires it.
5. Establish patch-scaled local dual estimates under shape regularity, and only then admissible overlapping hierarchical stability.
6. Use the new full 100/172 positive basis as an existence target, but redesign its many-to-one parent splits to achieve a stable local dual, meaningful control-site distribution, and robust second/third-order editing. The earlier 51/87 Hessian-and-radial-edge variant is a better-conditioned intermediate control-space experiment but does not reproduce P2.
7. Compare impulse, wave and checkerboard response against NURBS on the same parameter sites and controls, and study extraordinary-to-regular box-spline transitions separately.

**Stop/go condition:** If no positive transfer can be certified for an admissible fixed candidate fine basis, stop iterative coefficient tweaking and investigate a modified nested spline/dual structure or a mixed-degree extraordinary transition.


## 9. New full-space positive existence construction: 100/172 functions

### Complete Hermite-coordinate inventory

For an `n`-sector fan the globally assembled level-zero D–Y space has dimension `18n+10`, equivalently `10(n+1)` vertex jets (all derivatives of total order at most three) and `4(2n)` macro-edge midpoint jets. The existing v1.7 positive base has dimension `3(n+1)` and an invertible projection onto all vertex value/first-derivative data. The remaining nodal coordinates are:

- Seven second/third-order vertex jets per vertex: `7(n+1)` directions.
- All four edge-nodal directions on each of the `n` radial and `n` outer edges: `8n` directions.

There are `15n+7` such missing coordinates: **82** at valence 5 and **142** at valence 9. Together with the v1.7 base they form a linearly independent complete coordinate system for the published `18n+10`-dimensional global D–Y space, because the vertex value/first-derivative projection of the base is invertible.

### Exact local Bernstein domination

For each missing signed global nodal dual `H_k`, test the three positive v1.7 parents associated with the same center/outer vertex (center modes: center group; radial/outer edges: associated outer group). The exact-rational calculation in `full_space_positive.py` finds one parent `P` per mode such that every nonzero Bernstein coefficient of `H_k` lies at an index where the corresponding Bernstein coefficient of `P` is strictly positive. Hence

`r_k = min_{H_k[c] != 0} P[c]/|H_k[c]| > 0`.

Every missing mode passed this support and dominance test; no center reservoir was borrowed for an outer-associated mode. The smallest such ratios on the archived test fans are `3209/24000` (valence 5) and `1123/3600` (valence 9).

Group the missing modes by their assigned parent. Suppose parent `P` receives `m` signed modes `H_1,...,H_m`, with ratios `r_k`. Set `beta_k = r_k/(2m)`. Define

```
C_0 = (P - sum_k beta_k H_k) / (m+1)
C_k = (P + beta_k H_k) / (m+1),   k=1,...,m.
```

Since `|H_k[c]| <= P[c]/r_k` for all Bernstein indices `c`, we have

```
sum_k |beta_k H_k[c]| <= P[c]/2,
C_0[c] >= P[c]/(2(m+1)) >= 0,
C_k[c] >= (1-1/(2m))*P[c]/(m+1) >= 0.
```

Every new control is therefore Bernstein nonnegative **exactly**, with support contained in that of `P`. By direct cancellation `C_0+sum_k C_k=P`. The exact nonnegative transfer copies the old homogeneous parent data to every child. Each child inherits the same control-site coordinate, maintaining affine precision and rational surface geometry at insertion.

### Exact full rank and P5 completeness

The child span contains `P` by summation, and for each assigned mode

`H_k = ((m+1)C_k - P)/beta_k`.

Every `beta_k` is positive. Thus all missing Hermite nodal modes are recovered in the fine span. Combined with the invertible vertex-value/first-jet projection of the old basis, the new positive family spans the complete D–Y level-zero space. Its dimension and control count both equal `18n+10` (100 and 172). Since every global polynomial in `P5` belongs to the D–Y space, the full positive geometric family **exactly reproduces P5** using appropriate, generally signed coefficient vectors. Basis positivity does **not** require polynomial reproduction coefficients to be nonnegative.

### Recorded cost

| Full basis metric | Valence 5 | Valence 9 |
|---|---:|---:|
| Exact positive basis count | **100** | **172** |
| Rational Bernstein coefficients checked | **126,000** | **390,096** |
| Assigned missing jets | 82 | 142 |
| Maximum signed modes from one positive parent | 14 | 15 |
| Smallest perturbation `beta_k` | `43/4224` | `1123/108000` |
| Maximum outer support | 4 sectors | 4 sectors |
| Column-scaled Bernstein matrix condition, numerical | **16,612.70** | **23,846.02** |
| Partition, affine precision and positive fixed-mesh prolongation | Exact | Exact |
| True spatial refinement and level-uniform stability | Not established | Not established |

These condition numbers show that **existence alone does not solve numerical modeling quality**. The many-to-one parent splitting creates numerous related, nearly dependent functions and duplicated initial parameter sites. A user-facing geometric control net, high-quality curvature response, shape-regularity stability and spatial refinement remain absent. This is a certified finite fan **mathematical feasibility result**, not a generalized refinement algorithm.

### Intermediate Hessian experiment

`full_geometric_enrichment.py` also certifies a less aggressive positive split: every original vertex receives three independent second-order jet directions while the three radial-edge curvature modes remain positive and independent. This gives 51 functions at valence 5 and 87 at valence 9, with no support expansion. Nevertheless, their geometric spans still do not contain the quadratic monomials `x²`, `xy`, `y²` in the tested coefficient coordinates (reported least-squares relative residuals about `0.03–0.05` on fan 5 and `0.02–0.04` on fan 9). The missing mechanism is independent first-normal edge moments, originally coupled to vertex tangents in v1.7. This is one reason the **complete** Hermite inventory, rather than vertex Hessians alone, is required.

**Conclusion:** The crucial outstanding research theorem has narrowed from *whether a positive, local, P5-complete basis exists on a fixed extraordinary fan* to *whether such a basis can be **uniformly stable, appropriately scaled, and positively nested under actual spatial refinement**, with a practical control-site layout and compatible regular box-spline junctions*. The present exact instances answer only the former.

## Source papers and provenance

- O. Davydov and W. P. Yeo (2013), *Refinable C² Piecewise Quintic Polynomials on Powell–Sabin-12 Triangulations*, Journal of Computational and Applied Mathematics **240**, 62–73. Full manuscript: https://oleg-davydov.de/s25refinable.pdf . Source of genuine nested 42D local PS12 spaces and stable local MDS.
- H. Speleers (2010), *A normalized basis for quintic Powell–Sabin splines*, Computer Aided Geometric Design **27**, 438–457, https://doi.org/10.1016/j.cagd.2010.05.001 . Independent positive normalized basis; no transfer compatibility inferred.
- M. Knez and M. Šteblaj (2026/2027), *New approach for constructing edge B-spline-like basis functions for C¹ and C² splines over triangulations*, J. Comput. Appl. Math. **489**, 117914, https://doi.org/10.1016/j.cam.2026.117914 . Relevant geometry-based edge-basis approach, conditional on its own convexity assumptions.
- NPMS v1.2 exact two-macro edge-star certificates (`v12_reference/`), NPMS v1.3 positive edge extraction (`positive_physical_coefficients.json`, `v13_seam_reference/`), and NPMS v1.7 affine fan reference (`base_v17/`). These are project-specific restricted constructions, not published universal spline theorems.

## Reproduction

From the extracted directory, with `sympy`, `numpy`, and optionally `scipy` for old notebooks:

```sh
OPENBLAS_NUM_THREADS=1 python edge_enrichment.py
OPENBLAS_NUM_THREADS=1 python verify_p5_analytic.py
OPENBLAS_NUM_THREADS=1 python full_geometric_enrichment.py
OPENBLAS_NUM_THREADS=1 python full_space_positive.py
python verify_signed_edge_decomposition.py
python verify_scale_block.py
python verify_full_metadata.py  # structural quick-check; NOT equivalent to the exact coefficient rebuild
(cd v12_reference && OPENBLAS_NUM_THREADS=1 python verify_committed_certificates.py)
(cd v13_seam_reference && OPENBLAS_NUM_THREADS=1 python seam_world_exact.py)
sha256sum -c SHA256SUMS.txt
```

**Verification types:** exact rational/finite-dimensional identities and signs for construction; floating-point singular values are diagnostics; the isolated-block theorem uses a fixed-patch exact inverse and an affine-change-of-variables argument. No unsupported arbitrary-mesh claim is made.
