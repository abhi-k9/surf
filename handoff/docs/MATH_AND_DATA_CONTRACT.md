# Mathematical and data contract for the NPMS fixed-fan prototype

## 1. What is certified and what is not

The packaged full geometric bases for two rational planar fan domains have 100 (valence 5) and 172 (valence 9) basis functions and are exact-rationally constructed from an authentic 42-dimensional *local* Davydov–Yeo C² quintic Powell–Sabin-12 element, with global space dimension `18*n + 10`. The level-zero geometry basis is **nonnegative, partitioned, affine precise and P5 complete** on the specific domains. It is NOT a full positive geometric basis family for arbitrary meshes. Its sampled column-scaled coefficient condition numbers are approximately 16613 / 23846, so preserve diagnostics and avoid inverse-fitting or unqualified CAD-editing claims.

The 18/30 coarse-to-100/172 fine relation in this fixture is a **fixed-mesh enrichment**. It is not a D–Y uniform red mesh refinement. No verified level-independent lower stability bound, multiple extraordinary vertices, arbitrary topology or regular-box-spline transition is supplied.

The 42-dimensional local element is the analytic source; there is no reason to make production code invert 42x42 symbolic matrices at runtime. The values in `.bin` are the certified rational coefficients rounded to IEEE binary64 by the exporter.

## 2. Fan and PS12 geometry

The fan center is O=(0,0). Outer polygon vertex lists are in `fixtures/fan{5,9}_schema.json` and `macro_vertices_xy` of the fixture.

Macro sector i is the *ordered* triangle `T_i=[O,P_i,P_{i+1}]`, wrapping the last outer index. The reference element uses

```
v1 = (0, 3)  -> O
v2 = (0, 0)  -> P_i
v3 = (3, 0)  -> P_(i+1)
```

This is not an interchangeable permutation. Use supplied `microtri_vertices_xy[macro,cell,3,2]` to locate a point and compute barycentric coordinates in the microtriangle; do not recompute its geometry by assuming that cell 0 starts at the center. Each macro contains **12 PS12 microtriangles**. Each carries a quintic polynomial with **21 degree-five Bernstein ordinates** per global basis function.

With barycentric coordinates `(lambda0,lambda1,lambda2)` relative to the three vertices of that **microtriangle**, degree-d triangular Bernstein monomial is

`B_{ijk}^d(lambda) = d!/(i! j! k!) * lambda0^i * lambda1^j * lambda2^k`, with `i+j+k=d`.

The supplied Bernstein index order for d=5 is `for i=0..5: for j=0..(5-i): k=5-i-j`; see JSON `bernstein_order_ijk`. Never transpose or sort this order. For d=4,3 derivative temporary coefficient tables use the same nested-loop rule at that degree.

## 3. Runtime coefficient arrays

The **fine** basis binary data layout is `[function,macro,cell,bernstein_index]` (`float64`) containing every function's coefficient over every microcell of the two original fans. Functions with compact support have all-zero coefficient blocks outside that support; do not infer support from macro identity alone.

The binary header is `8s8I` in **little endian**: `"NPMSFAN1"` then `[version=1, degree=5, valence, macro_count, microcells_per_macro=12, bernstein_coefficients_per_cell=21, fine_function_count, coarse_parent_count]` as 8 uint32. Follow exactly `fixtures/fan{n}_binary_spec.json` for subsequent arrays/ordering and SHA-256. The code should reject unsupported version, truncated file, oversized lengths, NaNs, infinite values and corrupted counts *before* allocating large arrays. Do not native-cast a nonaligned header pointer. Treat arbitrary scene files as untrusted.

The `coarse_coefficients` in the same binary represent original valence-five 18-function / valence-nine 30-function carrier-free affine baseline. For each original coarse parent p and all coefficients `c`:

`coarse[p,c] = sum_{i: fine_parent_index[i]=p} fine[i,c]` (up to floating export rounding).

Hence homogeneous controls copy as `Q_fine[i] = Q_coarse[fine_parent_index[i]]`, preserving *any* positive-weight rational surface pointwise **on this unchanged PS12 triangulation**. The children inherit their parent's parameter-site XY coordinates in `fine_sites_uv`; these sites are heavily duplicated and sometimes outside the modeled polygon. There is no universal geometric meaning as intuitive cage vertices.

## 4. Basis derivatives in world parameter coordinates

Choose `lambda1=u, lambda2=v, lambda0=1-u-v`; for a microtriangle with parameter vertices `(t0,t1,t2)`, define `J=[t1-t0, t2-t0]` (columns), so `[x,y]^T=t0 + J [u,v]^T`.

For degree-five Bernstein coefficients `c_{ijk}`, the local barycentric-coordinate derivatives are exact finite differences:

`f_u = 5 sum_{i+j+k=4} (c_{i,j+1,k} - c_{i+1,j,k}) B_{ijk}^4`,

`f_v = 5 sum_{i+j+k=4} (c_{i,j,k+1} - c_{i+1,j,k}) B_{ijk}^4`,

`f_uu = 20 sum_{i+j+k=3} (c_{i,j+2,k} - 2 c_{i+1,j+1,k} + c_{i+2,j,k}) B_{ijk}^3`,

`f_uv = 20 sum_{i+j+k=3} (c_{i,j+1,k+1} - c_{i+1,j+1,k} - c_{i+1,j,k+1} + c_{i+2,j,k}) B_{ijk}^3`,

`f_vv = 20 sum_{i+j+k=3} (c_{i,j,k+2} - 2 c_{i+1,j,k+1} + c_{i+2,j,k}) B_{ijk}^3`.

Convert derivatives to physical parameter coordinates `(x,y)` by

`grad_xy = J^{-T} grad_uv`, `H_xy = J^{-T} H_uv J^{-1}`.

***Do not*** transform only first derivatives and keep second derivatives in cell local coordinates, nor compute physical seam C² by comparing untransformed derivatives across differently oriented cells. The affine mapping makes the second-derivative transformation above exact.

The independent Python implementation in `tools/reference_oracle.py` computes these by a separate source path from the archived mathematical construction.

## 5. Rational surface and derivatives

With 3D controls `P_i`, strictly positive weights `w_i`, positive basis `f_i(x,y)`:

`W = sum f_i w_i`; `A = sum f_i w_i P_i`; `S=A/W`.

Denote derivatives of A and W by subscripts `x,y,xx,xy,yy`. Use:

```
S_x  = (A_x - S W_x)/W
S_y  = (A_y - S W_y)/W
S_xx = (A_xx - S W_xx - 2 S_x W_x)/W
S_xy = (A_xy - S W_xy - S_x W_y - S_y W_x)/W
S_yy = (A_yy - S W_yy - 2 S_y W_y)/W
```

For regular surfaces, compute `n = cross(S_x,S_y) / norm(cross(S_x,S_y))`; if this normal length (relative to problem scale) is too small, return an explicit degenerate state. Fundamental form coefficients: `E=dot(S_x,S_x), F=dot(S_x,S_y), G=dot(S_y,S_y); e=dot(n,S_xx), f=dot(n,S_xy), g=dot(n,S_yy)`. Then `Gaussian K=(e*g-f*f)/(E*G-F*F)` and signed mean `H=(e*G-2*f*F+g*E)/(2*(E*G-F*F))`, with sign convention dependent on normal orientation. Principal curvatures are `H +/- sqrt(max(0,H*H-K))` only after checking metric determinant and numerical validity; near-degenerate or inconsistent values should be flagged, not silently corrected.

Basis `f_i` itself has partition `sum f_i=1`. Rational blending/influence for a selected control is `r_i = w_i f_i/W`; this also sums to 1 for positive weights. Under nonuniform weights, *the original stored control-site affine identity does not generally hold*.

## 6. Flat test and independently curved test

For the flat **unit-weight** preset set `P_i=(site_x,site_y,0), w_i=1`. The surface is exactly the parameter plane in the original rational construction:

`S=(x,y,0)`, `S_x=(1,0,0)`, `S_y=(0,1,0)`, `S_xx=S_xy=S_yy=(0,0,0)`.

This validates linear precision and all analytic derivative transforms. The curved/weighted preset assigns `w_i = 0.75+0.05*(i mod 9)` and `z_i=0.17*sin(0.31*i)+0.03*(parent_i mod 5)` with control X,Y retaining original sites. This is a **test field**, not a fair comparison against NURBS control cages; its 3D point and first/second derivative results are provided in `fixtures/*_golden_samples.json`.

Both fixtures include every PS12 internal seam; floating-point checks of 238 shared microedges with three points per edge attained maximum coefficient-basis derivative mismatch below `9e-15` across all values/gradients/Hessians. These are concrete sample residuals, not analytic continuity proofs; exact continuity is grounded in the research archive.

## 7. Numerical caveats

- The full 100/172 basis is strongly correlated. Avoid generic inverse solves without conditioning diagnostics and regularization; do not advertise robust control interpolation or independent editing quality as verified.
- On neighboring microcells, floating parameter point classification can choose either of two valid cells. Values and world first/second derivatives should agree to floating precision. Make classification tolerant of shared edges and deterministic for repeatability.
- The base function coefficients are rationally nonnegative in the certified research; after floating conversion tiny negative evaluation can be roundoff. Do not clip negative *meaningful* values, since it can break continuity and partition identities.
- Positive geometric basis alone does not guarantee immersion or absence of self-intersections for arbitrary 3D control positions.
- Signed analytic Hermite jet modes exist in the full D–Y space; only the specific geometric extraction furnished here is asserted nonnegative.
- Three control functions per original vertex in v1.7, many split children in the complete basis. These are not Greville-point cardinal interpolation functions; duplicated parameter sites are deliberate.
