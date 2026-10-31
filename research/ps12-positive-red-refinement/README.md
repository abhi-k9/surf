# NPMS: exact positive PS12 red subdivision (research milestone)

**Status:** Verified on the authentic quintic Powell–Sabin-12 reference split and numerically on the repository's existing fan5/fan9 fixtures. This is **stable refinement of a constrained Bézier coefficient representation**, NOT stable positive independent free-control refinement.

## Mathematical result

For one macrotriangle, canonical PS12 gives 12 microtriangles with 21 degree-five Bernstein coefficients each (252 coefficients). Uniform midpoint/red subdivision into four macrotriangles followed by PS12 subdivision produces 48 microtriangles and 1008 coefficients. Each fine microtriangle lies entirely within one coarse microtriangle (exact rational containment; affine invariant).

For each fine cell, let \(M_{kr}=\lambda_k(v_r^{fine})\), where \(\lambda_k\) are barycentrics relative to its containing coarse PS12 cell. Then \(M\ge0\), its columns sum to one, and the exact quintic restriction coefficients are

\[
R_{JI}=\frac{\binom{5}{I}}{\binom{5}{J}}
[\mu^J]\prod_{k=0}^2(\sum_{r=0}^2 M_{kr}\mu_r)^{I_k},\quad |I|=|J|=5.
\]

Concatenating all fine cells gives an exactly nonnegative, row-stochastic rational \(1008\times252\) matrix with **7518 strictly positive entries**. The canonical rational-row SHA-256 is `86b40ab48814183c45898f37a4cea0e7bbb4efe3df5b003c1c9d36b26b030092`.

For a fixed piecewise-quintic field and the concatenated (redundant, continuity-constrained) microcell Bézier coefficient arrays \(c^{(\ell)}\),

\[
\frac5{297}\lVert c^{(0)}\rVert_\infty\le
\lVert c^{(\ell)}\rVert_\infty\le\lVert c^{(0)}\rVert_\infty
\quad(\ell\ge0).
\]

The upper bound follows from the stochastic subdivision. The lower bound follows because the degree-five barycentric lattice collocation matrix \(A_{JI}=B_I^5(J/5)\) has exact inverse norm \(\lVert A^{-1}\rVert_\infty=297/5\). This is uniform in level and affine macrogeometry. It is **not** an inverse bound for a freely editable global spline basis.

For a positive-weight rational surface, apply \(R\) separately to the four homogeneous Bernstein coefficient fields. Each refined denominator coefficient is a positive convex combination of old denominator coefficients, and the resulting surface is *exactly unchanged*. Its fine Euclidean Bézier poles remain in convex hulls of the coarse poles. Global \(C^2\) continuity is inherited, because the underlying piecewise polynomial numerator and denominator are identical functions after subdivision; the original Davydov–Yeo nested-space theorem covers the additional supersmoothness rules.

## Evidence

- 48 fine microcells exactly contained in the 12 coarse cells; all 1008 exact row sums equal 1; 7518 nonnegative rational nonzero entries.
- Six independent source polynomials including four archived nonpolynomial 42D D–Y modes: 6048 exact Bernstein coefficient equalities and 576 exact pointwise equalities.
- All 42 D–Y source modes plus a quintic polynomial: **51,084 exact physical-coordinate C² seam-jet equalities** on the 66 interior fine-cell seams.
- Positive-weight rational numerator/denominator, including a nonpolynomial D–Y mode: exact preservation through three red levels (48, 192, 768 microtriangles).
- `surf` fan5 (100 functions): 504,000 refined coefficients, minimum 0, partition residual 2.00e−15, affine moment error ≤2.66e−15.
- `surf` fan9 (172 functions): 1,560,384 refined coefficients, minimum 0, partition residual 2.44e−15, affine moment error ≤2.39e−15.

All exact results use rational arithmetic. Fan fixture residuals are numerical and **do not** substitute for exact results. The original research package includes `PS12_Q5_RED_CSR_v1.json`, immutable fixture data, the independent verifiers and a SHA-256 manifest. It is attached in the project discussion.

## Safe implementation boundary

Do not expose 1008 fine coefficients as independent new spline handles: arbitrary changes violate the continuity constraints. Initial integration can optionally use the positive refinement for exact Bézier tessellation, geometry-preserving export, weight-bound propagation, and subdivision reference tests without changing any existing independent controls or fixture hashes.

The next **unsolved** step is to construct an independent, full-rank, local and positive fine D–Y basis \(E_f\) and a transfer \(T\ge0\) satisfying

\[
R E_c=E_fT,\qquad T\mathbf 1_{coarse}=\mathbf 1_{fine},
\]

then prove a uniform local-dual lower bound on overlapping hierarchical levels. Even a positive row-stochastic \(T\) proves only forward convex-hull/weight stability, not stable independent editing.

**Sources:** Davydov–Yeo (2013), https://oleg-davydov.de/s25refinable.pdf ; existing `PROVENANCE_AND_LIMITS.md` and certified `handoff/` research assets. Keep this milestone separate from production claims.
