# NPMS research agent — next proof obligation after positive edge integration

**Keep the current results unchanged.** New complete **100/172-function positive normalized P5-complete** bases now exist on two fixed fans, though they are ill-conditioned (~1.7×10⁴ and ~2.4×10⁴ column-scaled). The next work is to recover **stable local duals and genuine spatial refinement**, not repeat one-off positivity experiments. The smaller 33/57 and 51/87 bases offer conditioning/control tradeoffs.

## Primary question

For one admissible two-macro *interior* edge patch and its uniformly red-refined PS12 mesh, can we choose a complete nonnegative fine vertex/edge/detail extraction such that the coarse positive geometric family has an exact sparse **nonnegative** transfer into it, with no zero fine functions and a stable local dual bound?

1. Rebuild the full 42D reference and 8D edge-star exact coefficients; independently check the genuine D–Y constraints at macro interior split edges and **per-side** third-order edge-midpoint compatibility. Do not impose cross-macro C3 where the published space requires only C2 matching.
2. Choose a common rational geometry for original and red-refined triangles. Compute exact triangular Bernstein subdivision `D`, not approximate point matching.
3. Use the complete 100/172 positive normalized P5 basis as a **feasibility witness, not the proposed fine basis**: redesign local extraction to avoid 14–15 almost-duplicate children sharing one parent and copied parameter sites. Require a stable local dual, nonnegative Bernstein coefficients, true P5 completeness and useful control response. Keep signed fitting functionals distinct from positive surface controls.
4. After the fine basis is **fixed**, solve `D*E_coarse = E_fine*T` for `T>=0`. If infeasible, produce a rational Farkas certificate. Avoid the bilinear trap of optimizing E_fine and T simultaneously with an LP.
5. Test independence and shape-preserving homogeneous prolongation with arbitrary strictly positive rational weights. Enforce that fine basis functions have nontrivial supports; reject dead basis controls.
6. Measure support-scaled coefficient/dual norms and verify along at least 3 nested levels and overlapping admissible activation patterns. A forward-stochastic `T` is not an inverse stability theorem. Use the 7D reference isolated-block lower bound as a diagnostic target, not as a global proof.
7. Derive an extraordinary-to-regular quintic box-spline C2 transition with exact jet matching, if and only if the core positive spatial transfer survives gates 1-6.
8. Run identical impulse, checkerboard and curvature tests against genuine NURBS; do not compare different control-site conventions without disclosing them.

## Stop conditions

- Any negative Bernstein coefficient in a claimed positive basis.
- Any rank loss, zero function, incorrect C2/C3 convention, missing affine moment, or nonpositive rational weight.
- A refinement equality checked only at sampled points rather than exactly on a common Bernstein overlay.
- Any claimed uniform stability contradicted by an explicit shrinking-support or high-overlap test.

Produce code, fixture hashes, exact rational witness or exact infeasibility witness, numerical diagnostics labeled as such, and an updated research report. Keep old negative evidence and avoid unsupported "CAD ready" claims.
