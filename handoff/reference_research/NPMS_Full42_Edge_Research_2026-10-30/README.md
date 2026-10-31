# NPMS full 42D / positive edge-curvature integration

**Status:** exact positive **complete 100/172-dimensional geometric basis** existence certificates in the full 42D Davydov–Yeo quintic local architecture, plus better-conditioned restricted edge and Hessian enrichments, tested on two valence-5 and valence-9 planar fans. **Not** a uniformly stable spatial-refinement geometry kernel.

Read `RESEARCH_AUDIT.md` first; open `RESEARCH_EXPLORER.html` in a browser for a self-contained offline visual guide.

## Two levels of new construction

- Better-conditioned positive radial-edge basis **18→33** for valence 5; **30→57** for valence 9.
- An optional local Hessian-plus-edge basis gives **51/87** controls, independent second-order vertex jets but still incomplete quadratic reproduction.
- Exact full-dimensional positive normalized affine-precise bases **100/172**, spanning the entire level-zero D–Y spaces and therefore reproducing **P5**. However, the full bases have column-scaled coefficient condition numbers around **16,613/23,846**; they are existence certificates, not production bases.
- Every new edge mode supported on two macro sectors; assigned outer parent never exceeds four.
- Exact Bernstein positivity, partition of unity, affine-coordinate precision, independent edge-curvature controls, sparse nonnegative fixed-mesh transfer and rational geometry preservation.
- Child strength at least `10123/25600` on fan 5 and exactly `1/2` on fan 9.
- Full 42D PS12 element reproduces all 21 P5 monomials exactly; the **complete** positive 100/172 geometric bases inherit full P5 reproduction (coefficients may be signed), while the restricted 33/57 family does not.
- Archived eight-dimensional two-macro edge-star: seven independent positive modes plus one necessarily signed first-normal mode.
- Rigorously scale-independent **isolated/disjoint** positive edge-patch norm bound, but no overlapping multilevel stability proof.

## Reproduce

Python 3.11+; `numpy` and `sympy` required. The independently verified seam and positive fan code uses the bundled historical source files. Typical invocation:

```bash
python -m pip install -r requirements.txt
OPENBLAS_NUM_THREADS=1 python edge_enrichment.py
OPENBLAS_NUM_THREADS=1 python verify_p5_analytic.py
OPENBLAS_NUM_THREADS=1 python full_geometric_enrichment.py
OPENBLAS_NUM_THREADS=1 python full_space_positive.py
python verify_signed_edge_decomposition.py
python verify_scale_block.py
python verify_full_metadata.py   # fast structural check, does not repeat Bernstein extraction
(cd v12_reference && OPENBLAS_NUM_THREADS=1 python verify_committed_certificates.py)
(cd v13_seam_reference && OPENBLAS_NUM_THREADS=1 python seam_world_exact.py)
sha256sum -c SHA256SUMS.txt
```

`edge_enrichment.py`, `full_geometric_enrichment.py` and especially `full_space_positive.py` may require substantial SymPy rational computation because they rebuild the physical-coordinate 42D extraction on both fans. The full 172-control Bernstein certification may exceed short execution windows; use the fast metadata script for a separate limited check, not as a replacement for the coefficient verifier. No external downloads are required for validation. Research archives are untrusted code until examined; review scripts before execution.

## Structure

| Path | Purpose |
|---|---|
| `edge_enrichment.py` | Exact rational, outer-local positive edge-mode insertion and transfer check (33/57 controls) |
| `full_geometric_enrichment.py` | Exact positive Hessian+edge splits (51/87 controls); demonstrates remaining quadratic obstruction |
| `full_space_positive.py` | Exact full-dimensional, P5-complete positive basis existence construction (100/172 controls) |
| `full_positive_results.json` | Rational mode assignments, perturbation scales, conditioning of complete positive bases |
| `verify_full_metadata.py` | Fast independent rational assignment/scale/dimension checks; not a substitute for the full exact run |
| `verify_signed_edge_decomposition.py` | Exact 7-positive + 1-signed decomposition of the two-macro edge-star |
| `enrichment_results.json` | New coefficients, parent assignment and condition diagnostics |
| `base_v17/` | Archived 42D local reference and 18/30 affine-precise fan base |
| `positive_physical_coefficients.json` | Archived 3-per-radial-edge positive curvature modes |
| `verify_p5_analytic.py` | Exact local P5 analytic reproduction of all 21 monomials |
| `v12_reference/` | Independent 8D edge-star, 7 positive modes and fixed-patch certified bound |
| `verify_scale_block.py` | Affine-copy, disjoint-support block bound |
| `v13_seam_reference/` | Exact physical C2 seam verifier |
| `RESEARCH_AUDIT.md` | Full method, proofs, scope limitations, and next proof obligations |
| `AGENT_NEXT.md` | Precise next-experiment engineering and mathematical specification |

**No claim:** spatial red-refinement positivity, uniform conditioning of the full 100/172 positive bases or overlapping hierarchical bases, general manifolds, practical control placement, or NURBS curvature fairness. The restricted 33/57 positive basis remains polynomial-incomplete; the complete 100/172 basis is polynomial-complete but poorly conditioned.
