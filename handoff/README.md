# NPMS — Coding agent handoff: first practical 3D prototype

**Scope:** Implement an experimental, interactive C++ NPMS surface editor/evaluator for the **two authenticated fixed geometries** (valence 5 and valence 9), with headless tests and an independently generated numerical oracle. The package supplies loaded basis coefficients and scenes so the agent can build the app without any time-consuming symbolic reconstruction.

**Start here:** Send `agent_prompt/SEND_TO_CODING_AGENT.md` together with this entire ZIP to your coding agent. The prompt is self-contained and spells out all mandatory implementation gates and known mathematical limitations. Do not send only the prompt without the fixtures.

## Working assets

| File/directory | Contents |
|---|---|
| `agent_prompt/SEND_TO_CODING_AGENT.md` | Full implementation instruction to copy/send to a coding agent |
| `fixtures/fan5_full42_demo_v1.bin` | 100 complete fine basis functions; 18 original fixed-mesh parent controls; valence 5 |
| `fixtures/fan9_full42_demo_v1.bin` | 172 complete fine basis functions; 30 original parents; valence 9 |
| `fixtures/fan{5,9}_binary_spec.json` | Exact binary format, array shapes, SHA-256 and byte lengths |
| `fixtures/fan{5,9}_schema.json` | PS12 split geometry, ordered triangular Bernstein exponents, functions and parent mapping |
| `fixtures/fan{5,9}_full42_float64.npz` | Portable NumPy copy for independent numerical testing |
| `fixtures/fan{5,9}_golden_samples.json` | 120 + 216 independently evaluated test points with surface jets (through second derivatives) |
| `scenes/fan{5,9}_{plane,rational}.json` | Four deterministic 3D/weight scenes |
| `tools/reference_oracle.py` | Independent Python numerical Bernstein and rational jet evaluator (no SymPy) |
| `tools/check_handoff.py` | Fast numerical checks including every two-sided PS12 seam, affine precision, and binary hash |
| `tools/package_binary.py` | Regenerates C++-friendly binary and specs from `.npz` fixtures (no SymPy) |
| `tools/export_fixtures.py` | Recreates NPZ coefficients from authentic rational construction (slow symbolic reference only) |
| `docs/MATH_AND_DATA_CONTRACT.md` | Formal basis, derivatives, physical coordinate chain rule, rational shape and transfer semantics |
| `docs/ENGINEERING_TEST_MATRIX.md` | Mandatory testing, stop conditions, benchmark and screenshot expectations |
| `reference_research/NPMS_Full42_Davydov_Yeo_Positive_Research_2026-10-30.zip` | Untouched verified research archive as distributed before this handoff |
| `reference_research/NPMS_Full42_Edge_Research_2026-10-30/` | Extracted original exact Python, 42D local source, fixture metadata, history, audit and limitations |

## Fast independent verification

Requires Python 3.11+ and NumPy (SymPy only for historical exact/symbolic reconstruction). From this directory:

```sh
python -m pip install numpy
OPENBLAS_NUM_THREADS=1 python tools/check_handoff.py
OPENBLAS_NUM_THREADS=1 python tools/reference_oracle.py
python tools/package_binary.py
sha256sum -c MANIFEST_SHA256.txt
```

The first command may use whichever approved/offline dependency provisioning environment is available. The scripts never connect to network themselves. On Windows PowerShell set `$env:OPENBLAS_NUM_THREADS='1'` (or omit the environment variable); the Python scripts don't require Bash. Do not expect POSIX `sha256sum` to exist on Windows; `Get-FileHash -Algorithm SHA256` or Python's hashlib works.

Recorded checks performed on this handoff:

- Rebuilt complete floating-coefficient arrays from exact rational source for both original fans (a one-time symbolic computation).
- Verified in exact SymPy arithmetic that every coefficient is nonnegative; basis sums to one and fixed-mesh parent split identities hold. Compared affine moments inherited from original v1.7 basis.
- Exported f64 arrays and recovered all coefficients byte-for-byte from the C++ binary payload.
- Verified 336 independent point samples across microcells with both planar and positive-weight rational controls, through second derivatives; flat affine plane residual below `2e-9` in this oracle's assertions.
- Checked all 238 shared PS12 microedges (85+153), with three evaluation points per edge, for global function value, physical gradient and Hessian continuity; maximum observed discrepancy below `9e-15` in f64.
- Verified binary header sizes and SHA-256.

These are **engineering validations of two fixtures**. Mathematical proof-level claims and their scope are in the historical `RESEARCH_AUDIT.md`. Full valence-9 clean symbolic rebuild was previously time-limited, even though the original full rational run and this handoff's one-time fixture construction completed successfully. The new float golden tests must not be called independent exact-rational proofs.

## Mathematical and product limitations

- 42 dimensions per *local* macro-element, not 42 controls globally; fan global level-zero dimensions are 100/172.
- The original certified full positive family spans the entire global D–Y level-zero space on two **specific rational planar fan geometries**, but is poorly conditioned (numeric column-scaled coefficient condition ~16.6k/23.8k). Duplicate XY anchors are inevitable with this splitting construction.
- The 18/30-to-100/172 parent-child relationship is **fixed-mesh enrichment**, not spatial mesh refinement.
- No production-safe arbitrary topology, no refined-mesh positive transfer, no level-independent inverse stability, no certified regular-to-extraordinary chart join, no NURBS-quality curvature guarantee, no safe inverse fitting, and no self-intersection/immersion guarantee for arbitrary control edits.
- A control site may lie outside the polygon: parameter-site locations ensure affine precision, not directly intuitive/editable knot-based CAD control geometry.

## Paper/provenance references

- O. Davydov and W. P. Yeo, *Refinable C² Piecewise Quintic Polynomials on Powell–Sabin-12 Triangulations*, J. Comput. Appl. Math. 240 (2013), 62–73. https://oleg-davydov.de/s25refinable.pdf ; DOI: 10.1016/j.cam.2012.08.006 . This underpins the authentic nested spline-space construction; it does not provide the specific positive control geometry in this package.
- H. Speleers, *A normalized basis for quintic Powell–Sabin splines*, CAGD 27 (2010), 438–457; DOI: 10.1016/j.cagd.2010.05.001 . Independent normalized-basis result, not a proof of compatibility with this specific D–Y geometry extraction.
- Pixar OpenSubdiv: https://github.com/PixarAnimationStudios/OpenSubdiv — optional later benchmark for general-topology subdivision, NOT an NPMS core runtime dependency.

## Trust/security

Research source files are authored prototype code, not a security-audited dependency; review before running in a privileged environment. The legacy research archive contains Python `.pkl` files: never blindly unpickle external files. No executable dependency installation scripts, git pushes, live remote services or installers are needed to read the data. The binary is fixed little-endian with documented bounded layout; treat unknown inputs as untrusted. Validate any external dependency license and pin exact versions in the implementing repo.

## Deliverable expected from the coding agent

A cross-platform C++20 core plus interactive 3D viewer, independent headless sample/derivative tests, safe fixture loader and optional OBJ export, control editing/weights/undo, influence/curvature visualization, performance report and honest proof/limitation documentation. The agent should not attempt to deliver unrestricted manifold NPMS or verified multilevel refinement in this implementation project.
