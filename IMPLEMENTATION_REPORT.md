# NPMS prototype v0.1 — implementation and validation report

Source: the root of [abhi-k9/surf](https://github.com/abhi-k9/surf), branch `main`, checked out at `/workspace/surf`. The original standalone implementation was prepared and validated at `/workspace/npms`; historical evidence records retain that working directory. Repository integration checks are in `evidence/repository-tests.txt`. The uploaded archive and its extracted certified assets are preserved. The earlier standalone delivery archive `NPMS_Prototype_v0.1.zip` contains source, data, evidence, and an optional Linux x86-64 executable.

## Milestone results

**0 — Inventory/data contract: passed.** Both SHA-256-authenticated binary layouts match exact expected sizes/counts. Fan5 has 100 functions, 18 parents, 60 PS12 cells and nonnegative coefficient range [0, 0.38333333333333336]. Fan9 has 172 functions, 30 parents, 108 cells and range [0, 0.40833333333333333]. All original manifest entries still match. The reference-to-world vertex mapping and every ordered Bernstein triple are tested against the supplied schema. Inventories and support maps are in `evidence/fan*-inventory.json`.

**1 — Evaluator: passed.** A graphics-independent C++20 static library implements physical-coordinate analytical jets with double precision. For each barycentric coordinate it constructs a six-component derivative jet; polynomial multiplication differentiates each Bernstein monomial analytically through order two. This differs from the independent oracle's coefficient finite-difference identities. It handles zero barycentric coordinates without division by a coordinate. The quotient recurrence implements the supplied rational contract exactly. Weight scaling and centering homogeneous control positions reduce overflow and translation cancellation; these operations change neither basis functions nor geometry. No coefficient clipping, basis renormalization, seam averaging or finite-difference primary evaluator is used.

Every golden point is tested in both supplied scenes: 120 points on fan5 and 216 on fan9, for **672 surface jets / 12,096 scalar jet components**. Every individual basis value and all five derivatives are independently compared against the NPZ/Python oracle at those same points. Additional tests cover partition, affine moments, cell corners, actual physical mapping, finite-difference secondary diagnostics, every shared microedge at three points, fixed-mesh parent identities and rational control copying, malformed inputs, positive weights, regular/degenerate surfaces, and exact paraboloid/saddle principal curvature checks.

| Final maximum absolute error | Fan5 | Fan9 | Required bound |
|---|---:|---:|---:|
| Golden surface jets (both scenes) | 1.41917e-13 | 5.38458e-14 | 1e-9 |
| All individual basis jets vs independent NPZ oracle | 2.04281e-14 | 6.99441e-15 | 1e-9 |
| Unit-plane identity jets | 1.40828e-13 | 6.19181e-14 | 1e-9 |
| Partition residual | 1.77636e-15 | 1.99840e-15 | 1e-10 |
| Affine moment | 2.44249e-15 | 6.66134e-15 | 1e-10 |
| Two-sided basis value/gradient/Hessian | 3.90799e-14 | 1.52933e-14 | 1e-8 |
| Parent coefficient sum | 0 | 0 | 1e-12 |
| Homogeneous fixed-mesh rational transfer jets | 8.61533e-14 | 5.44009e-14 | 1e-9 |
| Secondary centered finite-difference derivative check | 8.10928e-9 | 7.10385e-9 | 1e-6 diagnostic |

**2 — Interactive application: passed on Linux.** A local Python-standard-library bridge runs persistent C++ evaluators. The dependency-free WebGL 2 UI supports both fans, all four provided scenes, an additional wave/checkerboard preset, orbit/pan/zoom, projected control picking, overlap cycling, parent/role/actual-support filters, XYZ/weight/group edits, undo/redo/reset, JSON save/reload and OBJ export. CPU tessellation subdivides the actual supplied PS12 cells and welds shared parameter points. Connectivity, boundary-edge counts, Euler characteristic, triangle counts and winding pass at levels 2/4/8/16. Rendering uses float32 display buffers; all numerical validation uses authoritative C++ binary64 results.

Browser tests verified editing changes scene state and the C++ regression separately verifies edits change surface jets. They also verified exact saved position/weight reload, undo/redo, grouped edits, invalid-weight rejection, control picking at a projected position, OBJ normals, and invalid curvature. Ten screenshots are in `evidence/`, with reproduction instructions in its README. The collapsed surface screenshot correctly shows one magenta point and an explicit notice rather than inventing a regular surface.

**3 — Diagnostics/performance: passed.** True rational selected-control influence, normal directions, signed mean/Gaussian curvature, partition residual, minimum basis value, macro/PS12/tessellation overlays, impulse and checkerboard experiments are available. Invalid normal/curvature is explicit in the API/UI. Mesh CPU time, logical buffer memory, draw submission percentiles, evaluation-batch distributions, process peak RSS, and observed browser frame cadence are measured. No arbitrary FPS threshold was used.

**4 — Optional comparator: deferred.** No NURBS/OpenSubdiv implementation or fairness/performance comparison is claimed. The next priority is validating more platforms and refining the editor without extending the certified mathematical scope.

## Actual test commands and outcomes

Executed on Debian 13, Linux x86-64, GCC 14.2.0, CMake 3.31.6, Python 3.12.14:

```sh
cmake --preset release
cmake --build --preset release --parallel 2
ctest --preset release -V
cmake --preset sanitize
cmake --build --preset sanitize --parallel 2
ctest --preset sanitize -V
OPENBLAS_NUM_THREADS=1 python handoff/tools/check_handoff.py
python tests/test_http.py build/release/npms_cli
python tests/test_viewer.py
python tools/benchmark.py
```

Final Release and sanitizer runs each passed **3/3 CTest suites**, including 208,350 C++ assertions, 15 malformed worker requests with recovery, full independent NPZ basis/surface comparisons, and the original manifest check. ASan, UBSan and LeakSanitizer reported no application errors. LeakSanitizer initially could not run under the managed sandbox's ptrace layer; the full test was then rerun outside that layer with detection enabled and passed. No leak check, mathematical tolerance or assertion was weakened. Browser tests initially needed CSP-compatible function predicates and an asynchronous file-read wait; these test-harness issues were corrected and the suite passed.

The separate HTTP test passes valid geometry, Host/Origin restrictions, named-path allowlist, JSON content type, malformed data, invalid weights and recovery. Core builds with warnings as errors. clang-format validation, JavaScript syntax checks and Python compilation pass. Source-to-archive clean-build evidence is recorded in `evidence/reproducible-build.txt`.

The Linux/Windows/macOS GitHub Actions workflow is included but has not been executed remotely. Windows, macOS, native desktop GPUs, accessibility across browsers, and long interactive sessions remain unverified. No test was skipped or expected-failed in the recorded Linux CTest runs. Optional tests are only registered when their prerequisites exist; this machine did have NumPy and executed them.

## Performance measurements

The following table is populated from the actual JSON benchmark records. Hardware: shared cloud AMD EPYC 9V74, five visible logical CPUs. Release uses GCC `-O3 -DNDEBUG` without fast-math; benchmark/evaluation/tessellation are serial. Timings are observations, not guaranteed budgets. Logical mesh bytes exclude fixture storage, allocator overhead and JSON/WebGL buffers.

| Measurement | Fan5 | Fan9 |
|---|---:|---:|
| 1,000 evaluations, median of 3 | 5.13 ms | 6.89 ms |
| 10,000 evaluations, median of 3 | 59.48 ms | 75.87 ms |
| 100,000 evaluations, median of 3 | 562.28 ms | 750.05 ms |
| PS12 mesh level 4 | 2.87 ms | 6.92 ms |
| PS12 mesh level 8 | 10.46 ms | 25.15 ms |
| PS12 mesh level 16 | 41.26 ms | 102.75 ms |
| PS12 mesh level 32 | 166.11 ms | 415.25 ms |


Detailed min/median/max for 1k/10k/100k batches are in `evidence/fan*-benchmark.json`. At level 32, meshes contain 61,440 / 110,592 triangles and about 8.24 / 14.84 MiB of logical geometry buffers. Peak benchmark-process RSS was about 15.2 / 29.7 MiB (the recorded `getrusage` field is the maximum child-process RSS seen by the harness to that point).

| Browser level 8 | Fan5 | Fan9 |
|---|---:|---:|
| Observed animated cadence | 42.3 FPS | 33.8 FPS |
| Frame cadence p50 / p95 | 21.6 / 36.0 ms | 27.8 / 42.0 ms |
| CPU submission p50 / p95 | 4.7 / 12.1 ms | 7.1 / 19.7 ms |
| Last edit-to-response observation | 61.7 ms | 89.3 ms |


Browser: Chromium 154.0.8037.57, WebGL 2 through ANGLE Vulkan **SwiftShader software rendering**, 1440×1100 viewport, level 8, 120 animated frames per fan. Frame cadence includes browser scheduling/compositing, while CPU submission measures the JS render call. An isolated GPU execution timer was not measured. Edit-to-response timings include C++ work, serialization and loopback transfer. These are not native-GPU or standalone workstation FPS claims.

## Review findings and remaining limitations

- Only the two authenticated supplied binaries are accepted. This is intentional. Supporting arbitrary modified fixtures would require a separately reviewed format/trust contract.
- Physical derivatives and curvature use the supplied parameter geometry; both differently oriented microcells and macroseams are checked. One-sided results are not averaged.
- The evaluator is double precision; very extreme finite coordinates/weight ratios can still overflow/underflow. It returns an explicit arithmetic/denominator error rather than clamping or publishing nonfinite geometry.
- Duplicate/correlated controls and conditioning ~16.6k/~23.8k remain properties of this full positive construction. The UI makes them navigable but does not establish intuitive shape fairness or robust inverse fitting.
- Editing uses picking plus numeric XYZ/weight inputs and group delta operations. There is no drag gizmo, adaptive tessellation, incremental basis cache or GPU spline evaluator. Level-32 JSON transfer and CPU draw-buffer rebuilding can be slow.
- No general topology, arbitrary valence, stable spatial refinement, regular extraordinary transition, self-intersection guarantee, trimming, fitting or production CAD claim is made. Fixed-mesh homogeneous transfer remains distinct from refinement.
- New evaluator/bridge helpers have regression and malformed-input coverage, but have not received an external security audit. The server is loopback-only and is not a public/multiuser service.
- The historical source papers and exact-rational scripts are provenance. This delivery did not rerun expensive exact proofs or establish a new theorem.

Next work, in priority order: execute the supplied Windows/macOS CI matrix and browser checks on native GPUs; add sustained interaction/fuzz coverage and a drag gizmo; reduce serialization and cache mesh/basis work after correctness checks; then consider a separate, genuinely quintic comparator with explicit effective-DOF/site mapping. Mathematical refinement/topology questions remain separate research tasks.

See [README.md](README.md) for complete Linux and Windows build/run commands, [PROVENANCE_AND_LIMITS.md](PROVENANCE_AND_LIMITS.md) for exact fixture hashes and research scope, and [evidence/README.md](evidence/README.md) for artifacts and screenshot reproduction.
