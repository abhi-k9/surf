# Engineering acceptance matrix and rejection rules

Run the fast fixture checks independently before any graphics work:

```bash
python -m pip install numpy
OPENBLAS_NUM_THREADS=1 python tools/check_handoff.py
OPENBLAS_NUM_THREADS=1 python tools/reference_oracle.py
```

Symbolic reference rebuilds require SymPy and can be expensive. See `reference_research/NPMS_Full42_Edge_Research_2026-10-30/README.md`. **Do not run them on every edit/commit.** The snapshot's exact validations and the binary fixture checks serve different purposes.

## Automated headless tests

| Area | Required test | Acceptance target |
|---|---|---|
| Loader | Parse both bin files with exact header, byte lengths, array sizes and SHA256 | Identical to `fixtures/*_binary_spec.json` |
| Invalid data | Truncated, wrong version/magic, pathological counts, overflow, malformed scene JSON, NaNs | Graceful error, no crash, no allocation bomb |
| Geometry mapping | Reference `(0,3),(0,0),(3,0)` -> center/current/next outer vertex; micro triangles | World parameter locations agree with fixture within 1e-13 |
| Bernstein index | Golden sample basis values at interior points; endpoints and microcell vertices | Within 1e-10; no index permutation |
| Basis positivity | 100/172 functions, 12 cells per sector | No meaningful negative values (roundoff tolerance 1e-12) |
| Unity | `sum B_i = 1` per sample | Maximum absolute residual < 1e-10 |
| Affine moment | `sum B_i a_i = (x,y)` for unit-weight original anchors | Max absolute residual < 1e-10 |
| Differential basis | Independently test all `B,B_x,B_y,B_xx,B_xy,B_yy` | Physical derivatives, correct chain rule; golden surface below |
| Planar surface | Flat plane preset across all golden samples | All 6 jet rows max abs error <= 1e-9 |
| Weighted surface | Rational preset across all golden samples | All 6 jet rows max abs error <= 1e-9 or 1e-8 relative for large-magnitude terms |
| Seams | Every shared PS12 microedge and macro radial seam (both sides), t=.2,.5,.8 | Value, gradient and Hessian maximum abs difference <= 1e-8 |
| Unit-weight plane curvature | Mean and Gaussian curvature of regular flat plane | Near zero (< 1e-8), valid normal orientation |
| Degeneracy | Every `P_i` identical or area collapse | Explicit `DegenerateSurface` flag; never NaN silently |
| Weights | 0, negative, NaN, ±Inf | Reject before evaluation; don't clamp |
| Parent transfer | Sum all children assigned to each parent equals stored coarse Bernstein coefficient on all microcells | Max coefficient error <= 1e-12 |
| Geometry preservation | Copy parent homogeneous controls to child controls; evaluate coarse/fine surfaces | All sample position/jets agree <= 1e-9 |
| Out-of-domain | Evaluate beyond outer fan polygon | Explicit outside-domain result; no extrapolated cell selection |
| Editing | Single control XYZ update; single weight update; undo/redo/reset/save/reload | Deterministic surface response, restored exact state |
| Tessellation | Vary per-microtriangle grid resolution (e.g. 2, 4, 8, 16) | Uniform shared-edge vertex sequences, no cracks, no T-junctions, consistent winding |
| UI | Orbit/pan/zoom, pick/filter/group, support heatmap, curvature invalid flag | Demonstrated on both fan scenes with screenshots |
| Memory/security | Read malformed/fuzzed file sizes and scene JSON under sanitizer | No memory errors or uncontrolled allocation |
| Performance | 1k/10k/100k evals, tessellation 4/8/16/32, both n=5/9, CPU and release build | Report measurements/hardware; *no arbitrary pass threshold* |

## Non-regression tests for repeated developer iterations

- A headless `test_fixture_contract` target should finish rapidly and run on every pull request.
- A `test_python_oracle` target should run on Linux in CI or locally if Python+NumPy available; it should not download dependencies at runtime.
- C++ tests must not call the same C++ function to compute both the result and its expected value. JSON golden data was generated independently in Python.
- Use numerical finite differences only as a **secondary diagnostic**; near microcell seams, centered finite differences can cross cells and require care. Compare analytically computed one-sided jets across edges instead.
- Validate at least one arbitrary **positive rational weight** scene and geometry preservation under copying `Q=(wP,w)`.
- A source fixture hash changing should fail the pipeline until reviewed and the reason recorded.
- `-Wall -Wextra -Wpedantic`/MSVC warning equivalents and sanitizer tests should run for the core; keep render-related GPU/driver tests separate.

## Suggested performance record

| Example | Data to capture |
|---|---|
| CPU single-point evaluations | Values and derivatives for 100/172 basis funcs, max/median/throughput |
| CPU tessellation | Number of triangles, build ms, bytes allocated, seam vertex dedup count |
| Viewer | Frame time p50/p95, CPU vs GPU, edit-to-redraw time, sample controls active |
| Memory | Resident peak, geometry buffer allocations, binary fixture residency |
| Configuration | OS, architecture, compiler/toolset version, Debug/Release, core count and CPU, graphics device/driver, tessellation density |

A regular NURBS or OpenSubdiv result is meaningful only with disclosed matched conditions (parameter domain, degree, control sites, effective degrees of freedom, comparable deformation magnitude). Exact geometric equivalence is not asserted.

## Explicit stop/go

1. **Stop at loader** if binary counts/order do not match; do not patch arrays by transposition until root-cause verified against schema.
2. **Stop at evaluator** if flat plane analytic jets fail; do not hide error using higher tessellation or normal recalculation.
3. **Stop at seams** if world C² matching fails; do not average derivatives across an incorrect seam.
4. **Stop at data** if any original certified fixture is modified; regenerate via the archived exact mathematical source and document changes, or roll back.
5. **Do not imply production readiness** if conditions, topology, refinement or curvature behavior remain unverified.

## Screenshot/demo expectations

Capture valence5 and valence9 in: (i) unit-weight flat state with macro/PS12 wireframe, (ii) weighted curved/edit state with analytic-normal shading, (iii) selected-control influence heatmap and grouped parent UI, (iv) valid and invalid curvature visualizations. Include a short reproducible procedure to obtain each state, not just pictures.
