# AGENT MISSION — Build the first usable NPMS surface editor/evaluator

**Status:** implementation request. You are the coding agent. **Start coding**, using the data and specifications in this attached handoff archive. Do not hand back a speculative design or a research proposal in place of a working program. It is fine to divide your work into independently testable commits, but deliver usable software and proof of the tests you ran. Do **not** attempt to solve the outstanding NPMS mathematical conjectures as a prerequisite for an MVP.

## A. Read these files before you change code

1. `README.md` — inventory, provenance, reproducibility and commands.
2. `docs/MATH_AND_DATA_CONTRACT.md` — exact meaning of the basis coefficient arrays, world-coordinate mapping, derivatives and rational surface.
3. `docs/ENGINEERING_TEST_MATRIX.md` — required deterministic and performance tests, quality gates and failure policies.
4. `fixtures/fan5_binary_spec.json`, `fixtures/fan9_binary_spec.json` — runtime binary layout, checksums.
5. `fixtures/fan5_schema.json` and `fixtures/fan9_schema.json` — function/microcell/Bernstein ordering.
6. `reference_research/NPMS_Full42_Edge_Research_2026-10-30/README.md`, `RESEARCH_AUDIT.md` — mathematical provenance/limitations. Read these before interpreting any ``P5 complete'' or ``positive'' claim.
7. The archived Python `tools/reference_oracle.py` and deterministic JSON golden samples; they are *validation references*, not a substitute for a working C++ evaluator.

**Important distinction:** We possess independently certified full *fixed-fan* 100/172-function, nonnegative, normalized, C² quintic bases on two specific planar parameter fans. They span their full level-zero Davydov–Yeo spaces and reproduce global P5 polynomials. They have bad measured coefficient conditioning (~16,613/~23,846), duplicate nominal control sites, no certified stable spatial refinement, no general topology, no verified regular extraordinary box-spline join, and no demonstrated NURBS-level editing fairness. The app must NOT suggest otherwise. The smaller restricted 33/57 basis has better conditioning but is *not* P5-complete; keep the two categories distinct.

## B. Product objective

Implement a **cross-platform, locally runnable interactive 3D NPMS demonstration** for Windows, Linux (including Arch/CachyOS/KDE), and preferably macOS. The priority is mathematical correctness, deterministic repeatability, usable editing and quality comparison—not elaborate rendering. No cloud services required.

### Required MVP behaviors

- Load either supplied fixed-fan fixture (valence 5 or 9), verify schema, counts, finite values, file size and SHA-256 where available, and clearly display demo-only scope.
- Evaluate positive piecewise-quintic basis functions from precomputed Bernstein coefficients for any valid 2D parameter point in the fan. Find the correct macrosector and Powell–Sabin microtriangle. Reject out-of-domain points; resolve boundary points deterministically.
- Compute basis values, **physical parameter derivatives through second order**, using a stable algorithm. Do not finite-difference derivatives as the primary evaluator. Implement CPU doubles, and enable headless use (no GUI dependency) as a library/test tool.
- Compute the **rational surface** `S(x,y) = sum_i B_i(x,y) w_i P_i / sum_i B_i(x,y) w_i`, with strictly positive weights. Evaluate its first/second derivatives, normal, principal/mean/Gaussian curvatures where the immersion is regular. Provide clear invalid/degenerate state when Jacobian rank collapses. Do not silently `NaN`/clamp or invent curvature.
- Tessellate the *actual 12 PS12 microtriangles per macrotriangle*, with consistent subdivision along shared edges and no cracks/T-junctions. Use analytic normals, and allow optional wireframe display of macro edges, PS12 internal edges and tessellation.
- Interactive orbit/pan/zoom, select a control, edit its XYZ or its individual positive weight, undo/redo, reset scene, save/reload scene, and visualize the resulting surface interactively. Control selection should show affected support region / basis influence field. Explicitly show duplicated initial control sites in the full basis; implement a practical grouping/filtering UI (by original parent, macro, role, support) so a user can navigate 100/172 controls.
- Load four provided scene presets (`scenes/fan{5,9}_{plane,rational}.json`) and validate that edits actually change the surface. The flat unit-weight plane is the primary identity test, not a visually interesting demo. Include at least one editable small-amplitude shape preset; preserve the exact loaded states for deterministic testing.
- Support basis diagnostic modes: partition residual, minimum basis value, influence heatmap for selected control, normal and curvature (with an explicit invalid-state color), macro/PS12 overlays and an optional simulated checkerboard edit pattern.
- A headless CLI or test executable should output deterministic surface sample values, derivatives, timing counters, numerical diagnostic summaries, and export an OBJ of the tessellated mesh (with normals if possible). GUI rendering quality should not be conflated with geometric accuracy.

**Recommended architecture:** C++20 core library, CMake presets, strictly optional viewer target using a maintained and license-compatible lightweight graphics/UI stack (e.g. GLFW + Dear ImGui + OpenGL; another well-supported equivalent is acceptable if you document it). Keep core evaluator entirely independent from graphics. The fixture loader should use standard C++ and fixed little-endian binary parsing (no NumPy runtime, no internet at run time). Prefer a small dependency surface and pinned/verified dependencies. Implement a Python-only test oracle separately from the C++ algorithms so the same bug cannot be shared accidentally. Include a headless smoke test on CI.

## C. Exact data you are given

The attached archive contains all required sample data. You do **not** need to run any expensive SymPy extraction to build the viewer:

- `fixtures/fan5_full42_demo_v1.bin`, `fixtures/fan9_full42_demo_v1.bin`: portable little-endian binary, with *both* the complete fine basis and its 18/30 v1.7 parent basis. Each fine basis has **100/172** functions with layout `[function, macro, cell, Bernstein]`, each macro has 12 microtriangles and 21 degree-5 coefficients per cell. Arrays for sites, parent groups, macro vertices and microtriangle vertices included. Follow `fixtures/*_binary_spec.json` EXACTLY; do not guess indexing.
- `fixtures/*_float64.npz`: the same fine coefficients in compressed NumPy format for tests. Do not make NPZ a runtime dependency in C++.
- `fixtures/*_schema.json`: exact PS12 reference vertex IDs, microcell ordering and degree-five Bernstein triple ordering.
- `fixtures/*_golden_samples.json`: **120** samples on fan5 and **216** samples on fan9. Each has target `surface_jet_6x3` for planar and weighted curved scenes, with rows `[S, S_x, S_y, S_xx, S_xy, S_yy]`, and basis partition/minimum checks. The oracle uses a separate Python implementation with explicit Bernstein derivative finite-difference identities and affine physical-coordinate chain rule.
- `scenes/`: deterministic plane and rational positive-weight control-position scenes.
- `tools/reference_oracle.py`, `tools/check_handoff.py`: run in an ordinary Python environment with NumPy; fixture-based independent checks, including *both* value and physical-coordinate C² seam continuity.
- `reference_research/`: the **original research ZIP and extracted exact-rational scripts**, `full_positive_results.json` mode assignments, the authentic 42D PS12 local basis coefficient matrix, edge-star certificates and signed/positive decomposition. Preserve these as read-only source/provenance. The historical code contains `.pkl` files: don't unpickle unknown code/data without inspecting provenance and trust boundary. Exact scripts are expensive; never require a full SymPy rebuild on every CI commit.

## D. Mathematical contract — implementation cannot change it

Let `f_i` denote *nonnegative piecewise quintic* basis functions. For the complete fine basis on the fixed fans:

`sum_i f_i(x,y) = 1`, `sum_i a_i f_i(x,y) = (x,y)` for the original stored parameter sites `a_i`. With positive weights `w_i`, rational geometry is

`W = sum f_i w_i`, `A = sum f_i w_i P_i`, `S = A/W`.

These satisfy partition and affine coordinate reproduction with **unit** weights and sites, not the blanket claim that arbitrary changed rational weights retain fixed affine-site interpolation. Homogeneous control transfers are separate from actual spatial refinement.

Map reference macro vertices **v1=(0,3), v2=(0,0), v3=(3,0)** respectively to **O=(0,0), P_i, P_(i+1)**; do not assume familiar vertex ordering. Use the exact supplied world-space microtriangle vertex coordinates to compute barycentric coordinates; the Bernstein triple `(i,j,k)` always applies to the corresponding `[vertex0,vertex1,vertex2]`. Derivative chain rule and quotient formulas must match `docs/MATH_AND_DATA_CONTRACT.md`.

Never silently renormalize already normalized basis functions. Do not clip negative coefficients or suppress seam errors. It is permissible to tolerate floating roundoff (~1e-12 relative) while diagnosing, but not to mask systematic errors. The precomputed coefficients were converted to `double` from certified rational inputs. Preserve coefficient ordering and exact certified archive intact.

## E. Phased delivery (complete each gate before expanding scope)

**Milestone 0: Inventory and fixture contract (mandatory before viewer).** Read all fixture specs, validate binary hashes, use a standalone headless loader, print dimensions, min/max coefficients, count of controls and grouped parent IDs. Provide a short `PROVENANCE_AND_LIMITS.md`. Do not auto-merge external repository changes; if the target repo is unknown, work in a new local source tree and document integration steps.

**Milestone 1: Evaluator correctness.** Implement Bernstein scalar and all first/second derivatives in physical fan coordinates. Demonstrate equality against every golden sample. Test positivity, sum-to-one, affine moments, continuity and rational jets on all supplied scenes. Establish a no-viewer core CI path.

**Milestone 2: Real interactive viewer.** CPU tessellation and analytic normals; orbit/pan/zoom, control picking/editing/weights, fast redraw, duplicate-site grouping/filtering, shape reset and persistence, visually inspect seams. Provide screenshots or a short screen recording and a build/run guide. Initial implementation may be CPU-driven; GPU basis evaluation is out of scope until correctness passes.

**Milestone 3: Diagnostics and reproducible performance benchmark.** Add selected-basis influence heatmap, macro/PS12 seams, curvature with degeneracy flags, regular edit/impulse/checkerboard modes, per-frame CPU profiling and memory measurements. Record results on valences 5/9 at multiple tessellation levels, including build type/hardware/platform. Do not claim a NURBS comparison is fair without matching control-site/effective DOF conventions and documenting that mapping.

**Milestone 4: Optional independent comparator.** Add a genuine regular tensor-product **quintic** B-spline/NURBS patch and/or OpenSubdiv baseline behind a separate module, only once milestone 1–3 have passed. Clearly distinguish regular chart vs extraordinary fan and record like-for-like sample/error criteria. No mandatory OpenSubdiv dependency to build NPMS core.

**Explicitly defer mathematical research tasks:** arbitrary manifold topology; arbitrary valence; C² regular-to-extraordinary patch joins; true D–Y red spatial refinement; nonnegative coarse/fine transfer under spatial subdivision; level-independent coefficient stability; high-order fitting; trimming/boolean/CAD interchange. Do not implement guessed versions of these or present the demo as solving them.

## F. Test and performance requirements

For each fan and fixture:

1. Binary and NPZ/oracle data shapes, finite-value and boundary checks; reject truncated/malicious files and `NaN`/`Inf`.
2. Basis positivity within roundoff and partition residual <= `1e-10` (double). Stored Bernstein coefficients are nonnegative; negative *nontrivial* evaluated values are a bug.
3. Unit-weight plane with original `a_i` XY sites: `S=(x,y,0)`, `S_x=(1,0,0)`, `S_y=(0,1,0)`, `S_xx=S_xy=S_yy=0` with max absolute error <= `1e-9` across supplied samples. Curvature zero where valid.
4. Rational nonuniform-weight scene: sample-wise surface position and all physical parameter derivatives <= `1e-9` absolute, with a relative threshold of `1e-8` for exceptional large-magnitude values. Compare against JSON golden tests; report worst case.
5. Across shared PS12 microcell seams and macrosector interfaces: compare one-sided positions, gradients and Hessians at multiple edge points. Reference Python reports largest physical Hessian mismatch < `9e-15` in double; a C++ bound of `1e-8` absolute is an initial engineering tolerance, not a proof.
6. Full 100/172 fine basis should numerically span its stored 18/30 coarse parent basis by grouping each control's parent ID. Check summed Bernstein coefficients across every microcell and rational surface invariance after exact **fixed-mesh** homogeneous control copying. Do **not** label this actual spatial refinement.
7. Weights must be finite and >0; negative/zero/NaN cause validation error. Degenerate parameter triangles or zero area must be rejected. Degenerate output surfaces should return a labeled invalid normal/curvature rather than divide by zero.
8. Implement consistency tests for Bernstein evaluation at cell corners, multiple points near boundaries, mixed coordinate transformations and finite differences used *only to validate* analytic first/second derivatives.
9. Run sanitizer-capable builds (ASan/UBSan where supported), warnings-as-errors in core, memory/runtime profiling and headless tests on the platforms CI supports. Keep CI economical: fast fixture checks per PR, heavy exact rational rebuild only on manual research command.
10. Provide useful benchmark output (ms for repeated evaluation batches and tessellation, basis function count, tessellation triangles, memory consumed, FPS with hardware/build info). No invented performance targets, and do not sacrifice correctness to satisfy an arbitrary FPS threshold.

See `docs/ENGINEERING_TEST_MATRIX.md` for the full gate checklist, including expected behaviors under malformed data.

## G. User interface and editing semantics

The positive full basis is an **existence construction**, not a conventional intuitive sparse CAD cage. Many fine controls inherit the *same* parametric anchor from the parent and are highly correlated. The agent must not pretend these are ordinary unique Greville sites. Display labels like `parent 07 / child outer_vertex 3 jet 6` from the supplied schema, with group editing and one-at-a-time controls. A selected control should show its *true normalized influence* `f_i w_i/W`, not approximate radial falloff. Avoid editing baselines by unconstrained adjustment of parameter anchors; those anchors are metadata for affine precision and aren't an inverse-geometry parameterization. The 3D control `P_i` can move freely, but weights remain strictly positive.

Treat separate control editing as an evaluation experiment, **not proof of intuitive shape design**. Provide undo/redo, visible curvature response, and a safe flat reset. Avoid claiming no self-intersections for arbitrary edits.

## H. Code quality, security, reproducibility

- C++20, RAII, no unsafe raw owning pointers, `-Wall -Wextra -Wpedantic` equivalents, consistent formatting/static analysis, address/undefined sanitizers where supported.
- CMake with **verified pinned** third-party dependencies; do not download and execute arbitrary scripts. Audit third-party license and keep source-to-binary build reproducible. Viewer deps optional; core and tests build offline.
- Validate all untrusted scene/fixture input, check integer overflow in buffer-size arithmetic and bound memory before allocating. Use little-endian loads explicitly (or reject other endianness with a meaningful error). Do not silently ignore unknown binary versions.
- No shell/network execution from scene files. Avoid unsafe Python pickle loading. Keep source fixtures unmodified; store user scenes separately.
- Clear separation: `libnpms` math, fixture loader, application GUI, benchmarks, tests, reference Python. Use immutable coefficient cache, independent mutable control scene state, and explicit validity flags.
- CI matrix with Linux and Windows; macOS if practical; compile headless core without graphics dependencies. Cache external deps conservatively; avoid superfluous repeated heavy CI.
- Honest documentation of numerical tolerances, limitations, demo topology and known poor basis conditioning.

## I. Required completion report

Finish with: source tree / branch or commit ID (if git configured), build commands for each tested platform, actual commands and results of automatic tests, where sample comparisons failed or passed, screenshots/video of valence 5 and 9 interactive editing, timings with hardware, main unresolved bugs, data/fixture hashes, and a next-work prioritized list. If a core gate fails, report it immediately and do not bypass it by weakening tests or changing the supplied mathematics without explicit approval.

**Success definition:** a real, testable local application showing correctly evaluated and editable 3D rational NPMS surfaces on the *two certified fan configurations*—not a general-topology CAD product and not merely a plan.
