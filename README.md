# NPMS · Beyond NURBS — interactive fixed-fan prototype

Working C++20 rational spline evaluator and offline interactive WebGL 2 surface editor for the two supplied certified fans. All geometry, analytical derivatives, normals, curvature, influence and CPU tessellation come from C++. The browser renders the computed mesh and edits scene state; it does not implement an alternative spline evaluator.

Read [PROVENANCE_AND_LIMITS.md](PROVENANCE_AND_LIMITS.md) and [IMPLEMENTATION_REPORT.md](IMPLEMENTATION_REPORT.md). This is an experimental modeling application for the fixed valence-5 and valence-9 fixtures, not a production CAD kernel.

## Build and run

Requirements: C++20 compiler (tested GCC 14.2), CMake 3.21+, and Python 3.11+ for the optional viewer bridge. Core library/CLI/tests have **no third-party dependencies**, no graphics dependency and no network requirement. NumPy is optional and only registers the independent oracle test if already installed. No dependency or fixture download occurs during configure, build, tests or runtime.

From this source directory, Linux/macOS:

```sh
cmake --preset release
cmake --build --preset release --parallel 2
ctest --preset release
python3 viewer/server.py --cli build/release/npms_cli
```

Open the loopback address printed by the server in your browser. Keep the terminal running. Ctrl+C stops the server and its C++ workers. WebGL 2 is required; Chrome/Chromium, Firefox or Edge with WebGL enabled should work. There is no public hosting or cloud service dependency.

Windows, from a Visual Studio 2022 Developer PowerShell with CMake/Python available:

```powershell
cmake --preset release
cmake --build --preset release --parallel 2
ctest --preset release
py viewer/server.py --cli build/release/Release/npms_cli.exe
```

With a single-configuration Windows generator, the executable is directly in `build/release/`. On Arch/CachyOS use the distribution's `base-devel`, `cmake` and `python` packages. These Windows/macOS instructions and their CI jobs are supplied but **were not executed here**. Linux x86-64 is the platform verified in this delivery.

The delivery archive also includes `bin/linux-x86_64/npms_cli`, built on Debian 13 / glibc 2.41. On a compatible Linux host, use `python3 viewer/server.py --cli bin/linux-x86_64/npms_cli`. Build from source on other hosts or older libc versions. The binary is an optional convenience, not a universal Linux binary.

A direct dependency-free core build is also possible:

```sh
c++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror -Iinclude src/core.cpp src/main.cpp -o npms_cli
```

## Editing and diagnostics

- Choose fan 5 or 9 and one of the four supplied plane/weighted scenes, or the additional small-amplitude wave/checkerboard experiment. Reset restores the selected preset. Original scene files are never overwritten.
- Choose a control in the explorer. Filter by original parent, role, or actual nonzero coefficient support sector. The metadata reports inherited duplicates. `Fit controls` reveals the large, repeated cage sites; click visible points to pick and click again to cycle overlaps. `Fit surface` returns to the modeled patch.
- Edit XYZ or a strictly positive weight, then Apply. The group checkbox applies the same position and weight delta to all children of that parent. The impulse button adds 0.5 to Z. Undo/redo restore exact scene snapshots (up to 100 operations). The original parameter anchors are never edited.
- Drag to orbit, Shift/right-drag to pan, scroll to zoom; Top view and Fit buttons reset framing.
- Toggle macro, actual PS12, or tessellation edges. Each microtriangle uses the same subdivision count, and shared parameter vertices are welded. Surface normals are analytical, not triangle-averaged.
- Influence shows the true rational value `f_i*w_i/W`, scaled by its **displayed sampled maximum**. Mean/Gaussian curvature use a symmetric nonlinear color scale with the sampled peak disclosed. Normal colors encode the analytic unit normal. Partition and minimum basis modes expose numerical residuals. Magenta means invalid normal/curvature; collapse intentionally tests this state.
- Save downloads scene JSON. Load validates its schema, fan, counts, finite coordinates and positive weights. Geometry and weights round-trip exactly; arbitrary metadata and fixture paths do not control filesystem access. OBJ exports current positions and analytical normals.

## Headless commands

```sh
build/release/npms_cli handoff 5 inventory
build/release/npms_cli handoff 9 sample handoff/scenes/fan9_rational.json
build/release/npms_cli handoff 5 obj handoff/scenes/fan5_rational.json 16
build/release/npms_cli handoff 9 benchmark
python3 tools/benchmark.py --cli build/release/npms_cli
```

`sample` emits one deterministic JSON record per supplied golden point, with rows `[S,Sx,Sy,Sxx,Sxy,Syy]`, curvature validity and diagnostics. `obj` writes `fan5.obj` or `fan9.obj` in the working directory. `benchmark` reports three-run evaluation distributions (1k, 10k, 100k samples), tessellation levels 4/8/16/32 and logical mesh-buffer bytes. Timings naturally vary; numerical output is deterministic for a given toolchain.

`worker` accepts one JSON request per stdin line, keeping the authenticated fixture resident. Example request: `{"xy":[0.2,0.2]}`. Add a `scene` object to override the weighted preset; alternatively send `{"resolution":8,"selected":0}` for a mesh. Unsupported/invalid input returns an error record; no shell, network or arbitrary path operation is exposed by worker requests.

## Validation

```sh
cmake --preset sanitize
cmake --build --preset sanitize --parallel 2
ctest --preset sanitize
OPENBLAS_NUM_THREADS=1 python3 handoff/tools/check_handoff.py
python3 tests/test_http.py build/release/npms_cli
```

ASan/UBSan/leak checking is supported with GCC/Clang. In this managed environment LeakSanitizer had to run outside the ptrace sandbox; all checks remained enabled. MSVC uses warnings as errors but the GNU sanitizer preset is not a Windows sanitizer configuration.

Optional browser tests require separately installed Playwright plus Chromium. With the local server on port 8765, run `python3 tests/test_viewer.py`. The test currently uses `/usr/bin/chromium`; adjust that **test harness** path for other systems. It records ten screenshots, edits and persistence, verifies projected control picking, and measures 120 animated frames on each fan. The application itself has no Playwright dependency.

**Do not run** `handoff/tools/reference_oracle.py` as a verifier: its main function regenerates golden files despite its historical docstring. `tests/test_oracle.py` imports its independent evaluator without invoking regeneration and checks all original hashes afterward. No pickle is loaded. Exact symbolic reconstruction is an explicitly separate research operation, not CI setup.

## Source layout and integration

- `include/npms/core.hpp`, `src/core.cpp`: reusable `npms` library; scene validation, basis/rational jets, curvature, tessellation, OBJ.
- `include/npms/json.hpp`, `sha256.hpp`: bounded JSON and SHA-256 helpers; no external library requirement.
- `src/main.cpp`: CLI/line-oriented worker and performance measurements.
- `viewer/`: optional Python standard-library loopback bridge and dependency-free HTML/CSS/WebGL client.
- `tests/`: numerical, malformed-input, independent oracle, HTTP and browser tests.
- `handoff/`: unmodified supplied assets and research provenance.
- `evidence/`: actual test output, inventories, screenshots, scenes, OBJ exports and benchmarks.

The source is maintained in [abhi-k9/surf](https://github.com/abhi-k9/surf) on `main`, at the repository root. The original implementation and recorded validation were prepared in `/workspace/npms`; the prototype has since been imported into `/workspace/surf` and its headless checks rerun there. Earlier evidence logs retain their actual original working directories. The complete original handoff is preserved under `handoff/`. This cloud task already has an isolated environment; additional Git worktrees are unnecessary.

The GitHub Actions configuration covers Linux, Windows and macOS headless builds, plus Linux sanitizer/NumPy regression. Those hosted jobs have **not** been run here. Review [THIRD_PARTY.md](THIRD_PARTY.md) for dependency and licensing scope.
