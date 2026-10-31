# Evidence index

- `release-tests.txt`: verbose CTest output for the final Release build; three suites passed.
- `sanitizer-tests.txt`: the same final code under ASan, UBSan and LeakSanitizer, with leak detection enabled outside the tracing sandbox; three suites passed.
- `handoff-validation.txt`: original independent Python handoff check, including 238 shared edges.
- `fan{5,9}-inventory.json`: authenticated binary dimensions, coefficient range, parent group sizes and actual support sectors.
- `fan{5,9}-benchmark.json`: 1k/10k/100k evaluation batches (three repeats) and PS12 tessellation at 4/8/16/32, with hardware/toolchain, memory and timings.
- `viewer-tests.txt`, `viewer-tests.json`: successful browser interaction report, Chromium/software-renderer identity, 120-frame profiles, and saved scene states.
- `http-tests.txt`: local bridge functional/security/recovery checks.
- `fan{5,9}-{plane,edited,influence,curvature,invalid}.png`: ten actual screenshots from the software, not illustrative renders.
- `fan{5,9}-saved-scene.json`: edited scenes saved and reloaded by the browser.
- `fan{5,9}-edited.obj`: browser exports of edited surfaces with C++ analytical normals.
- `fan{5,9}-cli.obj`: headless exports of the supplied weighted presets.

Screenshot reproduction: start the local server, run `python3 tests/test_viewer.py` with Playwright and `/usr/bin/chromium` installed. For each fan the test captures the unit plane with macro/PS12 overlays, changes control 0 Z and weight, captures shading, filters parent 1 for true influence, captures curvature, tests group editing/persistence/picking, and collapses geometry for the invalid-state screenshot. Files contain the actual states at capture time, not the final restored preset. The test then resets the fan.

Benchmark caveats: measurements come from this shared cloud Linux environment. Browser rendering is Chromium SwiftShader (software), not a workstation GPU. Frame cadence and CPU draw submission are distinct measurements; no isolated GPU timer result is claimed. Logical geometry-buffer bytes exclude allocator overhead, JSON serialization, renderer buffers and the immutable fixture arrays; process RSS includes the whole benchmark CLI. No frame-rate threshold or NURBS comparison is asserted.
