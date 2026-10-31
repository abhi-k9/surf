# Provenance and limits

The complete uploaded `NPMS_Prototype_Coding_Agent_Handoff_2026-10-30.zip` was extracted into `handoff/`. Its original files are unchanged. `handoff/MANIFEST_SHA256.txt` authenticates them against the uploaded manifest; this is integrity verification against the supplied source, not independent certification of its author. The original research ZIP remains inside `handoff/reference_research/`. Historical pickle files are preserved as data and are never loaded by this application or its tests.

The source of the mathematical construction is the supplied research audit, authentic local 42-dimensional Davydov–Yeo PS12 element, and supplied exact-rational construction scripts. The implementation consumes the precomputed binary64 coefficients, not reconstructed or modified coefficients. No expensive symbolic rebuild was executed for this implementation. Floating-point regression tests are numerical evidence, not new exact proofs.

| Fixture | Functions / original parents | PS12 cells | Binary bytes | SHA-256 |
|---|---:|---:|---:|---|
| fan5 | 100 / 18 | 60 | 1,194,888 | `3831b9c5aa314da4391522ea99e46ff340fe5ddaa00326ad9db39f53884ee68d` |
| fan9 | 172 / 30 | 108 | 3,674,664 | `a3db27d57bff39c4a2e49c268eb50ca91981bc1b4c6247876d3b444a38faad18` |

These are complete positive, normalized, affine-precise, P5-complete fixed-fan level-zero bases, according to the supplied certification. That claim applies only to the two supplied parameter geometries. The restricted historical 33/57-function family is different and is not P5-complete. The viewer uses only the full 100/172 family.

The original column-scaled coefficient condition numbers are about 16,613 and 23,846. Initial parameter sites are heavily duplicated and can lie outside the fan. They are affine-precision metadata, not unique Greville sites, interpolation points or an intuitive sparse cage. The UI groups children by their original parent, filters role and actual coefficient support, and cycles overlapping controls on repeated picks. Group position/weight editing applies a common delta to the selected parent's children; it is not a basis change.

Positive weights do not preserve affine-site interpolation under arbitrary weight changes. They also do not guarantee immersion, good curvature fairness or freedom from self-intersections. Normals and curvature carry a validity flag; invalid curvature is returned as JSON null and displayed in magenta. Completely collapsed geometry appears as a magenta point with an invalid-immersion notice because its triangles have no area.

The parent-to-child homogeneous copying test establishes rational invariance on the **same** PS12 mesh. This is fixed-mesh enrichment, not true spatial subdivision. No arbitrary topology/valence, red spatial-refinement stability, positive multilevel transfer, uniform conditioning, regular extraordinary C² chart join, robust fitting, trimming, booleans, CAD interchange, or production CAD readiness is claimed. An OBJ is a sampled triangle export only.

The optional regular quintic NURBS/OpenSubdiv comparator (handoff milestone 4) is deferred. No fairness, effective-DOF or performance comparison against NURBS is made.
