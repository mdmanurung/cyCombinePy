# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.2] - 2026-07-12

### Added
- Added `CITATION.cff` and citation documentation for citing cyCombinePy and
  the original cyCombine method paper.
- Added shared input validation for marker names, layers, finite marker
  matrices, and missing observation metadata across normalization, correction,
  evaluation, detection, and FCS I/O entry points.
- Added H5AD-safe correction reports, strict correction/confounding policies,
  and propagated report aggregation through `batch_correct`.
- Added deterministic exact-value tests, scientific/reproducibility checks,
  optional FlowSOM/inmoose integration guards, and an R-parity golden-test
  harness with a fixture-generation script.
- Added packaging, docs, CI, release, wheel, sdist, and optional-extras smoke
  checks for publication readiness.

### Changed
- Clarified README and docs scope/provenance wording for the AnnData workflow
  implemented by cyCombinePy.
- Tightened supported Python metadata to Python 3.10-3.12 and added the
  `py.typed` marker to the wheel.
- Stripped notebook execution outputs and documented notebooks as source-only
  examples.

## [0.1.1] - 2026-04-09

### Changed
- Version bump to 0.1.1; first release to PyPI and TestPyPI.

## [0.1.0.dev0] — Initial Python port

### Added
- Python port of the [cyCombine](https://github.com/biosurf/cyCombine) R package
  built natively on `AnnData`.
- Core batch-correction pipeline: `transform_asinh`, `normalize`, `create_som`,
  `correct_data`, and the end-to-end `batch_correct` orchestrator.
- ComBat correction via `inmoose.pycombat.pycombat_norm`.
- SOM clustering via `FlowSOM_Python`.
- FCS I/O via `pytometry` / `readfcs` (`cycombinepy.io.read_fcs_dir`).
- Evaluation utilities: `compute_emd` / `evaluate_emd`,
  `compute_mad` / `evaluate_mad`, and a `scib_metrics` wrapper.
- Batch-effect detection: `detect_batch_effect` and
  `detect_batch_effect_express`.
- Plotting helpers: `plot_density`, `plot_dimred`, `plot_emd_heatmap`.
- Unit tests covering the full pipeline (25 tests).

### Out of scope (not ported)
- Seurat / SingleCellExperiment wrappers
- Panel merging (`impute_across_panels`, `salvage_problematic`)
- `ComBat_seq` RNA-seq variant
- Alternative clustering backends (kohonen / FuseSOM / leiden / kmeans)
- `run_analysis` orchestrator
