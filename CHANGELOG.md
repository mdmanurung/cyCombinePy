# Changelog

All notable changes to this project are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
