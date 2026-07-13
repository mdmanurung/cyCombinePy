---
name: cycombinepy
description: Use when working with cyCombinePy for single-cell cytometry batch correction, including AnnData setup, FCS loading, arcsinh transformation, batch-effect diagnostics, FlowSOM clustering, per-cluster ComBat correction, correction reports, and EMD/MAD evaluation.
---

# cyCombinePy

Use this skill when a task involves the `cycombinepy` Python package or asks for
batch correction of single-cell cytometry data in AnnData.

## Workflow

1. Inspect the package version and available extras before assuming optional
   modules are installed.
2. Keep `adata.X` as the post-asinh, uncorrected matrix unless the user asks for
   a different convention.
3. Put normalized values used for clustering in
   `adata.layers["cycombine_normalized"]`.
4. Run correction with strict defaults for scientific work:
   `error_policy="raise"`, `confound_policy="raise"`, and
   `return_report=True`.
5. Inspect `adata.uns["cycombinepy_correction"]` before interpreting corrected
   results.
6. Evaluate batch effects before and after correction with EMD/MAD summaries and
   plots.

## Main APIs

- `transform_asinh`: arcsinh transform marker channels.
- `normalize`: batch-wise normalization for clustering.
- `create_som`: FlowSOM clustering.
- `correct_data`: per-cluster ComBat when cluster labels already exist.
- `batch_correct`: one-call normalization, clustering, and correction.
- `detect_batch_effect_express` and `detect_batch_effect`: pre-correction
  diagnostics.
- `compute_emd`, `evaluate_emd`, `compute_mad`, `evaluate_mad`: quantitative
  evaluation.
- `cycombinepy.io.read_fcs_dir`: FCS loading when the `io` extra is installed.
- `cycombinepy.plotting`: density, UMAP/PCA, and EMD heatmap helpers.

## Rules

- Validate that `adata.obs["batch"]` exists and has no missing values.
- Treat `sample`, `condition`, and `anchor` columns as optional metadata.
- Use `covar=` only for biological covariates that should be preserved.
- Use `anchor=` only for reference samples represented across batches.
- Do not overwrite `adata.X` with corrected values; use
  `adata.layers["cycombine_corrected"]`.
- If correction fails or covariates are confounded, report the exception and the
  correction report rather than silently relaxing policies.

For examples and edge cases, read `references/api_workflows.md`.

