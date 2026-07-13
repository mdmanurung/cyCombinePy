import json
from pathlib import Path


def _notebook_source(path: str) -> str:
    notebook = json.loads(Path(path).read_text())
    parts = []
    for cell in notebook["cells"]:
        source = cell.get("source", "")
        if isinstance(source, list):
            source = "".join(source)
        parts.append(source)
    return "\n".join(parts)


def _notebook_output_text(path: str) -> str:
    notebook = json.loads(Path(path).read_text())
    parts = []
    for cell in notebook["cells"]:
        for output in cell.get("outputs", []):
            if "text" in output:
                text = output["text"]
                parts.append("".join(text) if isinstance(text, list) else text)
            data = output.get("data", {})
            for mime_type, value in data.items():
                if mime_type.startswith("image/"):
                    continue
                parts.append("".join(value) if isinstance(value, list) else str(value))
    return "\n".join(parts)


def test_main_vignette_uses_current_report_and_layer_api():
    source = _notebook_source("docs/source/notebooks/cycombine.ipynb")

    assert "return_report=True" in source
    assert "cycombinepy_correction" in source
    assert "error_policy='raise'" in source
    assert "confound_policy='raise'" in source
    assert "NORMALIZED_LAYER = 'cycombine_normalized'" in source
    assert "layer=NORMALIZED_LAYER" in source
    assert "modular.X = adata.X.copy()" not in source


def test_detect_vignette_documents_current_validation_behavior():
    source = _notebook_source("docs/source/notebooks/detect_batch_effects.ipynb")

    assert "validate requested marker names" in source
    assert "`sample_key` is optional" in source
    assert "returned `mds` figure" in source


def test_citeseq_adt_vignette_documents_benchmark_workflow():
    source = _notebook_source("docs/source/notebooks/citeseq_adt_batch_correction.ipynb")

    assert "CITE-seq_pbmc_combined_preprocessed.h5mu" in source
    assert 'mdata.mod["prot"]' in source
    assert "cycombinepy.correct_data" in source
    assert "harmonypy.run_harmony" in source
    assert "scvi.model.TOTALVI.setup_mudata" in source
    assert "Harmony is therefore included only in the embedding benchmark" in source
    assert "totalVI log1p denoised" in source
    assert "benchmark_table.merge" not in source
    assert "display(expression_table)" in source
    assert "display(embedding_table)" in source


def test_citeseq_adt_vignette_uses_scib_metrics_and_has_rendered_outputs():
    path = "docs/source/notebooks/citeseq_adt_batch_correction.ipynb"
    notebook = json.loads(Path(path).read_text())
    source = _notebook_source(path)
    outputs = _notebook_output_text(path)

    assert "scib_metrics" in source
    assert "ilisi_knn" in source
    assert "silhouette_batch" in source
    assert any(cell.get("outputs") for cell in notebook["cells"] if cell["cell_type"] == "code")
    assert any(
        "image/png" in output.get("data", {})
        for cell in notebook["cells"]
        for output in cell.get("outputs", [])
    )
    assert "scvi-tools MuData:" in outputs
    assert "synthetic fallback" not in outputs
    assert "Harmony: ok" in outputs
    assert "totalVI: ok" in outputs
    assert "totalVI log1p denoised" in outputs
    assert "NaN" not in outputs


def test_vignettes_use_plain_academic_prose_without_em_dashes():
    em_dash = chr(0x2014)
    for path in [
        "docs/source/notebooks/cycombine.ipynb",
        "docs/source/notebooks/detect_batch_effects.ipynb",
        "docs/source/notebooks/citeseq_adt_batch_correction.ipynb",
    ]:
        assert em_dash not in _notebook_source(path)
