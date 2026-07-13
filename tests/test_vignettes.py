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


def test_vignettes_use_plain_academic_prose_without_em_dashes():
    em_dash = chr(0x2014)
    for path in [
        "docs/source/notebooks/cycombine.ipynb",
        "docs/source/notebooks/detect_batch_effects.ipynb",
    ]:
        assert em_dash not in _notebook_source(path)
