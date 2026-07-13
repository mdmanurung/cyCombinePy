"""Public API contract tests for cycombinepy."""

from __future__ import annotations


def test_top_level_all_exports_intentional_core_api():
    import cycombinepy as pc
    from cycombinepy.batch_correct import batch_correct
    from cycombinepy.correct import CombatCorrectionError, ConfoundedDesignError
    from cycombinepy.normalize import normalize

    exported = set(pc.__all__)

    assert {"batch_correct", "normalize"} <= exported
    assert {"CombatCorrectionError", "ConfoundedDesignError"} <= exported
    assert "io" not in exported
    assert "plotting" not in exported

    assert pc.batch_correct is batch_correct
    assert pc.normalize is normalize
    assert pc.CombatCorrectionError is CombatCorrectionError
    assert pc.ConfoundedDesignError is ConfoundedDesignError


def test_submodule_public_imports_work_without_top_level_exports():
    from cycombinepy import CombatCorrectionError, ConfoundedDesignError
    from cycombinepy.correct import CORRECTED_LAYER
    from cycombinepy.io import read_fcs_dir
    from cycombinepy.plotting import plot_density

    assert callable(read_fcs_dir)
    assert callable(plot_density)
    assert CORRECTED_LAYER == "cycombine_corrected"
    assert issubclass(CombatCorrectionError, RuntimeError)
    assert issubclass(ConfoundedDesignError, ValueError)
