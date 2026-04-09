"""cyCombinePy: Python port of cyCombine for batch correction of cytometry data."""

from cycombinepy._utils import check_confound, get_markers
from cycombinepy.batch_correct import batch_correct
from cycombinepy.cluster import create_som
from cycombinepy.correct import correct_data
from cycombinepy.detect import detect_batch_effect, detect_batch_effect_express
from cycombinepy.evaluate import (
    compute_emd,
    compute_mad,
    evaluate_emd,
    evaluate_mad,
)
from cycombinepy.normalize import normalize
from cycombinepy.preprocessing import transform_asinh

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "batch_correct",
    "check_confound",
    "compute_emd",
    "compute_mad",
    "correct_data",
    "create_som",
    "detect_batch_effect",
    "detect_batch_effect_express",
    "evaluate_emd",
    "evaluate_mad",
    "get_markers",
    "normalize",
    "transform_asinh",
]
