"""pycombine: Python port of cyCombine for batch correction of cytometry data."""

from pycombine._utils import check_confound, get_markers
from pycombine.batch_correct import batch_correct
from pycombine.cluster import create_som
from pycombine.correct import correct_data
from pycombine.detect import detect_batch_effect, detect_batch_effect_express
from pycombine.evaluate import (
    compute_emd,
    compute_mad,
    evaluate_emd,
    evaluate_mad,
)
from pycombine.normalize import normalize
from pycombine.preprocessing import transform_asinh

__version__ = "0.1.0.dev0"

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
