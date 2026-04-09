# Installation

## Requirements

`pycombine` supports Python 3.10 and newer. The core pipeline depends on
[`numpy`](https://numpy.org), [`pandas`](https://pandas.pydata.org),
[`scipy`](https://scipy.org), [`scikit-learn`](https://scikit-learn.org),
[`anndata`](https://anndata.readthedocs.io),
[`scanpy`](https://scanpy.readthedocs.io),
[`inmoose`](https://github.com/epigenelabs/inmoose) (ComBat),
[`flowsom`](https://github.com/saeyslab/FlowSOM_Python) (SOM clustering) and
[`formulaic`](https://matthewwardrop.github.io/formulaic/) (model-matrix
construction).

## Install from PyPI

```bash
pip install pycombine
```

This pulls in the core dependencies. The following optional extras enable
additional features:

| Extra         | Adds                                                                 |
| ------------- | -------------------------------------------------------------------- |
| `[io]`        | `pytometry` for reading raw FCS files via `pycombine.io.read_fcs_dir`|
| `[plotting]`  | `matplotlib` + `seaborn` for the `pycombine.plotting` helpers        |
| `[eval]`      | `scib-metrics` for `pycombine.evaluate.scib_metrics`                 |
| `[all]`       | Everything above                                                     |

Install with extras, e.g.:

```bash
pip install "pycombine[all]"
```

## Install from source (development)

Clone the repository and install in editable mode:

```bash
git clone https://github.com/mdmanurung/pyCombine.git
cd pyCombine
pip install -e ".[all,dev]"
```

The `[dev]` extra adds `pytest`, `pytest-cov`, `build`, and `twine` for the
release workflow. After installation you can run the test suite with:

```bash
pytest -q
```

## Building the documentation

The documentation site is built with [Sphinx](https://www.sphinx-doc.org/),
[`myst-nb`](https://myst-nb.readthedocs.io), and the
[`furo`](https://pradyunsg.me/furo/) theme. Install the `docs` extra and
build locally:

```bash
pip install -e ".[all,docs]"
sphinx-build -b html docs/source docs/build/html
```

The generated site lands in `docs/build/html/`. Open `index.html` in a
browser, or serve it locally with:

```bash
python -m http.server -d docs/build/html 8000
```

The tutorial notebooks under `docs/source/notebooks/` are committed with
their cell outputs, so rebuilding the docs does **not** require re-running
them. To refresh the executed outputs:

```bash
jupyter nbconvert --to notebook --execute --inplace \
    docs/source/notebooks/cycombine.ipynb \
    docs/source/notebooks/detect_batch_effects.ipynb
```

## Verifying the installation

```python
import pycombine as pc
print(pc.__version__)

# Functional smoke test: create a tiny AnnData and run the pipeline.
import anndata as ad
import numpy as np
import pandas as pd

rng = np.random.default_rng(0)
X = np.vstack([rng.normal(0, 1, (200, 5)),
               rng.normal(1, 1, (200, 5))])
obs = pd.DataFrame({"batch": ["a"] * 200 + ["b"] * 200})
obs.index = obs.index.astype(str)
adata = ad.AnnData(X=X, obs=obs)
adata.var_names = [f"CD{i}" for i in range(5)]

pc.batch_correct(adata, xdim=4, ydim=4, rlen=3, seed=0)
print("corrected layer:", adata.layers["cycombine_corrected"].shape)
```

If that prints a `(400, 5)` shape without raising, your installation is
working.
