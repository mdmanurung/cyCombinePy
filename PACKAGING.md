# Packaging and release

cycombinepy is packaged with [hatchling](https://hatch.pypa.io/) and published to
PyPI via GitHub Actions using [PyPI trusted publishing](https://docs.pypi.org/trusted-publishers/)
(OIDC, no long-lived API tokens).

Published artifacts support Python 3.10, 3.11, and 3.12.

## Local build

```bash
pip install build twine
python -m build
twine check dist/*
```

The build produces two artifacts in `dist/`:

- `cycombinepy-<version>.tar.gz`: source distribution (includes `src/`, `tests/`, `scripts/`,
  `inst/extdata/test.fcs`, `CITATION.cff`, `README.md`, `LICENSE`, `CHANGELOG.md`)
- `cycombinepy-<version>-py3-none-any.whl`: pure-Python wheel (includes only the
  `cycombinepy` package, the `cycombinepy/py.typed` marker, and license metadata)

Artifact checks should confirm that the wheel contains `cycombinepy/py.typed`
and that the sdist contains both `inst/extdata/test.fcs` and `CITATION.cff`.

The version is single-sourced from `src/cycombinepy/__init__.py` via
`[tool.hatch.version]`. Bump it there, not in `pyproject.toml`.

## Publication gates

Publication readiness is covered across the CI, release, and docs workflows:

- editable-install tests on Python 3.10, 3.11, and 3.12
- `python -m build` for the source distribution and wheel
- `twine check dist/*`
- wheel install smoke tests on Python 3.10, 3.11, and 3.12
- sdist install with `.[dev]` followed by `pytest -q -rs` on Python 3.12
- wheel install with `[all]` extras and optional dependency imports on Python 3.12
- documentation snippet tests and `sphinx-build -W` in the docs workflow

The sdist smoke intentionally runs the full test suite from the unpacked source
archive. That gate can be slower than the wheel smoke because it exercises the
release artifact as an install-and-test source tree.

## One-time GitHub / PyPI setup

### 1. Configure trusted publishers

Create a pending publisher on both sites pointing at this repo:

- **TestPyPI**: <https://test.pypi.org/manage/account/publishing/>
- **PyPI**:    <https://pypi.org/manage/account/publishing/>

Fields:

| Field | Value |
|---|---|
| PyPI Project Name | `cycombinepy` |
| Owner | `mdmanurung` |
| Repository name | `cyCombinePy` |
| Workflow name | `release.yml` |
| Environment name | `testpypi` (TestPyPI) / `pypi` (PyPI) |

### 2. Create matching GitHub environments

In the repo settings, open Environments and create two environments:

- `testpypi`
- `pypi`

Optionally protect the `pypi` environment with required reviewers so a tag push
blocks on manual approval before the real PyPI upload.

## Release flow

### Pre-flight

1. Ensure `main` is green on CI.
2. Bump the version in `src/cycombinepy/__init__.py` (e.g. `0.1.0.dev0` to
   `0.1.0`).
3. Update `CHANGELOG.md`: move items from `## [Unreleased]` into a new
   `## [<version>] - YYYY-MM-DD` section.
4. Check that `CITATION.cff` has the released version/date and is included in
   the sdist.
5. Commit: `git commit -am "Release v<version>"`.

### Dry run on TestPyPI

Two options:

**Option A: manual dispatch (no tag):**

```text
GitHub, Actions, Release, Run workflow, target: testpypi
```

Workflow builds the release artifacts, runs the artifact smoke gates, and uploads
to TestPyPI only. Manual dispatch cannot publish to real PyPI.

**Option B: pre-release tag:**

```bash
git tag v0.1.0rc1
git push origin v0.1.0rc1
```

A tag push always publishes to TestPyPI. It *also* attempts PyPI, so use only
release-candidate tags (`rcN`, `aN`, `bN`) that haven't been claimed on PyPI.

Verify install from TestPyPI:

```bash
python -m venv /tmp/pcv && source /tmp/pcv/bin/activate
pip install --index-url https://test.pypi.org/simple/ \
            --extra-index-url https://pypi.org/simple/ \
            cycombinepy==<version>
python -c "import cycombinepy; print(cycombinepy.__version__)"
```

### Real release

```bash
git tag v0.1.0
git push origin v0.1.0
```

The `release.yml` workflow:

1. Builds sdist + wheel
2. Runs `twine check`
3. Runs wheel, sdist, and extras smoke gates
4. Publishes to TestPyPI
5. Publishes to PyPI (gated by the `pypi` environment if you added a protection
   rule)

After the workflow finishes, create a GitHub Release from the tag with notes
pulled from `CHANGELOG.md`.

### Post-release

1. Bump `__version__` on `main` back to the next dev version (e.g.
   `0.1.3.dev0` after releasing `0.1.2`).
2. Commit and push.

## Troubleshooting

- **`twine check` fails on long description**: the README must be valid
  CommonMark; run `python -m readme_renderer README.md` to inspect.
- **`cycombinepy` name already taken on PyPI**: rename the project in
  `pyproject.toml` (e.g. `cycombinepy-cyto`) and update the trusted publisher.
- **OIDC publish fails with `not a trusted publisher`**: the workflow filename,
  environment name, or repo owner in the pending publisher does not match.
  Re-create the publisher with the exact values shown in the action log.
