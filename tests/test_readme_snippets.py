"""Executable README snippet tests."""

from __future__ import annotations

import re
from pathlib import Path


README = Path(__file__).resolve().parents[1] / "README.md"


def _extract_python_snippet(name: str) -> str:
    text = README.read_text(encoding="utf-8")
    start = f"<!-- cycombinepy-snippet:start {name} -->"
    end = f"<!-- cycombinepy-snippet:end {name} -->"

    assert start in text, f"README snippet start marker missing: {start}"
    assert end in text, f"README snippet end marker missing: {end}"

    marked = text.split(start, 1)[1].split(end, 1)[0]
    match = re.search(r"```python\n(?P<code>.*?)\n```", marked, re.DOTALL)
    assert match is not None, f"README snippet {name!r} must contain a Python block"
    return match.group("code")


def test_quickstart_synthetic_readme_snippet_executes():
    code = _extract_python_snippet("quickstart-synthetic")
    namespace = {"__name__": "__cycombinepy_readme_snippet__"}

    exec(compile(code, "README.md#quickstart-synthetic", "exec"), namespace)
