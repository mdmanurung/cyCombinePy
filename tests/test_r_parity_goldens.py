from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pytest


pytestmark = pytest.mark.parity

GOLDEN_ROOT = Path(__file__).parent / "golden"
MANIFESTS = sorted(GOLDEN_ROOT.glob("cycombine_r_*/manifest.json"))

EXPECTED_FIXTURE_FILES = {
    "input.csv",
    "normalize_scale.npz",
    "normalize_rank.npz",
    "emd.csv",
    "mad.csv",
    "corrected_fixed_labels.npz",
}
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def test_r_golden_manifest_schema_and_hashes() -> None:
    if not MANIFESTS:
        pytest.skip("R cyCombine golden fixtures are not generated yet")

    for manifest_path in MANIFESTS:
        manifest = _load_manifest(manifest_path)

        assert isinstance(manifest.get("schema_version"), int), manifest_path
        assert _non_empty_string(manifest.get("r_version")), manifest_path
        assert isinstance(manifest.get("package_versions"), dict), manifest_path
        assert manifest["package_versions"], manifest_path
        assert isinstance(manifest.get("random_seed"), int), manifest_path
        assert isinstance(manifest.get("provenance"), dict), manifest_path
        assert manifest["provenance"], manifest_path

        cycombine = manifest.get("cycombine")
        assert isinstance(cycombine, dict), manifest_path
        assert _non_empty_string(cycombine.get("version")), manifest_path
        assert _non_empty_string(cycombine.get("commit")), manifest_path

        hashes = manifest.get("sha256")
        assert isinstance(hashes, dict), manifest_path
        assert set(hashes) == EXPECTED_FIXTURE_FILES, manifest_path

        for fixture_name, expected_hash in hashes.items():
            assert SHA256_RE.fullmatch(expected_hash), fixture_name

            fixture_path = manifest_path.parent / fixture_name
            assert fixture_path.is_file(), fixture_name
            assert _sha256(fixture_path) == expected_hash


def _load_manifest(manifest_path: Path) -> dict[str, object]:
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    assert isinstance(manifest, dict)
    return manifest


def _non_empty_string(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
