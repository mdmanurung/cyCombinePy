from __future__ import annotations

from pathlib import Path

from cycombinepy._skills import install


def test_bundled_skill_contains_router_and_reference():
    skill_dir = install.bundled_skill_dir()

    assert (skill_dir / "SKILL.md").is_file()
    assert (skill_dir / "references" / "api_workflows.md").is_file()


def test_install_skill_copies_bundle_to_custom_destination(tmp_path):
    dest = tmp_path / "skills" / "cycombinepy"

    installed = install.install_skill("claude", dest=dest)

    assert installed == dest
    assert (dest / "SKILL.md").is_file()
    assert (dest / "references" / "api_workflows.md").is_file()


def test_install_skill_refuses_existing_destination_without_force(tmp_path):
    dest = tmp_path / "skills" / "cycombinepy"
    dest.mkdir(parents=True)

    try:
        install.install_skill("codex", dest=dest)
    except FileExistsError as exc:
        assert "--force" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected FileExistsError")


def test_main_print_path(capsys):
    code = install.main(["--print-path"])

    assert code == 0
    printed = Path(capsys.readouterr().out.strip())
    assert printed == install.bundled_skill_dir()

