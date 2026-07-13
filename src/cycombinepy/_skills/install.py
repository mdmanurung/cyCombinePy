"""Install the bundled cyCombinePy agent skill for Claude Code or Codex."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path


SKILL_NAME = "cycombinepy"
AGENTS = ("claude", "codex")


def bundled_skill_dir() -> Path:
    """Return the bundled skill directory inside the installed package."""
    return Path(__file__).resolve().parent / "data"


def default_dest(agent: str) -> Path:
    """Return the default personal skill directory for ``agent``."""
    if agent == "claude":
        return Path.home() / ".claude" / "skills" / SKILL_NAME
    if agent == "codex":
        codex_home = Path(
            # Respect CODEX_HOME when users keep Codex state outside ~/.codex.
            os.environ.get("CODEX_HOME", str(Path.home() / ".codex"))
        )
        return codex_home / "skills" / SKILL_NAME
    raise ValueError(f"Unknown agent: {agent!r}")


def install_skill(agent: str, dest: Path | None = None, force: bool = False) -> Path:
    """Copy the bundled skill into an agent's personal skill directory."""
    if agent not in AGENTS:
        raise ValueError(f"agent must be one of {AGENTS!r}; got {agent!r}")

    src = bundled_skill_dir()
    if not (src / "SKILL.md").is_file():
        raise FileNotFoundError(
            f"Bundled skill not found at {src}. The package may be missing "
            "its skill data."
        )

    if dest is None:
        dest = default_dest(agent)

    if dest.exists():
        if not force:
            raise FileExistsError(
                f"{dest} already exists. Re-run with --force to overwrite."
            )
        shutil.rmtree(dest)

    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(
        src,
        dest,
        ignore=shutil.ignore_patterns(".ipynb_checkpoints", "__pycache__"),
    )
    return dest


def _agents_from_arg(value: str) -> list[str]:
    if value == "all":
        return list(AGENTS)
    if value not in AGENTS:
        raise argparse.ArgumentTypeError(
            f"agent must be one of {', '.join((*AGENTS, 'all'))}"
        )
    return [value]


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point for ``cycombinepy-install-skills``."""
    parser = argparse.ArgumentParser(
        prog="cycombinepy-install-skills",
        description="Install the bundled cyCombinePy agent skill.",
    )
    parser.add_argument(
        "--agent",
        default="all",
        choices=(*AGENTS, "all"),
        help="Agent skill root to install into (default: all).",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=None,
        help=(
            "Destination directory. Only valid with --agent claude or "
            "--agent codex."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing skill installation.",
    )
    parser.add_argument(
        "--print-path",
        action="store_true",
        help="Print the bundled skill directory and exit without installing.",
    )
    args = parser.parse_args(argv)

    if args.print_path:
        print(bundled_skill_dir())
        return 0

    agents = _agents_from_arg(args.agent)
    if args.dest is not None and len(agents) > 1:
        print("error: --dest can only be used with one --agent", file=sys.stderr)
        return 2

    try:
        installed = [
            (agent, install_skill(agent, dest=args.dest, force=args.force))
            for agent in agents
        ]
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    for agent, dest in installed:
        print(f"Installed cyCombinePy skill for {agent} to {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
