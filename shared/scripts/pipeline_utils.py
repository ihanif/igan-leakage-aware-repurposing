"""
Shared utilities for the drug repurposing pipeline scripts.

All pipeline scripts import from here so that find_repo_root() and
load_config() stay in one place.
"""

from pathlib import Path

import yaml


def find_repo_root() -> Path:
    """Repository root, located from this file (shared/scripts/pipeline_utils.py).

    Does not depend on git, so it works from an archive download or when
    unpacked inside another git working tree.
    """
    return Path(__file__).resolve().parents[2]


def load_config(config_arg: str) -> tuple[dict, Path]:
    """Load a disease YAML config. Returns (config_dict, repo_root)."""
    repo_root = find_repo_root()
    path = Path(config_arg)
    if not path.is_absolute():
        path = repo_root / path
    with open(path) as f:
        return yaml.safe_load(f), repo_root
