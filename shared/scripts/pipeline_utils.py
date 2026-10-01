"""
Shared utilities for the drug repurposing pipeline scripts.

All pipeline scripts (01–10, cmap_analysis) import from here so that
find_repo_root() and load_config() stay in one place.
"""

import subprocess
from pathlib import Path

import yaml


def find_repo_root() -> Path:
    try:
        root = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], stderr=subprocess.DEVNULL
        ).decode().strip()
        return Path(root)
    except subprocess.CalledProcessError:
        return Path(__file__).resolve().parents[2]


def load_config(config_arg: str) -> tuple[dict, Path]:
    """Load a disease YAML config. Returns (config_dict, repo_root)."""
    repo_root = find_repo_root()
    path = Path(config_arg)
    if not path.is_absolute():
        path = repo_root / path
    with open(path) as f:
        return yaml.safe_load(f), repo_root
