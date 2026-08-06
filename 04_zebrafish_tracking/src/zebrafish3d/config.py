"""Experiment configuration loading and path resolution."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ExperimentConfig:
    """A JSON configuration whose paths are relative to ``project_root``."""

    path: Path
    project_root: Path
    data: dict[str, Any]

    def resolve(self, value: str | Path) -> Path:
        path = Path(value)
        return path if path.is_absolute() else self.project_root / path

    def view(self, name: str) -> dict[str, Any]:
        try:
            return self.data["views"][name]
        except KeyError as exc:
            raise KeyError(f"View {name!r} is not defined in {self.path}") from exc

    @property
    def expected_frame_count(self) -> int:
        counts = {
            name: int(view["last_frame"]) - int(view["first_frame"]) + 1
            for name, view in self.data["views"].items()
        }
        if len(set(counts.values())) != 1:
            raise ValueError(f"Aligned views have different frame counts: {counts}")
        count = next(iter(counts.values()))
        if count <= 0:
            raise ValueError(f"Frame ranges must be non-empty: {counts}")
        return count


def load_config(path: str | Path) -> ExperimentConfig:
    config_path = Path(path).expanduser().resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    for key in ("experiment_id", "views", "tracking", "reconstruction", "outputs"):
        if key not in data:
            raise ValueError(f"Missing required configuration key: {key}")
    if set(data["views"]) != {"left", "top"}:
        raise ValueError("The configuration must define exactly 'left' and 'top' views")

    root_value = data.get("project_root", ".")
    root = (config_path.parent / root_value).resolve()
    return ExperimentConfig(path=config_path, project_root=root, data=data)

