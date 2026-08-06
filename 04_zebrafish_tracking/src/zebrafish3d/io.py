"""Small, dependency-light helpers for frames and trajectory files."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

import numpy as np


IMAGE_SUFFIXES = {".jpg", ".jpeg", ".png"}


def sorted_frame_paths(directory: str | Path) -> list[Path]:
    directory = Path(directory)
    frames = [p for p in directory.iterdir() if p.suffix.lower() in IMAGE_SUFFIXES]
    try:
        return sorted(frames, key=lambda p: int(p.stem))
    except ValueError as exc:
        raise ValueError(
            f"Frame names in {directory} must have numeric stems, e.g. 000072.jpg"
        ) from exc


def mask_centroid(mask: np.ndarray) -> tuple[float, float, int]:
    """Return centroid ``(x, y)`` and area for a binary mask."""

    binary = np.asarray(mask).squeeze().astype(bool)
    if binary.ndim != 2:
        raise ValueError(f"Expected a 2D mask after squeeze, got shape {binary.shape}")
    rows, columns = np.nonzero(binary)
    if columns.size == 0:
        return float("nan"), float("nan"), 0
    return float(columns.mean()), float(rows.mean()), int(columns.size)


def load_track(path: str | Path) -> np.ndarray:
    """Load an ``(N, 2)`` center track from NPY or CSV."""

    path = Path(path)
    if path.suffix.lower() == ".npy":
        track = np.load(path, allow_pickle=False)
    elif path.suffix.lower() == ".csv":
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        track = np.asarray(
            [[float(row["x_px"]), float(row["y_px"])] for row in rows], dtype=float
        )
    else:
        raise ValueError(f"Unsupported trajectory format: {path.suffix}")

    track = np.asarray(track, dtype=float)
    if track.ndim != 2 or track.shape[1] != 2:
        raise ValueError(f"Expected an (N, 2) trajectory in {path}, got {track.shape}")
    track[np.all(track == -1, axis=1)] = np.nan
    return track


def remove_legacy_initial_duplicate(
    track: np.ndarray, expected_length: int
) -> tuple[np.ndarray, bool]:
    """Repair the duplicated frame zero produced by the original tracking script."""

    track = np.asarray(track, dtype=float)
    if len(track) == expected_length:
        return track, False
    if (
        len(track) == expected_length + 1
        and len(track) >= 2
        and np.allclose(track[0], track[1], equal_nan=True)
    ):
        return track[1:], True
    raise ValueError(
        f"Trajectory has {len(track)} samples; expected {expected_length}. "
        "Only the known duplicated-first-frame legacy case can be repaired automatically."
    )


def interpolate_missing(track: np.ndarray) -> np.ndarray:
    """Linearly fill NaNs, extending the nearest valid value at each edge."""

    result = np.asarray(track, dtype=float).copy()
    indices = np.arange(len(result))
    for column in range(result.shape[1]):
        valid = np.isfinite(result[:, column])
        if not valid.any():
            raise ValueError(f"Trajectory coordinate column {column} has no valid samples")
        result[:, column] = np.interp(
            indices, indices[valid], result[valid, column]
        )
    return result


def save_track(
    centers: np.ndarray,
    areas: np.ndarray,
    frame_paths: Iterable[Path],
    npy_path: str | Path,
    csv_path: str | Path,
) -> None:
    npy_path = Path(npy_path)
    csv_path = Path(csv_path)
    npy_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    np.save(npy_path, np.asarray(centers, dtype=float))

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "local_frame_index",
            "source_frame_index",
            "source_filename",
            "x_px",
            "y_px",
            "mask_area_px",
            "is_valid",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for local_index, (frame, center, area) in enumerate(
            zip(frame_paths, centers, areas, strict=True)
        ):
            valid = bool(np.all(np.isfinite(center)))
            writer.writerow(
                {
                    "local_frame_index": local_index,
                    "source_frame_index": int(frame.stem),
                    "source_filename": frame.name,
                    "x_px": f"{center[0]:.8f}" if valid else "nan",
                    "y_px": f"{center[1]:.8f}" if valid else "nan",
                    "mask_area_px": int(area),
                    "is_valid": int(valid),
                }
            )
