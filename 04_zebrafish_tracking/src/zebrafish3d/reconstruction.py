"""Fuse orthogonal image coordinates or triangulate calibrated camera rays."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

from .config import ExperimentConfig
from .io import interpolate_missing, load_track, remove_legacy_initial_duplicate


def orthographic_fusion(
    left_yz: np.ndarray,
    top_xy: np.ndarray,
    origin_px: np.ndarray,
    scale_units_per_px: np.ndarray,
    axis_sign: np.ndarray,
) -> np.ndarray:
    """Compatibility reconstruction: ``(X,Y,Z)=(top_x,top_y,left_y)``."""

    image_coordinates = np.column_stack((top_xy[:, 0], top_xy[:, 1], left_yz[:, 1]))
    return (image_coordinates - origin_px) * scale_units_per_px * axis_sign


def triangulate_dlt(
    top_points: np.ndarray,
    left_points: np.ndarray,
    top_projection: np.ndarray,
    left_projection: np.ndarray,
) -> np.ndarray:
    """Linear DLT triangulation from two ``3 x 4`` projection matrices."""

    top_projection = np.asarray(top_projection, dtype=float)
    left_projection = np.asarray(left_projection, dtype=float)
    if top_projection.shape != (3, 4) or left_projection.shape != (3, 4):
        raise ValueError("Each projection matrix must have shape (3, 4)")

    world_points = np.full((len(top_points), 3), np.nan, dtype=float)
    for index, ((top_u, top_v), (left_u, left_v)) in enumerate(
        zip(top_points, left_points, strict=True)
    ):
        if not np.all(np.isfinite([top_u, top_v, left_u, left_v])):
            continue
        system = np.vstack(
            (
                top_u * top_projection[2] - top_projection[0],
                top_v * top_projection[2] - top_projection[1],
                left_u * left_projection[2] - left_projection[0],
                left_v * left_projection[2] - left_projection[1],
            )
        )
        _, _, right_vectors = np.linalg.svd(system)
        homogeneous = right_vectors[-1]
        if np.isclose(homogeneous[3], 0.0):
            continue
        world_points[index] = homogeneous[:3] / homogeneous[3]
    return world_points


def reprojection_error(
    world_points: np.ndarray, image_points: np.ndarray, projection: np.ndarray
) -> np.ndarray:
    homogeneous = np.column_stack((world_points, np.ones(len(world_points))))
    projected = homogeneous @ np.asarray(projection, dtype=float).T
    predicted = projected[:, :2] / projected[:, 2:3]
    return np.linalg.norm(predicted - image_points, axis=1)


def _select_track_path(config: ExperimentConfig, view: dict[str, Any]) -> Path:
    current = config.resolve(view["track_npy"])
    if current.is_file():
        return current
    legacy_value = view.get("legacy_track_npy")
    if legacy_value:
        legacy = config.resolve(legacy_value)
        if legacy.is_file():
            return legacy
    raise FileNotFoundError(current)


def reconstruct_from_config(config: ExperimentConfig) -> tuple[Path, Path, list[str]]:
    expected = config.expected_frame_count
    warnings: list[str] = []
    tracks: dict[str, np.ndarray] = {}
    for view_name in ("left", "top"):
        path = _select_track_path(config, config.view(view_name))
        track, repaired = remove_legacy_initial_duplicate(load_track(path), expected)
        if repaired:
            warnings.append(f"Removed duplicated initial sample from legacy track: {path.name}")
        tracks[view_name] = track

    raw_left = tracks["left"].copy()
    raw_top = tracks["top"].copy()
    reconstruction = config.data["reconstruction"]
    if bool(reconstruction.get("interpolate_missing", True)):
        tracks = {name: interpolate_missing(track) for name, track in tracks.items()}

    method = reconstruction["method"]
    top_error = np.full(expected, np.nan)
    left_error = np.full(expected, np.nan)
    if method == "orthographic_fusion":
        xyz = orthographic_fusion(
            tracks["left"],
            tracks["top"],
            origin_px=np.asarray(reconstruction.get("origin_px", [0, 0, 0]), dtype=float),
            scale_units_per_px=np.asarray(
                reconstruction.get("scale_units_per_px", [1, 1, 1]), dtype=float
            ),
            axis_sign=np.asarray(reconstruction.get("axis_sign", [1, 1, 1]), dtype=float),
        )
        warnings.append(
            "Used orthographic coordinate fusion, not calibrated triangulation; "
            "the side-view horizontal coordinate is not used."
        )
    elif method == "dlt_triangulation":
        top_projection = np.asarray(reconstruction["top_projection_matrix"], dtype=float)
        left_projection = np.asarray(reconstruction["left_projection_matrix"], dtype=float)
        xyz = triangulate_dlt(
            tracks["top"], tracks["left"], top_projection, left_projection
        )
        top_error = reprojection_error(xyz, tracks["top"], top_projection)
        left_error = reprojection_error(xyz, tracks["left"], left_projection)
    else:
        raise ValueError(f"Unknown reconstruction method: {method}")

    output_csv = config.resolve(config.data["outputs"]["trajectory_csv"])
    output_npy = config.resolve(config.data["outputs"]["trajectory_npy"])
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_npy.parent.mkdir(parents=True, exist_ok=True)
    np.save(output_npy, xyz)

    paired_fps = float(reconstruction["paired_fps"])
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "frame_index",
            "time_s",
            "top_x_px",
            "top_y_px",
            "left_x_px",
            "left_y_px",
            "top_was_interpolated",
            "left_was_interpolated",
            "x",
            "y",
            "z",
            "top_reprojection_error_px",
            "left_reprojection_error_px",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index in range(expected):
            writer.writerow(
                {
                    "frame_index": index,
                    "time_s": f"{index / paired_fps:.8f}",
                    "top_x_px": f"{tracks['top'][index, 0]:.8f}",
                    "top_y_px": f"{tracks['top'][index, 1]:.8f}",
                    "left_x_px": f"{tracks['left'][index, 0]:.8f}",
                    "left_y_px": f"{tracks['left'][index, 1]:.8f}",
                    "top_was_interpolated": int(not np.all(np.isfinite(raw_top[index]))),
                    "left_was_interpolated": int(not np.all(np.isfinite(raw_left[index]))),
                    "x": f"{xyz[index, 0]:.8f}",
                    "y": f"{xyz[index, 1]:.8f}",
                    "z": f"{xyz[index, 2]:.8f}",
                    "top_reprojection_error_px": (
                        f"{top_error[index]:.8f}" if np.isfinite(top_error[index]) else "nan"
                    ),
                    "left_reprojection_error_px": (
                        f"{left_error[index]:.8f}" if np.isfinite(left_error[index]) else "nan"
                    ),
                }
            )
    return output_csv, output_npy, warnings

