"""Publication-oriented plotting without changing the reconstructed samples."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def plot_trajectory(csv_path: Path, output_paths: list[Path], units: str) -> None:
    try:
        import matplotlib
    except ImportError as exc:
        raise RuntimeError("Plotting requires matplotlib") from exc

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    xyz = np.asarray([[float(row[key]) for key in ("x", "y", "z")] for row in rows])
    times = np.asarray([float(row["time_s"]) for row in rows])
    if len(xyz) < 2:
        raise ValueError("At least two trajectory samples are required for plotting")

    points = xyz.reshape(-1, 1, 3)
    segments = np.concatenate((points[:-1], points[1:]), axis=1)
    figure = plt.figure(figsize=(6.2, 5.2), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    norm = matplotlib.colors.Normalize(vmin=float(times.min()), vmax=float(times.max()))
    collection = Line3DCollection(segments, cmap="viridis", norm=norm, linewidth=2.2)
    collection.set_array(times[:-1])
    axis.add_collection3d(collection)
    axis.scatter(*xyz[0], color="black", s=24, label="start", zorder=3)
    axis.scatter(*xyz[-1], color="red", s=24, label="end", zorder=3)

    minimum = np.nanmin(xyz, axis=0)
    maximum = np.nanmax(xyz, axis=0)
    padding = np.maximum((maximum - minimum) * 0.05, 1e-9)
    axis.set_xlim(minimum[0] - padding[0], maximum[0] + padding[0])
    axis.set_ylim(minimum[1] - padding[1], maximum[1] + padding[1])
    axis.set_zlim(minimum[2] - padding[2], maximum[2] + padding[2])
    axis.set_box_aspect(np.maximum(maximum - minimum, 1e-9))
    suffix = f" ({units})" if units else ""
    axis.set_xlabel(f"X{suffix}")
    axis.set_ylabel(f"Y{suffix}")
    axis.set_zlabel(f"Z{suffix}")
    axis.legend(frameon=False)
    colorbar = figure.colorbar(collection, ax=axis, pad=0.12, shrink=0.72)
    colorbar.set_label("Time (s)")

    for output_path in output_paths:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
