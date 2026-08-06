"""Command-line entry points for the complete reproducible workflow."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from .alignment import extract_aligned_views
from .config import ExperimentConfig, load_config
from .io import load_track, remove_legacy_initial_duplicate, sorted_frame_paths
from .plotting import plot_trajectory
from .reconstruction import reconstruct_from_config
from .tracking import track_view


def _validate(config: ExperimentConfig) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    try:
        expected = config.expected_frame_count
    except ValueError as exc:
        return [str(exc)], warnings

    for name in ("left", "top"):
        view = config.view(name)
        video = config.resolve(view["video"])
        frames_dir = config.resolve(view["frames_dir"])
        if not video.is_file():
            errors.append(f"Missing {name} video: {video}")
        if not frames_dir.is_dir():
            errors.append(f"Missing {name} aligned frame directory: {frames_dir}")
        else:
            frames = sorted_frame_paths(frames_dir)
            expected_indices = list(range(int(view["first_frame"]), int(view["last_frame"]) + 1))
            observed_indices = [int(path.stem) for path in frames]
            if observed_indices != expected_indices:
                errors.append(
                    f"{name} frame indices do not match the configured inclusive range"
                )

        track_candidates = [config.resolve(view["track_npy"])]
        if view.get("legacy_track_npy"):
            track_candidates.append(config.resolve(view["legacy_track_npy"]))
        existing = next((path for path in track_candidates if path.is_file()), None)
        if existing is None:
            warnings.append(f"No {name} trajectory exists yet")
        else:
            try:
                _, repaired = remove_legacy_initial_duplicate(load_track(existing), expected)
                if repaired:
                    warnings.append(f"{existing.name}: detected repairable duplicated first sample")
            except ValueError as exc:
                errors.append(f"{existing}: {exc}")

    checkpoint = config.resolve(config.data["tracking"]["checkpoint"])
    if not checkpoint.is_file():
        warnings.append(f"SAM 2 checkpoint has not been downloaded: {checkpoint}")
    if config.data["reconstruction"]["method"] == "orthographic_fusion":
        warnings.append("Reconstruction is in uncalibrated orthographic-fusion mode")
    return errors, warnings


def _plot_from_config(config: ExperimentConfig) -> None:
    outputs = config.data["outputs"]
    csv_path = config.resolve(outputs["trajectory_csv"])
    figure_paths = [config.resolve(value) for value in outputs["figures"]]
    units = str(config.data["reconstruction"].get("units", ""))
    plot_trajectory(csv_path, figure_paths, units)
    for path in figure_paths:
        print(f"Wrote {path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="zebrafish3d",
        description="Track a zebrafish in synchronized top/left videos and reconstruct a 3D path.",
    )
    parser.add_argument(
        "--config",
        default="configs/2025-03-19_17-22.json",
        help="Experiment JSON configuration (default: %(default)s)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract", help="Extract configured manual frame ranges")
    extract.add_argument("--overwrite", action="store_true", help="Rewrite existing JPG frames")

    track = subparsers.add_parser("track", help="Run SAM 2 tracking")
    track.add_argument("--view", choices=("left", "top", "all"), default="all")
    track.add_argument(
        "--point", nargs=2, type=float, metavar=("X", "Y"), help="Override one positive prompt"
    )

    subparsers.add_parser("reconstruct", help="Create the 3D trajectory CSV and NPY")
    subparsers.add_parser("plot", help="Plot an existing reconstructed trajectory")
    subparsers.add_parser("validate", help="Check data, frame ranges, tracks, and checkpoint")
    subparsers.add_parser("run", help="Run extract, tracking, reconstruction, and plotting")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        config = load_config(args.config)
        if args.command == "extract":
            extract_aligned_views(config, overwrite=args.overwrite)
            print(f"Validated/extracted {config.expected_frame_count} paired frames")
        elif args.command == "track":
            views = ("left", "top") if args.view == "all" else (args.view,)
            for view in views:
                npy_path, csv_path = track_view(config, view, args.point)
                print(f"Wrote {npy_path}")
                print(f"Wrote {csv_path}")
        elif args.command == "reconstruct":
            csv_path, npy_path, warnings = reconstruct_from_config(config)
            for warning in warnings:
                print(f"WARNING: {warning}", file=sys.stderr)
            print(f"Wrote {csv_path}")
            print(f"Wrote {npy_path}")
        elif args.command == "plot":
            _plot_from_config(config)
        elif args.command == "validate":
            errors, warnings = _validate(config)
            for warning in warnings:
                print(f"WARNING: {warning}")
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            if errors:
                return 1
            print(f"Validation passed: {config.expected_frame_count} paired frames")
        elif args.command == "run":
            extract_aligned_views(config)
            for view in ("left", "top"):
                track_view(config, view)
            _, _, warnings = reconstruct_from_config(config)
            for warning in warnings:
                print(f"WARNING: {warning}", file=sys.stderr)
            _plot_from_config(config)
        return 0
    except (FileNotFoundError, KeyError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

