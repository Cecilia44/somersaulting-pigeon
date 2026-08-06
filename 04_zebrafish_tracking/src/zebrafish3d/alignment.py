"""Extract manually synchronized frame ranges from the two source videos."""

from __future__ import annotations

from pathlib import Path

from .config import ExperimentConfig


def extract_frame_range(
    video_path: Path,
    output_dir: Path,
    first_frame: int,
    last_frame: int,
    overwrite: bool = False,
) -> list[Path]:
    """Extract a zero-based, inclusive frame range while retaining source indices."""

    if first_frame < 0 or last_frame < first_frame:
        raise ValueError(f"Invalid frame range {first_frame}..{last_frame}")
    if not video_path.is_file():
        raise FileNotFoundError(video_path)

    output_dir.mkdir(parents=True, exist_ok=True)
    expected = [output_dir / f"{index:06d}.jpg" for index in range(first_frame, last_frame + 1)]
    if not overwrite and all(path.is_file() for path in expected):
        return expected

    try:
        import cv2
    except ImportError as exc:
        raise RuntimeError("Frame extraction requires opencv-python") from exc

    capture = cv2.VideoCapture(str(video_path))
    if not capture.isOpened():
        raise RuntimeError(f"Could not open video: {video_path}")
    capture.set(cv2.CAP_PROP_POS_FRAMES, first_frame)
    try:
        for index, path in zip(range(first_frame, last_frame + 1), expected, strict=True):
            ok, frame = capture.read()
            if not ok:
                raise RuntimeError(f"Could not decode frame {index} from {video_path}")
            if overwrite or not path.exists():
                if not cv2.imwrite(str(path), frame, [cv2.IMWRITE_JPEG_QUALITY, 95]):
                    raise RuntimeError(f"Could not write frame: {path}")
    finally:
        capture.release()
    return expected


def extract_aligned_views(config: ExperimentConfig, overwrite: bool = False) -> None:
    counts: dict[str, int] = {}
    for name in ("left", "top"):
        view = config.view(name)
        frames = extract_frame_range(
            video_path=config.resolve(view["video"]),
            output_dir=config.resolve(view["frames_dir"]),
            first_frame=int(view["first_frame"]),
            last_frame=int(view["last_frame"]),
            overwrite=overwrite,
        )
        counts[name] = len(frames)
    if len(set(counts.values())) != 1:
        raise RuntimeError(f"Manual alignment produced unequal frame counts: {counts}")
