"""SAM 2 video segmentation and mask-centroid tracking."""

from __future__ import annotations

from contextlib import nullcontext
from pathlib import Path
from typing import Sequence

import numpy as np

from .config import ExperimentConfig
from .io import mask_centroid, save_track, sorted_frame_paths


def _interactive_point(frame_path: Path) -> list[list[float]]:
    import matplotlib.pyplot as plt
    from PIL import Image

    figure, axis = plt.subplots()
    axis.imshow(Image.open(frame_path))
    axis.set_title("Click the zebrafish once, then close the window")
    selected = plt.ginput(1, timeout=-1)
    plt.close(figure)
    if len(selected) != 1:
        raise RuntimeError("No prompt point was selected")
    return [[float(selected[0][0]), float(selected[0][1])]]


def _store_prediction(
    frame_index: int,
    object_ids: Sequence[int],
    mask_logits: Sequence,
    target_object_id: int,
    centers: np.ndarray,
    areas: np.ndarray,
) -> None:
    ids = [int(item) for item in object_ids]
    if target_object_id not in ids:
        return
    item_index = ids.index(target_object_id)
    mask = (mask_logits[item_index] > 0.0).detach().cpu().numpy()
    x, y, area = mask_centroid(mask)
    centers[frame_index] = (x, y)
    areas[frame_index] = area


def track_view(
    config: ExperimentConfig,
    view_name: str,
    point_override: Sequence[float] | None = None,
) -> tuple[Path, Path]:
    """Track one view and write one sample for each aligned frame."""

    try:
        import torch
        from sam2.build_sam import build_sam2_video_predictor
    except ImportError as exc:
        raise RuntimeError(
            "Tracking requires PyTorch and SAM 2. Install requirements.txt first."
        ) from exc

    view = config.view(view_name)
    tracking = config.data["tracking"]
    frame_dir = config.resolve(view["frames_dir"])
    frame_paths = sorted_frame_paths(frame_dir)
    if len(frame_paths) != config.expected_frame_count:
        raise ValueError(
            f"{view_name} has {len(frame_paths)} frames; expected {config.expected_frame_count}"
        )

    if point_override is not None:
        points = [[float(point_override[0]), float(point_override[1])]]
        labels = [1]
    else:
        prompt = view.get("prompt", {})
        points = prompt.get("points") or _interactive_point(frame_paths[0])
        labels = prompt.get("labels", [1] * len(points))
    if len(points) != len(labels):
        raise ValueError(f"Prompt points and labels differ in length for {view_name}")

    checkpoint = config.resolve(tracking["checkpoint"])
    if not checkpoint.is_file():
        raise FileNotFoundError(
            f"SAM 2 checkpoint not found: {checkpoint}. Run scripts/download_checkpoint.py"
        )
    requested_device = str(tracking.get("device", "auto"))
    if requested_device == "auto":
        if torch.cuda.is_available():
            device = "cuda"
        elif getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
            device = "mps"
        else:
            device = "cpu"
    else:
        device = requested_device

    predictor = build_sam2_video_predictor(
        tracking["model_config"], str(checkpoint), device=device
    )
    centers = np.full((len(frame_paths), 2), np.nan, dtype=float)
    areas = np.zeros(len(frame_paths), dtype=np.int64)
    object_id = int(tracking.get("object_id", 1))
    annotation_frame = int(tracking.get("annotation_frame", 0))
    point_array = np.asarray(points, dtype=np.float32)
    label_array = np.asarray(labels, dtype=np.int32)

    autocast_context = nullcontext()
    if str(device).startswith("cuda"):
        dtype_name = str(tracking.get("autocast_dtype", "bfloat16"))
        autocast_context = torch.autocast("cuda", dtype=getattr(torch, dtype_name))

    with torch.inference_mode(), autocast_context:
        state = predictor.init_state(video_path=str(frame_dir))
        predictor.reset_state(state)
        _, object_ids, mask_logits = predictor.add_new_points_or_box(
            inference_state=state,
            frame_idx=annotation_frame,
            obj_id=object_id,
            points=point_array,
            labels=label_array,
        )
        _store_prediction(
            annotation_frame, object_ids, mask_logits, object_id, centers, areas
        )
        for frame_index, object_ids, mask_logits in predictor.propagate_in_video(state):
            _store_prediction(
                int(frame_index), object_ids, mask_logits, object_id, centers, areas
            )

    npy_path = config.resolve(view["track_npy"])
    csv_path = config.resolve(view["track_csv"])
    save_track(centers, areas, frame_paths, npy_path, csv_path)
    return npy_path, csv_path

