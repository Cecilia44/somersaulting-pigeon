# Two-view zebrafish tracking and qualitative trajectory visualization

This repository turns manually synchronized top- and left-view zebrafish videos
into a qualitative three-coordinate center trajectory for visualizing turning
behavior. It records every manual choice in one JSON file and exposes the
complete workflow as command-line steps:

1. manually define the paired frame ranges;
2. extract the aligned frames;
3. prompt SAM 2 once per view and propagate the fish mask;
4. compute each mask's centroid;
5. fuse orthogonal image coordinates (or, for future calibrated data,
   triangulate camera rays); and
6. export auditable CSV/NPY data and a time-colored figure.

The included example is the recording `2025-03-19_17-22`. The original videos,
89 aligned frames per view, legacy center tracks, and legacy SVG are retained so
the cleaned pipeline can be checked against the earlier analysis.

## Repository layout

```text
configs/                 experiment choices, paths, prompts, and reconstruction
data/raw/                original left/top AVI files
data/aligned/            manually aligned frame ranges (source indices retained)
data/processed/          legacy and newly generated 2D centroid tracks
docs/METHODS.md          equations, assumptions, and reporting guidance
legacy/                  original scripts preserved verbatim for provenance
scripts/                 model-checkpoint downloader
src/zebrafish3d/         reusable pipeline implementation and CLI
tests/                   dependency-light numerical tests
outputs/                 generated 3D tables and figures
```

## Requirements

- Python 3.10–3.13 (Python 3.10 is recommended for the pinned environment)
- Git, because SAM 2 is installed from a pinned upstream commit
- A CUDA-capable NVIDIA GPU is strongly recommended for SAM 2.1 Hiera Large;
  extraction, reconstruction, tests, and plotting can run on CPU
- approximately 0.9 GB of disk space for the model checkpoint, plus the Python
  environment

The exact Python dependencies are pinned in `pyproject.toml` and installed by
`requirements.txt`: NumPy 2.2.4, Matplotlib 3.10.1, OpenCV 4.11.0.86, Pillow
11.1.0, PyTorch 2.6.0, torchvision 0.21.0, and SAM 2 at commit
`2b90b9f5ceec907a1c18123530e92e794ad901a4`.

Create an isolated environment from the repository root:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python scripts/download_checkpoint.py
```

For a machine-specific CUDA build, install the PyTorch 2.6.0/torchvision 0.21.0
wheels recommended for that CUDA runtime first, then run the requirements
command. Pip will keep the already satisfied versions.

## Quick start with the included data

Validate the recorded alignment and legacy results:

```bash
zebrafish3d --config configs/2025-03-19_17-22.json validate
```

Re-extract the manually selected ranges. Frame indices are zero-based and both
ends are inclusive. Existing complete frame sets are only validated unless
`--overwrite` is supplied.

```bash
zebrafish3d --config configs/2025-03-19_17-22.json extract
```

Run SAM 2 for both views. The original mouse clicks were not stored; with the
authors' approval, the saved configuration uses each legacy first-mask centroid
as the replacement positive prompt, so this command is non-interactive:

```bash
zebrafish3d --config configs/2025-03-19_17-22.json track --view all
```

To choose the fish manually for one view, remove that view's `prompt.points`
from the JSON or override it explicitly:

```bash
zebrafish3d --config configs/2025-03-19_17-22.json track --view top --point 1161.4 542.8
```

Reconstruct and plot:

```bash
zebrafish3d --config configs/2025-03-19_17-22.json reconstruct
zebrafish3d --config configs/2025-03-19_17-22.json plot
```

The reconstruction command prefers newly generated `*.centers.npy` tracks. If
they do not exist, it uses the supplied legacy arrays and safely removes their
known duplicated first sample. It never silently truncates unequal tracks.

## Included alignment and timing

| View | Source video | Reported video rate | Selected source frames | Paired frames |
|---|---|---:|---:|---:|
| Left | `Left_2025-03-19_17-22.avi` | 60 fps | 72–160 | 89 |
| Top | `Top_2025-03-19_17-22.avi` | about 57.58 fps average | 65–153 | 89 |

One visually corresponding frame was selected manually in the two videos, and
all subsequent frames were paired one-to-one from that point. The cleaned
extractor treats configured ranges as zero-based and inclusive; the supplied
JPGs are the authoritative aligned frames for this archived example. The
`time_s` column uses the original nominal rate of 60 paired samples per second.
Because the source rates differ and no clock-drift correction was performed,
that column is for visualization and not precision kinematic timing.

## Reconstruction caveat

No camera calibration or physical scale was acquired because the output is used
only to visualize turning behavior. The configured method therefore uses
`(X, Y, Z) = (top_x, top_y, left_y)` in pixel units. This is an orthographic
coordinate fusion, not metric stereo triangulation; the left-view horizontal
coordinate is not used and the three axis magnitudes are not comparable physical
distances.

The repository retains an optional DLT implementation for future experiments
that provide two `3 x 4` camera projection matrices. It is not used for the
included recording. See `docs/METHODS.md` for the exact distinction.

## Output files

- Per-view `*.centers.csv`: source frame, centroid in pixels, mask area, validity
- Per-view `*.centers.npy`: an `N x 2` floating-point center array; missing masks
  are `NaN`
- `trajectory_3d.csv`: original/interpolated 2D observations, interpolation
  flags, time, 3D coordinates, and (for calibrated DLT) reprojection errors
- `trajectory_3d.npy`: an `N x 3` floating-point trajectory
- `trajectory_3d.png` and `.svg`: unsmoothed, time-colored trajectory figures

## Tests

The numerical tests do not load SAM 2 or require a GPU:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

They cover missing-value interpolation, repair of the legacy duplicated-first-
frame bug, orthographic coordinate fusion, and calibrated DLT triangulation.

## Before public release

Add the final author list, repository license, citation metadata, and animal
protocol information. In the manuscript and figure legend, state that temporal
alignment was initialized by a manually selected corresponding frame and that
the displayed trajectory is an uncalibrated, qualitative pixel-coordinate
visualization of turning behavior. Do not report it as millimetre-scale 3D
reconstruction or use it for metric kinematics.
