# Computational method

## Scope of the analysis

The reconstructed path is intended only as a qualitative visualization of the
zebrafish's turning behavior. No camera calibration or physical length scale was
acquired. Accordingly, coordinates are camera-pixel coordinates and are not
used for metric position, distance, speed, or acceleration measurements.

## 1. Manual temporal alignment

One visually corresponding frame is selected manually in the two recordings.
All subsequent frames are then paired one-to-one from that synchronization
point. For the included recording, the left frame range is labelled 72–160 and
the top range 65–153, producing 89 paired observations. The cleaned extractor
defines these configured ranges as zero-based and inclusive; the supplied JPGs
are the authoritative aligned frames for the archived example.

This procedure does not estimate clock drift or correct dropped frames. It is
appropriate here because the paired sequence is used for qualitative turning
visualization. The `time_s` column uses the original nominal rate of 60 paired
samples per second and must not be treated as precision kinematic timing.

## 2. Segmentation and tracking

One positive point on the fish is supplied to SAM 2.1 Hiera Large on the first
aligned frame of each view. The original mouse clicks were not saved. With the
authors' approval, the configuration uses each legacy first-mask centroid as
the replacement positive prompt. SAM 2 propagates the object mask through all
frames. For each mask with area `A`, the tracked center is its pixel centroid:

`center_x = mean(mask_column_indices)`

`center_y = mean(mask_row_indices)`

Empty masks are represented by `NaN`. CSV output records the aligned frame,
original source-frame index, center, mask area, and validity flag. The code
writes by SAM 2's returned frame index, so every aligned frame has exactly one
row.

## 3. Missing observations

By default, missing x and y coordinates are interpolated independently and
linearly. Leading or trailing gaps use the nearest valid observation. The final
CSV explicitly marks which input observations were interpolated.

## 4. 3D reconstruction modes

### Current analysis: qualitative orthographic coordinate fusion

The experiment has no camera calibration or physical scale. For qualitative
visualization, the preserved method forms an image-coordinate trajectory as

`X = top_x`, `Y = top_y`, `Z = left_y`.

Configurable origins, axis signs, and units-per-pixel scales are applied after
this fusion. With the included configuration all scales are one, so coordinates
remain in pixels. Because the two cameras have different image dimensions,
axis magnitudes are not directly comparable physical distances. The left-view
horizontal coordinate is retained in the CSV for audit but is not used. This
result is described as a qualitative three-coordinate trajectory, not metric
stereo triangulation.

### Calibrated option: DLT triangulation

For a future calibrated experiment, if each camera's `3 x 4` projection matrix
is available, set `method` to
`dlt_triangulation` and provide `top_projection_matrix` and
`left_projection_matrix`. The implementation solves the four-equation DLT
system by SVD for every paired observation and reports per-view reprojection
errors. Output units follow the coordinate system used during calibration.

## 5. Visualization

Plots connect reconstructed observations directly and color line segments by
time. No spline or temporal smoothing is applied to the plotted coordinates.
This avoids introducing unreported overshoot into the displayed path.
