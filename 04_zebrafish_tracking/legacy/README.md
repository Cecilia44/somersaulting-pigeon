# Legacy scripts

These files are preserved verbatim for provenance. Use the package in `src/`
for reproducible analyses.

- `sam_track.py` hard-codes paths and appends frame zero twice: once immediately
  after prompting and once when SAM 2 propagation returns frame zero. The saved
  demo arrays therefore contain 90 points for 89 aligned frames.
- `plot_3d_traj.py` hard-codes another experiment ID, silently truncates unequal
  tracks, labels samples as seconds using an implicit 60 fps assumption, and
  applies a cubic spline only for display. Its reconstruction is coordinate
  fusion `(top_x, top_y, left_y)`, not calibrated stereo triangulation.

