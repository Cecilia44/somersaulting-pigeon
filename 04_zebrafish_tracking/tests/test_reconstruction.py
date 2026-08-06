from __future__ import annotations

import unittest

import numpy as np

from zebrafish3d.io import interpolate_missing, remove_legacy_initial_duplicate
from zebrafish3d.reconstruction import orthographic_fusion, triangulate_dlt


class TrajectoryIOTests(unittest.TestCase):
    def test_interpolate_missing_and_extend_edges(self) -> None:
        track = np.asarray([[np.nan, np.nan], [1.0, 2.0], [np.nan, np.nan], [3.0, 6.0]])
        actual = interpolate_missing(track)
        expected = np.asarray([[1.0, 2.0], [1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])
        np.testing.assert_allclose(actual, expected)

    def test_remove_known_legacy_duplicate(self) -> None:
        track = np.asarray([[1.0, 2.0], [1.0, 2.0], [3.0, 4.0]])
        repaired, changed = remove_legacy_initial_duplicate(track, expected_length=2)
        self.assertTrue(changed)
        np.testing.assert_allclose(repaired, [[1.0, 2.0], [3.0, 4.0]])


class ReconstructionTests(unittest.TestCase):
    def test_orthographic_fusion(self) -> None:
        top = np.asarray([[10.0, 20.0], [12.0, 24.0]])
        left = np.asarray([[99.0, 30.0], [98.0, 36.0]])
        actual = orthographic_fusion(
            left,
            top,
            origin_px=np.asarray([10.0, 20.0, 30.0]),
            scale_units_per_px=np.asarray([0.5, 0.5, 0.25]),
            axis_sign=np.asarray([1.0, -1.0, 1.0]),
        )
        np.testing.assert_allclose(actual, [[0.0, 0.0, 0.0], [1.0, -2.0, 1.5]])

    def test_dlt_recovers_known_point(self) -> None:
        top_projection = np.asarray(
            [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        )
        left_projection = np.asarray(
            [[1.0, 0.0, 0.0, -1.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]
        )
        world = np.asarray([2.0, 3.0, 4.0])
        top = np.asarray([[world[0] / world[2], world[1] / world[2]]])
        left = np.asarray([[(world[0] - 1.0) / world[2], world[1] / world[2]]])
        actual = triangulate_dlt(top, left, top_projection, left_projection)
        np.testing.assert_allclose(actual[0], world, atol=1e-10)


if __name__ == "__main__":
    unittest.main()

