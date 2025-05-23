from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from splipy.io import SPL

THIS_DIR = str(Path(__file__).parent)


class TestSPL(unittest.TestCase):
    def test_read(self):
        with SPL(THIS_DIR + "/geometries/ibeam.spl") as myfile:
            vol = myfile.read()[0]

        self.assertEqual(vol.pardim, 3)
        self.assertEqual(vol.dimension, 3)
        self.assertEqual(vol.order(), (5, 2, 5))
        self.assertEqual(vol.shape, (14, 2, 8))

        np.testing.assert_almost_equal(
            vol.knots("u", with_multiplicities=True),
            [
                0.0000000000000000e00,
                0.0000000000000000e00,
                0.0000000000000000e00,
                0.0000000000000000e00,
                0.0000000000000000e00,
                4.1666666666666674e-01,
                4.1666666666666674e-01,
                4.1666666666666674e-01,
                5.0000000000000011e-01,
                5.0000000000000011e-01,
                5.0000000000000011e-01,
                5.8333333333333337e-01,
                5.8333333333333337e-01,
                5.8333333333333337e-01,
                1.0000000000000000e00,
                1.0000000000000000e00,
                1.0000000000000000e00,
                1.0000000000000000e00,
                1.0000000000000000e00,
            ],
        )

        np.testing.assert_almost_equal(vol.knots("v", with_multiplicities=True), [0.0, 0.0, 1.0, 1.0])

        np.testing.assert_almost_equal(
            vol.knots("w", with_multiplicities=True),
            [
                0.0000000000000000e00,
                0.0000000000000000e00,
                0.0000000000000000e00,
                0.0000000000000000e00,
                0.0000000000000000e00,
                5.0000000000000000e-01,
                5.0000000000000000e-01,
                5.0000000000000000e-01,
                1.0000000000000000e00,
                1.0000000000000000e00,
                1.0000000000000000e00,
                1.0000000000000000e00,
                1.0000000000000000e00,
            ],
        )

        np.testing.assert_almost_equal(
            vol.controlpoints[:, 0, 0, 0],
            [
                2.5000000000000000e-01,
                2.5000000000000000e-01,
                2.5000000000000000e-01,
                2.3612428004820077e-01,
                1.0277514399035985e-01,
                1.0000000000000001e-01,
                1.0000000000000001e-01,
                1.0000000000000001e-01,
                1.0000000000000001e-01,
                1.0277514399035984e-01,
                2.3612428004820074e-01,
                2.5000000000000000e-01,
                2.5000000000000000e-01,
                2.5000000000000000e-01,
            ],
        )

        np.testing.assert_almost_equal(
            vol.controlpoints[:, 1, 0, 0],
            [
                -2.5000000000000000e-01,
                -2.5000000000000000e-01,
                -2.5000000000000000e-01,
                -2.3612428004820077e-01,
                -1.0277514399035985e-01,
                -1.0000000000000001e-01,
                -1.0000000000000001e-01,
                -1.0000000000000001e-01,
                -1.0000000000000001e-01,
                -1.0277514399035984e-01,
                -2.3612428004820074e-01,
                -2.5000000000000000e-01,
                -2.5000000000000000e-01,
                -2.5000000000000000e-01,
            ],
        )


if __name__ == "__main__":
    unittest.main()
