from __future__ import annotations

import unittest
from math import pi

import numpy as np

from splipy import Curve, Surface

try:
    import cv2

    has_image = True
except ImportError:
    has_image = False

if has_image:
    from splipy.utils.image import *


class TestImage(unittest.TestCase):
    @unittest.skipIf(not has_image, "Image module requires OpenCV 2")
    def test_curve(self):
        crv = image_curves("test/utils/disc.png")
        self.assertEqual(type(crv), list)
        crv = crv[0]
        self.assertEqual(type(crv), Curve)

        bb = crv.bounding_box()
        crv -= [(x[0] + x[1]) / 2 for x in bb]  # move circle center to the origin
        crv /= (bb[0][1] - bb[0][0]) / 2  # scale to have unit radius

        t = np.linspace(crv.start(0), crv.end(0), 113)
        x = crv(t)
        radius = np.sqrt(x[:, 0] ** 2 + x[:, 1] ** 2)

        self.assertTrue(np.allclose(radius, 1.0, atol=1e-2))
        self.assertAlmostEqual(crv.length(), 2 * pi, places=1)

    @unittest.skipIf(not has_image, "Image module requires OpenCV 2")
    def test_surface(self):
        disc = image_convex_surface("test/utils/disc.png")
        self.assertEqual(type(disc), Surface)

        bb = disc.bounding_box()
        disc -= [(x[0] + x[1]) / 2 for x in bb]  # move circle center to the origin
        disc /= (bb[0][1] - bb[0][0]) / 2  # scale to have unit radius

        self.assertAlmostEqual(disc.area(), pi, places=1)

    @unittest.skipIf(not has_image, "Image module requires OpenCV 2")
    def test_height_disc(self):
        hill = image_height("test/utils/disc.png", N=[30, 30], p=[3, 3])
        self.assertEqual(type(hill), Surface)

        # computes \int z(u,v) dxdv where x is an approximation to x={0 inside circle, 1 outside circle}, radius=0.48
        center = hill.center()  # computes integral
        area_under_graph = center[-1]  # look at z-component
        r = 0.484  # same image used for edge detector so need white area between circle edge and image edge
        expected_value = 1 - pi * r * r  # circle = 0, background = 1, unit domain [0,1]^2

        self.assertAlmostEqual(area_under_graph, expected_value, places=1)

    @unittest.skipIf(not has_image, "Image module requires OpenCV 2")
    def test_height_orient(self):
        surf = image_height("test/utils/gray_corners.png")  # rectangular input file
        self.assertEqual(type(surf), Surface)

        self.assertAlmostEqual(surf[0, 0, 0], 0)  # check x-coordinates
        self.assertAlmostEqual(surf[-1, 0, 0], 1)
        self.assertAlmostEqual(surf[0, -1, 0], 0)
        self.assertAlmostEqual(surf[-1, -1, 0], 1)

        self.assertAlmostEqual(surf[0, 0, 1], 0)  # check y-coordinates
        self.assertAlmostEqual(surf[-1, 0, 1], 0)
        self.assertAlmostEqual(surf[0, -1, 1], 1)
        self.assertAlmostEqual(surf[-1, -1, 1], 1)

        self.assertAlmostEqual(surf[0, 0, 2], 1.00, places=2)  # check z-coordinates (image gray-value)
        self.assertAlmostEqual(surf[-1, 0, 2], 160.0 / 255, places=2)
        self.assertAlmostEqual(surf[0, -1, 2], 0.00, places=2)
        self.assertAlmostEqual(surf[-1, -1, 2], 80.0 / 255, places=2)


if __name__ == "__main__":
    unittest.main()
