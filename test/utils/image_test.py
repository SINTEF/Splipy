# -*- coding: utf-8 -*-

from splipy import Curve, Surface
from math import pi
import numpy as np
import unittest
import sys

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
        crv = image_curves('splipy/utils/disc.png')
        self.assertEqual(type(crv), list)
        crv = crv[0]
        self.assertEqual(type(crv), Curve)

        bb = crv.bounding_box()
        crv -= [(x[0]+x[1])/2 for x in bb] # move circle center to the origin
        crv /= (bb[0][1]-bb[0][0])/2       # scale to have unit radius

        t = np.linspace(crv.start(0), crv.end(0), 113)
        x = crv(t)
        radius = np.sqrt(x[:,0]**2 + x[:,1]**2)

        self.assertTrue(np.allclose(radius, 1.0, atol=1e-2))
        self.assertAlmostEqual(crv.length(), 2*pi, places=1)

    @unittest.skipIf(not has_image, "Image module requires OpenCV 2")
    def test_surface(self):
        disc = image_convex_surface('splipy/utils/disc.png')
        self.assertEqual(type(disc), Surface)

        bb = disc.bounding_box()
        disc -= [(x[0]+x[1])/2 for x in bb] # move circle center to the origin
        disc /= (bb[0][1]-bb[0][0])/2       # scale to have unit radius

        self.assertAlmostEqual(disc.area(), pi, places=1)

if __name__ == '__main__':
    unittest.main()
