# -*- coding: utf-8 -*-

from splipy import Curve, Surface
from math import pi
import numpy as np
import unittest
import sys

if sys.version_info < (3,):
    from splipy.utils.Image import *

class TestImage(unittest.TestCase):

    @unittest.skipIf(sys.version_info >= (3,), "Image module not supported on Python 3")
    def test_curve(self):
        crv = image_curves('splipy/utils/disc.png')
        self.assertEqual(type(crv), Curve)

        bb = crv.bounding_box()
        crv -= [(x[0]+x[1])/2 for x in bb] # move circle center to the origin
        crv /= (bb[0][1]-bb[0][0])/2       # scale to have unit radius
        
        t = np.linspace(crv.start(0), crv.end(0), 113)
        x = crv(t)
        radius = np.sqrt(x[:,0]**2 + x[:,1]**2)

        self.assertTrue(np.allclose(radius, 1.0, atol=1e-2))
        self.assertAlmostEqual(crv.length(), 2*pi, places=1)

    @unittest.skipIf(sys.version_info >= (3,), "Image module not supported on Python 3")
    def test_surface(self):
        disc = image_convex_surface('splipy/utils/disc.png')
        self.assertEqual(type(disc), Surface)

        bb = disc.bounding_box()
        disc -= [(x[0]+x[1])/2 for x in bb] # move circle center to the origin
        disc /= (bb[0][1]-bb[0][0])/2       # scale to have unit radius

        self.assertAlmostEqual(disc.area(), pi, places=1)

if __name__ == '__main__':
    unittest.main()
