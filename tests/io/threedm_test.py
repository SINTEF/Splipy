# -*- coding: utf-8 -*-

try:
    from splipy.io import ThreeDM
    has_rhino = True
except ImportError:
    has_rhino = False

import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
import numpy as np
from numpy import sin, cos, pi
import os
import unittest

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class Test3DM(unittest.TestCase):

    @unittest.skipIf(not has_rhino, "3DM requires rhino3dm")
    def test_read_unit_curve(self):
        with ThreeDM(THIS_DIR + '/geometries/unit_curve.3dm') as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 2)
        self.assertEqual(len(crv),  2)
        self.assertTrue(np.allclose(crv[0], [0,0,0]))
        self.assertTrue(np.allclose(crv[1], [1,0,0]))

    @unittest.skipIf(not has_rhino, "3DM requires rhino3dm")
    def test_read_unit_square(self):
        with ThreeDM(THIS_DIR + '/geometries/unit_square.3dm') as myfile:
            srf = myfile.read()
        self.assertEqual(len(srf), 1)
        srf = srf[0]
        self.assertEqual(srf.order(0), 2)
        self.assertEqual(srf.order(1), 2)
        print(srf.shape)
        self.assertTrue(np.allclose(srf.shape, (2,2)))
        xrange, yrange, zrange = srf.bounding_box()
        self.assertTrue(np.allclose(xrange, [0,1]))
        self.assertTrue(np.allclose(yrange, [0,1]))
        self.assertTrue(np.allclose(zrange, [0,0]))

    @unittest.skipIf(not has_rhino, "3DM requires rhino3dm")
    def test_read_unit_cube(self):
        with ThreeDM(THIS_DIR + '/geometries/unit_cube.3dm') as myfile:
            srfs = myfile.read()
        self.assertEqual(len(srfs), 6)
        for s in srfs:
            self.assertAlmostEqual(s.area(), 1.0)

    @unittest.skipIf(not has_rhino, "3DM requires rhino3dm")
    def test_read_hexagon(self):
        with ThreeDM(THIS_DIR + '/geometries/unit_hexagon.3dm') as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 2)
        self.assertEqual(len(crv),  7)
        t = np.pi/3
        self.assertTrue(np.allclose(crv[0], [      1,      0,0]))
        self.assertTrue(np.allclose(crv[1], [ cos(t), sin(t),0]))
        self.assertTrue(np.allclose(crv[2], [-cos(t), sin(t),0]))
        self.assertTrue(np.allclose(crv[3], [     -1,      0,0]))
        self.assertTrue(np.allclose(crv[4], [-cos(t),-sin(t),0]))
        self.assertTrue(np.allclose(crv[5], [ cos(t),-sin(t),0]))

    @unittest.skipIf(not has_rhino, "3DM requires rhino3dm")
    def test_read_sphere(self):
        with ThreeDM(THIS_DIR + '/geometries/unit_sphere.3dm') as myfile:
            srf = myfile.read()
        self.assertEqual(len(srf), 1)
        srf = srf[0]
        self.assertEqual(srf.order(0), 3)
        self.assertEqual(srf.order(1), 3)
        self.assertAlmostEqual(srf.area(), 4*pi, places=2) # numerical integraion of area accurate to two digits

    @unittest.skipIf(not has_rhino, "3DM requires rhino3dm")
    def test_read_torus(self):
        with ThreeDM(THIS_DIR + '/geometries/torus_r1_R3.3dm') as myfile:
            srf = myfile.read()
        self.assertEqual(len(srf), 1)
        srf = srf[0]
        self.assertEqual(srf.order(0), 3)
        self.assertEqual(srf.order(1), 3)
        self.assertAlmostEqual(srf.area(), 4*pi*pi*3*1, places=2) # numerical integraion of area accurate to two digits


if __name__ == '__main__':
    unittest.main()
