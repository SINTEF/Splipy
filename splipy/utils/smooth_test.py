# -*- coding: utf-8 -*-

from splipy import *
from splipy.utils.smooth import *
import numpy as np
import unittest

class TestSmooth(unittest.TestCase):
    def test_surface(self):
        b1 = BSplineBasis(2)
        b2 = BSplineBasis(2)
        surf = Surface(b1, b2)
        surf.refine(3)

        # check unchanged for linear splines
        cp = surf.controlpoints.copy()
        smooth(surf)
        self.assertTrue(np.allclose(cp, surf.controlpoints))

        # check that something changed for quadratic splines
        surf.raise_order(1,1)
        surf.refine(1)
        cp = surf.controlpoints.copy()
        smooth(surf)
        self.assertFalse(np.allclose(cp, surf.controlpoints))

        # check that the average of {2,3,4} is 3 and {3,4,5} is 4
        surf[0:3,0:4,1] = [[2,3,4,5], [2,3,4,5], [2,4,3,5]]
        smooth(surf)
        self.assertAlmostEqual(surf[1,1,1], 3.0)
        self.assertAlmostEqual(surf[1,2,1], 4.0)

    def test_curve(self):
        crv = Curve(BSplineBasis(5), [[0.,1.], [1,1], [2,1], [3,4], [7,4]])
        smooth(crv)
        self.assertTrue(np.allclose(crv[:,:], np.array([[0,1], [1,1], [2,2], [4,3], [7,4]])))


if __name__ == '__main__':
    unittest.main()
