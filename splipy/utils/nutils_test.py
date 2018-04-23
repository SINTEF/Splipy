# -*- coding: utf-8 -*-

from splipy import Surface, Volume
from splipy.utils.nutils import *
import numpy as np
import unittest


class TestNutilsUtils(unittest.TestCase):

    def test_surface(self):
        surf = Surface()
        expected = np.array([[0,0], [0,1], [1,0], [1,1]])
        self.assertTrue(np.allclose(controlpoints(surf), expected))

    def test_volume(self):
        vol = Volume()
        expected = np.array([[0,0,0], [0,0,1], [0,1,0], [0,1,1],
                             [1,0,0], [1,0,1], [1,1,0], [1,1,1]])
        self.assertTrue(np.allclose(controlpoints(vol), expected))

    def test_multiplicity(self):
        surf = Surface()
        surf.refine(1)
        surf.raise_order(1)
        surf.refine(1)
        self.assertTrue(np.allclose(multiplicities(surf), [[3, 1, 2, 1, 3], [3, 1, 2, 1, 3]]))

if __name__ == '__main__':
    unittest.main()
