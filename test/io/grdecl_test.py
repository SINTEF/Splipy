# -*- coding: utf-8 -*-

from splipy.io import GRDECL
import unittest
import os
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestGRDECL(unittest.TestCase):

    def test_read(self):
        with GRDECL(THIS_DIR + '/test_geometries/EightCells.grdecl') as myfile:
            vol1, vol2, vol3 = myfile.read()
            print()
            print()
            print(vol1)
            print()
            print()
            print()
            print(vol2)
            print()
            print()
            print(vol3)
            print()
            print()

        self.assertEqual(vol1.pardim, 3)
        self.assertEqual(vol2.pardim, 3)
        self.assertEqual(vol3.pardim, 3)
        self.assertEqual(vol1.dimension, 3)
        # self.assertEqual(vol.order(), (5, 2, 5))
        # self.assertEqual(vol.shape, (14, 2, 8))

if __name__ == '__main__':
    unittest.main()
