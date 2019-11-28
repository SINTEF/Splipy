# -*- coding: utf-8 -*-

from splipy.io import GRDECL
import unittest
import os
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

class TestGRDECL(unittest.TestCase):

    def test_read(self):
        with GRDECL(THIS_DIR + '/geometries/EightCells.grdecl') as myfile:
            myfile.read()
            vol    = myfile.get_c0_mesh() # c0 mesh, compute average value at nodes
            mx_vol = myfile.get_mixed_cont_mesh() # c0 in (x,y)-directions, C^{-1} in z-direction
            dc_vol = myfile.get_cm1_mesh() # C^{-1} mesh in all directions

            # check volume output
            self.assertEqual(vol.pardim,    3)
            self.assertEqual(mx_vol.pardim, 3)
            self.assertEqual(dc_vol.pardim, 3)

            # in three spatial dimensions
            self.assertEqual(vol.dimension,    3)
            self.assertEqual(mx_vol.dimension, 3)
            self.assertEqual(dc_vol.dimension, 3)

            # all linear spline representations
            self.assertEqual(vol.order(),    (2, 2, 2))
            self.assertEqual(mx_vol.order(), (2, 2, 2))
            self.assertEqual(dc_vol.order(), (2, 2, 2))

            # different number of basis functions dending on continuity
            self.assertEqual(vol.shape,    (3,3,3))
            self.assertEqual(mx_vol.shape, (3,3,4))
            self.assertEqual(dc_vol.shape, (4,4,4))

if __name__ == '__main__':
    unittest.main()
