# -*- coding: utf-8 -*-

from math import pi
import unittest

import numpy as np

from splipy import TrimmedSurface, BSplineBasis
import splipy.curve_factory as cf


class TestSurface(unittest.TestCase):
    def test_constructor(self):
        trim_loop = [cf.line([.25,.25], [.75,.25]),
                     cf.line([.75,.25], [.75,.75]),
                     cf.line([.75,.75], [.25,.75]),
                     cf.line([.25,.75], [.25,.25])]
        # error with non-closed trimming loop
        with self.assertRaises(RuntimeError):
            s = TrimmedSurface(loops=[trim_loop[:3]])
        # error with trimming curves in 3D space
        with self.assertRaises(RuntimeError):
            loop_3d = [cf.line([.25,.25,1], [.75,.25,1]),
                       cf.line([.75,.25,1], [.75,.75,1]),
                       cf.line([.75,.75,1], [.25,.75,1]),
                       cf.line([.25,.75,1], [.25,.25,1])]
            s = TrimmedSurface(loops=[loop_3d])

        s = TrimmedSurface(loops=[trim_loop])
        self.assertEqual(s.dimension, 2)
        self.assertFalse(s.rational)


if __name__ == '__main__':
    unittest.main()
