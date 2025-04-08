from __future__ import annotations

import unittest

import splipy.curve_factory as cf
from splipy import TrimmedSurface


class TestSurface(unittest.TestCase):
    def test_constructor(self):
        trim_loop = [
            cf.line([0.25, 0.25], [0.75, 0.25]),
            cf.line([0.75, 0.25], [0.75, 0.75]),
            cf.line([0.75, 0.75], [0.25, 0.75]),
            cf.line([0.25, 0.75], [0.25, 0.25]),
        ]
        # error with non-closed trimming loop
        with self.assertRaises(RuntimeError):
            s = TrimmedSurface(loops=[trim_loop[:3]])
        # error with trimming curves in 3D space
        with self.assertRaises(RuntimeError):
            loop_3d = [
                cf.line([0.25, 0.25, 1], [0.75, 0.25, 1]),
                cf.line([0.75, 0.25, 1], [0.75, 0.75, 1]),
                cf.line([0.75, 0.75, 1], [0.25, 0.75, 1]),
                cf.line([0.25, 0.75, 1], [0.25, 0.25, 1]),
            ]
            s = TrimmedSurface(loops=[loop_3d])

        s = TrimmedSurface(loops=[trim_loop])
        self.assertEqual(s.dimension, 2)
        self.assertFalse(s.rational)


if __name__ == "__main__":
    unittest.main()
