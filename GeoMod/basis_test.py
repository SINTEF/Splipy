# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis
import numpy as np
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


class TestBasis(unittest.TestCase):
    def test_get_continuity(self):
        b = BSplineBasis(4, [0, 0, 0, 0, .3, 1, 1, 1.134, 1.134, 1.134, 2, 2, 2, 2])
        self.assertEqual(b.get_continuity(.3), 2)
        self.assertEqual(b.get_continuity(1), 1)
        self.assertEqual(b.get_continuity(1.134), 0)
        self.assertEqual(b.get_continuity(0), -1)
        self.assertEqual(b.get_continuity(2), -1)
        self.assertEqual(b.get_continuity(.4), np.inf)

    def test_errors(self):
        with self.assertRaises(ValueError):
            BSplineBasis(4, [1, 2, 3])
        with self.assertRaises(ValueError):
            BSplineBasis(4, [1, 2, 4, 8], periodic=1)
        with self.assertRaises(ValueError):
            BSplineBasis(4, [1, 2, 3, 7, 8], periodic=2)

    def test_greville(self):
        b = BSplineBasis(4, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3])
        self.assertAlmostEqual(b.greville(0), 0.0)
        self.assertAlmostEqual(b.greville(1), 1.0 / 3.0)
        self.assertAlmostEqual(b.greville(2), 1.0)
        self.assertAlmostEqual(b.greville(), [0.0, 1.0/3.0, 1.0, 2.0, 8.0/3.0, 3.0])

    def test_raise_order(self):
        with self.assertRaises(ValueError):
            BSplineBasis().raise_order(-1)

    def test_write_g2(self):
        buf = StringIO()
        BSplineBasis(3, [0,0,0,1,2,3,3,3]).write_g2(buf)
        self.assertEqual(buf.getvalue().strip(),
                         '5 3\n0.000000 0.000000 0.000000 1.000000 2.000000 '
                         '3.000000 3.000000 3.000000')

    def test_getitem(self):
        b = BSplineBasis(3, [0,0,0,1,2,2,2])
        self.assertEqual(b[0], 0.0)
        self.assertEqual(b[1], 0.0)
        self.assertEqual(b[3], 1.0)

    def test_repr(self):
        self.assertEqual(repr(BSplineBasis()), 'p=2, [ 0.  0.  1.  1.]')
        self.assertEqual(repr(BSplineBasis(periodic=0)), 'p=2, [ 0.  0.  1.  1.], C0-periodic')

if __name__ == '__main__':
    unittest.main()
