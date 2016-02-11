# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis, Volume
from GeoMod.SplineModel import *
import unittest
import numpy as np


class TestVertexDict(unittest.TestCase):

    def test(self):
        d = VertexDict(0, 1e-4)

        self.assertEqual(len(d), 0)
        self.assertEqual(d.internal, [])

        d[np.array([0, 0])] = 'a'
        self.assertEqual(d[np.array([0, 0])], 'a')
        self.assertEqual(d[np.array([9e-5, 9e-5])], 'a')

        with self.assertRaises(KeyError):
            d[np.array([1, 1])]
        with self.assertRaises(KeyError):
            d[np.array([1e-3, 1e-3])]
        with self.assertRaises(KeyError):
            d[np.array([0, 1e-3])]

        d[np.array([9e-5, -9e-5])] = 'b'
        self.assertEqual(d[np.array([0, 0])], 'b')

        self.assertEqual(d.items()[0][0][0], 0.0)
        self.assertEqual(d.items()[0][0][1], 0.0)
        self.assertEqual(d.items()[0][1], 'b')

        for k in d:
            self.assertEqual(k[0], 0.0)
            self.assertEqual(k[1], 0.0)

        del d[np.array([0, 0])]
        self.assertEqual(len(d), 0)
