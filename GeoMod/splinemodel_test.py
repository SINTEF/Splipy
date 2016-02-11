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


class TestOrientation(unittest.TestCase):

    def test_curves(self):
        c1 = Curve()
        c2 = Curve()

        tmp = Orientation.compute(c1)
        self.assertEqual(tmp.perm, (0,))
        self.assertEqual(tmp.flip, (False,))

        tmp = Orientation.compute(c1, c2)
        self.assertEqual(tmp.perm, (0,))
        self.assertEqual(tmp.flip, (False,))

        c1.reverse()
        tmp = Orientation.compute(c1, c2)
        self.assertEqual(tmp.perm, (0,))
        self.assertEqual(tmp.flip, (True,))

        c1 += (0.5, 0)
        with self.assertRaises(OrientationError):
            Orientation.compute(c1, c2)

        c1 -= (0.5, 0)
        c1.refine(2)
        with self.assertRaises(OrientationError):
            Orientation.compute(c1, c2)

        c2.refine(2)
        tmp = Orientation.compute(c1.controlpoints, c2.controlpoints)
        self.assertEqual(tmp.perm, (0,))
        self.assertEqual(tmp.flip, (True,))

        c2.set_dimension(3)
        with self.assertRaises(OrientationError):
            Orientation.compute(c1, c2)

    def test_surfaces(self):
        c = Curve()
        s1 = Surface()
        s2 = Surface()

        s1.refine(1, 0)
        s2.refine(1, 0)

        with self.assertRaises(OrientationError):
            Orientation.compute(c, s1)

        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (0, 1))
        self.assertTupleEqual(tmp.flip, (False, False))

        s1.reverse('u')
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (0, 1))
        self.assertTupleEqual(tmp.flip, (True, False))

        s2.reverse('v')
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (0, 1))
        self.assertTupleEqual(tmp.flip, (True, True))

        s2.reverse('u')
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (0, 1))
        self.assertTupleEqual(tmp.flip, (False, True))

        s2.reverse('v')
        s1.swap()
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (1, 0))
        self.assertTupleEqual(tmp.flip, (False, False))

        s1.reverse('u')
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (1, 0))
        self.assertTupleEqual(tmp.flip, (True, False))

        s2.reverse('v')
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (1, 0))
        self.assertTupleEqual(tmp.flip, (False, False))

        s2.reverse('u')
        tmp = Orientation.compute(s1, s2)
        self.assertTupleEqual(tmp.perm, (1, 0))
        self.assertTupleEqual(tmp.flip, (False, True))

    def test_compose(self):
        tmp = (Orientation((0, 1), (False, True)) *
               Orientation((0, 1), (True, True)))
        self.assertTupleEqual(tmp.perm, (0, 1))
        self.assertTupleEqual(tmp.flip, (True, False))

        tmp = (Orientation((0, 4, 1, 3, 2), (False, False, False, False, True)) *
               Orientation((4, 1, 3, 0, 2), (True, False, True, False, False)))
        self.assertTupleEqual(tmp.perm, (4, 2, 1, 0, 3))
        self.assertTupleEqual(tmp.flip, (True, False, False, False, False))

    def test_map(self):
        va = Volume()
        vb = Volume()
        vb.swap('u', 'w')
        vb.swap('u', 'v')
        vb.reverse('v')

        # va(u,v,w) corresponds to vb(w,u,1-v)
        ori = Orientation.compute(va, vb)

        self.assertTupleEqual(ori.perm, (2, 0, 1))
        self.assertTupleEqual(ori.perm_inv, (1, 2, 0))
        self.assertTupleEqual(ori.flip, (False, False, True))

        sections = [
            (( 0,  0,  0), ( 0,  0, -1), (), ()),
            ((-1,  0,  0), ( 0, -1, -1), (), ()),
            (( 0, -1,  0), ( 0,  0,  0), (), ()),
            ((-1, -1,  0), ( 0, -1,  0), (), ()),
            (( 0,  0, -1), (-1,  0, -1), (), ()),
            ((-1,  0, -1), (-1, -1, -1), (), ()),
            (( 0, -1, -1), (-1,  0,  0), (), ()),
            ((-1, -1, -1), (-1, -1,  0), (), ()),

            (( 0,  0, None), (None,  0, -1), (0,), (False,)),
            ((-1,  0, None), (None, -1, -1), (0,), (False,)),
            (( 0, -1, None), (None,  0,  0), (0,), (False,)),
            ((-1, -1, None), (None, -1,  0), (0,), (False,)),
            (( 0, None,  0), ( 0,  0, None), (0,), (True,)),
            ((-1, None,  0), ( 0, -1, None), (0,), (True,)),
            (( 0, None, -1), (-1,  0, None), (0,), (True,)),
            ((-1, None, -1), (-1, -1, None), (0,), (True,)),
            ((None,  0,  0), ( 0, None, -1), (0,), (False,)),
            ((None, -1,  0), ( 0, None,  0), (0,), (False,)),
            ((None,  0, -1), (-1, None, -1), (0,), (False,)),
            ((None, -1, -1), (-1, None,  0), (0,), (False,)),

            (( 0, None, None), (None,  0, None), (1, 0), (False, True)),
            ((-1, None, None), (None, -1, None), (1, 0), (False, True)),
            ((None,  0, None), (None, None, -1), (1, 0), (False, False)),
            ((None, -1, None), (None, None,  0), (1, 0), (False, False)),
            ((None, None,  0), ( 0, None, None), (0, 1), (False, True)),
            ((None, None, -1), (-1, None, None), (0, 1), (False, True)),

            ((None, None, None), (None, None, None), (2, 0, 1), (False, False, True)),
        ]

        for sec, mapped_sec, perm, flip in sections:
            view = ori.view_section(sec)
            self.assertTupleEqual(view.perm, perm)
            self.assertTupleEqual(view.flip, flip)
            self.assertTupleEqual(ori.map_section(sec), mapped_sec)
