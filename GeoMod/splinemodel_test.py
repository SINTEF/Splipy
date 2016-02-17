# -*- coding: utf-8 -*-

from GeoMod import *
from GeoMod.SplineModel import *
from GeoMod.ModState import state
import unittest
from math import pi
import numpy as np


class TestVertexDict(unittest.TestCase):

    def test(self):
        with state(controlpoint_absolute_tolerance=1e-4):
            d = VertexDict()

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


class TestObjectCatalogue(unittest.TestCase):

    def test_curves(self):
        curves = [Curve(),
                  Curve() + (1,0),
                  Curve() + (2,0),
                  Curve() + (1,1),
                  Curve() + (2,1)]
        curves[1].refine(2)
        curves[2].reverse()

        cat = ObjectCatalogue(1)
        for c in curves:
            cat(c)

        self.assertEqual(len(cat.nodes(1)), 5)
        self.assertEqual(len(cat.top_nodes()), 5)
        self.assertEqual(len(cat.nodes(0)), 7)

        vertex_nodes = [(cat(c.section(u=0, unwrap_points=False)).node,
                         cat(c.section(u=-1, unwrap_points=False)).node)
                        for c in curves]
        curve_nodes = [cat(c).node for c in curves]

        for cn, vn in zip(curve_nodes, vertex_nodes):
            self.assertIs(vn[0], cn.lower_nodes[0][0])
            self.assertIs(vn[1], cn.lower_nodes[0][1])

        self.assertSetEqual(vertex_nodes[0][0].higher_nodes[1], {curve_nodes[0]})
        self.assertSetEqual(vertex_nodes[0][1].higher_nodes[1], {curve_nodes[0], curve_nodes[1]})
        self.assertSetEqual(vertex_nodes[1][0].higher_nodes[1], {curve_nodes[0], curve_nodes[1]})
        self.assertSetEqual(vertex_nodes[1][1].higher_nodes[1], {curve_nodes[1], curve_nodes[2]})
        self.assertSetEqual(vertex_nodes[2][0].higher_nodes[1], {curve_nodes[2]})
        self.assertSetEqual(vertex_nodes[2][1].higher_nodes[1], {curve_nodes[1], curve_nodes[2]})
        self.assertSetEqual(vertex_nodes[3][0].higher_nodes[1], {curve_nodes[3]})
        self.assertSetEqual(vertex_nodes[3][1].higher_nodes[1], {curve_nodes[3], curve_nodes[4]})
        self.assertSetEqual(vertex_nodes[4][0].higher_nodes[1], {curve_nodes[3], curve_nodes[4]})
        self.assertSetEqual(vertex_nodes[4][1].higher_nodes[1], {curve_nodes[4]})

        c = cat(Curve() + (2,0))
        self.assertIs(c.node, curve_nodes[2])
        self.assertTupleEqual(c.orientation.perm, (0,))
        self.assertTupleEqual(c.orientation.flip, (True,))

        curve = Curve() + (1,0)
        curve.refine(1)
        self.assertIsNot(cat(curve).node, curve_nodes[1])

        curve = Curve() + (1,0)
        curve.refine(2)
        curve.reverse()
        c = cat(curve)
        self.assertIs(c.node, curve_nodes[1])
        self.assertTupleEqual(c.orientation.perm, (0,))
        self.assertTupleEqual(c.orientation.flip, (True,))

        self.assertEqual(len(cat.nodes(1)), 6)
        self.assertEqual(len(cat.top_nodes()), 6)
        self.assertEqual(len(cat.nodes(0)), 7)

    def test_surfaces(self):

        s1 = Surface()
        s2 = s1.clone()
        s2.rotate(pi/2)
        s2 += (2,0)
        s3 = s1 + (1,1)
        s3.reverse('v')
        s4 = s1 + (2,1)
        s4.reverse('u')
        s4.reverse('v')
        s4.swap()

        cat = ObjectCatalogue(2)
        cat(s1)
        cat(s2)
        cat(s3)
        cat(s4)

        self.assertIs(cat(s1).section(u=-1).node, cat(s2).section(v=-1).node)
        self.assertIs(cat(s2).section(u=-1).node, cat(s3).section(v=-1).node)
        self.assertIs(cat(s3).section(u=-1).node, cat(s4).section(v=-1).node)

        self.assertIs(cat(s1).section(-1, 0).node, cat(s2).section(0, -1).node)
        self.assertIs(cat(s1).section(-1, -1).node, cat(s2).section(-1, -1).node)
        self.assertIs(cat(s2).section(-1, -1).node, cat(s3).section(0, -1).node)
        self.assertIs(cat(s2).section(-1, 0).node, cat(s3).section(-1, -1).node)
        self.assertIs(cat(s3).section(-1, 0).node, cat(s4).section(0, -1).node)
        self.assertIs(cat(s3).section(-1, -1).node, cat(s4).section(-1, -1).node)

        self.assertTupleEqual(cat(s1).section(u=0).orientation.perm, (0,))
        self.assertTupleEqual(cat(s1).section(u=0).orientation.flip, (False,))
        self.assertTupleEqual(cat(s2).section(v=-1).orientation.perm, (0,))
        self.assertTupleEqual(cat(s2).section(v=-1).orientation.flip, (False,))
        self.assertTupleEqual(cat(s2).section(u=-1).orientation.perm, (0,))
        self.assertTupleEqual(cat(s2).section(u=-1).orientation.flip, (False,))
        self.assertTupleEqual(cat(s3).section(v=-1).orientation.perm, (0,))
        self.assertTupleEqual(cat(s3).section(v=-1).orientation.flip, (True,))
        self.assertTupleEqual(cat(s3).section(u=-1).orientation.perm, (0,))
        self.assertTupleEqual(cat(s3).section(u=-1).orientation.flip, (False,))
        self.assertTupleEqual(cat(s4).section(v=-1).orientation.perm, (0,))
        self.assertTupleEqual(cat(s4).section(v=-1).orientation.flip, (False,))

    def test_volumes(self):

        v1 = Volume()

        v2 = v1 + (1,0,0)
        v2.reverse('u')
        v2.reverse('v')

        v3 = v1 + (0,1,0)
        v3.reverse('w')
        v3.swap('u', 'w')

        v4 = v1 + (0,0,1)
        v4.swap('v', 'w')

        cat = ObjectCatalogue(3)
        cat(v1)
        cat(v2)
        cat(v3)
        cat(v4)

        self.assertIs(cat(v1).section(u=-1).node, cat(v2).section(u=-1).node)
        self.assertIs(cat(v1).section(v=-1).node, cat(v3).section(v=0).node)
        self.assertIs(cat(v1).section(w=-1).node, cat(v4).section(v=0).node)

        self.assertTupleEqual(cat(v2).section(u=-1).orientation.perm, (0, 1))
        self.assertTupleEqual(cat(v2).section(u=-1).orientation.flip, (True, False))
        self.assertTupleEqual(cat(v3).section(v=0).orientation.perm, (1, 0))
        self.assertTupleEqual(cat(v3).section(v=0).orientation.flip, (False, True))
        self.assertTupleEqual(cat(v4).section(v=0).orientation.perm, (0, 1))
        self.assertTupleEqual(cat(v4).section(v=0).orientation.flip, (False, False))

        self.assertIs(cat(v1).section(u=-1, v=0).node, cat(v2).section(u=-1, v=-1).node)
        self.assertIs(cat(v1).section(u=-1, v=-1).node, cat(v2).section(u=-1, v=0).node)
        self.assertIs(cat(v1).section(u=-1, w=0).node, cat(v2).section(u=-1, w=0).node)
        self.assertIs(cat(v1).section(u=-1, w=-1).node, cat(v2).section(u=-1, w=-1).node)
        self.assertIs(cat(v1).section(v=-1, u=0).node, cat(v3).section(v=0, w=0).node)
        self.assertIs(cat(v1).section(v=-1, w=-1).node, cat(v3).section(u=0, v=0).node)
        self.assertIs(cat(v1).section(u=0, w=-1).node, cat(v4).section(u=0, v=0).node)

        self.assertTupleEqual(cat(v1).section(u=-1, v=0).orientation.flip, (False,))
        self.assertTupleEqual(cat(v2).section(u=-1, v=-1).orientation.flip, (False,))
        self.assertTupleEqual(cat(v1).section(u=-1, w=0).orientation.flip, (False,))
        self.assertTupleEqual(cat(v2).section(u=-1, w=0).orientation.flip, (True,))
        self.assertTupleEqual(cat(v1).section(v=-1, u=0).orientation.flip, (False,))
        self.assertTupleEqual(cat(v3).section(v=0, w=0).orientation.flip, (True,))
        self.assertTupleEqual(cat(v3).section(u=0, v=0).orientation.flip, (False,))
        self.assertTupleEqual(cat(v4).section(u=0, v=0).orientation.flip, (False,))
