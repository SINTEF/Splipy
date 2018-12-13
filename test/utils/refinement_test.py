# -*- coding: utf-8 -*-

from splipy import *
from splipy.utils.refinement import *
import numpy as np
import unittest

class TestRefinement(unittest.TestCase):

    def test_geometric_curve(self):
        crv = Curve()
        geometric_refine(crv, 0.8, 10)

        self.assertEqual(len(crv), 12)
        k = crv.knots(0)
        for k0,k1,k2 in zip(k[:-2], k[1:-1], k[2:]):
            self.assertAlmostEqual((k2-k1)/(k1-k0), 0.8)

    def test_geometric_surface(self):
        surf = Surface()
        surf.raise_order(2,2)
        geometric_refine(surf, 0.7, 7, 0)
        geometric_refine(surf, 0.9, 8, 1)

        (knot1,knot2) = surf.knots()

        self.assertEqual(len(knot1), 9)
        self.assertEqual(len(knot2), 10)
        self.assertEqual(len(surf.knots(direction='u', with_multiplicities=True)), 15)
        self.assertEqual(len(surf.knots(direction='v', with_multiplicities=True)), 16)

        for k0,k1,k2 in zip(knot1[:-2], knot1[1:-1], knot1[2:]):
            self.assertAlmostEqual((k2-k1)/(k1-k0), 0.7)
        for k0,k1,k2 in zip(knot2[:-2], knot2[1:-1], knot2[2:]):
            self.assertAlmostEqual((k2-k1)/(k1-k0), 0.9)

    def test_geometric_volume(self):
        vol = Volume()
        vol.raise_order(1,2,3)
        geometric_refine(vol, 0.7, 4, 0)
        geometric_refine(vol, 0.9, 4, 1, reverse=True)
        geometric_refine(vol, 0.8, 2, 2)

        (knot1,knot2,knot3) = vol.knots()

        self.assertEqual(len(knot1), 6)
        self.assertEqual(len(knot2), 6)
        self.assertEqual(len(knot3), 4)

        for k0,k1,k2 in zip(knot1[:-2], knot1[1:-1], knot1[2:]):
            self.assertAlmostEqual((k2-k1)/(k1-k0), 0.7)
        for k0,k1,k2 in zip(knot2[:-2], knot2[1:-1], knot2[2:]):
            self.assertAlmostEqual((k1-k0)/(k2-k1), 0.9)
        for k0,k1,k2 in zip(knot3[:-2], knot3[1:-1], knot3[2:]):
            self.assertAlmostEqual((k2-k1)/(k1-k0), 0.8)

    def test_edge_surface(self):
        surf = Surface()
        surf.raise_order(2,2)
        edge_refine(surf, 10, 9, direction=0)

        knots = surf.knots(direction='u')
        self.assertEqual(len(knots), 11)

        # test that they form symmetric from each edge
        for low,high in zip(knots, knots[::-1]):
            self.assertAlmostEqual(low, 1-high)

    def test_edge_volume(self):
        vol = Volume()
        vol.raise_order(3,1,2)
        edge_refine(vol, 7, 4, direction=0)

        knots = vol.knots(direction='u')
        self.assertEqual(len(knots), 6)

        # test that they form symmetric from each edge
        for low,high in zip(knots, knots[::-1]):
            self.assertAlmostEqual(low, 1-high)

    def test_subdivide_surface(self):
        surf = Surface()
        surf.raise_order(2,2)
        surf.refine(30,30)
        pieces = subdivide([surf], 9)
        self.assertEqual(len(pieces), 100)

    def test_subdivide_volume(self):
        vol = Volume()
        vol.raise_order(2,1,3)
        vol.refine(10,19,12)
        pieces = subdivide([vol], (4,4,4))
        self.assertEqual(len(pieces), 125)


if __name__ == '__main__':
    unittest.main()
