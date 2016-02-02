# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis, Volume
import unittest
import numpy as np


class TestVolume(unittest.TestCase):
    def test_evaluate(self):
        # creating the identity mapping by different size for all directions
        vol = Volume(BSplineBasis(7), BSplineBasis(6), BSplineBasis(5))

        # call evaluation at a 2x3x4 grid of points
        u_val = np.linspace(0, 1, 2)
        v_val = np.linspace(0, 1, 3)
        w_val = np.linspace(0, 1, 4)
        value = vol(u_val, v_val, w_val)
        self.assertEqual(value.shape[0], 2)  # 2 u-evaluation points
        self.assertEqual(value.shape[1], 3)  # 3 v-evaluation points
        self.assertEqual(value.shape[2], 4)  # 4 w-evaluation points
        self.assertEqual(value.shape[3], 3)  # 3 dimensions (x,y,z)
        self.assertEqual(vol.order(), (7, 6, 5))
        for i, u in enumerate(u_val):
            for j, v in enumerate(v_val):
                for k, w in enumerate(w_val):
                    self.assertAlmostEqual(value[i, j, k, 0], u)  # identity map x=u
                    self.assertAlmostEqual(value[i, j, k, 1], v)  # identity map y=v
                    self.assertAlmostEqual(value[i, j, k, 2], w)  # identity map z=w

        # test errors and exceptions
        with self.assertRaises(ValueError):
            val = vol(-10, .5, .5)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = vol(+10, .3, .3)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = vol(.5, -10, .123)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = vol(.5, +10, .123)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = vol(.5, .2, +10)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = vol(.5, .2, -10)  # evalaute outside parametric domain

    def test_raise_order(self):
        # more or less random 3D volume with p=[2,2,1] and n=[4,3,2]
        controlpoints = [[0, 0, 0], [-1, 1, 0], [0, 2, 0], [1, -1, 0], [1, 0, 0], [1, 1, 0],
                         [2, 1, 0], [2, 2, 0], [2, 3, 0], [3, 0, 0], [4, 1, 0], [3, 2, 0],
                         [0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 2], [1, 0, 2], [1, 1, 2],
                         [2, 1, 2], [2, 2, 2], [2, 3, 2], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints)

        self.assertEqual(vol.order(), (3, 3, 2))
        evaluation_point1 = vol(0.23, 0.37, 0.44)  # pick some evaluation point (could be anything)

        vol.raise_order(1, 2, 4)

        self.assertEqual(vol.order(), (4, 5, 6))
        evaluation_point2 = vol(0.23, 0.37, 0.44)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # test a rational 3D volume
        controlpoints = [[0, 0, 1, 1], [-1, 1, .96, 1], [0, 2, 1, 1], [1, -1, 1, 1], [1, 0, .8, 1],
                         [1, 1, 1, 1], [2, 1, .89, 1], [2, 2, .9, 1], [2, 3, 1, 1], [3, 0, 1, 1],
                         [4, 1, 1, 1], [3, 2, 1, 1], [0, 0, 1, 2], [-1, 1, .7, 2], [0, 2, 1.3, 2],
                         [1, -1, 1, 2], [1, 0, .77, 2], [1, 1, 1, 2], [2, 1, .89, 1], [2, 2, .8, 4],
                         [2, 3, 1, 1], [3, 0, 1, 1], [4, 1, 1, 1], [3, 2, 1, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints, True)

        self.assertEqual(vol.order()[0], 3)
        self.assertEqual(vol.order()[1], 3)
        evaluation_point1 = vol(0.23, 0.37, 0.44)

        vol.raise_order(1, 2, 1)

        self.assertEqual(vol.order(), (4, 5, 3))
        evaluation_point2 = vol(0.23, 0.37, 0.44)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

    def test_insert_knot(self):
        # more or less random 3D volume with p=[2,2,1] and n=[4,3,2]
        controlpoints = [[0, 0, 0], [-1, 1, 0], [0, 2, 0], [1, -1, 0], [1, 0, 0], [1, 1, 0],
                         [2, 1, 0], [2, 2, 0], [2, 3, 0], [3, 0, 0], [4, 1, 0], [3, 2, 0],
                         [0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 2], [1, 0, 2], [1, 1, 2],
                         [2, 1, 2], [2, 2, 2], [2, 3, 2], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints)

        evaluation_point1 = vol(0.23, 0.37, 0.44)  # pick some evaluation point (could be anything)

        vol.insert_knot(0, .20)
        vol.insert_knot(0, .5)
        vol.insert_knot(0, .7)
        vol.insert_knot(1, .1)
        vol.insert_knot(1, 1.0 / 3)
        vol.insert_knot(2, .8)
        vol.insert_knot(2, .9)
        knot1, knot2, knot3 = vol.knots(with_multiplicities=True)
        self.assertEqual(len(knot1), 10)  # 7 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones
        self.assertEqual(len(knot3), 6)  # 4 to start with, 2 new ones
        self.assertEqual(vol.controlpoints.shape, (7, 5, 4, 3))

        evaluation_point2 = vol(0.23, 0.37, 0.44)

        # evaluation before and after insert_knot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # test a rational 3D volume
        controlpoints = [[0, 0, 1, 1], [-1, 1, .96, 1], [0, 2, 1, 1], [1, -1, 1, 1], [1, 0, .8, 1],
                         [1, 1, 1, 1], [2, 1, .89, 1], [2, 2, .9, 1], [2, 3, 1, 1], [3, 0, 1, 1],
                         [4, 1, 1, 1], [3, 2, 1, 1], [0, 0, 1, 2], [-1, 1, .7, 2], [0, 2, 1.3, 2],
                         [1, -1, 1, 2], [1, 0, .77, 2], [1, 1, 1, 2], [2, 1, .89, 1], [2, 2, .8, 4],
                         [2, 3, 1, 1], [3, 0, 1, 1], [4, 1, 1, 1], [3, 2, 1, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints, True)

        evaluation_point1 = vol(0.23, 0.37, 0.44)

        vol.insert_knot(0, [.20, .5, .7])
        vol.insert_knot(1, [.1, 1.0 / 3])
        vol.insert_knot(2, [.8, .9])
        knot1, knot2, knot3 = vol.knots(with_multiplicities=True)
        self.assertEqual(len(knot1), 10)  # 7 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones
        self.assertEqual(len(knot3), 6)  # 4 to start with, 2 new ones
        self.assertEqual(vol.controlpoints.shape, (7, 5, 4, 4))

        evaluation_point2 = vol(0.23, 0.37, 0.44)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

    def test_force_rational(self):
        # more or less random 3D volume with p=[3,2,1] and n=[4,3,2]
        controlpoints = [[0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 2], [1, 0, 2], [1, 1, 2],
                         [2, 1, 2], [2, 2, 2], [2, 3, 2], [3, 0, 0], [4, 1, 0], [3, 2, 0],
                         [0, 0, 3], [-1, 1, 3], [0, 2, 3], [1, -1, 5], [1, 0, 5], [1, 1, 5],
                         [2, 1, 4], [2, 2, 4], [2, 3, 4], [3, 0, 2], [4, 1, 2], [3, 2, 2]]
        basis1 = BSplineBasis(4, [0, 0, 0, 0, 2, 2, 2, 2])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints)

        evaluation_point1 = vol(0.23, .66, .32)
        control_point1 = vol[0]
        vol.force_rational()
        evaluation_point2 = vol(0.23, .66, .32)
        control_point2 = vol[0]
        # ensure that volume has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        # ensure that we include rational weights of 1
        self.assertEqual(len(control_point1), 3)
        self.assertEqual(len(control_point2), 4)
        self.assertEqual(control_point2[3], 1)
        self.assertEqual(vol.rational, True)

    def test_swap_parametrization(self):
        # more or less random 3D volume with p=[3,2,1] and n=[4,3,2]
        controlpoints = [[0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 2], [1, 0, 2], [1, 1, 2],
                         [2, 1, 2], [2, 2, 2], [2, 3, 2], [3, 0, 0], [4, 1, 0], [3, 2, 0],
                         [0, 0, 3], [-1, 1, 3], [0, 2, 3], [1, -1, 5], [1, 0, 5], [1, 1, 5],
                         [2, 1, 4], [2, 2, 4], [2, 3, 4], [3, 0, 2], [4, 1, 2], [3, 2, 2]]
        basis1 = BSplineBasis(4, [0, 0, 0, 0, 2, 2, 2, 2])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints)

        evaluation_point1 = vol(0.23, .56, .12)
        control_point1 = vol[1]  # this is control point i=(1,0,0), when n=(4,3,2)
        self.assertEqual(vol.order(), (4, 3, 2))
        vol.swap_parametrization(0, 1)
        evaluation_point2 = vol(0.56, .23, .12)
        control_point2 = vol[3]  # this is control point i=(0,1,0), when n=(3,4,2)
        self.assertEqual(vol.order(), (3, 4, 2))

        # ensure that volume has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # check that the control points have re-ordered themselves
        self.assertEqual(control_point1[0], control_point2[0])
        self.assertEqual(control_point1[1], control_point2[1])
        self.assertEqual(control_point1[2], control_point2[2])

        vol.swap_parametrization(1, 2)
        evaluation_point3 = vol(.56, .12, .23)
        control_point3 = vol[6]  # this is control point i=(0,0,1), when n=(3,2,4)
        self.assertEqual(vol.order(), (3, 2, 4))

        # ensure that volume has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point3[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point3[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point3[2])

        # check that the control points have re-ordered themselves
        self.assertEqual(control_point1[0], control_point3[0])
        self.assertEqual(control_point1[1], control_point3[1])
        self.assertEqual(control_point1[2], control_point3[2])

    def test_split(self):
        # more or less random 3D volume with p=[3,2,1] and n=[4,3,2]
        controlpoints = [[0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 2], [1, 0, 2], [1, 1, 2],
                         [2, 1, 2], [2, 2, 2], [2, 3, 2], [3, 0, 0], [4, 1, 0], [3, 2, 0],
                         [0, 0, 3], [-1, 1, 3], [0, 2, 3], [1, -1, 5], [1, 0, 5], [1, 1, 5],
                         [2, 1, 4], [2, 2, 4], [2, 3, 4], [3, 0, 2], [4, 1, 2], [3, 2, 2]]
        basis1 = BSplineBasis(4, [0, 0, 0, 0, 2, 2, 2, 2])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis3 = BSplineBasis(2, [0, 0, 1, 1])
        vol = Volume(basis1, basis2, basis3, controlpoints)
        split_u_vol = vol.split(0, [.1, .2, .3, .4, .5, .6])
        split_v_vol = vol.split(1, .1)
        split_w_vol = vol.split(2, [.4, .5, .6])

        self.assertEqual(len(split_u_vol), 7)
        self.assertEqual(len(split_v_vol), 2)
        self.assertEqual(len(split_w_vol), 4)

        # check that the u-vector is properly split
        self.assertAlmostEqual(split_u_vol[0].start()[0], 0.0)
        self.assertAlmostEqual(split_u_vol[0].end()[0], 0.1)
        self.assertAlmostEqual(split_u_vol[1].start()[0], 0.1)
        self.assertAlmostEqual(split_u_vol[1].end()[0], 0.2)
        self.assertAlmostEqual(split_u_vol[2].start()[0], 0.2)
        self.assertAlmostEqual(split_u_vol[2].end()[0], 0.3)
        self.assertAlmostEqual(split_u_vol[6].start()[0], 0.6)
        self.assertAlmostEqual(split_u_vol[6].end()[0], 2.0)
        # check that the other vectors remain unchanged
        self.assertAlmostEqual(split_u_vol[2].start()[1], 0.0)
        self.assertAlmostEqual(split_u_vol[2].end()[1], 1.0)
        self.assertAlmostEqual(split_u_vol[2].start()[2], 0.0)
        self.assertAlmostEqual(split_u_vol[2].end()[2], 1.0)
        # check that the v-vector is properly split
        self.assertAlmostEqual(split_v_vol[0].start()[1], 0.0)
        self.assertAlmostEqual(split_v_vol[0].end()[1], 0.1)
        self.assertAlmostEqual(split_v_vol[1].start()[1], 0.1)
        self.assertAlmostEqual(split_v_vol[1].end()[1], 1.0)
        # check that the others remain unchanged
        self.assertAlmostEqual(split_v_vol[1].start()[0], 0.0)
        self.assertAlmostEqual(split_v_vol[1].end()[0], 2.0)
        self.assertAlmostEqual(split_v_vol[1].start()[2], 0.0)
        self.assertAlmostEqual(split_v_vol[1].end()[2], 1.0)
        # check that the w-vector is properly split
        self.assertAlmostEqual(split_w_vol[1].start()[2], 0.4)
        self.assertAlmostEqual(split_w_vol[1].end()[2], 0.5)
        self.assertAlmostEqual(split_w_vol[2].start()[2], 0.5)
        self.assertAlmostEqual(split_w_vol[2].end()[2], 0.6)
        # check that the others remain unchanged
        self.assertAlmostEqual(split_w_vol[1].start()[0], 0.0)
        self.assertAlmostEqual(split_w_vol[1].end()[0], 2.0)
        self.assertAlmostEqual(split_w_vol[1].start()[1], 0.0)
        self.assertAlmostEqual(split_w_vol[1].end()[1], 1.0)

        # check that evaluations remain unchanged
        pt1 = vol(0.23, 0.12, 0.3)

        self.assertAlmostEqual(split_u_vol[2].evaluate(.23, .12, .3)[0], pt1[0])
        self.assertAlmostEqual(split_u_vol[2].evaluate(.23, .12, .3)[1], pt1[1])
        self.assertAlmostEqual(split_u_vol[2].evaluate(.23, .12, .3)[2], pt1[2])

        self.assertAlmostEqual(split_v_vol[1].evaluate(.23, .12, .3)[0], pt1[0])
        self.assertAlmostEqual(split_v_vol[1].evaluate(.23, .12, .3)[1], pt1[1])
        self.assertAlmostEqual(split_v_vol[1].evaluate(.23, .12, .3)[2], pt1[2])

        self.assertAlmostEqual(split_w_vol[0].evaluate(.23, .12, .3)[0], pt1[0])
        self.assertAlmostEqual(split_w_vol[0].evaluate(.23, .12, .3)[1], pt1[1])
        self.assertAlmostEqual(split_w_vol[0].evaluate(.23, .12, .3)[2], pt1[2])

    def test_reparam(self):
        # identity mapping, control points generated from knot vector
        basis1 = BSplineBasis(4, [2,2,2,2,3,6,7,7,7,7])
        basis2 = BSplineBasis(3, [-3,-3,-3,20,30,31,31,31])
        basis3 = BSplineBasis(5, [0,0,0,0,0,8,8,8,8,8])
        vol = Volume(basis1, basis2, basis3)

        self.assertAlmostEqual(vol.start(0),  2)
        self.assertAlmostEqual(vol.end(0),    7)
        self.assertAlmostEqual(vol.start(1), -3)
        self.assertAlmostEqual(vol.end(1),   31)
        self.assertAlmostEqual(vol.start(2),  0)
        self.assertAlmostEqual(vol.end(2),    8)

        vol.reparam((4,10), (0,9), (2,3))
        self.assertAlmostEqual(vol.start(0),  4)
        self.assertAlmostEqual(vol.end(0),   10)
        self.assertAlmostEqual(vol.start(1),  0)
        self.assertAlmostEqual(vol.end(1),    9)
        self.assertAlmostEqual(vol.start(2),  2)
        self.assertAlmostEqual(vol.end(2),    3)

        vol.reparam((5,11), direction=0)
        self.assertAlmostEqual(vol.start(0),  5)
        self.assertAlmostEqual(vol.end(0),   11)
        self.assertAlmostEqual(vol.start(1),  0)
        self.assertAlmostEqual(vol.end(1),    9)
        self.assertAlmostEqual(vol.start(2),  2)
        self.assertAlmostEqual(vol.end(2),    3)

        vol.reparam((5,11), direction=1)
        self.assertAlmostEqual(vol.start(0),  5)
        self.assertAlmostEqual(vol.end(0),   11)
        self.assertAlmostEqual(vol.start(1),  5)
        self.assertAlmostEqual(vol.end(1),   11)
        self.assertAlmostEqual(vol.start(2),  2)
        self.assertAlmostEqual(vol.end(2),    3)

        vol.reparam((5,11), direction=2)
        self.assertAlmostEqual(vol.start(0),  5)
        self.assertAlmostEqual(vol.end(0),   11)
        self.assertAlmostEqual(vol.start(1),  5)
        self.assertAlmostEqual(vol.end(1),   11)
        self.assertAlmostEqual(vol.start(2),  5)
        self.assertAlmostEqual(vol.end(2),   11)

        vol.reparam((-9,9))
        self.assertAlmostEqual(vol.start(0), -9)
        self.assertAlmostEqual(vol.end(0),    9)
        self.assertAlmostEqual(vol.start(1),  0)
        self.assertAlmostEqual(vol.end(1),    1)
        self.assertAlmostEqual(vol.start(2),  0)
        self.assertAlmostEqual(vol.end(2),    1)

        vol.reparam()
        self.assertAlmostEqual(vol.start(0),  0)
        self.assertAlmostEqual(vol.end(0),    1)
        self.assertAlmostEqual(vol.start(1),  0)
        self.assertAlmostEqual(vol.end(1),    1)
        self.assertAlmostEqual(vol.start(2),  0)
        self.assertAlmostEqual(vol.end(2),    1)

        vol.reparam((4,10), (0,9), (2,7))
        vol.reparam(direction=1)
        self.assertAlmostEqual(vol.start(0),  4)
        self.assertAlmostEqual(vol.end(0),   10)
        self.assertAlmostEqual(vol.start(1),  0)
        self.assertAlmostEqual(vol.end(1),    1)
        self.assertAlmostEqual(vol.start(2),  2)
        self.assertAlmostEqual(vol.end(2),    7)

if __name__ == '__main__':
    unittest.main()
