# -*- coding: utf-8 -*-

from splipy import BSplineBasis, Volume
import splipy.volume_factory as VolumeFactory
import unittest
import numpy as np
from math import pi, sqrt


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

    def test_evaluate_nontensor(self):
        vol = Volume(BSplineBasis(7), BSplineBasis(7), BSplineBasis(5))

        u_val = [0, 0.1, 0.9, 0.3]
        v_val = [0.2, 0.3, 0.9, 0.4]
        w_val = [0.3, 0.5, 0.5, 0.0]
        value = vol(u_val, v_val, w_val, tensor=False)

        self.assertEqual(value.shape[0], 4)
        self.assertEqual(value.shape[1], 3)

        for i, (u, v, w) in enumerate(zip(u_val, v_val, w_val)):
            self.assertAlmostEqual(value[i, 0], u)  # identity map x=u
            self.assertAlmostEqual(value[i, 1], v)  # identity map y=v
            self.assertAlmostEqual(value[i, 2], w)  # identity map z=w

    def test_indexing(self):
        v = Volume()

        self.assertEqual(v[0][0], 0.0)
        self.assertEqual(v[0][1], 0.0)
        self.assertEqual(v[0][2], 0.0)
        self.assertEqual(v[1][0], 1.0)
        self.assertEqual(v[1][1], 0.0)
        self.assertEqual(v[1][2], 0.0)
        self.assertEqual(v[-1][0], 1.0)
        self.assertEqual(v[-1][1], 1.0)
        self.assertEqual(v[-1][2], 1.0)

        self.assertEqual(v[:][0,0], 0.0)
        self.assertEqual(v[:][1,0], 1.0)
        self.assertEqual(v[:][1,1], 0.0)
        self.assertEqual(v[:][2,1], 1.0)
        self.assertEqual(v[:][2,2], 0.0)
        self.assertEqual(v[:][5,0], 1.0)
        self.assertEqual(v[:][5,1], 0.0)
        self.assertEqual(v[:][5,2], 1.0)

        self.assertEqual(v[0,0,0][0], 0.0)
        self.assertEqual(v[0,0,0][1], 0.0)
        self.assertEqual(v[0,0,0][2], 0.0)
        self.assertEqual(v[0,1,0][0], 0.0)
        self.assertEqual(v[0,1,0][1], 1.0)
        self.assertEqual(v[0,1,0][2], 0.0)
        self.assertEqual(v[0,:,1][0,0], 0.0)
        self.assertEqual(v[0,:,1][0,1], 0.0)
        self.assertEqual(v[0,:,1][0,2], 1.0)
        self.assertEqual(v[0,:,1][1,0], 0.0)
        self.assertEqual(v[0,:,1][1,1], 1.0)
        self.assertEqual(v[0,:,1][1,2], 1.0)

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

    def test_lower_order(self):
        b = BSplineBasis(4, [0,0,0,0,.2, .3, .3, .6, .9, 1,1,1,1])
        t = b.greville()
        Y, X, Z = np.meshgrid(t,t,t)
        cp = np.zeros((len(t), len(t), len(t), 3))
        cp[...,0] = X*(1-Y)
        cp[...,1] = X*X
        cp[...,2] = Z**2 + 2

        vol  = VolumeFactory.interpolate(cp, [b,b,b])
        vol2 = vol.lower_order(1) # still in space, vol2 is *also* exact
        u = np.linspace(0,1, 5)
        v = np.linspace(0,1, 6)
        w = np.linspace(0,1, 7)
        self.assertTrue(np.allclose( vol(u,v,w), vol2(u,v,w) ))
        self.assertTupleEqual(vol.order(), (4,4,4))
        self.assertTupleEqual(vol2.order(), (3,3,3))


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

        vol.insert_knot(.20,        0)
        vol.insert_knot( .5,       'u')
        vol.insert_knot( .7,        0)
        vol.insert_knot( .1,        1)
        vol.insert_knot( 1.0 / 3,   1)
        vol.insert_knot( .8,        2)
        vol.insert_knot( .9,       'W')
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

        vol.insert_knot([.20, .5, .7], 0)
        vol.insert_knot([.1, 1.0 / 3], 1)
        vol.insert_knot([.8, .9], 2)
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

    def test_swap(self):
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
        vol.swap(0, 1)
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

        vol.swap(1, 2)
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
        split_u_vol = vol.split([.1, .2, .3, .4, .5, .6], 0)
        split_v_vol = vol.split(.1, 1)
        split_w_vol = vol.split([.4, .5, .6], 2)

        self.assertEqual(len(split_u_vol), 7)
        self.assertEqual(len(split_v_vol), 2)
        self.assertEqual(len(split_w_vol), 4)

        # check that the u-vector is properly split
        self.assertAlmostEqual(split_u_vol[0].start(0), 0.0)
        self.assertAlmostEqual(split_u_vol[0].end(0),   0.1)
        self.assertAlmostEqual(split_u_vol[1].start(0), 0.1)
        self.assertAlmostEqual(split_u_vol[1].end(0),   0.2)
        self.assertAlmostEqual(split_u_vol[2].start(0), 0.2)
        self.assertAlmostEqual(split_u_vol[2].end(0),   0.3)
        self.assertAlmostEqual(split_u_vol[6].start(0), 0.6)
        self.assertAlmostEqual(split_u_vol[6].end(0),   2.0)
        # check that the other vectors remain unchanged
        self.assertAlmostEqual(split_u_vol[2].start(1), 0.0)
        self.assertAlmostEqual(split_u_vol[2].end(1),   1.0)
        self.assertAlmostEqual(split_u_vol[2].start(2), 0.0)
        self.assertAlmostEqual(split_u_vol[2].end(2),   1.0)
        # check that the v-vector is properly split
        self.assertAlmostEqual(split_v_vol[0].start(1), 0.0)
        self.assertAlmostEqual(split_v_vol[0].end(1),   0.1)
        self.assertAlmostEqual(split_v_vol[1].start(1), 0.1)
        self.assertAlmostEqual(split_v_vol[1].end(1),   1.0)
        # check that the others remain unchanged
        self.assertAlmostEqual(split_v_vol[1].start(0), 0.0)
        self.assertAlmostEqual(split_v_vol[1].end(0),   2.0)
        self.assertAlmostEqual(split_v_vol[1].start(2), 0.0)
        self.assertAlmostEqual(split_v_vol[1].end(2),   1.0)
        # check that the w-vector is properly split
        self.assertAlmostEqual(split_w_vol[1].start(2), 0.4)
        self.assertAlmostEqual(split_w_vol[1].end(2),   0.5)
        self.assertAlmostEqual(split_w_vol[2].start(2), 0.5)
        self.assertAlmostEqual(split_w_vol[2].end(2),   0.6)
        # check that the others remain unchanged
        self.assertAlmostEqual(split_w_vol[1].start(0), 0.0)
        self.assertAlmostEqual(split_w_vol[1].end(0),   2.0)
        self.assertAlmostEqual(split_w_vol[1].start(1), 0.0)
        self.assertAlmostEqual(split_w_vol[1].end(1),   1.0)

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

    def test_reverse(self):
        # identity mapping, control points generated from knot vector
        basis1 = BSplineBasis(4, [2,2,2,2,3,6,12,12,12,12])
        basis2 = BSplineBasis(3, [-3,-3,-3,20,30,33,33,33])
        basis3 = BSplineBasis(5, [0,0,0,0,0,8,8,8,8,8])
        vol = Volume(basis1, basis2, basis3)

        u = np.linspace( 2,12,5)
        v = np.linspace(-3,33,5)
        w = np.linspace( 0, 8,5)

        pt = vol(u,v,w)

        vol.reverse('v')
        pt2 = vol(u,v[::-1],w)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)
        self.assertAlmostEqual(vol.start('v'),  -3)
        self.assertAlmostEqual(vol.end('v'),    33)

        vol.reverse(2)
        pt2 = vol(u,v[::-1],w[::-1])
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)
        self.assertAlmostEqual(vol.start('w'),   0)
        self.assertAlmostEqual(vol.end('w'),     8)

    def test_make_identical(self):
        basis1 = BSplineBasis(4, [-1,-1,0,0,1,1,2,2], periodic=1)
        basis2 = BSplineBasis(3, [-1,0,0,1,1,2],      periodic=0)
        basis3 = BSplineBasis(2)
        vol1 = Volume()
        vol2 = Volume(basis1, basis2, basis3)
        vol1.refine(1)
        Volume.make_splines_identical(vol1,vol2)

        for v in (vol1, vol2):
            self.assertEqual(v.periodic(0), False)
            self.assertEqual(v.periodic(1), False)
            self.assertEqual(v.periodic(2), False)

            self.assertEqual(v.order(), (4,3,2))
            self.assertAlmostEqual(len(v.knots(0, True)), 11)
            self.assertAlmostEqual(len(v.knots(1, True)), 8)
            self.assertAlmostEqual(len(v.knots(2, True)), 5)

    def test_bounding_box(self):
        vol = Volume()
        bb = vol.bounding_box()
        self.assertAlmostEqual(bb[0][0], 0 )
        self.assertAlmostEqual(bb[0][1], 1 )
        self.assertAlmostEqual(bb[1][0], 0 )
        self.assertAlmostEqual(bb[1][1], 1 )
        self.assertAlmostEqual(bb[2][0], 0 )
        self.assertAlmostEqual(bb[2][1], 1 )

        vol.refine(2)
        vol.rotate(pi/4, [1,0,0])
        vol += (1,0,1)
        bb = vol.bounding_box()
        self.assertAlmostEqual(bb[0][0], 1 )
        self.assertAlmostEqual(bb[0][1], 2 )
        self.assertAlmostEqual(bb[1][0], -sqrt(2)/2 )
        self.assertAlmostEqual(bb[1][1],  sqrt(2)/2 )
        self.assertAlmostEqual(bb[2][0], 1 )
        self.assertAlmostEqual(bb[2][1], 1+sqrt(2) )

    def test_controlpoint_access(self):
        v = Volume()
        v.refine(1)
        self.assertAlmostEqual(v[0,0,0, 0] ,  0)
        self.assertAlmostEqual(v[0,0,0][0] ,  0)
        self.assertAlmostEqual(v[0,1,0][0] ,  0)
        self.assertAlmostEqual(v[0,1,0][1] , .5)
        self.assertAlmostEqual(v[0,1,0, 1] , .5)
        self.assertAlmostEqual(v[0,0,2][2] ,  1)
        self.assertAlmostEqual(v[4][0]     , .5)
        self.assertAlmostEqual(v[4][1]     , .5)
        self.assertAlmostEqual(v[4][2]     ,  0)
        self.assertAlmostEqual(v[13][0]    , .5)
        self.assertAlmostEqual(v[13][1]    , .5)
        self.assertAlmostEqual(v[13][2]    , .5)
        self.assertAlmostEqual(v[14][0]    ,  1)
        self.assertAlmostEqual(v[14][1]    , .5)
        self.assertAlmostEqual(v[14][2]    , .5)

        v[0]        = [.1, .1, .1]
        v[1,1,1]    = [.6, .6, .6]
        v[1,0,0][0] = .4
        self.assertAlmostEqual(v[0,0,0][0], .1)
        self.assertAlmostEqual(v[0,0,0][1], .1)
        self.assertAlmostEqual(v[0,0,0][2], .1)
        self.assertAlmostEqual(v[13][0]   , .6)
        self.assertAlmostEqual(v[13][1]   , .6)
        self.assertAlmostEqual(v[13][2]   , .6)
        self.assertAlmostEqual(v[1][0]    , .4)
        self.assertAlmostEqual(v[1][1]    ,  0)
        self.assertAlmostEqual(v[1][2]    ,  0)

        v[:,0,0]         = 13
        v[0,1,0][:]      = 12
        # v[0,2,1]         = [9,8,7]
        # v[0,2,0]         = [9,8,7]
        v[0,2,1::-1]   = [9,8,7]
        v[1,2,1::-1,:] = [[6,5,4],[3,2,1]]
        self.assertAlmostEqual(v[1,0,0,0], 13)
        self.assertAlmostEqual(v[1,0,0,1], 13)
        self.assertAlmostEqual(v[2,0,0,2], 13)
        self.assertAlmostEqual(v[0,1,0,0], 12)
        self.assertAlmostEqual(v[0,1,0,2], 12)
        self.assertAlmostEqual(v[0,2,2,1],  1)
        self.assertAlmostEqual(v[0,2,1,0],  9)
        self.assertAlmostEqual(v[0,2,1,1],  8)
        self.assertAlmostEqual(v[0,2,0,2],  7)
        self.assertAlmostEqual(v[1,2,1,0],  6)
        self.assertAlmostEqual(v[1,2,1,1],  5)
        self.assertAlmostEqual(v[1,2,1,2],  4)
        self.assertAlmostEqual(v[1,2,0,0],  3)
        self.assertAlmostEqual(v[1,2,0,1],  2)
        self.assertAlmostEqual(v[1,2,0,2],  1)

    def test_volume(self):
        v = Volume()
        self.assertAlmostEqual(v.volume(), 1.0)
        v -= (.5, .5, 0)
        v[:,:,1,0:2] = 0.0 # squeeze top together, creating a pyramid
        self.assertAlmostEqual(v.volume(), 1.0/3)

    def test_operators(self):
        v = Volume()
        v.raise_order(1,1,2)
        v.refine(3,2,1)

        # test translation operator
        v2 = v + [1,0,0]
        v3 = [1,0,0] + v
        v += [1,0,0]
        self.assertTrue(np.allclose(v2.controlpoints, v3.controlpoints))
        self.assertTrue(np.allclose(v.controlpoints,  v3.controlpoints))

        # test scaling operator
        v2 = v * 3
        v3 = 3 * v
        v *= 3
        self.assertTrue(np.allclose(v2.controlpoints, v3.controlpoints))
        self.assertTrue(np.allclose(v.controlpoints,  v3.controlpoints))

if __name__ == '__main__':
    unittest.main()
