# -*- coding: utf-8 -*-

from splipy import SplineObject, Curve, BSplineBasis
import splipy.curve_factory as CurveFactory
from math import sqrt, pi
import numpy as np
import unittest

class TestCurve(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        crv = Curve(controlpoints=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        val = crv(0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(crv.dimension, 3)

        # test 2D constructor
        crv2 = Curve()
        val = crv(0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(crv2.dimension, 2)

    def test_evaluate(self):
        # create the mapping
        # x(t) = 2t + 1
        # y(t) = 2t(1-t)
        # z(t) = 0
        controlpoints = [[1, 0, 0], [2, 1, 0], [3, 0, 0]]
        crv = Curve(BSplineBasis(3), controlpoints)

        # startpoint evaluation
        val = crv(0.0)
        self.assertAlmostEqual(val[0], 1.0)
        self.assertAlmostEqual(val[1], 0.0)
        self.assertAlmostEqual(val[2], 0.0)

        # inner evaluation
        val = crv(0.4)
        self.assertAlmostEqual(val[0], 1.8)
        self.assertAlmostEqual(val[1], 0.48)
        self.assertAlmostEqual(val[2], 0.0)

        # endpoint evaluation
        val = crv(1.0)
        self.assertAlmostEqual(val[0], 3.0)
        self.assertAlmostEqual(val[1], 0.0)
        self.assertAlmostEqual(val[2], 0.0)

        # test evaluation at multiple points
        val = crv([0.0, 0.4, 0.8, 1.0])
        self.assertEqual(len(val.shape), 2)  # return matrix
        self.assertEqual(val.shape[0], 4)  # 4 evaluation points
        self.assertEqual(val.shape[1], 3)  # (x,y,z) results
        self.assertAlmostEqual(val[0, 0], 1.0)
        self.assertAlmostEqual(val[0, 1], 0.0)
        self.assertAlmostEqual(val[0, 2], 0.0)  # startpt evaluation
        self.assertAlmostEqual(val[1, 0], 1.8)
        self.assertAlmostEqual(val[1, 1], 0.48)
        self.assertAlmostEqual(val[1, 2], 0.0)  # inner evaluation
        self.assertAlmostEqual(val[3, 0], 3.0)
        self.assertAlmostEqual(val[3, 1], 0.0)
        self.assertAlmostEqual(val[3, 2], 0.0)  # endpt evaluation

        # test errors and exceptions
        with self.assertRaises(ValueError):
            val = crv(-10)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = crv(+10)  # evalaute outside parametric domain

    def test_flip_parametrization(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0, 0, 0], [1, 1, 0], [2, -1, 0], [3, 0, 0]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, .3, 1, 1, 1]), controlpoints)

        p1 = crv(0.23)
        crv.reverse()
        p2 = crv(0.77)
        self.assertAlmostEqual(p1[0], p2[0])
        self.assertAlmostEqual(p1[1], p2[1])
        self.assertAlmostEqual(p1[2], p2[2])

    def test_force_rational(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0, 0, 0], [1, 1, 0], [2, -1, 0], [3, 0, 0]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, .3, 1, 1, 1]), controlpoints)

        evaluation_point1 = crv(0.23)
        control_point1 = crv[0]
        crv.force_rational()
        evaluation_point2 = crv(0.23)
        control_point2 = crv[0]
        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        # ensure that we include rational weights of 1
        self.assertEqual(len(control_point1), 3)
        self.assertEqual(len(control_point2), 4)
        self.assertEqual(control_point2[3], 1)
        self.assertEqual(crv.rational, True)

    def test_raise_order(self):
        # non-uniform knot vector of a squiggly quadratic n=5 curve in 3D
        controlpoints = [[0, 0, 0], [1, 1, 1], [2, -1, 0], [3, 0, -1], [0, 0, -5]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, .3, .4, 1, 1, 1]), controlpoints)

        evaluation_point1 = crv(0.37)
        crv.raise_order(2)
        evaluation_point2 = crv(0.37)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # ensure that curve has the right order
        self.assertEqual(crv.order(0), 5)

        # check integer type for argument
        with self.assertRaises(TypeError):
            crv.raise_order(0.5)

        # check logic error for negative argument
        with self.assertRaises(Exception):
            crv.raise_order(-1)

    def test_lower_order(self):
        basis = BSplineBasis(4, [0,0,0,0,.2, .3, .3, .6, .9, 1,1,1,1])
        interp_pts = basis.greville()
        x = [[t*(1-t), t**2] for t in interp_pts]

        crv = CurveFactory.interpolate(x, basis) # function in space, exact representation kept

        crv2 = crv.lower_order(1) # still in space, crv2 is *also* exact

        t = np.linspace(0,1, 13)
        self.assertTrue(np.allclose( crv(t), crv2(t)) )
        self.assertEqual(crv.order(0),  4)
        self.assertEqual(crv2.order(0), 3)
        self.assertEqual(crv.continuity(0.3),  1)
        self.assertEqual(crv.continuity(0.6),  2)
        self.assertEqual(crv2.continuity(0.3), 1)
        self.assertEqual(crv2.continuity(0.6), 1)


    def test_insert_knot(self):
        # non-uniform knot vector of a squiggly quadratic n=5 curve in 3D
        controlpoints = [[0, 0, 0], [1, 1, 1], [2, -1, 0], [3, 0, -1], [0, 0, -5]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, .3, .4, 1, 1, 1]), controlpoints)

        evaluation_point1 = crv(0.37)
        crv.insert_knot(.2)
        crv.insert_knot(.3)
        crv.insert_knot(.9)
        evaluation_point2 = crv(0.37)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        self.assertEqual(len(crv.knots(0, with_multiplicities=True)), 11)

        # test knot insertion on single knot span
        crv = Curve(BSplineBasis(5), [[0, 0, 0], [1, 1, 1], [2, -1, 0], [3, 0, -1], [0, 0, -5]])
        evaluation_point1 = crv(0.27)
        crv.insert_knot(.2)
        evaluation_point2 = crv(0.27)
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # test knot insertion on first of two-knot span
        crv = Curve(
            BSplineBasis(3, [0, 0, 0, .5, 1, 1, 1]), [[0, 0, 0], [2, -1, 0], [3, 0, -1], [0, 0, -5]
                                                      ])
        evaluation_point1 = crv(0.27)
        crv.insert_knot(.2)
        evaluation_point2 = crv(0.27)
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # test knot insertion on last of two-knot span
        crv = Curve(
            BSplineBasis(3, [0, 0, 0, .5, 1, 1, 1]), [[0, 0, 0], [2, -1, 0], [3, 0, -1], [0, 0, -5]
                                                      ])
        evaluation_point1 = crv(0.27)
        crv.insert_knot(.9)
        evaluation_point2 = crv(0.27)
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # test knot insertion down to C0 basis
        crv = Curve(BSplineBasis(3), [[0, 0, 0], [2, -1, 0], [0, 0, -5]])
        evaluation_point1 = crv(0.27)
        crv.insert_knot(.4)
        crv.insert_knot(.4)
        evaluation_point2 = crv(0.27)
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # test rational curves, here a perfect circle represented as n=9, p=2-curve
        s = 1.0 / sqrt(2)
        controlpoints = [[1, 0, 1], [s, s, s], [0, 1, 1], [-s, s, s], [-1, 0, 1], [-s, -s, s],
                         [0, -1, 1], [s, -s, s], [1, 0, 1]]
        crv = Curve(BSplineBasis(3, [-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5]), controlpoints, True)
        evaluation_point1 = crv(0.37)
        crv.insert_knot(.2)
        crv.insert_knot(.3)
        crv.insert_knot(.9)
        evaluation_point2 = crv(0.37)
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertEqual(len(crv.knots(0, with_multiplicities=True)), 15)

        # test errors and exceptions
        with self.assertRaises(TypeError):
            crv.insert_knot(1, 2, 3)  # too many arguments
        with self.assertRaises(ValueError):
            crv.insert_knot(1, 2)  # direction=2 is illegal for curves
        with self.assertRaises(TypeError):
            crv.insert_knot()  # too few arguments
        with self.assertRaises(ValueError):
            crv.insert_knot(-0.2)  # Outside-domain error
        with self.assertRaises(ValueError):
            crv.insert_knot(4.4)  # Outside-domain error

    def test_reparam(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0, 0, 0], [1, 1, 0], [2, -1, 0], [3, 0, 0]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, 1.32, 3, 3, 3]), controlpoints)

        # get some info on the initial curve
        knots1 = crv.knots(0)
        evaluation_point1 = crv(1.20)
        self.assertEqual(knots1[0], 0)
        self.assertEqual(knots1[-1], 3)

        # reparametrize
        crv.reparam((6.0, 9.0))

        # get some info on the reparametrized curve
        knots2 = crv.knots(0)
        evaluation_point2 = crv(7.20)
        self.assertEqual(knots2[0], 6)
        self.assertEqual(knots2[-1], 9)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # normalize, i.e. set domain to [0,1]
        crv.reparam()

        # get some info on the normalized curve
        knots3 = crv.knots(0)
        evaluation_point3 = crv(0.40)
        self.assertEqual(knots3[0], 0)
        self.assertEqual(knots3[-1], 1)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point3[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point3[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point3[2])

        # test errors and exceptions
        with self.assertRaises(ValueError):
            crv.reparam((9, 3))
        with self.assertRaises(TypeError):
            crv.reparam(("one", "two"))

    def test_split(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0, 0, 1], [1, 1, 0], [2, -1, 0], [3, 0, 0]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, .7, 1, 1, 1]), controlpoints)

        # get some info on the initial curve
        evaluation_point1 = crv(0.50)
        evaluation_point2 = crv(0.70)
        evaluation_point3 = crv(0.33)

        # split curves away from knot
        new_curves_050 = crv.split(0.50)
        self.assertEqual(len(new_curves_050), 2)
        self.assertEqual(
            len(new_curves_050[0].knots(0, with_multiplicities=True)), 6)  # open knot vector [0,0,0,.5,.5,.5]
        self.assertEqual(
            len(new_curves_050[1].knots(0, with_multiplicities=True)), 7)  # open knot vector [.5,.5,.5,.7,1,1,1]

        # split curves at existing knot
        new_curves_070 = crv.split(0.70)
        self.assertEqual(len(new_curves_070), 2)
        self.assertEqual(
            len(new_curves_070[0].knots(0, with_multiplicities=True)), 6)  # open knot vector [0,0,0,.7,.7,.7]
        self.assertEqual(
            len(new_curves_070[1].knots(0, with_multiplicities=True)), 6)  # open knot vector [.7,.7,.7,1,1,1]

        # split curves multiple points
        new_curves_all = crv.split([0.50, 0.70])
        self.assertEqual(len(new_curves_all), 3)
        self.assertEqual(
            len(new_curves_all[0].knots(0, with_multiplicities=True)), 6)  # open knot vector [0,0,0,.5,.5,.5]
        self.assertEqual(
            len(new_curves_all[1].knots(0, with_multiplicities=True)), 6)  # open knot vector [.5,.5,.5,.7,.7,.7]
        self.assertEqual(
            len(new_curves_all[2].knots(0, with_multiplicities=True)), 6)  # open knot vector [.7,.7,.7,1,1,1]

        # compare all curves which exist at parametric point 0.5
        for c in new_curves_050 + [new_curves_070[0]] + new_curves_all[0:2]:
            new_curve_evaluation = c(0.50)
            self.assertAlmostEqual(evaluation_point1[0], new_curve_evaluation[0])
            self.assertAlmostEqual(evaluation_point1[1], new_curve_evaluation[1])
            self.assertAlmostEqual(evaluation_point1[2], new_curve_evaluation[2])

            # compare all curves which exist at parametric point 0.33
        for c in [new_curves_050[0]] + [new_curves_070[0]] + [new_curves_all[0]]:
            new_curve_evaluation = c(0.33)
            self.assertAlmostEqual(evaluation_point3[0], new_curve_evaluation[0])
            self.assertAlmostEqual(evaluation_point3[1], new_curve_evaluation[1])
            self.assertAlmostEqual(evaluation_point3[2], new_curve_evaluation[2])

        # compare all curves which exist at parametric point 0.7
        for c in [new_curves_050[1]] + new_curves_070 + new_curves_all[1:3]:
            new_curve_evaluation = c(0.70)
            self.assertAlmostEqual(evaluation_point2[0], new_curve_evaluation[0])
            self.assertAlmostEqual(evaluation_point2[1], new_curve_evaluation[1])
            self.assertAlmostEqual(evaluation_point2[2], new_curve_evaluation[2])

            # test errors and exceptions
        with self.assertRaises(TypeError):
            crv.split(.1, .2, .3)  # too many arguments
        with self.assertRaises(Exception):
            crv.split(-0.2)  # GoTools returns error on outside-domain errors
        with self.assertRaises(Exception):
            crv.split(1.4)  # GoTools returns error on outside-domain errors

    def test_periodic_split(self):
        # non-uniform rational knot vector of a periodic cubic n=7 curve
        controlpoints = [[1, 0, 1], [1, 1, .7], [0, 1, .89], [-1, 1, 0.5], [-1, 0, 1], [-1,-.5,1], [1, -.5, 1]]
        basis = BSplineBasis(4, [-3, -2, -1, 0, 1, 2, 2.5, 4, 5, 6, 7, 8, 9, 9.5], 2)
        crv = Curve(basis, controlpoints, rational=True)
        crv2 = crv.split(1)   # split at knot value
        crv3 = crv.split(6.5) # split outside existing knot

        self.assertEqual(len(crv),   7)
        self.assertEqual(len(crv2), 10)
        self.assertEqual(len(crv3), 11)

        self.assertEqual(crv.periodic(),  True)
        self.assertEqual(crv2.periodic(), False)
        self.assertEqual(crv3.periodic(), False)

        t = np.linspace(6.5, 8, 13) # domain where all parameter values are the same
        pt  = crv( t)
        pt2 = crv2(t)
        pt3 = crv3(t)

        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(pt-pt3), 0.0)

    def test_center(self):
        # create the geometric mapping x(t) = t, y(t) = t^3,  t=[0,1]
        cp = [[0, 0], [1.0/3, 0], [2.0/3, 0], [1.0, 1]]
        basis = BSplineBasis(4)
        crv = Curve(basis, cp)
        center = crv.center()

        self.assertAlmostEqual(center[0], .5)
        self.assertAlmostEqual(center[1], .25)

    def test_derivative(self):
        # testing the parametrization x(t) = [.5*x^3 / ((1-x)^3+.5*x^3), 0]
        cp = [[0,0,1], [0,0,0], [0,0,0], [.5,0,.5]]
        crv = Curve(BSplineBasis(4), cp, rational=True)
        def expect_derivative(x):
            return 6*(1-x)**2*x**2/(x**3 - 6*x**2 + 6*x - 2)**2
        def expect_derivative_2(x):
            return -12*x*(x**5 - 3*x**4 + 2*x**3 + 4*x**2 - 6*x + 2)/(x**3 - 6*x**2 + 6*x - 2)**3

        # insert a few more knots to spice things up
        crv.insert_knot([.2, .71])

        self.assertAlmostEqual(crv.derivative(0.32, 1)[0], expect_derivative(0.32))
        self.assertAlmostEqual(crv.derivative(0.32, 1)[1], 0)
        self.assertAlmostEqual(crv.derivative(0.71, 1)[0], expect_derivative(0.71))

        self.assertAlmostEqual(crv.derivative(0.22, 2)[0], expect_derivative_2(0.22))
        self.assertAlmostEqual(crv.derivative(0.22, 2)[1], 0)
        self.assertAlmostEqual(crv.derivative(0.86, 2)[0], expect_derivative_2(0.86))

    def test_tangent_and_normal(self):
        crv = CurveFactory.circle()
        crv.set_dimension(3)
        t = np.linspace(crv.start(0), crv.end(0), 13)
        X = crv(t)
        T = crv.tangent(t)
        B = crv.binormal(t)
        N = crv.normal(t)

        # check correct size
        self.assertEqual(len(B.shape), 2) # returns a matrix
        self.assertEqual(B.shape[0],  13)
        self.assertEqual(B.shape[1],   3)

        # check all evaluation points
        for (x,t,b,n) in zip(X,T,B,N):
            # for the circle we have x*t=0 and x=-n
            self.assertAlmostEqual(np.dot(x,t), 0.0)
            self.assertTrue(np.allclose( x, -n) )
            # for all curves we have that t,b,n are orthogonal
            self.assertAlmostEqual(np.dot(t,b), 0.0)
            self.assertAlmostEqual(np.dot(b,n), 0.0)
            self.assertAlmostEqual(np.dot(n,t), 0.0)
            # for planar (2D) curves, we have b=[0,0,1]
            self.assertTrue(np.allclose( b, [0,0,1]) )

        # check that evaluations work for single-valued input
        t = crv.tangent(.23)
        b = crv.binormal(.23)
        n = crv.normal(.23)
        self.assertEqual(len(t.shape), 1) # is a vector (not matrix)
        self.assertEqual(len(b.shape), 1)
        self.assertEqual(len(n.shape), 1)
        self.assertAlmostEqual(np.linalg.norm(t), 1.0) # should be normalized
        self.assertAlmostEqual(np.linalg.norm(b), 1.0)
        self.assertAlmostEqual(np.linalg.norm(n), 1.0)

    def test_append(self):
        crv  = Curve(BSplineBasis(3), [[0,0], [1,0], [0,1]])
        crv2 = Curve(BSplineBasis(4), [[0,1,0], [0,1,1], [0,2,1], [0,2,2]])
        crv2.insert_knot(0.5)

        crv3 = crv.clone()
        crv3.append(crv2)

        expected_knot = [0,0,0,0,1,1,1,1.5,2,2,2,2]
        self.assertEqual(crv3.order(direction=0),  4)
        self.assertEqual(crv3.rational, False)
        self.assertAlmostEqual(np.linalg.norm(crv3.knots(0, True)-expected_knot), 0.0)

        t = np.linspace(0,1,11)
        pt        = np.zeros((11,3))
        pt[:,:-1] = crv(t)
        pt2       = crv3(t)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)
        pt  = crv2(t)
        pt2 = crv3(t+1.0)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

    def test_length(self):
        crv = Curve()
        self.assertAlmostEqual(crv.length(), 1.0)
        crv = Curve(BSplineBasis(2, [-1,-1,1,2,3,3]), [[0,0,0], [1,0,0], [1,0,3],[1,10,3]])
        self.assertAlmostEqual(crv.length(), 14.0)

    def test_make_periodic(self):
        my_cps = np.array([[0, -1], [1, 0], [0, 1], [-1, 0], [0, -1]], dtype=np.float)

        crv = Curve(BSplineBasis(2, [0, 0, 1, 2, 3, 4, 4]), my_cps, rational=False)
        crv = crv.make_periodic(0)
        cps = [[0, -1], [1, 0], [0, 1], [-1, 0]]
        self.assertAlmostEqual(np.linalg.norm(crv.controlpoints - cps), 0.0)

        crv = Curve(BSplineBasis(2, [0, 0, 1, 2, 3, 4, 4]), my_cps, rational=False)
        crv = crv.make_periodic(0)
        cps = [[0, -1], [1, 0], [0, 1], [-1, 0]]
        self.assertAlmostEqual(np.linalg.norm(crv.controlpoints - cps), 0.0)

        crv = Curve(BSplineBasis(3, [0, 0, 0, 1, 2, 3, 3, 3]), my_cps, rational=False)
        crv = crv.make_periodic(0)
        cps = [[0, -1], [1, 0], [0, 1], [-1, 0]]
        self.assertAlmostEqual(np.linalg.norm(crv.controlpoints - cps), 0.0)

        crv = Curve(BSplineBasis(3, [0, 0, 0, 1, 2, 3, 3, 3]), my_cps, rational=False)
        crv = crv.make_periodic(1)
        cps = [[-1, 0], [1, 0], [0, 1]]
        self.assertAlmostEqual(np.linalg.norm(crv.controlpoints - cps), 0.0)

    def test_make_periodic_reconstruct(self):
        orig = Curve(
            BSplineBasis(4, [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7], 2),
            [[1, 1], [2, 2], [3, 3], [4, 4]],
            rational=True,
        )
        recons = orig.split(0).make_periodic(2)

        self.assertAlmostEqual(np.linalg.norm(orig.controlpoints - recons.controlpoints), 0.0)
        self.assertAlmostEqual(np.linalg.norm(orig.bases[0].knots - recons.bases[0].knots), 0.0)

    def test_curvature(self):
        # linear curves have zero curvature
        crv = Curve()
        self.assertAlmostEqual(crv.curvature(.3), 0.0)
        # test multiple evaluation points
        t = np.linspace(0,1, 10)
        k = crv.curvature(t)
        self.assertTrue(np.allclose(k, 0.0))

        # test circle
        crv = CurveFactory.circle(r=3) + [1,1]
        t = np.linspace(0,2*pi, 10)
        k = crv.curvature(t)
        self.assertTrue(np.allclose(k, 1.0/3.0)) # circles: k = 1/r

        # test 3D (np.cross has different behaviour in 2D/3D)
        crv.set_dimension(3)
        k = crv.curvature(t)
        self.assertTrue(np.allclose(k, 1.0/3.0)) # circles: k = 1/r

    def test_torsion(self):
        # planar curves have zero torsion
        controlpoints = [[1, 0, 1], [1, 1, .7], [0, 1, .89], [-1, 1, 0.5], [-1, 0, 1], [-1,-.5,1], [1, -.5, 1]]
        basis = BSplineBasis(4, [-3, -2, -1, 0, 1, 2, 2.5, 4, 5, 6, 7, 8, 9, 9.5], 2)
        crv = Curve(basis, controlpoints, rational=True)
        self.assertAlmostEqual(crv.torsion(.3), 0.0)

        # test multiple evaluation points
        t = np.linspace(0,1, 10)
        k = crv.torsion(t)
        self.assertEqual(len(k.shape), 1)    # is vector
        self.assertEqual(k.shape[0],  10)    # length 10 vector
        self.assertTrue(np.allclose(k, 0.0)) # zero torsion for planar curves

        # test helix: [a*cos(t), a*sin(t), b*t]
        t = np.linspace(0, 6*pi, 300)
        a = 3.0
        b = 2.0
        x = np.array([a*np.cos(t), a*np.sin(t), b*t])
        crv = CurveFactory.cubic_curve(x.T, t=t)
        t = np.linspace(0, 6*pi, 10)
        k = crv.torsion(t)

        # this is a helix approximation, hence atol=1e-3
        self.assertTrue(np.allclose(k, b/(a**2+b**2), atol=1e-3)) # helix have const. torsion


if __name__ == '__main__':
    unittest.main()
