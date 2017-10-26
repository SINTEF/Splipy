# -*- coding: utf-8 -*-

from splipy import Surface, BSplineBasis
import splipy.surface_factory as SurfaceFactory
from math import pi
import numpy as np
import unittest

class TestSurface(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        cp = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]]
        surf = Surface(controlpoints=cp)
        val = surf(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf[0]), 3)

        # test 2D constructor
        cp = [[0, 0], [1, 0], [0, 1], [1, 1]]
        surf2 = Surface(controlpoints=cp)
        val = surf2(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf2[0]), 2)

        # test rational 2D constructor
        cp = [[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
        surf3 = Surface(controlpoints=cp, rational=True)
        val = surf3(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf3[0]), 3)

        # test rational 3D constructor
        cp = [[0, 0, 0, 1], [1, 0, 0, 1], [0, 1, 0, 1], [1, 1, 0, 1]]
        surf4 = Surface(controlpoints=cp, rational=True)
        val = surf4(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf4[0]), 4)

        # test constructor with single basis
        b = BSplineBasis(4)
        surf = Surface(b,b)
        surf.insert_knot(.3, 'u') # change one, but not the other
        self.assertEqual(len(surf.knots('u')), 3)
        self.assertEqual(len(surf.knots('v')), 2)

        # TODO: Include a default constructor specifying nothing, or just polynomial degrees, or just knot vectors.
        #       This should create identity mappings

        # test errors and exceptions
        controlpoints = [[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]]
        with self.assertRaises(ValueError):
            basis1 = BSplineBasis(2, [1, 1, 0, 0])
            basis2 = BSplineBasis(2, [0, 0, 1, 1])
            surf = Surface(basis1, basis2, controlpoints)  # illegal knot vector
        with self.assertRaises(ValueError):
            basis1 = BSplineBasis(2, [0, 0, .5, 1, 1])
            basis2 = BSplineBasis(2, [0, 0, 1, 1])
            surf = Surface(basis1, basis2, controlpoints)  # too few controlpoints
        # TODO: Create fail tests for rational surfaces with weights equal to zero
        #       Create fail tests for providing too few control points
        #       Create fail tests for providing too many control points

    def test_evaluate(self):
        # knot vector [t_1, t_2, ... t_{n+p+1}]
        # polynomial degree p (order-1)
        # n basis functions N_i(t), for i=1...n
        # the power basis {1,t,t^2,t^3,...} can be expressed as:
        # 1     = sum         N_i(t)
        # t     = sum ts_i  * N_i(t)
        # t^2   = sum t2s_i * N_i(t)
        # ts_i  = sum_{j=i+1}^{i+p}   t_j / p
        # t2s_i = sum_{j=i+1}^{i+p-1} sum_{k=j+1}^{i+p} t_j*t_k / (p 2)
        # (p 2) = binomial coefficent

        # creating the mapping:
        #   x(u,v) = u^2*v + u(1-v)
        #   y(u,v) = v
        controlpoints = [[0, 0], [1.0 / 4, 0], [3.0 / 4, 0], [.75, 0], [0, 1], [0, 1], [.5, 1], [1,
                                                                                                 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .5, 1, 1, 1])
        basis2 = BSplineBasis(2, [0, 0, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)

        # call evaluation at a 5x4 grid of points
        val = surf([0, .2, .5, .6, 1], [0, .2, .4, 1])
        self.assertEqual(len(val.shape), 3)  # result should be wrapped in 3-index tensor
        self.assertEqual(val.shape[0], 5)  # 5 evaluation points in u-direction
        self.assertEqual(val.shape[1], 4)  # 4 evaluation points in v-direction
        self.assertEqual(val.shape[2], 2)  # 2 coordinates (x,y)

        # check evaluation at (0,0)
        self.assertAlmostEqual(val[0][0][0], 0.0)
        self.assertAlmostEqual(val[0][0][1], 0.0)
        # check evaluation at (.2,0)
        self.assertAlmostEqual(val[1][0][0], 0.2)
        self.assertAlmostEqual(val[1][0][1], 0.0)
        # check evaluation at (.2,.2)
        self.assertAlmostEqual(val[1][1][0], 0.168)
        self.assertAlmostEqual(val[1][1][1], 0.2)
        # check evaluation at (.5,.4)
        self.assertAlmostEqual(val[2][2][0], 0.4)
        self.assertAlmostEqual(val[2][2][1], 0.4)
        # check evaluation at (.6,1)
        self.assertAlmostEqual(val[3][3][0], 0.36)
        self.assertAlmostEqual(val[3][3][1], 1)

        # test errors and exceptions
        with self.assertRaises(ValueError):
            val = surf(-10, .5)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf(+10, .3)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf(.5, -10)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf(.5, +10)  # evalaute outside parametric domain

    def test_derivative(self):
        # knot vector [t_1, t_2, ... t_{n+p+1}]
        # polynomial degree p (order-1)
        # n basis functions N_i(t), for i=1...n
        # the power basis {1,t,t^2,t^3,...} can be expressed as:
        # 1     = sum         N_i(t)
        # t     = sum ts_i  * N_i(t)
        # t^2   = sum t2s_i * N_i(t)
        # ts_i  = sum_{j=i+1}^{i+p}   t_j / p
        # t2s_i = sum_{j=i+1}^{i+p-1} sum_{k=j+1}^{i+p} t_j*t_k / (p 2)
        # (p 2) = binomial coefficent

        # creating the mapping:
        #   x(u,v) = u^2*v + u(1-v)
        #   y(u,v) = v
        controlpoints = [[0, 0], [1.0 / 4, 0], [3.0 / 4, 0], [.75, 0], [0, 1], [0, 1], [.5, 1], [1,
                                                                                                 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .5, 1, 1, 1])
        basis2 = BSplineBasis(2, [0, 0, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)

        # call evaluation at a 5x4 grid of points
        val = surf.derivative([0, .2, .5, .6, 1], [0, .2, .4, 1], d=(1, 0))
        self.assertEqual(len(val.shape), 3)  # result should be wrapped in 3-index tensor
        self.assertEqual(val.shape[0], 5)  # 5 evaluation points in u-direction
        self.assertEqual(val.shape[1], 4)  # 4 evaluation points in v-direction
        self.assertEqual(val.shape[2], 2)  # 2 coordinates (x,y)

        self.assertAlmostEqual(surf.derivative(.2, .2, d=(1, 0))[0], .88)  # dx/du=2uv+(1-v)
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(1, 0))[1], 0)  # dy/du=0
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(0, 1))[0], -.16)  # dx/dv=u^2-u
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(0, 1))[1], 1)  # dy/dv=1
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(1, 1))[0], -.60)  # d2x/dudv=2u-1
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(2, 0))[0], 0.40)  # d2x/dudu=2v
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(3, 0))[0], 0.00)  # d3x/du3=0
        self.assertAlmostEqual(surf.derivative(.2, .2, d=(0, 2))[0], 0.00)  # d2y/dv2=0

        # test errors and exceptions
        with self.assertRaises(ValueError):
            val = surf.derivative(-10, .5)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf.derivative(+10, .3)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf.derivative(.5, -10)  # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf.derivative(.5, +10)  # evalaute outside parametric domain

    def test_raise_order(self):
        # more or less random 2D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0, 0], [-1, 1], [0, 2], [1, -1], [1, 0], [1, 1], [2, 1], [2, 2], [2, 3],
                         [3, 0], [4, 1], [3, 2]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)

        self.assertEqual(surf.order()[0], 3)
        self.assertEqual(surf.order()[1], 3)
        evaluation_point1 = surf(0.23, 0.37)  # pick some evaluation point (could be anything)

        surf.raise_order(1, 2)

        self.assertEqual(surf.order()[0], 4)
        self.assertEqual(surf.order()[1], 5)
        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

        # test a rational 2D surface
        controlpoints = [[0, 0, 1], [-1, 1, .96], [0, 2, 1], [1, -1, 1], [1, 0, .8], [1, 1, 1],
                         [2, 1, .89], [2, 2, .9], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints, True)

        self.assertEqual(surf.order()[0], 3)
        self.assertEqual(surf.order()[1], 3)
        evaluation_point1 = surf(0.23, 0.37)

        surf.raise_order(1, 2)

        self.assertEqual(surf.order()[0], 4)
        self.assertEqual(surf.order()[1], 5)
        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

    def test_insert_knot(self):
        # more or less random 2D surface with p=[3,2] and n=[4,3]
        controlpoints = [[0, 0], [-1, 1], [0, 2], [1, -1], [1, 0], [1, 1], [2, 1], [2, 2], [2, 3],
                         [3, 0], [4, 1], [3, 2]]
        basis1 = BSplineBasis(4, [0, 0, 0, 0, 2, 2, 2, 2])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)

        # pick some evaluation point (could be anything)
        evaluation_point1 = surf(0.23, 0.37)

        surf.insert_knot(.22,     0)
        surf.insert_knot(.5,      0)
        surf.insert_knot(.7,      0)
        surf.insert_knot(.1,      1)
        surf.insert_knot(1.0 / 3, 1)
        knot1, knot2 = surf.knots(with_multiplicities=True)
        self.assertEqual(len(knot1), 11)  # 8 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones

        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after InsertKnot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

        # test a rational 2D surface
        controlpoints = [[0, 0, 1], [-1, 1, .96], [0, 2, 1], [1, -1, 1], [1, 0, .8], [1, 1, 1],
                         [2, 1, .89], [2, 2, .9], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints, True)

        evaluation_point1 = surf(0.23, 0.37)

        surf.insert_knot(.22,     0)
        surf.insert_knot(.5,      0)
        surf.insert_knot(.7,      0)
        surf.insert_knot(.1,      1)
        surf.insert_knot(1.0 / 3, 1)
        knot1, knot2 = surf.knots(with_multiplicities=True)
        self.assertEqual(len(knot1), 10)  # 7 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones

        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after InsertKnot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

        # test errors and exceptions
        with self.assertRaises(TypeError):
            surf.insert_knot(1, 2, 3)  # too many arguments
        with self.assertRaises(ValueError):
            surf.insert_knot("tree-fiddy", .5)  # wrong argument type
        with self.assertRaises(ValueError):
            surf.insert_knot(0, -0.2)  # Outside-domain error
        with self.assertRaises(ValueError):
            surf.insert_knot(1, 1.4)  # Outside-domain error

    def test_force_rational(self):
        # more or less random 3D surface with p=[3,2] and n=[4,3]
        controlpoints = [[0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 1], [1, 0, 1], [1, 1, 1],
                         [2, 1, 1], [2, 2, 1], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(4, [0, 0, 0, 0, 2, 2, 2, 2])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)

        evaluation_point1 = surf(0.23, .66)
        control_point1 = surf[0]
        surf.force_rational()
        evaluation_point2 = surf(0.23, .66)
        control_point2 = surf[0]
        # ensure that surface has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        # ensure that we include rational weights of 1
        self.assertEqual(len(control_point1), 3)
        self.assertEqual(len(control_point2), 4)
        self.assertEqual(control_point2[3], 1)
        self.assertEqual(surf.rational, True)

    def test_swap(self):
        # more or less random 3D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 1], [1, 0, .5], [1, 1, 1],
                         [2, 1, 1], [2, 2, .5], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, .64, 2, 2, 2])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)

        evaluation_point1 = surf(0.23, .56)
        control_point1 = surf[1]  # this is control point i=(1,0), when n=(4,3)
        surf.swap()
        evaluation_point2 = surf(0.56, .23)
        control_point2 = surf[3]  # this is control point i=(0,1), when n=(3,4)

        # ensure that surface has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # check that the control points have re-ordered themselves
        self.assertEqual(control_point1[0], control_point2[0])
        self.assertEqual(control_point1[1], control_point2[1])
        self.assertEqual(control_point1[2], control_point2[2])

    def test_split(self):
        # test a rational 2D surface
        controlpoints = [[0, 0, 1], [-1, 1, .96], [0, 2, 1], [1, -1, 1], [1, 0, .8], [1, 1, 1],
                         [2, 1, .89], [2, 2, .9], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [1, 1, 1, 1.4, 5, 5, 5])
        basis2 = BSplineBasis(3, [2, 2, 2, 7, 7, 7])
        surf = Surface(basis1, basis2, controlpoints, True)

        split_u_surf = surf.split([1.1, 1.6, 4], 0)
        split_v_surf = surf.split(3.1, 1)

        self.assertEqual(len(split_u_surf), 4)
        self.assertEqual(len(split_v_surf), 2)

        # check that the u-vector is properly split
        self.assertAlmostEqual(split_u_surf[0].start()[0], 1.0)
        self.assertAlmostEqual(split_u_surf[0].end()[0], 1.1)
        self.assertAlmostEqual(split_u_surf[1].start()[0], 1.1)
        self.assertAlmostEqual(split_u_surf[1].end()[0], 1.6)
        self.assertAlmostEqual(split_u_surf[2].start()[0], 1.6)
        self.assertAlmostEqual(split_u_surf[2].end()[0], 4.0)
        self.assertAlmostEqual(split_u_surf[3].start()[0], 4.0)
        self.assertAlmostEqual(split_u_surf[3].end()[0], 5.0)
        # check that the v-vectors remain unchanged
        self.assertAlmostEqual(split_u_surf[2].start()[1], 2.0)
        self.assertAlmostEqual(split_u_surf[2].end()[1], 7.0)
        # check that the v-vector is properly split
        self.assertAlmostEqual(split_v_surf[0].start()[1], 2.0)
        self.assertAlmostEqual(split_v_surf[0].end()[1], 3.1)
        self.assertAlmostEqual(split_v_surf[1].start()[1], 3.1)
        self.assertAlmostEqual(split_v_surf[1].end()[1], 7.0)
        # check that the u-vector remain unchanged
        self.assertAlmostEqual(split_v_surf[1].start()[0], 1.0)
        self.assertAlmostEqual(split_v_surf[1].end()[0], 5.0)

        # check that evaluations remain unchanged
        pt1 = surf(3.23, 2.12)

        self.assertAlmostEqual(split_u_surf[2].evaluate(3.23, 2.12)[0], pt1[0])
        self.assertAlmostEqual(split_u_surf[2].evaluate(3.23, 2.12)[1], pt1[1])

        self.assertAlmostEqual(split_v_surf[0].evaluate(3.23, 2.12)[0], pt1[0])
        self.assertAlmostEqual(split_v_surf[0].evaluate(3.23, 2.12)[1], pt1[1])

    def test_reparam(self):
        # identity mapping, control points generated from knot vector
        basis1 = BSplineBasis(4, [2,2,2,2,3,6,7,7,7,7])
        basis2 = BSplineBasis(3, [-3,-3,-3,20,30,31,31,31])
        surf = Surface(basis1, basis2)

        self.assertAlmostEqual(surf.start(0),  2)
        self.assertAlmostEqual(surf.end(0),    7)
        self.assertAlmostEqual(surf.start(1), -3)
        self.assertAlmostEqual(surf.end(1),   31)

        surf.reparam((4,10), (0,9))
        self.assertAlmostEqual(surf.start(0),  4)
        self.assertAlmostEqual(surf.end(0),   10)
        self.assertAlmostEqual(surf.start(1),  0)
        self.assertAlmostEqual(surf.end(1),    9)

        surf.reparam((5,11), direction=0)
        self.assertAlmostEqual(surf.start(0),  5)
        self.assertAlmostEqual(surf.end(0),   11)
        self.assertAlmostEqual(surf.start(1),  0)
        self.assertAlmostEqual(surf.end(1),    9)

        surf.reparam((5,11), direction='v')
        self.assertAlmostEqual(surf.start(0),  5)
        self.assertAlmostEqual(surf.end(0),   11)
        self.assertAlmostEqual(surf.start(1),  5)
        self.assertAlmostEqual(surf.end(1),   11)

        surf.reparam((-9,9))
        self.assertAlmostEqual(surf.start(0), -9)
        self.assertAlmostEqual(surf.end(0),    9)
        self.assertAlmostEqual(surf.start(1),  0)
        self.assertAlmostEqual(surf.end(1),    1)

        surf.reparam()
        self.assertAlmostEqual(surf.start(0),  0)
        self.assertAlmostEqual(surf.end(0),    1)
        self.assertAlmostEqual(surf.start(1),  0)
        self.assertAlmostEqual(surf.end(1),    1)

        surf.reparam((4,10), (0,9))
        surf.reparam(direction=1)
        self.assertAlmostEqual(surf.start(0),  4)
        self.assertAlmostEqual(surf.end(0),   10)
        self.assertAlmostEqual(surf.start(1),  0)
        self.assertAlmostEqual(surf.end(1),    1)

    def test_periodic_split(self):
        # create a double-periodic spline on the knot vector [-1,0,0,1,1,2,2,3,3,4,4,5]*pi/2
        surf = SurfaceFactory.torus()

        surf2 = surf.split( pi/2, 0) # split on existing knot
        surf3 = surf.split( 1.23, 1) # split between knots
        surf4 = surf2.split(pi,   1) # split both periodicities

        # check periodicity tags
        self.assertEqual(surf.periodic(0), True)
        self.assertEqual(surf.periodic(1), True)
        self.assertEqual(surf2.periodic(0), False)
        self.assertEqual(surf2.periodic(1), True)
        self.assertEqual(surf3.periodic(0), True)
        self.assertEqual(surf3.periodic(1), False)
        self.assertEqual(surf4.periodic(0), False)
        self.assertEqual(surf4.periodic(1), False)

        # check parametric domain boundaries
        self.assertAlmostEqual(surf2.start(0),   pi/2)
        self.assertAlmostEqual(surf2.end(0),   5*pi/2)
        self.assertAlmostEqual(surf3.start(0),      0)
        self.assertAlmostEqual(surf3.end(0),     2*pi)
        self.assertAlmostEqual(surf3.start(1),   1.23)
        self.assertAlmostEqual(surf3.end(1),   1.23+2*pi)

        # check knot vector lengths
        self.assertEqual(len(surf2.knots(0, True)), 12)
        self.assertEqual(len(surf2.knots(1, True)), 12)
        self.assertEqual(len(surf3.knots(0, True)), 12)
        self.assertEqual(len(surf3.knots(1, True)), 14)
        self.assertEqual(len(surf4.knots(0, True)), 12)
        self.assertEqual(len(surf4.knots(1, True)), 12)

        # check that evaluation is unchanged over a 9x9 grid of shared parametric coordinates
        u   = np.linspace(pi/2, 2*pi, 9)
        v   = np.linspace(pi,   2*pi, 9)
        pt  = surf( u,v)
        pt2 = surf2(u,v)
        pt3 = surf3(u,v)
        pt4 = surf4(u,v)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)
        self.assertAlmostEqual(np.linalg.norm(pt-pt3), 0.0)
        self.assertAlmostEqual(np.linalg.norm(pt-pt4), 0.0)

    def test_make_identical(self):
        basis1 = BSplineBasis(4, [-1,-1,0,0,1,2,9,9,10,10,11,12], periodic=1)
        basis2 = BSplineBasis(2, [0,0,.25,1,1])
        cp    = [[0,0,0,1], [1,0,0,1.1], [2,0,0,1], [0,1,0,.7], [1,1,0,.8], [2,1,0,1], [0,2,0,1], [1,2,0,1.2], [2,2,0,1]]
        surf1 = Surface(BSplineBasis(3), BSplineBasis(3), cp, True) # rational 3D
        surf2 = Surface(basis1, basis2)                             # periodic 2D
        surf1.insert_knot([0.25, .5, .75], direction='v')

        Surface.make_splines_identical(surf1,surf2)

        for s in (surf1, surf2):
            self.assertEqual(s.periodic(0), False)
            self.assertEqual(s.periodic(1), False)
            self.assertEqual(s.dimension,   3)
            self.assertEqual(s.rational, True)
            self.assertEqual(s.order(), (4,3))

            self.assertAlmostEqual(len(s.knots(0, True)), 12)
            self.assertAlmostEqual(len(s.knots(1, True)), 10)

            self.assertEqual(s.bases[1].continuity(.25), 0)
            self.assertEqual(s.bases[1].continuity(.75), 1)
            self.assertEqual(s.bases[0].continuity(.2), 2)
            self.assertEqual(s.bases[0].continuity(.9), 1)

    def test_center(self):
        # make an ellipse at (2,1)
        surf = SurfaceFactory.disc(3)
        surf.scale((3,1))
        surf += (2,1)
        center = surf.center()
        self.assertEqual(len(center), 2)
        self.assertAlmostEqual(center[0], 2.0)
        self.assertAlmostEqual(center[1], 1.0)

    def test_reverse(self):
        basis1 = BSplineBasis(4, [2,2,2,2,3,6,7,7,7,7])
        basis2 = BSplineBasis(3, [-3,-3,-3,20,30,31,31,31])
        surf  = Surface(basis1, basis2)
        surf2 = Surface(basis1, basis2)
        surf3 = Surface(basis1, basis2)

        surf2.reverse('v')
        surf3.reverse('u')

        for i in range(6):
            # loop over surf forward, and surf2 backward (in 'v'-direction)
            for (cp1, cp2) in zip(surf[i,:,:], surf2[i,::-1,:]):
                self.assertAlmostEqual(cp1[0], cp2[0])
                self.assertAlmostEqual(cp1[1], cp2[1])

        for j in range(5):
            # loop over surf forward, and surf3 backward (in 'u'-direction)
            for (cp1, cp2) in zip(surf[:,j,:], surf3[::-1,j,:]):
                self.assertAlmostEqual(cp1[0], cp2[0])
                self.assertAlmostEqual(cp1[1], cp2[1])


    def test_repr(self):
        self.assertEqual(repr(Surface()), 'p=2, [ 0.  0.  1.  1.]\n'
                                          'p=2, [ 0.  0.  1.  1.]\n'
                                          '[ 0.  0.]\n'
                                          '[ 1.  0.]\n'
                                          '[ 0.  1.]\n'
                                          '[ 1.  1.]\n')

    def test_edges(self):
        (umin, umax, vmin, vmax) = Surface().edges()
        # check controlpoints
        self.assertAlmostEqual(umin[0,0], 0)
        self.assertAlmostEqual(umin[0,1], 0)
        self.assertAlmostEqual(umin[1,0], 0)
        self.assertAlmostEqual(umin[1,1], 1)

        self.assertAlmostEqual(umax[0,0], 1)
        self.assertAlmostEqual(umax[0,1], 0)
        self.assertAlmostEqual(umax[1,0], 1)
        self.assertAlmostEqual(umax[1,1], 1)

        self.assertAlmostEqual(vmin[0,0], 0)
        self.assertAlmostEqual(vmin[0,1], 0)
        self.assertAlmostEqual(vmin[1,0], 1)
        self.assertAlmostEqual(vmin[1,1], 0)

        # check a slightly more general surface
        cp = [[0,0], [.5, -.5], [1, 0],
              [-.6,1], [1, 1], [2, 1.4],
              [0,2], [.8, 3], [2, 2.4]]
        surf = Surface(BSplineBasis(3), BSplineBasis(3), cp)
        edg  = surf.edges()
        u   = np.linspace(0,1, 9)
        v   = np.linspace(0,1, 9)

        pt  = surf(0,v)
        pt2 = edg[0](v)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        pt  = surf(1,v)
        pt2 = edg[1](v)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        pt  = surf(u,0).reshape(9,2)
        pt2 = edg[2](u)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        pt  = surf(u,1).reshape(9,2)
        pt2 = edg[3](u)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

    def test_normal(self):
        surf = SurfaceFactory.sphere(1)
        surf.swap()
        u    = np.linspace(surf.start(0),surf.end(0), 9)
        v    = np.linspace(surf.start(1),surf.end(1), 9)
        xpts = surf(u,v)
        npts = surf.normal(u,v)

        # check that the normal is pointing out of the unit ball on a 9x9 evaluation grid
        for (i,j) in zip(xpts, npts):
            for (x,n) in zip(i,j):
                self.assertAlmostEqual(n[0], x[0])
                self.assertAlmostEqual(n[1], x[1])
                self.assertAlmostEqual(n[2], x[2])

        # check single value input
        n = surf.normal(0,0)
        self.assertEqual(len(n), 3)

        # test 2D surface
        s = Surface()
        n = s.normal(.5,.5)
        self.assertEqual(len(n), 3)
        self.assertAlmostEqual(n[0], 0.0)
        self.assertAlmostEqual(n[1], 0.0)
        self.assertAlmostEqual(n[2], 1.0)
        n = s.normal([.25, .5], [.1,.2,.3,.4,.5,.6,.7,.8,.9])
        self.assertEqual(n.shape[0], 2)
        self.assertEqual(n.shape[1], 9)
        self.assertEqual(n.shape[2], 3)
        self.assertAlmostEqual(n[1,4,0], 0.0)
        self.assertAlmostEqual(n[1,4,1], 0.0)
        self.assertAlmostEqual(n[1,4,2], 1.0)

        # test errors
        s = Surface(BSplineBasis(3), BSplineBasis(3), [[0]]*9) # 1D-surface
        with self.assertRaises(RuntimeError):
            s.normal(.5, .5)

    def test_area(self):
        s = Surface()
        self.assertAlmostEqual(s.area(), 1.0)
        s = Surface(controlpoints=[[-2,0,1], [2,0,1], [-1,1,1], [1,1,1]])
        self.assertAlmostEqual(s.area(), 3.0)

    def test_const_par_crv(self):
        # more or less random 2D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0, 0], [-1, 1], [0, 2], [1, -1], [1, 0], [1, 1], [2, 1], [2, 2], [2, 3],
                         [3, 0], [4, 1], [3, 2]]
        basis1 = BSplineBasis(3, [0, 0, 0, .4, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        surf = Surface(basis1, basis2, controlpoints)
        surf.refine(1)

        # try general knot in u-direction
        crv = surf.const_par_curve(0.3, 'u')
        v = np.linspace(0,1,13)
        self.assertTrue(np.allclose(surf(0.3, v), crv(v)))

        # try existing knot in u-direction
        crv = surf.const_par_curve(0.4, 'u')
        v = np.linspace(0,1,13)
        self.assertTrue(np.allclose(surf(0.4, v), crv(v)))

        # try general knot in v-direction
        crv = surf.const_par_curve(0.3, 'v')
        u = np.linspace(0,1,13)
        self.assertTrue(np.allclose(surf(u, 0.3).reshape(13,2), crv(u)))

        # try start-point
        crv = surf.const_par_curve(0.0, 'v')
        u = np.linspace(0,1,13)
        self.assertTrue(np.allclose(surf(u, 0.0).reshape(13,2), crv(u)))

        # try end-point
        crv = surf.const_par_curve(1.0, 'v')
        u = np.linspace(0,1,13)
        self.assertTrue(np.allclose(surf(u, 1.0).reshape(13,2), crv(u)))




if __name__ == '__main__':
    unittest.main()
