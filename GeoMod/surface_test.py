from Surface import *
import unittest

class TestSurface(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        cp = [[0,0,0], [1,0,0], [0,1,0], [1,1,0]]
        surf = Surface(controlpoints=cp)
        val = surf(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf[0]), 3)

        # test 2D constructor
        cp = [[0,0], [1,0], [0,1], [1,1]]
        surf2 = Surface(controlpoints=cp)
        val  = surf2(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf2[0]), 2)

        # test rational 2D constructor
        cp = [[0,0,1], [1,0,1], [0,1,1], [1,1,1]]
        surf3 = Surface(controlpoints=cp, rational=True)
        val = surf3(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf3[0]), 3)

        # test rational 3D constructor
        cp = [[0,0,0,1], [1,0,0,1], [0,1,0,1], [1,1,0,1]]
        surf4 = Surface(controlpoints=cp, rational=True)
        val = surf4(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf4[0]), 4)

        # TODO: Include a default constructor specifying nothing, or just polynomial degrees, or just knot vectors.
        #       This should create identity mappings

        # test errors and exceptions
        controlpoints = [[0,0,1], [1,0,1], [0,1,1], [1,1,1]]
        with self.assertRaises(ValueError):
            basis1 = BSplineBasis(2, [1,1,0,0])
            basis2 = BSplineBasis(2, [0,0,1,1])
            surf = Surface(basis1, basis2, controlpoints)        # illegal knot vector
        with self.assertRaises(ValueError):
            basis1 = BSplineBasis(2, [0,0,.5,1,1])
            basis2 = BSplineBasis(2, [0,0,1,1])
            surf = Surface(basis1,basis2, controlpoints)         # too few controlpoints
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
        controlpoints = [[0,0],[1.0/4,0],[3.0/4,0],[.75,0],  [0,1],[0,1],[.5,1],[1,1]]
        basis1 = BSplineBasis(3, [0,0,0,.5,1,1,1])
        basis2 = BSplineBasis(2, [0,0,1,1])
        surf = Surface(basis1, basis2, controlpoints)

        # call evaluation at a 5x4 grid of points
        val = surf([0,.2,.5,.6,1], [0,.2,.4,1])
        self.assertEqual(len(val.shape), 3) # result should be wrapped in 3-index tensor
        self.assertEqual(val.shape[0],   5) # 5 evaluation points in u-direction
        self.assertEqual(val.shape[1],   4) # 4 evaluation points in v-direction
        self.assertEqual(val.shape[2],   2) # 2 coordinates (x,y)

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
            val = surf(-10, .5) # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf(+10, .3) # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf(.5, -10) # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf(.5, +10) # evalaute outside parametric domain

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
        controlpoints = [[0,0],[1.0/4,0],[3.0/4,0],[.75,0],  [0,1],[0,1],[.5,1],[1,1]]
        basis1 = BSplineBasis(3, [0,0,0,.5,1,1,1])
        basis2 = BSplineBasis(2, [0,0,1,1])
        surf = Surface(basis1, basis2, controlpoints)

        # call evaluation at a 5x4 grid of points
        val = surf.evaluate_derivative([0,.2,.5,.6,1], [0,.2,.4,1], 1,0)
        self.assertEqual(len(val.shape), 3) # result should be wrapped in 3-index tensor
        self.assertEqual(val.shape[0],   5) # 5 evaluation points in u-direction
        self.assertEqual(val.shape[1],   4) # 4 evaluation points in v-direction
        self.assertEqual(val.shape[2],   2) # 2 coordinates (x,y)

        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,1,0)[0],  .88) # dx/du=2uv+(1-v)
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,1,0)[1],    0) # dy/du=0
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,0,1)[0], -.16) # dx/dv=u^2-u
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,0,1)[1],    1) # dy/dv=1
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,1,1)[0], -.60) # d2x/dudv=2u-1
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,2,0)[0], 0.40) # d2x/dudu=2v
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,3,0)[0], 0.00) # d3x/du3=0
        self.assertAlmostEqual(surf.evaluate_derivative(.2,.2,0,2)[0], 0.00) # d2y/dv2=0

        # test errors and exceptions
        with self.assertRaises(ValueError):
            val = surf.evaluate_derivative(-10, .5) # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf.evaluate_derivative(+10, .3) # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf.evaluate_derivative(.5, -10) # evalaute outside parametric domain
        with self.assertRaises(ValueError):
            val = surf.evaluate_derivative(.5, +10) # evalaute outside parametric domain


    def test_raise_order(self):
        # more or less random 2D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0,0],[-1,1],[0,2],  [1,-1],[1,0],[1,1],  [2,1],[2,2],[2,3],  [3,0],[4,1],[3,2]]
        basis1 = BSplineBasis(3, [0,0,0,.4,1,1,1])
        basis2 = BSplineBasis(3, [0,0,0,1,1,1])
        surf = Surface(basis1, basis2, controlpoints)

        self.assertEqual(surf.get_order()[0], 3)
        self.assertEqual(surf.get_order()[1], 3)
        evaluation_point1 = surf(0.23, 0.37) # pick some evaluation point (could be anything)

        surf.raise_order(1,2)

        self.assertEqual(surf.get_order()[0], 4)
        self.assertEqual(surf.get_order()[1], 5)
        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])


        # test a rational 2D surface 
        controlpoints = [[0,0,1],[-1,1,.96],[0,2,1],  [1,-1,1],[1,0,.8],[1,1,1],  [2,1,.89],[2,2,.9],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        basis1 = BSplineBasis(3, [0,0,0,.4,1,1,1])
        basis2 = BSplineBasis(3, [0,0,0,1,1,1])
        surf = Surface(basis1, basis2, controlpoints, True)

        self.assertEqual(surf.get_order()[0], 3)
        self.assertEqual(surf.get_order()[1], 3)
        evaluation_point1 = surf(0.23, 0.37)

        surf.raise_order(1,2)

        self.assertEqual(surf.get_order()[0], 4)
        self.assertEqual(surf.get_order()[1], 5)
        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

    def test_insert_knot(self):
        # more or less random 2D surface with p=[3,2] and n=[4,3]
        controlpoints = [[0,0],[-1,1],[0,2],  [1,-1],[1,0],[1,1],  [2,1],[2,2],[2,3],  [3,0],[4,1],[3,2]]
        basis1 = BSplineBasis(4, [0,0,0,0,2,2,2,2])
        basis2 = BSplineBasis(3, [0,0,0,1,1,1])
        surf = Surface(basis1, basis2, controlpoints)

        # pick some evaluation point (could be anything)
        evaluation_point1 = surf(0.23, 0.37)

        surf.insert_knot(0, .22)
        surf.insert_knot(0, .5)
        surf.insert_knot(0, .7)
        surf.insert_knot(1, .1)
        surf.insert_knot(1, 1.0/3)
        knot1, knot2 = surf.get_knots(True)
        self.assertEqual(len(knot1), 11) # 8 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones

        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after InsertKnot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])


        # test a rational 2D surface 
        controlpoints = [[0,0,1],[-1,1,.96],[0,2,1],  [1,-1,1],[1,0,.8],[1,1,1],  [2,1,.89],[2,2,.9],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        basis1 = BSplineBasis(3, [0,0,0,.4,1,1,1])
        basis2 = BSplineBasis(3, [0,0,0,1,1,1])
        surf = Surface(basis1, basis2, controlpoints, True)

        evaluation_point1 = surf(0.23, 0.37)

        surf.insert_knot(0, .22)
        surf.insert_knot(0, .5)
        surf.insert_knot(0, .7)
        surf.insert_knot(1, .1)
        surf.insert_knot(1, 1.0/3)
        knot1, knot2 = surf.get_knots(True)
        self.assertEqual(len(knot1), 10) # 7 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones

        evaluation_point2 = surf(0.23, 0.37)

        # evaluation before and after InsertKnot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

        # test errors and exceptions
        with self.assertRaises(TypeError):
            surf.insert_knot(1, 2, 3)          # too many arguments
        with self.assertRaises(TypeError):
            surf.insert_knot(1)                # too few arguments
        with self.assertRaises(ValueError):
            surf.insert_knot("tree-fiddy", .5) # wrong argument type
        with self.assertRaises(ValueError):
            surf.insert_knot(0, -0.2)          # Outside-domain error
        with self.assertRaises(ValueError):
            surf.insert_knot(1, 1.4)           # Outside-domain error

    def test_force_rational(self):
        # more or less random 3D surface with p=[3,2] and n=[4,3]
        controlpoints = [[0,0,1],[-1,1,1],[0,2,1],  [1,-1,1],[1,0,1],[1,1,1],  [2,1,1],[2,2,1],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        basis1 = BSplineBasis(4, [0,0,0,0,2,2,2,2])
        basis2 = BSplineBasis(3, [0,0,0,1,1,1])
        surf = Surface(basis1, basis2, controlpoints)

        evaluation_point1 = surf(0.23, .66)
        control_point1    = surf[0]
        surf.force_rational()
        evaluation_point2 = surf(0.23, .66)
        control_point2    = surf[0]
        # ensure that surface has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        # ensure that we include rational weights of 1
        self.assertEqual(len(control_point1), 3)
        self.assertEqual(len(control_point2), 4)
        self.assertEqual(control_point2[3], 1)
        self.assertEqual(surf.rational, True)

    def test_swap_parametrization(self):
        # more or less random 3D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0,0,1],[-1,1,1],[0,2,1],  [1,-1,1],[1,0,.5],[1,1,1],  [2,1,1],[2,2,.5],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        basis1 = BSplineBasis(3, [0,0,0,.64,2,2,2])
        basis2 = BSplineBasis(3, [0,0,0,1,1,1])
        surf = Surface(basis1, basis2, controlpoints)

        evaluation_point1 = surf(0.23, .56)
        control_point1    = surf[1]  # this is control point i=(1,0), when n=(4,3)
        surf.swap_parametrization()
        evaluation_point2 = surf(0.56, .23)
        control_point2    = surf[3]  # this is control point i=(0,1), when n=(3,4)

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
        controlpoints = [[0,0,1],[-1,1,.96],[0,2,1],  [1,-1,1],[1,0,.8],[1,1,1],  [2,1,.89],[2,2,.9],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        basis1 = BSplineBasis(3, [1,1,1, 1.4,5,5,5])
        basis2 = BSplineBasis(3, [2,2,2,7,7,7])
        surf = Surface(basis1, basis2, controlpoints, True)

        split_u_surf = surf.split(0, [1.1, 1.6, 4])
        split_v_surf = surf.split(1, 3.1)

        self.assertEqual(len(split_u_surf), 4)
        self.assertEqual(len(split_v_surf), 2)
        
        # check that the u-vector is properly split
        self.assertAlmostEqual(split_u_surf[0].start()[0], 1.0)
        self.assertAlmostEqual(split_u_surf[0].end()[0],   1.1)
        self.assertAlmostEqual(split_u_surf[1].start()[0], 1.1)
        self.assertAlmostEqual(split_u_surf[1].end()[0],   1.6)
        self.assertAlmostEqual(split_u_surf[2].start()[0], 1.6)
        self.assertAlmostEqual(split_u_surf[2].end()[0],   4.0)
        self.assertAlmostEqual(split_u_surf[3].start()[0], 4.0)
        self.assertAlmostEqual(split_u_surf[3].end()[0],   5.0)
        # check that the v-vectors remain unchanged
        self.assertAlmostEqual(split_u_surf[2].start()[1], 2.0)
        self.assertAlmostEqual(split_u_surf[2].end()[1],   7.0)
        # check that the v-vector is properly split
        self.assertAlmostEqual(split_v_surf[0].start()[1], 2.0)
        self.assertAlmostEqual(split_v_surf[0].end()[1],   3.1)
        self.assertAlmostEqual(split_v_surf[1].start()[1], 3.1)
        self.assertAlmostEqual(split_v_surf[1].end()[1],   7.0)
        # check that the u-vector remain unchanged
        self.assertAlmostEqual(split_v_surf[1].start()[0], 1.0)
        self.assertAlmostEqual(split_v_surf[1].end()[0],   5.0)

        # check that evaluations remain unchanged
        pt1 = surf(3.23, 2.12)

        self.assertAlmostEqual(split_u_surf[2].evaluate(3.23,2.12)[0], pt1[0])
        self.assertAlmostEqual(split_u_surf[2].evaluate(3.23,2.12)[1], pt1[1])

        self.assertAlmostEqual(split_v_surf[0].evaluate(3.23,2.12)[0], pt1[0])
        self.assertAlmostEqual(split_v_surf[0].evaluate(3.23,2.12)[1], pt1[1])





if __name__ == '__main__':
    unittest.main()
