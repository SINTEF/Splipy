from Surface import *
import unittest

class TestSurface(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        controlpoints = [[0,0,0], [1,0,0], [0,1,0], [1,1,0]]
        surf = Surface(2,2, [0,0,1,1], [0,0,1,1], controlpoints)
        val = surf.Evaluate(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf[0]), 3)

        # test 2D constructor
        controlpoints = [[0,0], [1,0], [0,1], [1,1]]
        surf2 = Surface(2,2, [0,0,1,1], [0,0,1,1], controlpoints)
        val  = surf2.Evaluate(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf2[0]), 2)

        # test rational 2D constructor
        controlpoints = [[0,0,1], [1,0,1], [0,1,1], [1,1,1]]
        surf3 = Surface(2,2, [0,0,1,1], [0,0,1,1], controlpoints, True)
        val = surf3.Evaluate(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf3[0]), 2)

        # test rational 3D constructor
        controlpoints = [[0,0,0,1], [1,0,0,1], [0,1,0,1], [1,1,0,1]]
        surf4 = Surface(2,2, [0,0,1,1], [0,0,1,1], controlpoints, True)
        val = surf4.Evaluate(0.5, 0.5)
        self.assertEqual(val[0], 0.5)
        self.assertEqual(len(surf4[0]), 3)

        # TODO: Include a default constructor specifying nothing, or just polynomial degrees, or just knot vectors.
        #       This should create identity mappings

        # test errors and exceptions
        controlpoints = [[0,0,1], [1,0,1], [0,1,1], [1,1,1]]
        with self.assertRaises(ValueError):
            surf = Surface(2,2,[1,1,0,0], [0,0,1,1], controlpoints)            # illegal knot vector
        with self.assertRaises(ValueError):
            surf = Surface(2,2,[0,0,.5,1,1], [0,0,1,1], controlpoints)         # too few controlpoints
        # TODO: Create fail tests for rational surfaces with weights equal to zero
        #       Create fail tests for providing too few control points
        #       Create fail tests for providing too many control points

    def test_evaluate(self):
        # knot vector [t_1, t_2, ... t_{n+p+1}]
        # polynomial degree p
        # n basis functions N_i(t), for i=1...n
        # the power basis {1,t,t^2,t^3,...} can be expressed as:
        # 1     = sum         N_i(t)
        # t     = sum ts_i  * N_i(t)
        # t^2   = sum t2s_i * N_i(t)
        # ts_i  = sum_{j=i+1}^{i+p}   t_j / p
        # t2s_i = sum_{j=i+1}^{i+p-1} sum_{k=j+1}^{i+p} t_j*t_k / (p 2)
        # (p 2) = binomial coefficent


        # creating the mapping:
        #   x(u,v) = u^2*v^2
        #   y(u,v) = u*v
        controlpoints = [[0,0],[0,0],[0,0],  [0,0],[0,.25],[0,.5],  [0,0], [0,.5], [1,1]]
        surf = Surface(3,3, [0,0,0,1,1,1], [0,0,0,1,1,1], controlpoints)

        # startpoint evaluation
        val = surf.Evaluate(0,0)     
        self.assertAlmostEqual(val[0], 0.0)
        self.assertAlmostEqual(val[1], 0.0)

        # startpoint with derivatives
        val = surf.Evaluate(0,0,  1)  
        self.assertAlmostEqual(val[0][0], 0.0)
        self.assertAlmostEqual(val[0][1], 0.0)
        self.assertAlmostEqual(val[1][0], 0.0) # dx/du = 2*u*v^2
        self.assertAlmostEqual(val[1][1], 0.0) # dy/du = v
        self.assertAlmostEqual(val[2][0], 0.0) # dx/dv = 2*u^2*v
        self.assertAlmostEqual(val[2][1], 0.0) # dy/dv = u      

        # midpoint with derivatives
        val = surf.Evaluate(0.5, 0.5,   1)  
        self.assertAlmostEqual(val[0][0], 0.0625)
        self.assertAlmostEqual(val[0][1], 0.25)
        self.assertAlmostEqual(val[1][0], 0.25)  # dx/du = 2*u*v^2
        self.assertAlmostEqual(val[1][1], 0.5)   # dy/du = v
        self.assertAlmostEqual(val[2][0], 0.25)  # dx/dv = 2*u^2*v
        self.assertAlmostEqual(val[2][1], 0.5)   # dy/dv = u      

        # mid top edge with derivatives (sensitive to left-evaluation)
        val = surf.Evaluate(0.5, 1.0,   1)  
        self.assertAlmostEqual(val[0][0], 0.25)
        self.assertAlmostEqual(val[0][1], 0.5)
        self.assertAlmostEqual(val[1][0], 1.0) # dx/du = 2*u*v^2
        self.assertAlmostEqual(val[1][1], 1.0) # dy/du = v
        self.assertAlmostEqual(val[2][0], 0.5) # dx/dv = 2*u^2*v
        self.assertAlmostEqual(val[2][1], 0.5) # dy/dv = u      

        # mid right edge with derivatives (sensitive to left-evaluation)
        val = surf.Evaluate(1.0, 0.5,   1)  
        self.assertAlmostEqual(val[0][0], 0.25)
        self.assertAlmostEqual(val[0][1], 0.5)
        self.assertAlmostEqual(val[1][0], 0.5) # dx/du = 2*u*v^2
        self.assertAlmostEqual(val[1][1], 0.5) # dy/du = v
        self.assertAlmostEqual(val[2][0], 1.0) # dx/dv = 2*u^2*v
        self.assertAlmostEqual(val[2][1], 1.0) # dy/dv = u      

        # top right corner with derivatives (sensitive to left-evaluation)
        val = surf.Evaluate(1.0, 1.0,   1)  
        self.assertAlmostEqual(val[0][0], 1.0)
        self.assertAlmostEqual(val[0][1], 1.0)
        self.assertAlmostEqual(val[1][0], 2.0) # dx/du = 2*u*v^2
        self.assertAlmostEqual(val[1][1], 1.0) # dy/du = v
        self.assertAlmostEqual(val[2][0], 2.0) # dx/dv = 2*u^2*v
        self.assertAlmostEqual(val[2][1], 1.0) # dy/dv = u      

        # second derivatives at midpoint 
        val = surf.Evaluate(0.5, 0.5,   2)  
        self.assertAlmostEqual(len(val),  6)
        self.assertAlmostEqual(val[3][0],  0.5) # d2x/du2  = 2*v^2
        self.assertAlmostEqual(val[3][1],  0.0) # d2y/du2  = 0
        self.assertAlmostEqual(val[4][0],  1.0) # d2x/dudv = 4*u*v
        self.assertAlmostEqual(val[4][1],  1.0) # d2y/dudv = 1
        self.assertAlmostEqual(val[5][0],  0.5) # d2x/dv2  = 2*u^2
        self.assertAlmostEqual(val[5][1],  0.0) # d2y/dv2  = 0

        # third derivatives
        val = surf.Evaluate(0.5, 0.5,   3)  
        self.assertAlmostEqual(len(val),  10)
        self.assertAlmostEqual(val[6][0],  0.0) # d3x/du3   = 0
        self.assertAlmostEqual(val[6][1],  0.0) # d3y/du3   = 0
        self.assertAlmostEqual(val[7][0],  2.0) # d3x/du2dv = 4*v
        self.assertAlmostEqual(val[7][1],  0.0) # d3y/du2dv = 0
        self.assertAlmostEqual(val[8][0],  2.0) # d3x/dudv2 = 4*u
        self.assertAlmostEqual(val[9][0],  0.0) # d3x/dv3   = 0

        # fourth derivatives
        val = surf.Evaluate(0.5, 0.5,   4)  
        self.assertAlmostEqual(len(val),  15)
        self.assertAlmostEqual(val[11][0],  0.0) # d4x/du3dv  = 0
        self.assertAlmostEqual(val[12][0],  4.0) # d4x/du2dv2 = 4
        self.assertAlmostEqual(val[13][0],  0.0) # d4x/dudv3  = 0

        # check integer type for derivative
        with self.assertRaises(TypeError):
            val = surf.Evaluate(0.5,0.5,  1.5) 

        # GoTools throws exception for negative derivatives
        with self.assertRaises(Exception):
            val = surf.Evaluate(0.5,0.5,  -1) 

    def test_raise_order(self):
        # more or less random 2D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0,0],[-1,1],[0,2],  [1,-1],[1,0],[1,1],  [2,1],[2,2],[2,3],  [3,0],[4,1],[3,2]]
        surf = Surface(3,3, [0,0,0,.4,1,1,1], [0,0,0,1,1,1], controlpoints)

        self.assertEqual(surf.GetOrder()[0], 3)
        self.assertEqual(surf.GetOrder()[1], 3)
        evaluation_point1 = surf.Evaluate(0.23, 0.37) # pick some evaluation point (could be anything)

        surf.RaiseOrder(1,2)

        self.assertEqual(surf.GetOrder()[0], 4)
        self.assertEqual(surf.GetOrder()[1], 5)
        evaluation_point2 = surf.Evaluate(0.23, 0.37)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])


        # test a rational 2D surface 
        controlpoints = [[0,0,1],[-1,1,.96],[0,2,1],  [1,-1,1],[1,0,.8],[1,1,1],  [2,1,.89],[2,2,.9],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        surf = Surface(3,3, [0,0,0,.4,1,1,1], [0,0,0,1,1,1], controlpoints, True)

        self.assertEqual(surf.GetOrder()[0], 3)
        self.assertEqual(surf.GetOrder()[1], 3)
        evaluation_point1 = surf.Evaluate(0.23, 0.37)

        surf.RaiseOrder(1,2)

        self.assertEqual(surf.GetOrder()[0], 4)
        self.assertEqual(surf.GetOrder()[1], 5)
        evaluation_point2 = surf.Evaluate(0.23, 0.37)

        # evaluation before and after RaiseOrder should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

    def test_insert_knot(self):
        # more or less random 2D surface with p=[3,2] and n=[4,3]
        controlpoints = [[0,0],[-1,1],[0,2],  [1,-1],[1,0],[1,1],  [2,1],[2,2],[2,3],  [3,0],[4,1],[3,2]]
        surf = Surface(4,3, [0,0,0,0,2,2,2,2], [0,0,0,1,1,1], controlpoints)

        # pick some evaluation point (could be anything)
        evaluation_point1 = surf.Evaluate(0.23, 0.37)

        surf.InsertKnot(0, .22)
        surf.InsertKnot(0, .5)
        surf.InsertKnot(0, .7)
        surf.InsertKnot(1, .1)
        surf.InsertKnot(1, 1.0/3)
        knot1, knot2 = surf.GetKnots(True)
        self.assertEqual(len(knot1), 11) # 8 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones

        evaluation_point2 = surf.Evaluate(0.23, 0.37)

        # evaluation before and after InsertKnot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])


        # test a rational 2D surface 
        controlpoints = [[0,0,1],[-1,1,.96],[0,2,1],  [1,-1,1],[1,0,.8],[1,1,1],  [2,1,.89],[2,2,.9],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        surf = Surface(3,3, [0,0,0,.4,1,1,1], [0,0,0,1,1,1], controlpoints, True)

        evaluation_point1 = surf.Evaluate(0.23, 0.37)

        surf.InsertKnot(0, .22)
        surf.InsertKnot(0, .5)
        surf.InsertKnot(0, .7)
        surf.InsertKnot(1, .1)
        surf.InsertKnot(1, 1.0/3)
        knot1, knot2 = surf.GetKnots(True)
        self.assertEqual(len(knot1), 10) # 7 to start with, 3 new ones
        self.assertEqual(len(knot2), 8)  # 6 to start with, 2 new ones

        evaluation_point2 = surf.Evaluate(0.23, 0.37)

        # evaluation before and after InsertKnot should remain unchanged
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

        # test errors and exceptions
        with self.assertRaises(TypeError):
            surf.InsertKnot(1, 2, 3)          # too many arguments
        with self.assertRaises(TypeError):
            surf.InsertKnot(1)                # too few arguments
        with self.assertRaises(TypeError):
            surf.InsertKnot("tree-fiddy", .5) # wrong argument type
        with self.assertRaises(ValueError):
            surf.InsertKnot(0, -0.2)          # Outside-domain error
        with self.assertRaises(ValueError):
            surf.InsertKnot(1, 1.4)           # Outside-domain error

    def test_force_rational(self):
        # more or less random 3D surface with p=[3,2] and n=[4,3]
        controlpoints = [[0,0,1],[-1,1,1],[0,2,1],  [1,-1,1],[1,0,1],[1,1,1],  [2,1,1],[2,2,1],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        surf = Surface(4,3, [0,0,0,0,2,2,2,2], [0,0,0,1,1,1], controlpoints)

        evaluation_point1 = surf.Evaluate(0.23, .66)
        control_point1    = surf[0]
        surf.ForceRational()
        evaluation_point2 = surf.Evaluate(0.23, .66)
        control_point2    = surf[0]
        # ensure that surface has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        # ensure that we include rational weights of 1
        self.assertEqual(len(control_point1), 3)
        self.assertEqual(len(control_point2), 3)
        self.assertEqual(surf.rational, True)

    def test_swap_parametrization(self):
        # more or less random 3D surface with p=[2,2] and n=[4,3]
        controlpoints = [[0,0,1],[-1,1,1],[0,2,1],  [1,-1,1],[1,0,.5],[1,1,1],  [2,1,1],[2,2,.5],[2,3,1],  [3,0,1],[4,1,1],[3,2,1]]
        surf = Surface(3,3, [0,0,0,.64,2,2,2], [0,0,0,1,1,1], controlpoints)

        evaluation_point1 = surf.Evaluate(0.23, .56)
        control_point1    = surf[1]  # this is control point i=(1,0), when n=(4,3)
        surf.SwapParametrization()
        evaluation_point2 = surf.Evaluate(0.56, .23)
        control_point2    = surf[3]  # this is control point i=(0,1), when n=(3,4)

        # ensure that surface has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # check that the control points have re-ordered themselves
        self.assertEqual(control_point1[0], control_point2[0])
        self.assertEqual(control_point1[1], control_point2[1])
        self.assertEqual(control_point1[2], control_point2[2])

if __name__ == '__main__':
    unittest.main()
