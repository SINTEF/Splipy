from GoTools import Curve
from math import sqrt
import unittest

class TestCurve(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        crv = Curve(2, [0, 0, 1, 1], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        val = crv.Evaluate(0.5)
        self.assertEqual(val[0], 0.5)

        # test 2D constructor
        crv2 = Curve(2, [0, 0, 1, 1], [[0.0, 0.0], [1.0, 0.0]])
        val  = crv.Evaluate(0.5)
        self.assertEqual(val[0], 0.5)

    def test_evaluate(self):
        controlpoints = [[0,0,0],  [1,1,0],  [2,0,0]]
        crv = Curve(3, [0, 0, 0, 1, 1, 1], controlpoints)

        # startpoint evaluation
        val = crv.Evaluate(0.0)     
        self.assertEqual(val[0], 0.0)
        self.assertEqual(val[1], 0.0)

        # startpoint with derivatives
        val = crv.Evaluate(0.0, 1)  
        self.assertEqual(val[0][0], 0.0)
        self.assertEqual(val[0][1], 0.0)
        self.assertEqual(val[1][0], 2.0)
        self.assertEqual(val[1][1], 2.0)

        # midpoint with derivatives
        val = crv.Evaluate(0.5, 1)  
        self.assertEqual(val[0][0], 1.0)
        self.assertEqual(val[0][1], 0.5)
        self.assertEqual(val[1][0], 2.0)
        self.assertEqual(val[1][1], 0.0)

        # end with derivatives (sensitive to left-evaluation)
        val = crv.Evaluate(1.0, 1)  
        self.assertEqual(val[0][0],  2.0)
        self.assertEqual(val[0][1],  0.0)
        self.assertEqual(val[1][0],  2.0)
        self.assertEqual(val[1][1], -2.0)

        # second derivatives at midpoint
        val = crv.Evaluate(0.5, 2)  
        self.assertEqual(val[2][0],  0.0)
        self.assertEqual(val[2][1], -4.0)

        # third derivatives should all vanish everywhere
        val = crv.Evaluate(0.5, 3)  
        self.assertEqual(val[3][0],  0.0)
        self.assertEqual(val[3][1],  0.0)

        # check integer type for derivative
        with self.assertRaises(TypeError):
            val = crv.Evaluate(0.5, 1.5) 

        # GoTools throws exception for negative derivatives
        with self.assertRaises(Exception):
            val = crv.Evaluate(0.5, -1) 

    def test_flip_parametrization(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0,0,0],  [1,1,0],  [2,-1,0],  [3,0,0]]
        crv = Curve(3, [0, 0, 0, .3, 1, 1, 1], controlpoints)

        p1 = crv.Evaluate(0.23)
        crv.FlipParametrization()
        p2 = crv.Evaluate(0.77)
        self.assertAlmostEqual(p1[0], p2[0])
        self.assertAlmostEqual(p1[1], p2[1])
        self.assertAlmostEqual(p1[2], p2[2])

    def test_force_rational(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0,0,0],  [1,1,0],  [2,-1,0],  [3,0,0]]
        crv = Curve(3, [0, 0, 0, .3, 1, 1, 1], controlpoints)

        evaluation_point1 = crv.Evaluate(0.23)
        control_point1    = crv[0]
        crv2              = crv.ForceRational()
        evaluation_point2 = crv2.Evaluate(0.23)
        control_point2    = crv2[0]
        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])
        # ensure that we include rational weights of 1
        self.assertEqual(len(control_point1), 3)
        self.assertEqual(len(control_point2), 4)
        self.assertEqual(control_point2[3], 1)

    def test_raise_order(self):
        # non-uniform knot vector of a squiggly quadratic n=5 curve in 3D
        controlpoints = [[0,0,0],  [1,1,1],  [2,-1,0],  [3,0,-1], [0,0,-5]]
        crv = Curve(3, [0, 0, 0, .3, .4, 1, 1, 1], controlpoints)

        crv
        evaluation_point1 = crv.Evaluate(0.37)
        crv2              = crv.RaiseOrder(2)
        evaluation_point2 = crv2.Evaluate(0.37)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # ensure that curve has the right order
        self.assertEqual(crv.GetOrder(),  5)
        self.assertEqual(crv2.GetOrder(), 5)

        # check integer type for argument
        with self.assertRaises(TypeError):
            crv.RaiseOrder(0.5);

        # check logic error for negative argument (gotools cast this error)
        with self.assertRaises(Exception):
            crv.RaiseOrder(-1);

    def test_insert_knot(self):
        # non-uniform knot vector of a squiggly quadratic n=5 curve in 3D
        controlpoints = [[0,0,0],  [1,1,1],  [2,-1,0],  [3,0,-1], [0,0,-5]]
        crv = Curve(3, [0, 0, 0, .3, .4, 1, 1, 1], controlpoints)

        evaluation_point1 = crv.Evaluate(0.37)
        crv.InsertKnot(.2)
        crv.InsertKnot(.3)
        crv.InsertKnot(.9)
        evaluation_point2 = crv.Evaluate(0.37)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # ensure that curve has the knot length
        self.assertEqual(len(crv.GetKnots(True)),  11)

        # test rational curves, here a perfect circle represented as n=9, p=2-curve
        s = 1.0/sqrt(2)
        controlpoints = [[1,0,1], [s,s,s], [0,1,1], [-s,s,s], [-1,0,1], [-s,-s,s], [0,-1,1], [s,-s,s], [1,0,1]]
        crv = Curve(3, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4], controlpoints, True)

        evaluation_point1 = crv.Evaluate(0.37)
        crv.InsertKnot(.2)
        crv.InsertKnot(.3)
        crv.InsertKnot(.9)
        evaluation_point2 = crv.Evaluate(0.37)

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])

        # ensure that curve has the knot length
        self.assertEqual(len(crv.GetKnots(True)),  15)


        # test errors and exceptions
        with self.assertRaises(TypeError):
            crv.InsertKnot(1, 2)          # too many arguments
        with self.assertRaises(TypeError):
            crv.InsertKnot()              # too few arguments
        with self.assertRaises(TypeError):
            crv.InsertKnot("tree-fiddy")  # wrong argument type
        with self.assertRaises(ValueError):
            crv.InsertKnot(-0.2)          # Outside-domain error
        with self.assertRaises(ValueError):
            crv.InsertKnot( 4.4)          # Outside-domain error


    def test_reparametrize(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0,0,0],  [1,1,0],  [2,-1,0],  [3,0,0]]
        crv = Curve(3, [0, 0, 0, 1.32, 3, 3, 3], controlpoints)

        # get some info on the initial curve
        knots1            = crv.GetKnots()
        evaluation_point1 = crv.Evaluate(1.20)
        self.assertEqual(knots1[0],  0);
        self.assertEqual(knots1[-1], 3);

        # reparametrize
        crv.ReParametrize(6.0, 9.0);

        # get some info on the reparametrized curve
        knots2            = crv.GetKnots()
        evaluation_point2 = crv.Evaluate(7.20)
        self.assertEqual(knots2[0],  6);
        self.assertEqual(knots2[-1], 9);

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point2[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point2[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point2[2])

        # normalize, i.e. set domain to [0,1]
        crv.ReParametrize()

        # get some info on the normalized curve
        knots3            = crv.GetKnots()
        evaluation_point3 = crv.Evaluate(0.40)
        self.assertEqual(knots3[0],  0);
        self.assertEqual(knots3[-1], 1);

        # ensure that curve has not chcanged, by comparing evaluation of it
        self.assertAlmostEqual(evaluation_point1[0], evaluation_point3[0])
        self.assertAlmostEqual(evaluation_point1[1], evaluation_point3[1])
        self.assertAlmostEqual(evaluation_point1[2], evaluation_point3[2])

        # test errors and exceptions
        with self.assertRaises(ValueError):
            crv.ReParametrize(9, 3)
        with self.assertRaises(TypeError):
            crv.ReParametrize("one", "two")

    def test_split(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0,0,1],  [1,1,0],  [2,-1,0],  [3,0,0]]
        crv = Curve(3, [0, 0, 0, 0.7, 1, 1, 1], controlpoints)

        # get some info on the initial curve
        evaluation_point1 = crv.Evaluate(0.50)
        evaluation_point2 = crv.Evaluate(0.70)
        evaluation_point3 = crv.Evaluate(0.33)

        # split curves away from knot
        new_curves_050 = crv.Split(0.50)
        self.assertEqual(len(new_curves_050), 2)
        self.assertEqual(len(new_curves_050[0].GetKnots(True)), 6) # open knot vector [0,0,0,.5,.5,.5]
        self.assertEqual(len(new_curves_050[1].GetKnots(True)), 7) # open knot vector [.5,.5,.5,.7,1,1,1]

        # split curves at existing knot
        new_curves_070 = crv.Split(0.70)
        self.assertEqual(len(new_curves_070), 2)
        self.assertEqual(len(new_curves_070[0].GetKnots(True)), 6) # open knot vector [0,0,0,.7,.7,.7]
        self.assertEqual(len(new_curves_070[1].GetKnots(True)), 6) # open knot vector [.7,.7,.7,1,1,1]

        # split curves multiple points
        new_curves_all = crv.Split([0.50, 0.70])
        self.assertEqual(len(new_curves_all), 3)
        self.assertEqual(len(new_curves_all[0].GetKnots(True)), 6) # open knot vector [0,0,0,.5,.5,.5]
        self.assertEqual(len(new_curves_all[1].GetKnots(True)), 6) # open knot vector [.5,.5,.5,.7,.7,.7]
        self.assertEqual(len(new_curves_all[2].GetKnots(True)), 6) # open knot vector [.7,.7,.7,1,1,1]
        
        # compare all curves which exist at parametric point 0.5 and 0.33
        for c in new_curves_050 + [new_curves_070[0]] + new_curves_all[0:2]:
            new_curve_evaluation = c.Evaluate(0.50)
            self.assertAlmostEqual(evaluation_point1[0], new_curve_evaluation[0])
            self.assertAlmostEqual(evaluation_point1[1], new_curve_evaluation[1])
            self.assertAlmostEqual(evaluation_point1[2], new_curve_evaluation[2])
            new_curve_evaluation = c.Evaluate(0.33)
            self.assertAlmostEqual(evaluation_point3[0], new_curve_evaluation[0])
            self.assertAlmostEqual(evaluation_point3[1], new_curve_evaluation[1])
            self.assertAlmostEqual(evaluation_point3[2], new_curve_evaluation[2])

        # compare all curves which exist at parametric point 0.7
        for c in [new_curves_050[1]] + new_curves_070 + new_curves_all[1:3]:
            new_curve_evaluation = c.Evaluate(0.70)
            self.assertAlmostEqual(evaluation_point2[0], new_curve_evaluation[0])
            self.assertAlmostEqual(evaluation_point2[1], new_curve_evaluation[1])
            self.assertAlmostEqual(evaluation_point2[2], new_curve_evaluation[2])
            
        # test errors and exceptions
        with self.assertRaises(TypeError):
            crv.Split(.1, .2, .3)    # too many arguments
        with self.assertRaises(ValueError):
            crv.Split("tree-fiddy")  # wrong argument type
        with self.assertRaises(Exception):
            crv.Split(-0.2)          # GoTools returns error on outside-domain errors
        with self.assertRaises(Exception):
            crv.Split( 1.4)          # GoTools returns error on outside-domain errors


if __name__ == '__main__':
    unittest.main()
