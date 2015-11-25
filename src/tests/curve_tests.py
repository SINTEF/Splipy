from GoTools import Curve
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
        control_point1    = crv[0]
        crv2              = crv.RaiseOrder(2)
        evaluation_point2 = crv2.Evaluate(0.37)
        control_point2    = crv2[0]

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


if __name__ == '__main__':
    unittest.main()
