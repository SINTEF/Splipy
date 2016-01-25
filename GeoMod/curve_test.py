from GeoMod import Curve, BSplineBasis
from math import sqrt
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
        crv.flip_parametrization()
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

        crv
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

        # check logic error for negative argument (gotools cast this error)
        with self.assertRaises(Exception):
            crv.raise_order(-1)

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
        crv = Curve(BSplineBasis(3, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]), controlpoints, True)
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
            crv.insert_knot(1, 2)  # too many arguments
        with self.assertRaises(TypeError):
            crv.insert_knot()  # too few arguments
        with self.assertRaises(ValueError):
            crv.insert_knot(-0.2)  # Outside-domain error
        with self.assertRaises(ValueError):
            crv.insert_knot(4.4)  # Outside-domain error

    def test_reparametrize(self):
        # non-uniform knot vector of a squiggly quadratic n=4 curve
        controlpoints = [[0, 0, 0], [1, 1, 0], [2, -1, 0], [3, 0, 0]]
        crv = Curve(BSplineBasis(3, [0, 0, 0, 1.32, 3, 3, 3]), controlpoints)

        # get some info on the initial curve
        knots1 = crv.knots(0)
        evaluation_point1 = crv(1.20)
        self.assertEqual(knots1[0], 0)
        self.assertEqual(knots1[-1], 3)

        # reparametrize
        crv.reparametrize((6.0, 9.0))

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
        crv.reparametrize()

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
            crv.reparametrize((9, 3))
        with self.assertRaises(TypeError):
            crv.reparametrize(("one", "two"))

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


if __name__ == '__main__':
    unittest.main()
