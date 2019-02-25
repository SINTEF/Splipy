# -*- coding: utf-8 -*-

from splipy import BSplineBasis, Curve
import splipy.curve_factory as CurveFactory
from math import pi, sqrt, cos, sin
import numpy as np
from numpy.linalg import norm
import unittest
try:
    import nutils
    has_nutils = True
except ImportError:
    has_nutils = False

class TestCurveFactory(unittest.TestCase):
    def test_line(self):

        # 2D line
        c = CurveFactory.line([1, 1], [2, 0])
        self.assertEqual(c.order(0), 2)  # linear discretization
        self.assertEqual(len(c), 2)  # two control points
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c[0][0], 1)
        self.assertEqual(c[0][1], 1)
        self.assertEqual(c[-1][0], 2)
        self.assertEqual(c[-1][1], 0)

        # 3D line
        c = CurveFactory.line([1, 2, 3], [8, 7, 6])
        self.assertEqual(c.order(0), 2)  # linear discretization
        self.assertEqual(len(c), 2)  # two control points
        self.assertEqual(c.dimension, 3)

    def test_n_gon(self):
        # test default 5 side n-gon
        c = CurveFactory.n_gon()
        self.assertEqual(len(c), 5)
        self.assertEqual(len(c.knots(0, True)), 8)
        self.assertEqual(c.order(0), 2)
        # evaluate at second corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.end(0) / 5.0)[0], cos(2 * pi / 5))
        self.assertAlmostEqual(c.evaluate(c.end(0) / 5.0)[1], sin(2 * pi / 5))
        # evaluate at fourh corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.end(0) / 5.0 * 4)[0], cos(2 * pi / 5 * 4))
        self.assertAlmostEqual(c.evaluate(c.end(0) / 5.0 * 4)[1], sin(2 * pi / 5 * 4))

        # test a radius 3 septagon
        c = CurveFactory.n_gon(n=7, r=3)
        self.assertEqual(len(c), 7)
        # evaluate at third corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.end(0) / 7.0)[0], 3 * cos(2 * pi / 7))
        self.assertAlmostEqual(c.evaluate(c.end(0) / 7.0)[1], 3 * sin(2 * pi / 7))

        # test errors
        with self.assertRaises(ValueError):
            c = CurveFactory.n_gon(r=-2.5)
        with self.assertRaises(ValueError):
            c = CurveFactory.n_gon(n=1)

    def test_polygon(self):
        pts = [[1,0], [1,1], [0,1], [0,2], [6,2]]
        c = CurveFactory.polygon(pts)
        expected_knots = [0,0,1,2,3,9,9]
        actual_knots   = c.knots(0,True)

        self.assertEqual(len(c),      5)
        self.assertEqual(c.order(0),  2)
        self.assertEqual(c.dimension, 2)
        self.assertAlmostEqual(norm(expected_knots - actual_knots), 0.0)

        c = CurveFactory.polygon([0,0], [1,0], [0,1], [-1,0], relative=True)
        self.assertEqual(len(c), 4)
        self.assertAlmostEqual(c[2][0], 1)
        self.assertAlmostEqual(c[2][1], 1)
        self.assertAlmostEqual(c[3][0], 0)
        self.assertAlmostEqual(c[3][1], 1)

    def test_circle(self):

        # unit circle of radius 1
        c = CurveFactory.circle()
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt), 1.0)  # check radius=1
        # test evaluation at key points
        pt = c.evaluate(0)  # evaluate at 0
        self.assertAlmostEqual(pt[0], 1.0)
        self.assertAlmostEqual(pt[1], 0.0)
        pt = c.evaluate(c.end(0) / 4)  # evaluate at pi/2
        self.assertAlmostEqual(pt[0], 0.0)
        self.assertAlmostEqual(pt[1], 1.0)
        pt = c.evaluate(c.end(0) / 2)  # evaluate at pi
        self.assertAlmostEqual(pt[0], -1.0)
        self.assertAlmostEqual(pt[1], 0.0)
        pt = c.evaluate(c.end(0) * 2)  # evaluate outside domain (test periodic)
        self.assertAlmostEqual(pt[0], 1.0)
        self.assertAlmostEqual(pt[1], 0.0)

        # circle of radius different from 1
        c = CurveFactory.circle(3)
        # test evaluation at 25 points for radius=3, outside domain
        t = np.linspace(c.start(0) - 3, c.end(0) + 2, 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt, 2), 3.0)  # check radius=3
        self.assertAlmostEqual(c.length(), 6*pi, places=3)

        # circle not at origin
        c = CurveFactory.circle(1, center=(1,0,0), normal=(1,1,1))
        # test evaluation at 25 points
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt-[1,0,0], 2), 1.0)  # check radius=1
            self.assertAlmostEqual(pt[0]+pt[1]+pt[2] - 1, 0.0) # in plane x+y+z=1
        self.assertAlmostEqual(c.length(), 2*pi, places=3)

        # test alt circle
        c = CurveFactory.circle(r=3, type='p4C1')
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        self.assertTrue(np.allclose(x[:,0]**2 + x[:,1]**2, 3.0**2))
        for k in c.knots(0):
            self.assertEqual(c.continuity(k), 1)

        # test circle with different xaxis
        c = CurveFactory.circle(xaxis=(1,1,0))
        t = np.linspace(c.start(0), c.end(0), 5)
        x = c.evaluate(t)
        self.assertTrue(np.allclose(x[:,0]**2 + x[:,1]**2, 1.0**2))
        self.assertTrue(np.allclose(x[0,:], [ 1/sqrt(2), 1/sqrt(2)]))
        self.assertTrue(np.allclose(x[1,:], [-1/sqrt(2), 1/sqrt(2)]))

        # test using all parameters (in x+y+z=1 plane)
        c = CurveFactory.circle(r=sqrt(2), normal=(1,1,1), center=(1,0,0), xaxis=(-1,1,0))
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        self.assertTrue(np.allclose((x[:,0]-1)**2 + x[:,1]**2 + x[:,2]**2, 2.0))
        self.assertTrue(np.allclose(c(  0 ), [ 0,1,0]))
        self.assertTrue(np.allclose(c( pi ), [2,-1,0]))

        # test errors and exceptions
        with self.assertRaises(ValueError):
            c = CurveFactory.circle(-2.5)  # negative radius

    def test_circle_segment(self):

        # basic circle segment
        c = CurveFactory.circle_segment(pi * 0.9)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt, 2), 1.0)  # check radius=1

        # radius 7 circle segment
        c = CurveFactory.circle_segment(pi * 1.87, 7)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=7
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt, 2), 7.0)  # check radius=7

        # negative theta
        c = CurveFactory.circle_segment(-pi/2)
        self.assertEqual(c.rational, True)
        self.assertTrue(np.allclose(c(0),     [1,0]))
        self.assertTrue(np.allclose(c(-pi/4), [1/sqrt(2),-1/sqrt(2)]))
        self.assertTrue(np.allclose(c(-pi/2), [0,-1]))

        # boundary case with one knot span circle segment
        c = CurveFactory.circle_segment(2 * pi / 3)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        self.assertEqual(len(c.knots(0)), 2)
        self.assertFalse(c.periodic(0))
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt, 2), 1.0)  # check radius=1

        # boundary case with full circle
        c = CurveFactory.circle_segment(2 * pi)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        self.assertEqual(len(c.knots(0)), 5)
        self.assertTrue(c.periodic(0))
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt, 2), 1.0)  # check radius=1

        # test errors and exceptions
        with self.assertRaises(ValueError):
            c = CurveFactory.circle_segment(3 * pi)  # outside domain
        with self.assertRaises(ValueError):
            c = CurveFactory.circle_segment(-3 * pi)  # outside domain
        with self.assertRaises(ValueError):
            c = CurveFactory.circle_segment(pi, -2)  # negative radius

    def test_circle_segment_from_three_points(self):
        # quarter circle (xy-plane)
        c = CurveFactory.circle_segment_from_three_points([1,0], [1.0/sqrt(2), 1.0/sqrt(2)], [0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        self.assertEqual(c.dimension, 2)
        self.assertTrue(c.rational)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt), 1.0)  # check radius=1
        self.assertTrue(np.allclose(c[0],  [1,0,1]))         # check endpoints
        self.assertTrue(np.allclose(c[-1], [0,1,1]))

        # quarter circle (x=y plane)
        c = CurveFactory.circle_segment_from_three_points([1.0/sqrt(2), 1.0/sqrt(2),0], [.5, .5, 1/sqrt(2)], [0,0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        self.assertEqual(c.dimension, 3)
        self.assertTrue(c.rational)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt), 1.0)  # check radius=1
        self.assertTrue(np.allclose(x[:,0], x[:,1]))         # check x=y plane
        self.assertTrue(np.allclose(c[-1], [0,0,1,1]))       # check endpoints

        # one-eight circle ([1,-1,0] normal, center in (0,0,0) )
        c = CurveFactory.circle_segment_from_three_points([.5, .5, 1/sqrt(2)], [2**(-3.0/2), 2**(-3.0/2), sqrt(3)/2], [0,0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt), 1.0)  # check radius=1
        self.assertTrue(np.allclose(x[:,0], x[:,1]))         # check x=y plane
        self.assertTrue(np.allclose(c[-1], [0,0,1,1]))       # check endpoints

        # one-eight circle ([1,-1,0] normal, center in (1,0,0))
        c = CurveFactory.circle_segment_from_three_points([1.5, .5, 1/sqrt(2)], [1+2**(-3.0/2), 2**(-3.0/2), sqrt(3)/2], [1,0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt-np.array([1,0,0])), 1.0)  # check radius=1
        self.assertTrue(np.allclose(x[:,0]-1, x[:,1]))       # check (x-1)=y plane
        self.assertTrue(np.allclose(c[-1], [1,0,1,1]))       # check endpoints

        # theta > pi
        c = CurveFactory.circle_segment_from_three_points([4,0], [0,4], [0,-4])
        self.assertTrue(np.allclose(c[0],  [4, 0, 1])) # check startpoint
        self.assertTrue(np.allclose(c[-1], [0,-4, 1])) # check endpoint
        t = (c.start(0)+c.end(0))/2.0
        self.assertTrue(np.allclose(c(t), [-2*sqrt(2),2*sqrt(2)])) # midpoint evaluation
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt), 4.0)  # check radius=4

        # theta > pi in x=3/4 y plane (r=5), center=(4,3,0)
        c = CurveFactory.circle_segment_from_three_points([4,3,-5], [8,6,0], [0,0,0])
        self.assertTrue(np.allclose(c[0],  [4, 3,-5, 1])) # check startpoint
        self.assertTrue(np.allclose(c[-1], [0, 0, 0, 1])) # check endpoint
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(norm(pt-np.array([4,3,0])), 5.0)  # check radius=5
        # check plane normal
        self.assertTrue(np.allclose(c.binormal(0),  [ 3.0/5,-4.0/5, 0]))
        self.assertTrue(np.allclose(c.binormal(.5), [ 3.0/5,-4.0/5, 0])) 


    def test_cubic_curve(self):
        t = np.linspace(0,1,80)  # interpolation points
        s = np.linspace(0,1,150) # evaluation points
        x = np.zeros((80,2))     # physical points

        # test PERIODIC curves
        x[:,0] = 16*t*t*(1-t)*(1-t)
        x[:,1] = 1-x[:,0]
        crv = CurveFactory.cubic_curve(x, CurveFactory.Boundary.PERIODIC, t)
        y = crv(s)
        # exact solution is t^4, cubic approximation gives t^3, needs atol
        self.assertTrue(np.allclose(y[:,0],   16*s*s*(1-s)*(1-s), atol=1e-2))
        self.assertTrue(np.allclose(y[:,1], 1-16*s*s*(1-s)*(1-s), atol=1e-2))

        # test FREE boundary type
        x[:,0] = 12*t*t*t         + 2*t
        x[:,1] = 12*t*t*t - 3*t*t 
        crv = CurveFactory.cubic_curve(x, CurveFactory.Boundary.FREE, t=t)
        y = crv(s)
        self.assertTrue(np.allclose(y[:,0], 12*s*s*s         + 2*s))
        self.assertTrue(np.allclose(y[:,1], 12*s*s*s - 3*s*s      ))

        # test NATURAL boundary type (x''(t)=0 at the boundary)
        x[:,0] = t
        x[:,1] = 3
        crv = CurveFactory.cubic_curve(x, CurveFactory.Boundary.NATURAL, t=t)
        y = crv(s)
        self.assertTrue(np.allclose(y[:,0],   s))
        self.assertTrue(np.allclose(y[:,1], 3.0))

        # test TANGENT boundary type (x'(t)=g at both endpoints)
        x[:,0] =   t*t +   t - 1
        x[:,1] = 3*t*t - 3*t + 1
        dx     = [[2*t[ 0] + 1, 6*t[ 0]-3],
                  [2*t[-1] + 1, 6*t[-1]-3]]
        crv = CurveFactory.cubic_curve(x, CurveFactory.Boundary.TANGENT, t=t, tangents=dx)
        y = crv(s)
        self.assertTrue(np.allclose(y[:,0],   s*s +   s - 1))
        self.assertTrue(np.allclose(y[:,1], 3*s*s - 3*s + 1))

        # test HERMITE boundary type (x'(t)=g(t) at all interior points)
        x[:,0] =   t*t +   t - 1
        x[:,1] = 3*t*t - 3*t + 1
        dx     = np.vstack([[2*t + 1], [6*t-3]]).T
        crv = CurveFactory.cubic_curve(x, CurveFactory.Boundary.HERMITE, t=t, tangents=dx)
        y = crv(s)
        self.assertTrue(np.allclose(y[:,0],   s*s +   s - 1))
        self.assertTrue(np.allclose(y[:,1], 3*s*s - 3*s + 1))

    def test_bezier(self):
        crv = CurveFactory.bezier([[0,0], [0,1], [1,1], [1,0], [2,0], [2,1],[1,1]])
        self.assertEqual(len(crv.knots(0)), 3)
        self.assertTrue(np.allclose(crv(0), [0,0]))
        t = crv.tangent(0)
        self.assertTrue(np.allclose(t/norm(t), [0,1]))
        t = crv.tangent(1, above=False)
        self.assertTrue(np.allclose(t/norm(t), [0,-1]))
        t = crv.tangent(1, above=True)
        self.assertTrue(np.allclose(t/norm(t), [1,0]))
        self.assertTrue(np.allclose(crv(1), [1,0]))
        self.assertTrue(crv.order(0), 4)

        # test the exact same curve, only with relative keyword
        crv = CurveFactory.bezier([[0,0], [0,1], [1,0], [0,-1], [1,0], [0,1],[1,0]], relative=True)
        self.assertEqual(len(crv.knots(0)), 3)
        self.assertTrue(np.allclose(crv(0), [0,0]))
        t = crv.tangent(0)
        self.assertTrue(np.allclose(t/norm(t), [0,1]))
        t = crv.tangent(1, above=False)
        self.assertTrue(np.allclose(t/norm(t), [0,-1]))
        t = crv.tangent(1, above=True)
        self.assertTrue(np.allclose(t/norm(t), [1,0]))
        self.assertTrue(np.allclose(crv(1), [1,0]))
        self.assertTrue(crv.order(0), 4)

    def test_ellipse(self):
        # test (x/1)^2 + (y/5)^2 = 1
        c = CurveFactory.ellipse(1,5)
        t = np.linspace(c.start(0), c.end(0), 25)
        for pt in c.evaluate(t):
            x,y = pt
            self.assertAlmostEqual((x/1)**2 + (y/5)**2, 1)

        # test (x-.3/1)^2 + (y-6/5)^2 = 1
        c = CurveFactory.ellipse(1,5, center=(.3, 6))
        t = np.linspace(c.start(0), c.end(0), 25)
        for pt in c.evaluate(t):
            x,y = pt
            self.assertAlmostEqual(((x-.3)/1)**2 + ((y-6)/5)**2, 1)

        # test ellipse along x=y axis
        c = CurveFactory.ellipse(1,2, xaxis=(1,1))
        t = np.linspace(c.start(0), c.end(0), 25)
        for pt in c.evaluate(t):
            x,y = pt
            s = 1/sqrt(2)
            self.assertAlmostEqual(((s*x + s*y)/1)**2 + ((s*x - s*y)/2)**2, 1)

        # test ellipse in 3D
        c = CurveFactory.ellipse(1,2, normal=(0,1,0), xaxis=(1,0,1))
        t = np.linspace(c.start(0), c.end(0), 25)
        for pt in c.evaluate(t):
            x,y,z = pt
            s = 1/sqrt(2)
            self.assertAlmostEqual(((s*x + s*z)/1)**2 + ((s*x - s*z)/2)**2, 1)

if __name__ == '__main__':
    unittest.main()
