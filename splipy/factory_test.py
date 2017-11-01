# -*- coding: utf-8 -*-

from splipy import BSplineBasis, Curve, Surface, Volume
import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
import splipy.volume_factory as VolumeFactory
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
        self.assertAlmostEqual(np.linalg.norm(expected_knots - actual_knots), 0.0)

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
            self.assertAlmostEqual(np.linalg.norm(pt), 1.0)  # check radius=1
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
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 3.0)  # check radius=3
        self.assertAlmostEqual(c.length(), 6*pi, places=3)

        # circle not at origin
        c = CurveFactory.circle(1, center=(1,0,0), normal=(1,1,1))
        # test evaluation at 25 points
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt-[1,0,0], 2), 1.0)  # check radius=1
            self.assertAlmostEqual(pt[0]+pt[1]+pt[2] - 1, 0.0) # in plane x+y+z=1
        self.assertAlmostEqual(c.length(), 2*pi, places=3)

        # test alt circle
        c = CurveFactory.circle(r=3, type='p4C1')
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        self.assertTrue(np.allclose(x[:,0]**2 + x[:,1]**2, 3.0**2))
        for k in c.knots(0):
            self.assertEqual(c.continuity(k), 1)


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
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 1.0)  # check radius=1

        # radius 7 circle segment
        c = CurveFactory.circle_segment(pi * 1.87, 7)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=7
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 7.0)  # check radius=7

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
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 1.0)  # check radius=1

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
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 1.0)  # check radius=1

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
            self.assertAlmostEqual(np.linalg.norm(pt), 1.0)  # check radius=1
        self.assertTrue(np.allclose(c[0],  [1,0,1]))         # check endpoints
        self.assertTrue(np.allclose(c[-1], [0,1,1]))

        # quarter circle (x=y plane)
        c = CurveFactory.circle_segment_from_three_points([1.0/sqrt(2), 1.0/sqrt(2),0], [.5, .5, 1/sqrt(2)], [0,0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        self.assertEqual(c.dimension, 3)
        self.assertTrue(c.rational)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt), 1.0)  # check radius=1
        self.assertTrue(np.allclose(x[:,0], x[:,1]))         # check x=y plane
        self.assertTrue(np.allclose(c[-1], [0,0,1,1]))       # check endpoints

        # one-eight circle ([1,-1,0] normal, center in (0,0,0) )
        c = CurveFactory.circle_segment_from_three_points([.5, .5, 1/sqrt(2)], [2**(-3.0/2), 2**(-3.0/2), sqrt(3)/2], [0,0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt), 1.0)  # check radius=1
        self.assertTrue(np.allclose(x[:,0], x[:,1]))         # check x=y plane
        self.assertTrue(np.allclose(c[-1], [0,0,1,1]))       # check endpoints

        # one-eight circle ([1,-1,0] normal, center in (1,0,0))
        c = CurveFactory.circle_segment_from_three_points([1.5, .5, 1/sqrt(2)], [1+2**(-3.0/2), 2**(-3.0/2), sqrt(3)/2], [1,0,1])
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt-np.array([1,0,0])), 1.0)  # check radius=1
        self.assertTrue(np.allclose(x[:,0]-1, x[:,1]))       # check (x-1)=y plane
        self.assertTrue(np.allclose(c[-1], [1,0,1,1]))       # check endpoints

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


class TestSurfaceFactory(unittest.TestCase):
    def test_square(self):
        surf = SurfaceFactory.square((4, 5))
        self.assertEqual(surf.dimension, 2)
        self.assertEqual(surf.rational, False)
        self.assertEqual(surf.order(), (2, 2))

    def test_curve_interpolation(self):
        basis = BSplineBasis(4, [0, 0, 0, 0, .3, .9, 1, 1, 1, 1])
        t = np.array(basis.greville())
        # create the mapping (x,y,z)=(t^2, 1-t, t^3+2*t)
        x_pts = np.zeros((len(t), 3))
        x_pts[:, 0] = t * t
        x_pts[:, 1] = 1 - t
        x_pts[:, 2] = t * t * t + 2 * t
        crv = CurveFactory.interpolate(x_pts, basis)
        self.assertEqual(crv.order(0), 4)
        self.assertAlmostEqual(crv(.4)[0], .4**2)  # x=t^2
        self.assertAlmostEqual(crv(.4)[1], 1 - .4)  # y=1-t
        self.assertAlmostEqual(crv(.4)[2], .4**3 + 2 * .4)  # z=t^3+2t
        self.assertAlmostEqual(crv(.5)[0], .5**2)  # x=t^2
        self.assertAlmostEqual(crv(.5)[1], 1 - .5)  # y=1-t
        self.assertAlmostEqual(crv(.5)[2], .5**3 + 2 * .5)  # z=t^3+2t

    def test_disc(self):
        # radial disc
        surf = SurfaceFactory.disc()
        x = surf.evaluate([0, 1], [0, pi / 4, pi / 2, pi])
        self.assertAlmostEqual(x[0][0][0], 0)
        self.assertAlmostEqual(x[0][0][1], 0)
        self.assertAlmostEqual(x[1][0][0], 1)
        self.assertAlmostEqual(x[1][0][1], 0)
        self.assertAlmostEqual(x[1][1][0], 1 / sqrt(2))
        self.assertAlmostEqual(x[1][1][1], 1 / sqrt(2))
        self.assertAlmostEqual(x[1][2][0], 0)
        self.assertAlmostEqual(x[1][2][1], 1)
        self.assertAlmostEqual(surf.area(), pi, places=3)

        # radial disc of size different from 1
        surf = SurfaceFactory.disc(4)
        # test evaluation at 25 points for radius=4
        v = np.linspace(surf.start('v'), surf.end('v'),25)
        u = np.linspace(surf.start('u'), surf.end('u'), 5)
        x = surf.evaluate(u, v)
        for circles, i in zip(x, range(5)):
            for pt in circles:
                self.assertAlmostEqual(norm(pt, 2), 4.0 * i / 4)  # check radius
        self.assertAlmostEqual(surf.area(), 4*4*pi, places=3)

        # square disc
        surf = SurfaceFactory.disc(3, type='square')
        # evaluate on all 4 edges, 5 pts on each edge
        u = np.linspace(0, 1, 5)
        v = np.linspace(0, 1, 5)
        x = surf.evaluate(u, v)
        for pt in np.array(x[0, :, :]):  # umin edge
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 3.0)  # check radius
        for pt in np.array(x[-1, :, :]):  # umax edge
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 3.0)  # check radius
        for pt in np.array(x[:, 0, :]):  # vmin edge
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 3.0)  # check radius
        for pt in np.array(x[:, -1, :]):  # vmax edge
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 3.0)  # check radius
        self.assertAlmostEqual(surf.area(), 3*3*pi, places=2)

    def test_revolve(self):
        # square torus
        square = CurveFactory.n_gon(4)
        square.rotate(pi / 2, (1, 0, 0))
        square.translate((2, 0, 0))  # in xz-plane with corners at (3,0),(2,1),(1,0),(2,-1)
        surf = SurfaceFactory.revolve(square)
        surf.reparam()  # set parametric space to (0,1)^2
        v = np.linspace(0, 1, 13)
        x = surf.evaluate(0, v)  # outer ring evaluation u=0
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 3.0)  # check radius=3
        x = surf.evaluate(.25, v)  # top ring evaluation u=.25
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(pt[0] * pt[0] + pt[1] * pt[1], 2 * 2)  # check radius=2
            self.assertAlmostEqual(pt[2], 1)  # check height=1
        x = surf.evaluate(.375, v)  # mid inner ring evaluation u=.375
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(pt[0] * pt[0] + pt[1] * pt[1], 1.5 * 1.5)  # check radius=1.5
            self.assertAlmostEqual(pt[2], .5)  # check height=0.5

        # incomplete revolve
        c    = CurveFactory.line([1,0], [0,1], relative=True)
        surf = SurfaceFactory.revolve(c, theta=4.2222, axis=[0,1,0])
        surf.reparam()
        u = np.linspace(0,1,7)
        v = np.linspace(0,1,7)
        x = surf(u,v)
        for uPt in x:
            for pt in uPt:
                self.assertAlmostEqual(pt[0]**2 + pt[2]**2, 1.0) # radius 1 from y-axis


    def test_surface_torus(self):
        # default torus
        surf = SurfaceFactory.torus(1, 3)
        start = surf.start()
        end = surf.end()
        # check a 13 evaluation points of v-evaluation (around the z-axis)
        v = np.linspace(start[1], end[1], 13)
        # check minor-circle u=0 (outmost ring)
        x = surf.evaluate(0, v)
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(pt[0] * pt[0] + pt[1] * pt[1], 4 * 4)  # check radius=4
            self.assertAlmostEqual(pt[2], 0)  # check height=0
        # check minor-circle u=pi (innermost ring)
        x = surf.evaluate(pi, v)
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(pt[0] * pt[0] + pt[1] * pt[1], 2 * 2)  # check radius=2
            self.assertAlmostEqual(pt[2], 0)  # check height=0
        # check minor-circle u=pi/2 (topmost ring)
        x = surf.evaluate(pi / 2, v)
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(pt[0] * pt[0] + pt[1] * pt[1], 3 * 3)  # check radius=3
            self.assertAlmostEqual(pt[2], 1)  # check height=1
        # check minor-circle u=3*pi/2 (mid-evaluation)
        x = surf.evaluate(3 * pi / 4, v)
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(pt[0] * pt[0] + pt[1] * pt[1],
                                   (3 - 1.0 / sqrt(2))**2)  # check radius=3-1/sqrt(2)
            self.assertAlmostEqual(pt[2], 1.0 / sqrt(2))  # check height=1/sqrt(2)

    def test_sphere(self):
        # unit ball
        surf = SurfaceFactory.sphere()
        # test 7x7 grid for radius = 1
        for u in np.linspace(surf.start()[0], surf.end()[0], 7):
            for v in np.linspace(surf.start()[1], surf.end()[1], 7):
                self.assertAlmostEqual(np.linalg.norm(surf(u, v), 2), 1.0)
        self.assertAlmostEqual(surf.area(), 4*pi, places=3)

    def test_volume_sphere(self):
        # unit ball
        ball = VolumeFactory.sphere(type='square')
        # test 7x7 grid for radius = 1
        for surf in ball.faces():
            for u in np.linspace(surf.start(0), surf.end(0), 7):
                for v in np.linspace(surf.start(1), surf.end(1), 7):
                    self.assertAlmostEqual(np.linalg.norm(surf(u, v), 2), 1.0)
            self.assertAlmostEqual(surf.area(), 4*pi/6, places=3)

        # unit ball
        ball = VolumeFactory.sphere(type='radial')
        # test 7x7 grid for radius = 1
        for u in np.linspace(surf.start(0), surf.end(0), 7):
            for v in np.linspace(surf.start(1), surf.end(1), 7):
                self.assertAlmostEqual(np.linalg.norm(ball(u, v, 0), 2), 1.0) # w=0 is outer shell
                self.assertAlmostEqual(np.linalg.norm(ball(u, v, 1), 2), 0.0) # w=1 is the degenerate core
        self.assertAlmostEqual(ball.faces()[4].area(), 4*pi, places=3)

    def test_cylinder_surface(self):
        # unit cylinder
        surf = SurfaceFactory.cylinder()
        # test 7x7 grid for xy-radius = 1 and v=z
        for u in np.linspace(surf.start()[0], surf.end()[0], 7):
            for v in np.linspace(surf.start()[1], surf.end()[1], 7):
                x = surf(u, v)
                self.assertAlmostEqual(np.linalg.norm(x[:2], 2), 1.0)  # (x,y) coordinates to z-axis
                self.assertAlmostEqual(x[2], v)  # z coordinate should be linear
        self.assertAlmostEqual(surf.area(), 2*pi, places=3)

        # cylinder with parameters
        surf = SurfaceFactory.cylinder(r=2, h=5, center=(0,0,1), axis=(1,0,0))
        for u in np.linspace(surf.start()[0], surf.end()[0], 7):
            for v in np.linspace(surf.start()[1], surf.end()[1], 7):
                x = surf(u, v)
                self.assertAlmostEqual(x[1]**2+(x[2]-1)**2, 2.0**2) # distance to (z-1)=y=0
                self.assertAlmostEqual(x[0], v*5)                   # x coordinate should be linear
        self.assertAlmostEqual(surf.area(), 2*2*pi*5, places=3)

    def test_edge_curves(self):
        # create an arrow-like 2D geometry with the pointy end at (-1,1) towards up and left
        # mixes rational and non-rational curves with different parametrization spaces
        c1 = CurveFactory.circle_segment(pi / 2)
        c2 = Curve(BSplineBasis(2, [0, 0, 1, 2, 2]), [[0, 1], [-1, 1], [-1, 0]])
        c3 = CurveFactory.circle_segment(pi / 2)
        c3.rotate(pi)
        c4 = Curve(BSplineBasis(2), [[0, -1], [1, 0]])

        surf = SurfaceFactory.edge_curves(c1, c2, c3, c4)

        # srf spits out parametric space (0,1)^2, so we sync these up to input curves
        c3.reverse()
        c4.reverse()
        c1.reparam()
        c2.reparam()
        c3.reparam()
        c4.reparam()

        for u in np.linspace(0, 1, 7):
            self.assertAlmostEqual(surf(u, 0)[0], c1(u)[0])  # x-coord, bottom crv
            self.assertAlmostEqual(surf(u, 0)[1], c1(u)[1])  # y-coord, bottom crv
        for u in np.linspace(0, 1, 7):
            self.assertAlmostEqual(surf(u, 1)[0], c3(u)[0])  # x-coord, top crv
            self.assertAlmostEqual(surf(u, 1)[1], c3(u)[1])  # y-coord, top crv
        for v in np.linspace(0, 1, 7):
            self.assertAlmostEqual(surf(0, v)[0], c4(v)[0])  # x-coord, left crv
            self.assertAlmostEqual(surf(0, v)[1], c4(v)[1])  # y-coord, left crv
        for v in np.linspace(0, 1, 7):
            self.assertAlmostEqual(surf(1, v)[0], c2(v)[0])  # x-coord, right crv
            self.assertAlmostEqual(surf(1, v)[1], c2(v)[1])  # y-coord, right crv

        # add a case where opposing sites have mis-matching rationality
        crvs = Surface().edges() # returned in order umin, umax, vmin, vmax
        crvs[0].force_rational()
        crvs[1].reverse()
        crvs[2].reverse()
        # input curves should be clockwise oriented closed loop
        srf = SurfaceFactory.edge_curves(crvs[0], crvs[3], crvs[1], crvs[2])
        crvs[1].reverse()
        u = np.linspace(0,1,7)
        self.assertTrue(np.allclose(srf(u,0).reshape((7,2)), crvs[0](u)))
        self.assertTrue(np.allclose(srf(u,1).reshape((7,2)), crvs[1](u)))

        # test self-organizing curve ordering when they are not sequential
        srf = SurfaceFactory.edge_curves(crvs[0], crvs[2].reverse(), crvs[3], crvs[1])
        u = np.linspace(0,1,7)
        self.assertTrue(np.allclose(srf(u,0).reshape((7,2)), crvs[0](u)))
        self.assertTrue(np.allclose(srf(u,1).reshape((7,2)), crvs[1](u)))

        # test error handling
        with self.assertRaises(ValueError):
            srf = SurfaceFactory.edge_curves(crvs + (Curve(),)) # 5 input curves

    @unittest.skipIf(not has_nutils, "EdgeCurves with poisson solver requires nutils")
    def test_edgecurves(self):
        # create an arrow-like 2D geometry with the pointy end at (-1,1) towards up and left
        # rebuild to avoid rational representations
        c1 = CurveFactory.circle_segment(pi / 2).rebuild(3,11)
        c2 = Curve(BSplineBasis(2, [0, 0, 1, 2, 2]), [[0, 1], [-1, 1], [-1, 0]])
        c3 = CurveFactory.circle_segment(pi / 2).rebuild(3,11)
        c3.rotate(pi)
        c4 = Curve(BSplineBasis(2), [[0, -1], [1, 0]]).rebuild(3,10)
        c4 = c4.rebuild(4,11)

        surf = SurfaceFactory.edge_curves([c1, c2, c3, c4], type='poisson')
        
        # check right dimensions of output
        self.assertEqual(surf.shape[0], 11) # 11 controlpoints in the circle segment
        self.assertEqual(surf.shape[1], 13) # 11 controlpoints in c4, +2 for C0-knot in c1
        self.assertEqual(surf.order(0), 3)
        self.assertEqual(surf.order(1), 4)

        # check symmetry: all interior points lie on the y=-x line
        u = np.linspace(0,1,7)
        pts = surf(u,0.5)
        for x in pts[:,0,:]:
            self.assertAlmostEqual(x[0], -x[1])

        # check that c1 edge conforms to surface edge
        u = np.linspace(0,1,7)
        c1.reparam()
        pts_surf = surf(u,0.0)
        pts_c1   = c1(u)
        for (xs,xc) in zip(pts_surf[:,0,:], pts_c1):
            self.assertTrue(np.allclose(xs, xc))

        # check that c2 edge conforms to surface edge
        v = np.linspace(0,1,7)
        c2.reparam()
        pts_surf = surf(1.0,v)
        pts_c2   = c2(v)
        for (xs,xc) in zip(pts_surf[:,0,:], pts_c2):
            self.assertTrue(np.allclose(xs, xc))

    def test_thicken(self):
        c = Curve()                       # 2D curve from (0,0) to (1,0)
        s = SurfaceFactory.thicken(c, .5) # extend to y=[-.5, .5]
        self.assertTupleEqual(s.order(), (2,2))
        self.assertTupleEqual(s.start(), (0,0))
        self.assertTupleEqual(s.end(),   (1,1))
        self.assertTupleEqual(s.bounding_box()[0], (0.0,1.0))
        self.assertTupleEqual(s.bounding_box()[1], (-.5, .5))

        # test a case with vanishing velocity. x'(t)=0, y'(t)=0 for t=0
        c = Curve(BSplineBasis(3), [[0,0],[0,0],[1,0]]) # x(t)=t^2, y(t)=0
        s = SurfaceFactory.thicken(c, .5)
        self.assertTupleEqual(s.order(), (3,2))
        self.assertTupleEqual(s.start(), (0,0))
        self.assertTupleEqual(s.end(),   (1,1))
        self.assertTupleEqual(s.bounding_box()[0], (0.0,1.0))
        self.assertTupleEqual(s.bounding_box()[1], (-.5, .5))

        def myThickness(t):
            return t**2
        c = Curve(BSplineBasis(3))
        s = SurfaceFactory.thicken(c, myThickness)
        self.assertTupleEqual(s.order(), (3,2))
        self.assertTupleEqual(s.start(), (0,0))
        self.assertTupleEqual(s.end(),   (1,1))
        self.assertTupleEqual(s.bounding_box()[0], ( 0.0, 1.0))
        self.assertTupleEqual(s.bounding_box()[1], (-1.0, 1.0))

    def test_interpolate(self):
        t = np.linspace(0, 1, 7)
        V,U = np.meshgrid(t,t)
        x = np.zeros((7,7,3))
        x[:,:,0] = U*U + V*V
        x[:,:,1] = U - V
        x[:,:,2] = U*U*V*V*V - U*V + 3
        b1 = BSplineBasis(3, [0,0,0,.33,.66,.7, .8,1,1,1])
        b2 = BSplineBasis(4, [0,0,0,0,.1, .2, .5,1,1,1,1])

        surf = SurfaceFactory.interpolate(x, [b1,b2], [t,t])
        t = np.linspace(0, 1, 13)
        V,U = np.meshgrid(t,t)
        x = surf(t,t)

        self.assertTrue(np.allclose(x[:,:,0], U*U + V*V))
        self.assertTrue(np.allclose(x[:,:,1], U   - V))
        self.assertTrue(np.allclose(x[:,:,2], U*U*V*V*V - U*V + 3))

    def test_l2(self):
        t = np.linspace(0, 1, 110)
        V,U = np.meshgrid(t,t)
        x = np.zeros((110,110,3))
        x[:,:,0] = U*U + V*V
        x[:,:,1] = U - V
        x[:,:,2] = U*U*V*V*V - U*V + 3
        b1 = BSplineBasis(3, [0,0,0,.33,.66,.7, .8,1,1,1])
        b2 = BSplineBasis(4, [0,0,0,0,.1, .2, .5,1,1,1,1])

        surf = SurfaceFactory.least_square_fit(x, [b1,b2], [t,t])
        t = np.linspace(0, 1, 13)
        V,U = np.meshgrid(t,t)
        x = surf(t,t)

        self.assertTrue(np.allclose(x[:,:,0], U*U + V*V))
        self.assertTrue(np.allclose(x[:,:,1], U   - V))
        self.assertTrue(np.allclose(x[:,:,2], U*U*V*V*V - U*V + 3))




class TestVolumeFactory(unittest.TestCase):

    def test_cylinder_volume(self):
        # unit cylinder
        vol = VolumeFactory.cylinder()
        # test 5x5x5 grid for xy-radius = w and v=z
        for u in np.linspace(vol.start()[0], vol.end()[0], 5):
            for v in np.linspace(vol.start()[1], vol.end()[1], 5):
                for w in np.linspace(vol.start()[2], vol.end()[2], 5):
                    x = vol(u, v, w)
                    self.assertAlmostEqual(
                        np.linalg.norm(x[:2], 2), u)  # (x,y) coordinates to z-axis
                    self.assertAlmostEqual(x[2], w)  # z coordinate should be linear
        self.assertAlmostEqual(vol.volume(), pi, places=3)

        # cylinder with parameters
        vol = VolumeFactory.cylinder(r=2, h=5, center=(0,0,1), axis=(1,0,0))
        for u in np.linspace(vol.start()[0], vol.end()[0], 5):
            for v in np.linspace(vol.start()[1], vol.end()[1], 5):
                for w in np.linspace(vol.start()[2], vol.end()[2], 5):
                    x = vol(u, v, w)
                    self.assertAlmostEqual(x[1]**2+(x[2]-1)**2, u**2)
                    self.assertAlmostEqual(x[0], w*5)
        self.assertAlmostEqual(vol.volume(), 2*2*pi*5, places=3)

    def test_edge_surfaces(self):
        # test 3D surface vs 2D rational surface

        # more or less random 3D surface with p=[2,2] and n=[3,4]
        controlpoints = [[0, 0, 1], [-1, 1, 1], [0, 2, 1], [1, -1, 1], [1, 0, .5], [1, 1, 1],
                         [2, 1, 1], [2, 2, .5], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis2 = BSplineBasis(3, [0, 0, 0, .64, 2, 2, 2])
        top = Surface(basis1, basis2, controlpoints)

        # more or less random 2D rational surface with p=[1,2] and n=[3,4]
        controlpoints = [[0, 0, 1], [-1, 1, .96], [0, 2, 1], [1, -1, 1], [1, 0, .8], [1, 1, 1],
                         [2, 1, .89], [2, 2, .9], [2, 3, 1], [3, 0, 1], [4, 1, 1], [3, 2, 1]]
        basis1 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis2 = BSplineBasis(2, [0, 0, .4, .44, 1, 1])
        bottom = Surface(basis1, basis2, controlpoints, True)

        vol = VolumeFactory.edge_surfaces(bottom, top)

        # set parametric domain to [0,1]^2 for easier comparison
        top.reparam()
        bottom.reparam()

        # verify on 7x7x2 evaluation grid
        for u in np.linspace(0, 1, 7):
            for v in np.linspace(0, 1, 7):
                for w in np.linspace(0, 1, 2):  # rational basis, not linear in w-direction
                    self.assertAlmostEqual(
                        vol(u, v, w)[0], bottom(u, v)[0] *
                        (1 - w) + top(u, v)[0] * w)  # x-coordinate
                    self.assertAlmostEqual(
                        vol(u, v, w)[1], bottom(u, v)[1] *
                        (1 - w) + top(u, v)[1] * w)  # y-coordinate
                    self.assertAlmostEqual(
                        vol(u, v, w)[2], 0 * (1 - w) + top(u, v)[2] * w)  # z-coordinate

        # test 3D surface vs 2D surface
        controlpoints = [[0, 0], [-1, 1], [0, 2], [1, -1], [1, 0], [1, 1], [2, 1], [2, 2], [2, 3],
                         [3, 0], [4, 1], [3, 2]]
        basis1 = BSplineBasis(3, [0, 0, 0, 1, 1, 1])
        basis2 = BSplineBasis(2, [0, 0, .4, .44, 1, 1])
        bottom = Surface(basis1, basis2, controlpoints)  # non-rational!

        vol = VolumeFactory.edge_surfaces(bottom, top)  # also non-rational!

        # verify on 5x5x7 evaluation grid
        for u in np.linspace(0, 1, 5):
            for v in np.linspace(0, 1, 5):
                for w in np.linspace(0, 1, 7):  # include inner evaluation points
                    self.assertAlmostEqual(
                        vol(u, v, w)[0], bottom(u, v)[0] *
                        (1 - w) + top(u, v)[0] * w)  # x-coordinate
                    self.assertAlmostEqual(
                        vol(u, v, w)[1], bottom(u, v)[1] *
                        (1 - w) + top(u, v)[1] * w)  # y-coordinate
                    self.assertAlmostEqual(
                        vol(u, v, w)[2], 0 * (1 - w) + top(u, v)[2] * w)  # z-coordinate


    def test_edge_surfaces_six_sides(self):
        # create the unit cube
        vol = Volume()
        vol.raise_order(2,2,2)
        vol.refine(3)

        edges = vol.faces()

        # edge_surface should give back the same unit cube
        vol2 = VolumeFactory.edge_surfaces(edges)

        # check discretization
        self.assertEqual(vol2.order(0), 4)
        self.assertEqual(vol2.order(1), 4)
        self.assertEqual(vol2.order(2), 4)

        self.assertEqual(len(vol2.knots(0)), 5) # [0,.25,.5,.75,1]
        self.assertEqual(len(vol2.knots(1)), 5)
        self.assertEqual(len(vol2.knots(2)), 5)

        # check a 5x5x5 evaluation grid
        u = np.linspace(0,1,5)
        v = np.linspace(0,1,5)
        w = np.linspace(0,1,5)
        pt  = vol( u,v,w)
        pt2 = vol2(u,v,w)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

    def test_surface_loft(self):
        crv1 = Curve(BSplineBasis(3, range(11), 1), [[1,-1], [1,0], [1,1], [-1,1], [-1,0], [-1,-1]])
        crv2 = CurveFactory.circle(2) + (0,0,1)
        crv3 = Curve(BSplineBasis(4, range(11), 2), [[1,-1,2], [1,1,2], [-1,1,2], [-1,-1,2]])
        crv4 = CurveFactory.circle(2) + (0,0,3)
        surf = SurfaceFactory.loft(crv1, crv2, crv3, crv4)

        crv1.set_dimension(3) # for convenience when evaluating
        t = np.linspace( 0, 1, 13)

        u = np.linspace(crv1.start(0), crv1.end(0), 13)
        pt  = crv1(u)
        pt2 = surf(t,0).reshape(13,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        u = np.linspace(crv2.start(0), crv2.end(0), 13)
        pt  = crv2(u)
        pt2 = surf(t,1).reshape(13,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        u = np.linspace(crv3.start(0), crv3.end(0), 13)
        pt  = crv3(u)
        pt2 = surf(t,2).reshape(13,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        u = np.linspace(crv4.start(0), crv4.end(0), 13)
        pt  = crv4(u)
        pt2 = surf(t,3).reshape(13,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

    def test_volume_loft(self):
        crv1 = Curve(BSplineBasis(3, range(11), 1), [[1,-1], [1,0], [1,1], [-1,1], [-1,0], [-1,-1]])
        crv2 = CurveFactory.circle(2) + (0,0,1)
        crv3 = Curve(BSplineBasis(4, range(11), 2), [[1,-1,2], [1,1,2], [-1,1,2], [-1,-1,2]])
        crv4 = CurveFactory.circle(2) + (0,0,3)

        surf = []
        for c in [crv1, crv2, crv3, crv4]:
            c2 = c.clone()
            c2.project('z')
            surf.append(SurfaceFactory.edge_curves(c, c2))

        vol = VolumeFactory.loft(surf)

        surf[0].set_dimension(3) # for convenience when evaluating
        t = np.linspace( 0, 1, 9)
        s = np.linspace( 0, 1, 9)

        u = np.linspace(surf[0].start(0), surf[0].end(0), 9)
        v = np.linspace(surf[0].start(1), surf[0].end(1), 9)
        pt  = surf[0](u,v)
        pt2 = vol(s,t,0).reshape(9,9,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        u = np.linspace(surf[1].start(0), surf[1].end(0), 9)
        u = np.linspace(surf[1].start(1), surf[1].end(1), 9)
        pt  = surf[1](u,v)
        pt2 = vol(s,t,1).reshape(9,9,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        u = np.linspace(surf[2].start(0), surf[2].end(0), 9)
        v = np.linspace(surf[2].start(1), surf[2].end(1), 9)
        pt  = surf[2](u,v)
        pt2 = vol(s,t,2).reshape(9,9,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

        u = np.linspace(surf[3].start(0), surf[3].end(0), 9)
        v = np.linspace(surf[3].start(1), surf[3].end(1), 9)
        pt  = surf[3](u,v)
        pt2 = vol(s,t,3).reshape(9,9,3)
        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)

    def test_revolve(self):
        # square torus
        square = Surface() + (1,0)
        square.rotate(pi / 2, (1, 0, 0)) # in xz-plane with corners at (1,0),(2,0),(2,1),(1,1)

        vol = VolumeFactory.revolve(square)
        vol.reparam()  # set parametric space to (0,1)^3
        u = np.linspace(0, 1, 7)
        v = np.linspace(0, 1, 7)
        w = np.linspace(0, 1, 7)
        x = vol.evaluate(u,v,w)
        V,U,W = np.meshgrid(v,u,w)
        R = np.sqrt(x[:,:,:,0]**2 + x[:,:,:,1]**2)

        self.assertEqual(np.allclose(R, U+1), True)
        self.assertEqual(np.allclose(x[:,:,:,2], V), True)
        self.assertAlmostEqual(vol.volume(), 2*pi*1.5, places=3)

        # test incomplete reolve
        vol = VolumeFactory.revolve(square, theta=pi/3)
        vol.reparam()  # set parametric space to (0,1)^3
        u = np.linspace(0, 1, 7)
        v = np.linspace(0, 1, 7)
        w = np.linspace(0, 1, 7)
        x = vol.evaluate(u,v,w)
        V,U,W = np.meshgrid(v,u,w)
        R = np.sqrt(x[:,:,:,0]**2 + x[:,:,:,1]**2)

        self.assertEqual(np.allclose(R, U+1), True)
        self.assertEqual(np.allclose(x[:,:,:,2], V), True)
        self.assertAlmostEqual(vol.volume(), 2*pi*1.5/6, places=3)
        self.assertTrue(np.all(x >= 0)) # completely contained in first octant

    def test_interpolate(self):
        t = np.linspace(0, 1, 5)
        V,U,W = np.meshgrid(t,t,t)
        x = np.zeros((5,5,5,3))
        x[:,:,:,0] = U*U + V*V
        x[:,:,:,1] = U*U*W
        x[:,:,:,2] = W*W*W + U*V
        b1 = BSplineBasis(3, [0,0,0,.33,.66,1,1,1])
        b2 = BSplineBasis(4, [0,0,0,0,.5,1,1,1,1])
        b3 = BSplineBasis(4, [0,0,0,0,.2,1,1,1,1])

        vol = VolumeFactory.interpolate(x, [b1,b2,b3], [t,t,t])
        t = np.linspace(0, 1, 7)
        V,U,W = np.meshgrid(t,t,t)
        x = vol(t,t,t)

        self.assertTrue(np.allclose(x[:,:,:,0], U*U + V*V))
        self.assertTrue(np.allclose(x[:,:,:,1], U*U*W))
        self.assertTrue(np.allclose(x[:,:,:,2], W*W*W + U*V))

    def test_l2(self):
        t = np.linspace(0, 1, 80)
        V,U,W = np.meshgrid(t,t,t)
        x = np.zeros((80,80,80,3))
        x[:,:,:,0] = U*U + V*V
        x[:,:,:,1] = U*U*W
        x[:,:,:,2] = W*W*W + U*V
        b1 = BSplineBasis(3, [0,0,0,.33,.66,1,1,1])
        b2 = BSplineBasis(4, [0,0,0,0,.5,1,1,1,1])
        b3 = BSplineBasis(4, [0,0,0,0,.2,1,1,1,1])

        vol = VolumeFactory.least_square_fit(x, [b1,b2,b3], [t,t,t])
        t = np.linspace(0, 1, 7)
        V,U,W = np.meshgrid(t,t,t)
        x = vol(t,t,t)

        self.assertTrue(np.allclose(x[:,:,:,0], U*U + V*V))
        self.assertTrue(np.allclose(x[:,:,:,1], U*U*W))
        self.assertTrue(np.allclose(x[:,:,:,2], W*W*W + U*V))

if __name__ == '__main__':
    unittest.main()
