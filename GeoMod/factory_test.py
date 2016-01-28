# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis, Curve, Surface
import GeoMod.CurveFactory as CurveFactory
import GeoMod.SurfaceFactory as SurfaceFactory
import GeoMod.VolumeFactory as VolumeFactory
from math import pi, sqrt, cos, sin
import numpy as np
import unittest


class TestFactory(unittest.TestCase):
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
        self.assertEqual(len(c.knots(0)), 8)
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
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(0), c.end(0), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 1.0)  # check radius=1

        # boundary case with full circle
        c = CurveFactory.circle_segment(2 * pi)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        self.assertEqual(len(c.knots(0)), 4)
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

        # radial disc of size different from 1
        surf = SurfaceFactory.disc(4)
        # test evaluation at 25 points for radius=4
        v = np.linspace(0, 2 * pi, 25)
        u = 1
        x = surf.evaluate(u, v)
        for pt in np.array(x[0, :, :]):
            self.assertAlmostEqual(np.linalg.norm(pt, 2), 4.0)  # check radius

        # square disc
        surf = SurfaceFactory.disc(3, 'square')
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

    def test_revolve(self):
        # square torus
        square = CurveFactory.n_gon(4)
        square.rotate(pi / 2, (1, 0, 0))
        square.translate((2, 0, 0))  # in xz-plane with corners at (3,0),(2,1),(1,0),(2,-1)
        surf = SurfaceFactory.revolve(square)
        surf.reparametrize()  # set parametric space to (0,1)^2
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

    def test_cylinder_surface(self):
        # unit cylinder
        surf = SurfaceFactory.cylinder()
        # test 7x7 grid for xy-radius = 1 and v=z
        for u in np.linspace(surf.start()[0], surf.end()[0], 7):
            for v in np.linspace(surf.start()[1], surf.end()[1], 7):
                x = surf(u, v)
                self.assertAlmostEqual(np.linalg.norm(x[:2], 2), 1.0)  # (x,y) coordinates to z-axis
                self.assertAlmostEqual(x[2], v)  # z coordinate should be linear

    def test_edge_curves(self):
        # create an arrow-like 2D geometry with the pointy end at (-1,1) towards up and left
        # mixes rational and non-rational curves with different parametrization spaces
        c1 = CurveFactory.circle_segment(pi / 2)
        c2 = Curve(BSplineBasis(2, [0, 0, 1, 2, 2]), [[0, 1], [-1, 1], [-1, 0]])
        c3 = CurveFactory.circle_segment(pi / 2)
        c3.rotate(pi)
        c4 = Curve(BSplineBasis(2), [[0, -1], [1, 0]])

        surf = SurfaceFactory.edge_curves([c1, c2, c3, c4])

        # srf spits out parametric space (0,1)^2, so we sync these up to input curves
        c3.flip_parametrization()
        c4.flip_parametrization()
        c1.reparametrize()
        c2.reparametrize()
        c3.reparametrize()
        c4.reparametrize()

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

    def test_cylinder_volume(self):
        # unit cylinder
        vol = VolumeFactory.cylinder()
        # test 5x5x5 grid for xy-radius = w and v=z
        for u in np.linspace(vol.start()[0], vol.end()[0], 5):
            for v in np.linspace(vol.start()[1], vol.end()[1], 5):
                for w in np.linspace(vol.start()[2], vol.end()[2], 5):
                    x = vol(u, v, w)
                    self.assertAlmostEqual(
                        np.linalg.norm(x[:2], 2), w)  # (x,y) coordinates to z-axis
                    self.assertAlmostEqual(x[2], v)  # z coordinate should be linear

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

        vol = VolumeFactory.edge_surfaces([bottom, top])

        # set parametric domain to [0,1]^2 for easier comparison
        top.reparametrize()
        bottom.reparametrize()

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

        vol = VolumeFactory.edge_surfaces([bottom, top])  # also non-rational!

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


if __name__ == '__main__':
    unittest.main()
