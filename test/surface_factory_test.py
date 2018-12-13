# -*- coding: utf-8 -*-

from splipy import BSplineBasis, Curve, Surface
import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
from math import pi, sqrt
import numpy as np
from numpy.linalg import norm
import unittest
try:
    import nutils
    has_nutils = True
except ImportError:
    has_nutils = False

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
        self.assertTrue(np.allclose(x[0,0], [0,0]))
        self.assertTrue(np.allclose(x[1,0], [1,0]))
        self.assertTrue(np.allclose(x[1,1], [1/sqrt(2),1/sqrt(2)]))
        self.assertTrue(np.allclose(x[1,2], [0,1]))
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

        # test xaxis
        surf = SurfaceFactory.disc(r=3, xaxis=(0,1,0))
        u = surf.end('u')
        self.assertTrue(np.allclose(surf(u,0), [0,3]))

        # test normal
        surf = SurfaceFactory.disc(r=2, normal=(1,1,0), xaxis=(0,0,1))
        u = surf.end('u')
        self.assertTrue(np.allclose(surf(u,0), [0,0,2]))

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

        # test xaxis
        surf = SurfaceFactory.cylinder(r=sqrt(2), xaxis=(1,1,0))
        self.assertTrue(np.allclose(surf(0,0), [1,1,0]))

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
    def test_edge_curves_poisson(self):
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

    @unittest.skipIf(not has_nutils, "EdgeCurves with elasticity solver requires nutils")
    def test_edge_curves_elasticity(self):
        # create an arrow-like 2D geometry with the pointy end at (-1,1) towards up and left
        # rebuild to avoid rational representations
        c1 = CurveFactory.circle_segment(pi / 2).rebuild(3,11)
        c2 = Curve(BSplineBasis(2, [0, 0, 1, 2, 2]), [[0, 1], [-1, 1], [-1, 0]])
        c3 = CurveFactory.circle_segment(pi / 2).rebuild(3,11)
        c3.rotate(pi)
        c4 = Curve(BSplineBasis(2), [[0, -1], [1, 0]]).rebuild(3,10)
        c4 = c4.rebuild(4,11)

        surf = SurfaceFactory.edge_curves([c1, c2, c3, c4], type='elasticity')

        # check right dimensions of output
        self.assertEqual(surf.shape[0], 11) # 11 controlpoints in the circle segment
        self.assertEqual(surf.shape[1], 13) # 11 controlpoints in c4, +2 for C0-knot in c1
        self.assertEqual(surf.order(0), 3)
        self.assertEqual(surf.order(1), 4)

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

    @unittest.skipIf(not has_nutils, "EdgeCurves with finitestrain solver requires nutils")
    def test_edge_curves_finitestrain(self):
        # create an arrow-like 2D geometry with the pointy end at (-1,1) towards up and left
        # rebuild to avoid rational representations
        c1 = CurveFactory.circle_segment(pi / 2).rebuild(3,11)
        c2 = Curve(BSplineBasis(2, [0, 0, 1, 2, 2]), [[0, 1], [-1, 1], [-1, 0]])
        c3 = CurveFactory.circle_segment(pi / 2).rebuild(3,11)
        c3.rotate(pi)
        c4 = Curve(BSplineBasis(2), [[0, -1], [1, 0]]).rebuild(3,10)
        c4 = c4.rebuild(4,11)

        surf = SurfaceFactory.edge_curves([c1, c2, c3, c4], type='finitestrain')

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

    @unittest.skipIf(not has_nutils, "EdgeCurves with finitestrain solver requires nutils")
    def test_edge_curves_finitestrain_lshape(self):
        # Create an L-shape geometry with an interior 270-degree angle at the origin (u=.5, v=1)
        c1 = CurveFactory.polygon([[-1, 1], [-1,-1], [1,-1]])
        c2 = CurveFactory.polygon([[ 1,-1], [ 1, 0]])
        c3 = CurveFactory.polygon([[ 1, 0], [ 0, 0], [0, 1]])
        c4 = CurveFactory.polygon([[ 0, 1], [-1, 1]])
        c1.refine(2).raise_order(1)
        c2.refine(2).raise_order(1)
        surf = SurfaceFactory.edge_curves([c1, c2, c3, c4], type='finitestrain')
        # the quickest way to check for self-intersecting geometry here is that
        # the normal is pointing the wrong way: down z-axis instead of up
        surf.reparam().set_dimension(3)
        self.assertTrue(surf.normal(0.5, 0.98)[2] > 0.0)
        # also check that no controlpoints leak away into the first quadrant
        self.assertFalse(np.any(np.logical_and(surf[:,:,0] > 0, surf[:,:,1] > 0)))

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

        # test 3D geometry
        c = Curve()
        c.set_dimension(3)
        s = SurfaceFactory.thicken(c, 1) # cylinder along x-axis with h=1, r=1
        for u in np.linspace(s.start(0), s.end(0), 5):
            for v in np.linspace(s.start(1), s.end(1), 5):
                x = s(u, v)
                self.assertAlmostEqual(x[1]**2+x[2]**2, 1.0**2) # distance to x-axis
                self.assertAlmostEqual(x[0], u)                 # x coordinate should be linear

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

if __name__ == '__main__':
    unittest.main()
