# -*- coding: utf-8 -*-

from splipy import BSplineBasis, Curve, Surface, Volume
import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
import splipy.volume_factory as VolumeFactory
from math import pi, sqrt
import numpy as np
from numpy.linalg import norm
import unittest
try:
    import nutils
    has_nutils = True
except ImportError:
    has_nutils = False

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

        # test xaxis
        vol = VolumeFactory.cylinder(r=sqrt(2), xaxis=(1,1,0))
        u = vol.end('u')
        self.assertTrue(np.allclose(vol(u,0,0), [1,1,0]))

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

        # test axis revolve
        vol = VolumeFactory.revolve(Surface()+(1,1), theta=pi/3, axis=(1,0,0))
        vol.reparam()  # set parametric space to (0,1)^3
        u = np.linspace(0, 1, 7)
        v = np.linspace(0, 1, 7)
        w = np.linspace(0, 1, 7)
        x = vol.evaluate(u,v,w)
        V,U,W = np.meshgrid(v,u,w)
        R = np.sqrt(x[:,:,:,1]**2 + x[:,:,:,2]**2)

        self.assertEqual(np.allclose(R, V+1), True)
        self.assertEqual(np.allclose(x[:,:,:,0], U+1), True)
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

if __name__ == '__main__':
    unittest.main()
