from BSplineBasis import *
import CurveFactory
import SurfaceFactory
import VolumeFactory
import numpy as np
import unittest

class TestFactory(unittest.TestCase):

    def test_line(self):

        ### 2D line
        c = CurveFactory.line([1,1], [2,0])
        self.assertEqual(c.get_order(), 2) # linear discretization
        self.assertEqual(len(c), 2)        # two control points
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c[0][0], 1)
        self.assertEqual(c[0][1], 1)
        self.assertEqual(c[-1][0], 2)
        self.assertEqual(c[-1][1], 0)

        ### 3D line
        c = CurveFactory.line([1,2,3], [8,7,6])
        self.assertEqual(c.get_order(), 2) # linear discretization
        self.assertEqual(len(c), 2)        # two control points
        self.assertEqual(c.dimension, 3)

    def test_n_gon(self):
        pi = np.pi
        ### test default 5 side n-gon
        c = CurveFactory.n_gon()
        self.assertEqual(len(c), 5)
        self.assertEqual(len(c.get_knots()), 6)
        self.assertEqual(c.get_order(), 2)
        # evaluate at second corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0)[0], np.cos(2*pi/5))
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0)[1], np.sin(2*pi/5))
        # evaluate at fourh corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0*4)[0], np.cos(2*pi/5*4))
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0*4)[1], np.sin(2*pi/5*4))

        ### test a radius 3 septagon
        c = CurveFactory.n_gon(n=7, r=3)
        self.assertEqual(len(c), 7)
        # evaluate at third corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.stop()/7.0)[0], 3*np.cos(2*pi/7))
        self.assertAlmostEqual(c.evaluate(c.stop()/7.0)[1], 3*np.sin(2*pi/7))

    def test_circle(self):

        ### unit circle of radius 1
        c = CurveFactory.circle()
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt), 1.0) # check radius=1
        # test evaluation at key points
        pt = c.evaluate(0)         # evaluate at 0
        self.assertAlmostEqual(pt[0], 1.0)
        self.assertAlmostEqual(pt[1], 0.0)
        pt = c.evaluate(c.end()/4) # evaluate at pi/2
        self.assertAlmostEqual(pt[0], 0.0)
        self.assertAlmostEqual(pt[1], 1.0)
        pt = c.evaluate(c.end()/2) # evaluate at pi
        self.assertAlmostEqual(pt[0], -1.0)
        self.assertAlmostEqual(pt[1],  0.0)
        pt = c.evaluate(c.end()*2) # evaluate outside domain (test periodic)
        self.assertAlmostEqual(pt[0], 1.0)
        self.assertAlmostEqual(pt[1], 0.0)
        
        ### circle of radius different from 1
        c = CurveFactory.circle(3)
        # test evaluation at 25 points for radius=3, outside domain
        t = np.linspace(c.start()-3, c.end()+2, 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius=3

        ### test errors and exceptions
        with self.assertRaises(ValueError):
            c = CurveFactory.circle(-2.5) # negative radius

    def test_circle_segment(self):
        pi = np.pi

        ### basic circle segment
        c = CurveFactory.circle_segment(np.pi*0.9)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 1.0) # check radius=1

        ### radius 7 circle segment
        c = CurveFactory.circle_segment(pi*1.87, 7)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=7
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 7.0) # check radius=7

        ### boundary case with one knot span circle segment
        c = CurveFactory.circle_segment(2*pi/3)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        self.assertEqual(len(c.get_knots()), 2)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 1.0) # check radius=1

        ### boundary case with full circle
        c = CurveFactory.circle_segment(2*pi)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        self.assertEqual(len(c.get_knots()), 4)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 1.0) # check radius=1

        ### test errors and exceptions
        with self.assertRaises(ValueError):
            c = CurveFactory.circle_segment(3*pi)  # outside domain
        with self.assertRaises(ValueError):
            c = CurveFactory.circle_segment(-3*pi) # outside domain
        with self.assertRaises(ValueError):
            c = CurveFactory.circle_segment(pi, -2) # negative radius

    def test_square(self):
        surf = SurfaceFactory.square((4,5))
        self.assertEqual(surf.dimension, 2)
        self.assertEqual(surf.rational, False)
        self.assertEqual(surf.get_order(), (2,2))

    def test_curve_interpolation(self):
        basis = BSplineBasis(4, [0,0,0,0,.3,.9,1,1,1,1])
        t = np.array(basis.greville())
        # create the mapping (x,y,z)=(t^2, 1-t, t^3+2*t)
        x_pts = np.zeros((len(t),3))
        x_pts[:,0] = t*t
        x_pts[:,1] = 1-t
        x_pts[:,2] = t*t*t + 2*t
        crv = CurveFactory.interpolate(x_pts, basis)
        self.assertEqual(crv.get_order(), 4)
        self.assertAlmostEqual(crv(.4)[0], .4**2)      # x=t^2
        self.assertAlmostEqual(crv(.4)[1], 1-.4)       # y=1-t
        self.assertAlmostEqual(crv(.4)[2], .4**3+2*.4) # z=t^3+2t
        self.assertAlmostEqual(crv(.5)[0], .5**2)      # x=t^2
        self.assertAlmostEqual(crv(.5)[1], 1-.5)       # y=1-t
        self.assertAlmostEqual(crv(.5)[2], .5**3+2*.5) # z=t^3+2t

    def test_disc(self):
        pi = np.pi
        ### radial disc
        surf = SurfaceFactory.disc()
        x = surf.evaluate([0,1], [0,pi/4,pi/2,pi])
        self.assertAlmostEqual(x[0][0][0], 0)
        self.assertAlmostEqual(x[0][0][1], 0)
        self.assertAlmostEqual(x[1][0][0], 1)
        self.assertAlmostEqual(x[1][0][1], 0)
        self.assertAlmostEqual(x[1][1][0], 1/np.sqrt(2))
        self.assertAlmostEqual(x[1][1][1], 1/np.sqrt(2))
        self.assertAlmostEqual(x[1][2][0], 0)
        self.assertAlmostEqual(x[1][2][1], 1)

        ### radial disc of size different from 1
        surf = SurfaceFactory.disc(4)
        # test evaluation at 25 points for radius=4
        v = np.linspace(0,2*pi, 25)
        u = 1
        x = surf.evaluate(u,v)
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 4.0) # check radius

        ### square disc
        surf = SurfaceFactory.disc(3, 'square')
        # evaluate on all 4 edges, 5 pts on each edge
        u = np.linspace(0,1, 5)
        v = np.linspace(0,1, 5)
        x = surf.evaluate(u,v)
        for pt in np.array(x[0,:,:]): # umin edge
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius
        for pt in np.array(x[-1,:,:]): # umax edge
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius
        for pt in np.array(x[:,0,:]): # vmin edge
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius
        for pt in np.array(x[:,-1,:]): # vmax edge
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius

    def test_revolve(self):
        pi = np.pi
        ### square torus
        square = CurveFactory.n_gon(4)
        square.rotate(pi/2, (1,0,0))
        square.translate((2,0,0)) # in xz-plane with corners at (3,0),(2,1),(1,0),(2,-1)
        surf = SurfaceFactory.revolve(square)
        surf.reparametrize() # set parametric space to (0,1)^2
        v = np.linspace(0,1,13)
        x = surf.evaluate(0,v) # outer ring evaluation u=0
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius=3
        x = surf.evaluate(.25,v) # top ring evaluation u=.25
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(pt[0]*pt[0]+pt[1]*pt[1], 2*2) # check radius=2
            self.assertAlmostEqual(pt[2], 1)                     # check height=1
        x = surf.evaluate(.375,v) # mid inner ring evaluation u=.375
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(pt[0]*pt[0]+pt[1]*pt[1], 1.5*1.5) # check radius=1.5
            self.assertAlmostEqual(pt[2], .5)                        # check height=0.5

    def test_surface_torus(self):
        pi = np.pi
        ### default torus
        surf = SurfaceFactory.torus(1,3)
        start = surf.start()
        end   = surf.end()
        # check a 13 evaluation points of v-evaluation (around the z-axis)
        v = np.linspace(start[1],end[1],13)
        # check minor-circle u=0 (outmost ring)
        x = surf.evaluate(0,v)
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(pt[0]*pt[0]+pt[1]*pt[1], 4*4) # check radius=4
            self.assertAlmostEqual(pt[2], 0)                     # check height=0
        # check minor-circle u=pi (innermost ring)
        x = surf.evaluate(pi,v)
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(pt[0]*pt[0]+pt[1]*pt[1], 2*2) # check radius=2
            self.assertAlmostEqual(pt[2], 0)                     # check height=0
        # check minor-circle u=pi/2 (topmost ring)
        x = surf.evaluate(pi/2,v)
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(pt[0]*pt[0]+pt[1]*pt[1], 3*3) # check radius=3
            self.assertAlmostEqual(pt[2], 1)                     # check height=1
        # check minor-circle u=3*pi/2 (mid-evaluation)
        x = surf.evaluate(3*pi/4,v)
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(pt[0]*pt[0]+pt[1]*pt[1], (3-1.0/np.sqrt(2))**2) # check radius=3-1/sqrt(2)
            self.assertAlmostEqual(pt[2], 1.0/np.sqrt(2))             # check height=1/sqrt(2)

    def test_sphere(self):
        pi = np.pi
        ### unit ball
        surf = SurfaceFactory.sphere()
        # test 7x7 grid for radius = 1
        for u in np.linspace(surf.start()[0], surf.end()[0], 7):
            for v in np.linspace(surf.start()[1], surf.end()[1], 7):
                self.assertAlmostEqual(np.linalg.norm(surf(u,v),2), 1.0)

    def test_cylinder_surface(self):
        pi = np.pi
        ### unit cylinder
        surf = SurfaceFactory.cylinder()
        # test 7x7 grid for xy-radius = 1 and v=z
        for u in np.linspace(surf.start()[0], surf.end()[0], 7):
            for v in np.linspace(surf.start()[1], surf.end()[1], 7):
                x = surf(u,v)
                self.assertAlmostEqual(np.linalg.norm(x[:2],2), 1.0) # (x,y) coordinates to z-axis
                self.assertAlmostEqual(x[2],  v)                     # z coordinate should be linear

    def test_cylinder_volume(self):
        pi = np.pi
        ### unit cylinder
        vol = VolumeFactory.cylinder()
        # test 5x5x5 grid for xy-radius = w and v=z
        for u in np.linspace(vol.start()[0], vol.end()[0], 5):
            for v in np.linspace(vol.start()[1], vol.end()[1], 5):
                for w in np.linspace(vol.start()[2], vol.end()[2], 5):
                    x = vol(u,v,w)
                    self.assertAlmostEqual(np.linalg.norm(x[:2],2), w) # (x,y) coordinates to z-axis
                    self.assertAlmostEqual(x[2],  v)                   # z coordinate should be linear

if __name__ == '__main__':
    unittest.main()
