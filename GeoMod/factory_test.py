from Factory import *
import unittest

class TestFactory(unittest.TestCase):

    def test_line(self):

        ### 2D line
        c = line([1,1], [2,0])
        self.assertEqual(c.get_order(), 2) # linear discretization
        self.assertEqual(len(c), 2)        # two control points
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c[0][0], 1)
        self.assertEqual(c[0][1], 1)
        self.assertEqual(c[-1][0], 2)
        self.assertEqual(c[-1][1], 0)

        ### 3D line
        c = line([1,2,3], [8,7,6])
        self.assertEqual(c.get_order(), 2) # linear discretization
        self.assertEqual(len(c), 2)        # two control points
        self.assertEqual(c.dimension, 3)

    def test_n_gon(self):
        ### test default 5 side n-gon
        c = n_gon()
        self.assertEqual(len(c), 5)
        self.assertEqual(len(c.get_knots()), 6)
        self.assertEqual(c.get_order(), 2)
        # evaluate at second corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0)[0], cos(2*pi/5))
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0)[1], sin(2*pi/5))
        # evaluate at fourh corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0*4)[0], cos(2*pi/5*4))
        self.assertAlmostEqual(c.evaluate(c.stop()/5.0*4)[1], sin(2*pi/5*4))

        ### test a radius 3 septagon
        c = n_gon(n=7, r=3)
        self.assertEqual(len(c), 7)
        # evaluate at third corner (clockwise from (1,0) )
        self.assertAlmostEqual(c.evaluate(c.stop()/7.0)[0], 3*cos(2*pi/7))
        self.assertAlmostEqual(c.evaluate(c.stop()/7.0)[1], 3*sin(2*pi/7))

    def test_circle(self):

        ### unit circle of radius 1
        c = circle()
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
        c = circle(3)
        # test evaluation at 25 points for radius=3, outside domain
        t = np.linspace(c.start()-3, c.end()+2, 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 3.0) # check radius=3

        ### test errors and exceptions
        with self.assertRaises(ValueError):
            c = circle(-2.5) # negative radius

    def test_circle_segment(self):

        ### basic circle segment
        c = circle_segment(pi*0.9)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 1.0) # check radius=1

        ### radius 7 circle segment
        c = circle_segment(pi*1.87, 7)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        # test evaluation at 25 points for radius=7
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 7.0) # check radius=7

        ### boundary case with one knot span circle segment
        c = circle_segment(2*pi/3)
        self.assertEqual(c.dimension, 2)
        self.assertEqual(c.rational, True)
        self.assertEqual(len(c.get_knots()), 2)
        # test evaluation at 25 points for radius=1
        t = np.linspace(c.start(), c.end(), 25)
        x = c.evaluate(t)
        for pt in np.array(x):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 1.0) # check radius=1

        ### boundary case with full circle
        c = circle_segment(2*pi)
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
            c = circle_segment(3*pi)  # outside domain
        with self.assertRaises(ValueError):
            c = circle_segment(-3*pi) # outside domain
        with self.assertRaises(ValueError):
            c = circle_segment(pi, -2) # negative radius

    def test_square(self):
        surf = square(4,5)
        self.assertEqual(surf.dimension, 2)
        self.assertEqual(surf.rational, False)
        self.assertEqual(surf.get_order(), (2,2))

    def test_disc(self):
        ### radial disc
        surf = disc()
        x = surf.evaluate([0,1], [0,pi/4,pi/2,pi])
        self.assertAlmostEqual(x[0][0][0], 0)
        self.assertAlmostEqual(x[0][0][1], 0)
        self.assertAlmostEqual(x[1][0][0], 1)
        self.assertAlmostEqual(x[1][0][1], 0)
        self.assertAlmostEqual(x[1][1][0], 1/sqrt(2))
        self.assertAlmostEqual(x[1][1][1], 1/sqrt(2))
        self.assertAlmostEqual(x[1][2][0], 0)
        self.assertAlmostEqual(x[1][2][1], 1)

        ### radial disc of size different from 1
        surf = disc(4)
        # test evaluation at 25 points for radius=4
        v = np.linspace(0,2*pi, 25)
        u = 1
        x = surf.evaluate(u,v)
        for pt in np.array(x[0,:,:]):
            self.assertAlmostEqual(np.linalg.norm(pt,2), 4.0) # check radius

        ### square disc
        surf = disc(3, 'square')
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



if __name__ == '__main__':
    unittest.main()
