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


if __name__ == '__main__':
    unittest.main()
