# -*- coding: utf-8 -*-

from splipy import BSplineBasis
import numpy as np
import unittest


class TestBasis(unittest.TestCase):
    def test_get_continuity(self):
        b = BSplineBasis(4, [0, 0, 0, 0, .3, 1, 1, 1.134, 1.134, 1.134, 2, 2, 2, 2])
        self.assertEqual(b.continuity(.3), 2)
        self.assertEqual(b.continuity(1), 1)
        self.assertEqual(b.continuity(1.134), 0)
        self.assertEqual(b.continuity(0), -1)
        self.assertEqual(b.continuity(2), -1)
        self.assertEqual(b.continuity(.4), np.inf)

    def test_errors(self):
        with self.assertRaises(ValueError):
            BSplineBasis(4, [1, 2, 3])
        with self.assertRaises(ValueError):
            BSplineBasis(4, [1, 2, 4, 8], periodic=1)
        with self.assertRaises(ValueError):
            BSplineBasis(4, [1, 2, 3, 7, 8], periodic=2)
        with self.assertRaises(ValueError):
            BSplineBasis(3, [0,0,0,3,2,2,2])
        with self.assertRaises(ValueError):
            BSplineBasis(3, [0,1,2,3,5,6], periodic=1)
        with self.assertRaises(ValueError):
            BSplineBasis(-1, [0,0,1,1])

    def test_num_functions(self):
        b = BSplineBasis(4, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3])
        self.assertEqual(b.num_functions(), 6)
        b = BSplineBasis(4, [-2, -1, 0, 0, 1, 2, 3, 3, 4, 5], periodic=1)
        self.assertEqual(b.num_functions(), 4)

    def test_start_end(self):
        b = BSplineBasis(4, [-2, -1, 0, 0, 1, 2, 3, 3, 4, 5], periodic=1)
        self.assertAlmostEqual(b.start(), 0)
        self.assertAlmostEqual(b.end(),   3)
        b = BSplineBasis(3, [0, 1, 1, 2, 3, 4, 4, 6])
        self.assertAlmostEqual(b.start(), 1)
        self.assertAlmostEqual(b.end(),   4)

    def test_reverse(self):
        b = BSplineBasis(3, [0,0,0,1,3,3,6,10,15,21,21,21])
        b.reverse()
        expect = [0,0,0,6,11,15,18,18,20,21,21,21]
        self.assertAlmostEqual(np.linalg.norm(b.knots - expect), 0)

    def test_reparam(self):
        b = BSplineBasis(5, [-1, 0, 1, 1, 1, 3, 9, 10, 11, 11, 11, 13, 19], periodic=1)
        b.reparam(11, 21)
        expect = [9, 10, 11, 11, 11, 13, 19, 20, 21, 21, 21, 23, 29]
        self.assertAlmostEqual(np.linalg.norm(b.knots - expect), 0)
        b.reparam()
        expect = [-.2, -.1, 0, 0, 0, .2, .8, .9, 1.0, 1.0, 1.0, 1.2, 1.8]
        self.assertAlmostEqual(np.linalg.norm(b.knots - expect), 0)

    def test_greville(self):
        b = BSplineBasis(4, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3])
        self.assertAlmostEqual(b.greville(0), 0.0)
        self.assertAlmostEqual(b.greville(1), 1.0 / 3.0)
        self.assertAlmostEqual(b.greville(2), 1.0)
        self.assertAlmostEqual(b.greville(), [0.0, 1.0/3.0, 1.0, 2.0, 8.0/3.0, 3.0])

    def test_raise_order(self):
        # test normal knot vector
        b  = BSplineBasis(4, [0, 0, 0, 0, 1, 2, 3, 3, 3, 3])
        b2 = b.raise_order(2)
        expect =  [0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3]
        self.assertAlmostEqual(np.linalg.norm(b2.knots - expect), 0)
        self.assertEqual(b2.order, 6)

        # test periodic knot vector
        b = BSplineBasis(5, [-1, 0, 1, 1, 1, 3, 9, 10, 11, 11, 11, 13, 19], periodic=1)
        b2 = b.raise_order(1)
        expect =  [0, 0, 1, 1, 1, 1, 3, 3, 9, 9, 10, 10, 11, 11, 11, 11, 13, 13]
        self.assertAlmostEqual(np.linalg.norm(b2.knots - expect), 0)
        self.assertEqual(b2.order,    6)
        self.assertEqual(b2.periodic, 1)

        with self.assertRaises(ValueError):
            BSplineBasis().raise_order(-1)

    def test_roll(self):
        b = BSplineBasis(3, [-1, 0, 0, 2, 3, 4, 4, 6], periodic=0)
        b.roll(3)
        expect = [0, 2, 3, 4, 4, 6, 7, 8]
        self.assertAlmostEqual(np.linalg.norm(b.knots - expect), 0)
        b.roll(2)
        expect = [2, 3, 4, 4, 6, 7, 8, 8]
        self.assertAlmostEqual(np.linalg.norm(b.knots - expect), 0)

        with self.assertRaises(IndexError):
            b.roll(19)

    def test_getitem(self):
        b = BSplineBasis(3, [0,0,0,1,2,2,2])
        self.assertEqual(b[0], 0.0)
        self.assertEqual(b[1], 0.0)
        self.assertEqual(b[3], 1.0)

    def test_repr(self):
        self.assertEqual(repr(BSplineBasis()), 'p=2, [ 0.  0.  1.  1.]')
        self.assertEqual(repr(BSplineBasis(periodic=0)), 'p=2, [-1.  0.  1.  2.], C0-periodic')

    def test_roll(self):
        b = BSplineBasis(4, [-2, -1, -1, 0, 2, 4, 6.5, 7, 8, 8, 9, 11, 13, 15.5], periodic=2)
        b.roll(3)
        self.assertEqual(len(b.knots), 14)
        self.assertAlmostEqual(b.knots[0], 0)
        self.assertAlmostEqual(b.knots[1], 2)
        self.assertAlmostEqual(b.knots[2], 4)
        self.assertAlmostEqual(b.knots[3], 6.5)
        self.assertAlmostEqual(b.knots[4], 7)
        self.assertAlmostEqual(b.knots[5], 8)
        self.assertAlmostEqual(b.knots[6], 8)
        self.assertAlmostEqual(b.knots[7], 9)
        self.assertAlmostEqual(b.knots[8], 11)
        self.assertAlmostEqual(b.knots[9], 13)
        self.assertAlmostEqual(b.knots[10], 15.5)
        self.assertAlmostEqual(b.knots[11], 16)
        self.assertAlmostEqual(b.knots[12], 17)
        self.assertAlmostEqual(b.knots[13], 17)

    def test_integrate(self):
        # create the linear functions x(t) = [1-t, t] on t=[0,1]
        b = BSplineBasis()
        self.assertAlmostEqual(b.integrate(0,1)[0], 0.5)
        self.assertAlmostEqual(b.integrate(0,1)[1], 0.5)
        self.assertAlmostEqual(b.integrate(.25,.5)[0], 5.0/32)
        self.assertAlmostEqual(b.integrate(.25,.5)[1], 3.0/32)

        # create the quadratic functions x(t) = [(1-t)^2, 2t(1-t), t^2] on t=[0,1]
        b = BSplineBasis(3)
        self.assertAlmostEqual(b.integrate(0,1)[0], 1.0/3)
        self.assertAlmostEqual(b.integrate(0,1)[1], 1.0/3)
        self.assertAlmostEqual(b.integrate(0,1)[2], 1.0/3)
        self.assertAlmostEqual(b.integrate(.25,.5)[0], 19.0/192)
        self.assertAlmostEqual(b.integrate(.25,.5)[1], 11.0/96)
        self.assertAlmostEqual(b.integrate(.25,.5)[2], 7.0/192)

        # create periodic quadratic functions on [0,3]. This is 3 functions, which are all
        # translated versions of the one below:
        #        | 1/2 t^2          t=[0,1]
        # N[3] = { -t^2 + 3t - 3/2  t=[1,2]
        #        | 1/2 (3-t)^2      t=[2,3]
        b = BSplineBasis(3, [-2,-1,0,1,2,3,4,5], periodic=1)
        self.assertEqual(len(b.integrate(0,3)), 3) # returns 3 functions
        self.assertAlmostEqual(b.integrate(0,3)[0], 1)
        self.assertAlmostEqual(b.integrate(0,3)[1], 1)
        self.assertAlmostEqual(b.integrate(0,3)[2], 1)
        self.assertAlmostEqual(b.integrate(0,1)[0], 1.0/6)
        self.assertAlmostEqual(b.integrate(0,1)[1], 2.0/3)
        self.assertAlmostEqual(b.integrate(0,1)[2], 1.0/6)
        self.assertAlmostEqual(b.integrate(0,2)[0], 2.0/6)
        self.assertAlmostEqual(b.integrate(0,2)[1], 5.0/6)
        self.assertAlmostEqual(b.integrate(0,2)[2], 5.0/6)
    
    def test_matches(self):
        b1 = BSplineBasis(3, [0,0,0,1,2,3,4,4,4])
        b2 = BSplineBasis(3, [1,1,1,2,3,4,5,5,5])
        b3 = BSplineBasis(4, [1,1,1,2,3,4,5,5,5])
        b4 = BSplineBasis(4, [2,2,2,4,6,8,10,10,10])
        b5 = BSplineBasis(3, [-1,0,0,1,2,3,4,4,5], periodic=0)
        b6 = BSplineBasis(3, [-1,0,0,1,2,3,4,4,5])
        b7 = BSplineBasis(3, [0,0,0,2,3,3,3])
        b8 = BSplineBasis(3, [5,5,5,7,11,11,11])

        self.assertEqual(b1.matches(b2), True)
        self.assertEqual(b1.matches(b2, reverse=True), True)
        self.assertEqual(b2.matches(b1), True)
        self.assertEqual(b1.matches(b3), False)
        self.assertEqual(b2.matches(b3), False)
        self.assertEqual(b1.matches(b4), False)
        self.assertEqual(b4.matches(b4), True)
        self.assertEqual(b5.matches(b6), False)
        self.assertEqual(b7.matches(b8, reverse=True), True)
        self.assertEqual(b8.matches(b7, reverse=True), True)

    def test_make_periodic(self):
        b = BSplineBasis(4, [1,1,1,1,2,3,4,5,5,5,5])

        c = b.make_periodic(0)
        self.assertEqual(c.periodic, 0)
        self.assertEqual(c.order, 4)
        self.assertAlmostEqual(c.knots[0], 0)
        self.assertAlmostEqual(c.knots[1], 1)
        self.assertAlmostEqual(c.knots[2], 1)
        self.assertAlmostEqual(c.knots[3], 1)
        self.assertAlmostEqual(c.knots[4], 2)
        self.assertAlmostEqual(c.knots[5], 3)
        self.assertAlmostEqual(c.knots[6], 4)
        self.assertAlmostEqual(c.knots[7], 5)
        self.assertAlmostEqual(c.knots[8], 5)
        self.assertAlmostEqual(c.knots[9], 5)
        self.assertAlmostEqual(c.knots[10], 6)

        c = b.make_periodic(1)
        self.assertEqual(c.periodic, 1)
        self.assertEqual(c.order, 4)
        self.assertAlmostEqual(c.knots[0], -1)
        self.assertAlmostEqual(c.knots[1], 0)
        self.assertAlmostEqual(c.knots[2], 1)
        self.assertAlmostEqual(c.knots[3], 1)
        self.assertAlmostEqual(c.knots[4], 2)
        self.assertAlmostEqual(c.knots[5], 3)
        self.assertAlmostEqual(c.knots[6], 4)
        self.assertAlmostEqual(c.knots[7], 5)
        self.assertAlmostEqual(c.knots[8], 5)
        self.assertAlmostEqual(c.knots[9], 6)
        self.assertAlmostEqual(c.knots[10], 7)

        c = b.make_periodic(2)
        self.assertEqual(c.periodic, 2)
        self.assertEqual(c.order, 4)
        self.assertAlmostEqual(c.knots[0], -2)
        self.assertAlmostEqual(c.knots[1], -1)
        self.assertAlmostEqual(c.knots[2], 0)
        self.assertAlmostEqual(c.knots[3], 1)
        self.assertAlmostEqual(c.knots[4], 2)
        self.assertAlmostEqual(c.knots[5], 3)
        self.assertAlmostEqual(c.knots[6], 4)
        self.assertAlmostEqual(c.knots[7], 5)
        self.assertAlmostEqual(c.knots[8], 6)
        self.assertAlmostEqual(c.knots[9], 7)
        self.assertAlmostEqual(c.knots[10], 8)

if __name__ == '__main__':
    unittest.main()
