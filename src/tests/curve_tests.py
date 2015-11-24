from GoTools import Curve
import unittest

class TestCurve(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        crv = Curve(2, [0, 0, 1, 1], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        val = crv.Evaluate(0.5)
        self.assertEqual(val[0], 0.5)

        # test 2D constructor
        crv2 = Curve(2, [0, 0, 1, 1], [[0.0, 0.0], [1.0, 0.0]])
        val  = crv.Evaluate(0.5)
        self.assertEqual(val[0], 0.5)

    def test_evaluate(self):
        controlpoints = [[0,0,0],  [1,1,0],  [2,0,0]]
        crv = Curve(3, [0, 0, 0, 1, 1, 1], controlpoints)

        # startpoint evaluation
        val = crv.Evaluate(0.0)     
        self.assertEqual(val[0], 0.0)
        self.assertEqual(val[1], 0.0)

        # startpoint with derivatives
        val = crv.Evaluate(0.0, 1)  
        self.assertEqual(val[0][0], 0.0)
        self.assertEqual(val[0][1], 0.0)
        self.assertEqual(val[1][0], 2.0)
        self.assertEqual(val[1][1], 2.0)

        # midpoint with derivatives
        val = crv.Evaluate(0.5, 1)  
        self.assertEqual(val[0][0], 1.0)
        self.assertEqual(val[0][1], 0.5)
        self.assertEqual(val[1][0], 2.0)
        self.assertEqual(val[1][1], 0.0)

        # end with derivatives (sensitive to left-evaluation)
        val = crv.Evaluate(1.0, 1)  
        self.assertEqual(val[0][0],  2.0)
        self.assertEqual(val[0][1],  0.0)
        self.assertEqual(val[1][0],  2.0)
        self.assertEqual(val[1][1], -2.0)

        # second derivatives at midpoint
        val = crv.Evaluate(0.5, 2)  
        self.assertEqual(val[2][0],  0.0)
        self.assertEqual(val[2][1], -4.0)

        # third derivatives should all vanish everywhere
        val = crv.Evaluate(0.5, 3)  
        self.assertEqual(val[3][0],  0.0)
        self.assertEqual(val[3][1],  0.0)

        # check integer type for derivative
        with self.assertRaises(TypeError):
            val = crv.Evaluate(0.5, 1.5) 

        # GoTools throws exception for negative derivatives
        with self.assertRaises(Exception):
            val = crv.Evaluate(0.5, -1) 

if __name__ == '__main__':
    unittest.main()
