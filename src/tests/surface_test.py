from GoTools import Surface
import unittest

class TestSurface(unittest.TestCase):
    def test_constructor(self):
        # test 3D constructor
        controlpoints = [[0,0,0], [1,0,0], [0,1,0], [1,1,0]];
        surf = Surface(2,2, [0,0,1,1], [0,0,1,1], controlpoints);
        val = surf.Evaluate(0.5, 0.5)
        self.assertEqual(val[0], 0.5)

        # test 2D constructor
        controlpoints = [[0,0], [1,0], [0,1], [1,1]];
        surf2 = Surface(2,2, [0,0,1,1], [0,0,1,1], controlpoints);
        val  = surf2.Evaluate(0.5, 0.5)
        self.assertEqual(val[0], 0.5)

    def test_evaluate(self):
        # unit square, identity mapping, p=2
        controlpoints = [[0,0], [.5,0],[1,0],  [0,.5], [.5,.5],[1,.5],  [0,1], [.5,1], [1,1]];
        surf = Surface(3,3, [0,0,0,1,1,1], [0,0,0,1,1,1], controlpoints)

        # startpoint evaluation
        val = surf.Evaluate(0,0)     
        self.assertEqual(val[0], 0.0)
        self.assertEqual(val[1], 0.0)

        # startpoint with derivatives
        val = surf.Evaluate(0,0,  1)  
        self.assertEqual(val[0][0], 0.0)
        self.assertEqual(val[0][1], 0.0)
        self.assertEqual(val[1][0], 1.0)
        self.assertEqual(val[1][1], 0.0)
        self.assertEqual(val[2][0], 0.0)
        self.assertEqual(val[2][1], 1.0)

        # midpoint with derivatives
        val = surf.Evaluate(0.5, 0.5,   2)  
        self.assertEqual(val[0][0], 0.5)
        self.assertEqual(val[0][1], 0.5)
        self.assertEqual(val[1][0], 1.0)
        self.assertEqual(val[1][1], 0.0)
        self.assertEqual(val[2][0], 0.0)
        self.assertEqual(val[2][1], 1.0)

        # mid top edge with derivatives (sensitive to left-evaluation)
        val = surf.Evaluate(0.5, 1.0,   1)  
        self.assertEqual(val[0][0], 0.5)
        self.assertEqual(val[0][1], 1.0)
        self.assertEqual(val[1][0], 1.0)
        self.assertEqual(val[1][0], 1.0)
        self.assertEqual(val[1][1], 0.0)
        self.assertEqual(val[2][0], 0.0)
        self.assertEqual(val[2][1], 1.0)

        # mid right edge with derivatives (sensitive to left-evaluation)
        val = surf.Evaluate(1.0, 0.5,   1)  
        self.assertEqual(val[0][0], 1.0)
        self.assertEqual(val[0][1], 0.5)
        self.assertEqual(val[1][0], 1.0)
        self.assertEqual(val[1][1], 0.0)
        self.assertEqual(val[2][0], 0.0)
        self.assertEqual(val[2][1], 1.0)

        # top right corner with derivatives (sensitive to left-evaluation)
        val = surf.Evaluate(1.0, 1.0,   1)  
        self.assertEqual(val[0][0], 1.0)
        self.assertEqual(val[0][1], 1.0)
        self.assertEqual(val[1][0], 1.0)
        self.assertEqual(val[1][1], 0.0)
        self.assertEqual(val[2][0], 0.0)
        self.assertEqual(val[2][1], 1.0)

        # second derivatives at midpoint (is actually zero due to identity mapping)
        val = surf.Evaluate(0.5, 0.5,   2)  
        self.assertEqual(val[3][0],  0.0)
        self.assertEqual(val[3][1],  0.0)
        self.assertEqual(val[4][0],  0.0)
        self.assertEqual(val[4][1],  0.0)
        self.assertEqual(val[5][0],  0.0)
        self.assertEqual(val[5][1],  0.0)

        # third derivatives should all vanish everywhere
        val = surf.Evaluate(0.5, 0.5,   3)  
        self.assertEqual(len(val),  10)
        self.assertEqual(val[6][0],  0.0)
        self.assertEqual(val[6][1],  0.0)
        self.assertEqual(val[7][0],  0.0)
        self.assertEqual(val[7][1],  0.0)
        self.assertEqual(val[8][0],  0.0)
        self.assertEqual(val[9][0],  0.0)

        # check integer type for derivative
        with self.assertRaises(TypeError):
            val = surf.Evaluate(0.5,0.5,  1.5) 

        # GoTools throws exception for negative derivatives
        with self.assertRaises(Exception):
            val = surf.Evaluate(0.5,0.5,  -1) 

if __name__ == '__main__':
    unittest.main()
