from BSplineBasis import *
import numpy as np
import unittest

class TestBasis(unittest.TestCase):
    def test_get_continuity(self):
        b = BSplineBasis(4, [0,0,0,0,.3,1,1,1.134,1.134,1.134, 2,2,2,2])
        self.assertEqual(b.get_continuity( .3  ), 2)
        self.assertEqual(b.get_continuity(1    ), 1)
        self.assertEqual(b.get_continuity(1.134), 0)
        self.assertEqual(b.get_continuity(0    ), -1)
        self.assertEqual(b.get_continuity(2    ), -1)
        self.assertEqual(b.get_continuity(.4   ), np.inf)


if __name__ == '__main__':
    unittest.main()
