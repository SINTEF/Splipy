from GoTools import Curve
import unittest

class TestCurve(unittest.TestCase):
    def test_constructor(self):
        crv = Curve(2, [0, 0, 1, 1], [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        val = crv.Evaluate(0.5)
        self.assertEqual(val[0], 0.5)

if __name__ == '__main__':
    unittest.main()
