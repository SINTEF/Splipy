from __future__ import annotations

import unittest
from pathlib import Path

import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory
from splipy.io import SVG

THIS_DIR = str(Path(__file__).parent)


class TestSVG(unittest.TestCase):
    def test_read_linear(self):
        with SVG(THIS_DIR + "/geometries/linear_curve.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 2)
        self.assertEqual(len(crv), 2)

    def test_read_cubic(self):
        with SVG(THIS_DIR + "/geometries/cubic_curve.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 4)
        self.assertEqual(len(crv), 7)

    def test_read_cubic_three_knot_span(self):
        with SVG(THIS_DIR + "/geometries/cubic_three_knot_span_curve.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 4)
        self.assertEqual(len(crv), 10)

    def test_read_multiple_curves(self):
        with SVG(THIS_DIR + "/geometries/three_curves.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 3)
        self.assertEqual(crv[0].order(0), 4)
        self.assertEqual(len(crv[0]), 10)
        self.assertEqual(crv[1].order(0), 4)
        self.assertEqual(len(crv[1]), 10)
        self.assertEqual(crv[2].order(0), 2)
        self.assertEqual(len(crv[2]), 4)

    def test_write_curve(self):
        crv = CurveFactory.polygon([[0, 0], [1, 0], [1, 1], [0, 1]])
        with SVG("output.svg") as myfile:
            myfile.write(crv)
        self.assertTrue(Path("output.svg").is_file())
        Path("output.svg").unlink()

    def test_write_surface(self):
        srf = SurfaceFactory.disc(type="square").rebuild([4, 4], [7, 7])
        with SVG("output.svg") as myfile:
            myfile.write(srf)
        self.assertTrue(Path("output.svg").is_file())
        Path("output.svg").unlink()

    def test_read_github(self):
        with SVG(THIS_DIR + "/geometries/github.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 4)
        self.assertEqual(len(crv), 76)

    def test_read_california(self):
        with SVG(THIS_DIR + "/geometries/california.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 1)
        crv = crv[0]
        self.assertEqual(crv.order(0), 4)
        self.assertEqual(len(crv), 307)

    def test_read_isengard(self):
        with SVG(THIS_DIR + "/geometries/isengard.svg") as myfile:
            crv = myfile.read()
        self.assertEqual(len(crv), 9)
        for i in range(9):
            bb = crv[i].bounding_box()
            self.assertGreater(bb[0][0], -1)
            self.assertLess(bb[0][1], 201)
            self.assertGreater(bb[1][0], -1)
            self.assertLess(bb[1][1], 221)


if __name__ == "__main__":
    unittest.main()
