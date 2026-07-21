from __future__ import annotations

import unittest
from math import pi

import numpy as np
from scipy.integrate import quad

import splipy.curve_factory as cf
import splipy.surface_factory as sf
import splipy.volume_factory as vf
from splipy.curve_factory import Boundary


def relerr(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


class TestQuadratureConvergence(unittest.TestCase):
    """Length/area/volume/error integrate square roots of polynomials, which a fixed
    order+1 Gauss rule cannot resolve (#188). These tolerances are tighter than that
    rule reaches, so they fail unless the quadrature is refined until it converges.
    """

    def test_curve_length_multispan(self):
        # multi-span cubic with strong per-span speed variation; oracle = adaptive quad
        crv = cf.cubic_curve(
            np.array([[0, 0], [0.2, 2.5], [0.4, 0.1], [5, 0], [9.6, 0.1], [9.8, 2.5], [10, 0]], dtype=float),
            boundary=Boundary.NATURAL,
        )

        def speed(t):
            return float(np.linalg.norm(np.asarray(crv.derivative(t))))

        knots = np.unique(crv.knots(0))
        oracle = sum(quad(speed, a, b, epsabs=1e-13, epsrel=1e-13)[0] for a, b in zip(knots[:-1], knots[1:]))
        self.assertLess(relerr(crv.length(), oracle), 1e-8)

    def test_curve_length_rational_circle(self):
        # rational (NURBS) curve: speed is a root of a rational function
        self.assertAlmostEqual(cf.circle(r=1).length(), 2 * pi, places=9)
        self.assertAlmostEqual(cf.circle(r=3).length(), 2 * pi * 3, places=8)

    def test_surface_area_nurbs_sphere(self):
        self.assertAlmostEqual(sf.sphere(r=1).area(), 4 * pi, places=9)
        self.assertAlmostEqual(sf.sphere(r=2).area(), 4 * pi * 4, places=8)

    def test_volume_rational_solids(self):
        self.assertAlmostEqual(vf.sphere(r=1).volume(), 4 / 3 * pi, places=9)
        self.assertAlmostEqual(vf.cylinder(r=1, h=2).volume(), 2 * pi, places=9)

    def test_volume_nonrational_stays_exact(self):
        # polynomial integrand: order+1 is already exact and must remain so after refining
        cube = vf.cube(size=2)
        self.assertAlmostEqual(cube.volume(), 8.0, places=12)

    def test_curve_error_l2(self):
        # L2 error against a non-polynomial target; oracle = adaptive quad of |x_h - x|^2
        crv = cf.circle(r=2)

        def target(t):
            a = np.asarray(t, float)
            return np.stack([2 * np.cos(a), 2 * np.sin(a)], axis=-1)

        err2, _ = crv.error(target)

        def integrand(t):
            diff = np.asarray(crv(t)) - np.array([2 * np.cos(t), 2 * np.sin(t)])
            return float(np.dot(diff, diff))

        knots = np.unique(crv.knots(0))
        oracle = sum(
            quad(integrand, a, b, epsabs=1e-13, epsrel=1e-13)[0] for a, b in zip(knots[:-1], knots[1:])
        )
        self.assertLess(relerr(float(np.sum(err2)), oracle), 1e-8)


if __name__ == "__main__":
    unittest.main()
