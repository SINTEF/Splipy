from __future__ import annotations

import datetime
import subprocess
from pathlib import Path

import numpy as np

from . import general

# Get random generator
rng = np.random.default_rng()


def diff_curve():
    return """
        u    = np.linspace(crv.start(0), crv.end(0), 11)
        du   = crv.derivative(u)
        du2  = crv2(u)
"""


def diff_surface():
    return """
        u    = np.linspace(surf.start(0), surf.end(0), 7)
        v    = np.linspace(surf.start(1), surf.end(1), 7)
        du   = surf.derivative(u,v, d=(1,0))
        du2  = surf2(u,v)
        dv   = surf.derivative(u,v, d=(0,1))
        dv2  = surf3(u,v)

"""


def diff_volume():
    return """
        u    = np.linspace(vol.start(0), vol.end(0), 5)
        v    = np.linspace(vol.start(1), vol.end(1), 5)
        w    = np.linspace(vol.start(2), vol.end(2), 5)
        du   = vol.derivative(u,v,w, d=(1,0,0))
        du2  = vol2(u,v,w)
        dv   = vol.derivative(u,v,w, d=(0,1,0))
        dv2  = vol3(u,v,w)
        dw   = vol.derivative(u,v,w, d=(0,0,1))
        dw2  = vol4(u,v,w)
"""


def differentiate(p):
    obj_name = ["", "crv", "surf", "vol"]
    pardim = len(p)
    result = ""
    for i in range(pardim):
        result += f"        {obj_name[pardim]}{i + 2} = {obj_name[pardim]}.get_derivative_spline({i})\n"
        expected_order = list(p)
        expected_order[i] = p[i] - 1
        for j in range(pardim):
            result += (
                f"        self.assertEqual({obj_name[pardim]}{i + 2}"
                f".order(direction={j}), {expected_order[j]})\n"
            )
    return result


# dump large sets of control points
np.set_printoptions(threshold=np.inf, linewidth=100)
np.set_printoptions(threshold=np.inf)
# open file and write headers
f = Path("derivative_spline_test.py").open("w")
f.write("# --- Automatic generated test file  ---\n")
f.write("# Generator    : " + Path(__file__).name + "\n")
f.write("# Date         : " + str(datetime.date.today()) + "\n")
f.write("# Git revision : " + subprocess.check_output(["git", "rev-parse", "HEAD"]).decode())
f.write("""
import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis
from math import sqrt
import unittest


class TestDerivativeSpline(unittest.TestCase):
""")

evaluate = [None, diff_curve, diff_surface, diff_volume]
for baseP in [2, 5]:
    for dim in [2, 3]:
        for rational in [False]:
            for periodic in [-1, 0, 1, 2]:
                for pardim in range(1, 4):
                    p = np.array([baseP] * pardim) + rng.integers(0, 3, pardim)
                    if periodic >= p[0] - 1:
                        continue
                    if dim < 3 and pardim > 2:
                        continue
                    n = p + 3
                    n += rng.integers(-2, 1, pardim)
                    n[0] = np.maximum(2 * (p[0] + periodic), n[0])
                    cp = general.gen_controlpoints(n, dim, rational, periodic)
                    knot = [
                        general.gen_knot(n[i], p[i], (i == 0) * (periodic + 1) - 1) for i in range(pardim)
                    ]
                    f.write("    def test_" + general.get_name(n, p, dim, rational, periodic) + "(self):\n")
                    cp_str = repr(cp).replace("\n", "") if cp.shape[0] > 30 else repr(cp)
                    f.write("        controlpoints = np." + cp_str + "\n")
                    general.write_basis(f, p, knot, periodic)
                    general.write_object_creation(f, rational, pardim)
                    f.write(differentiate(p))
                    f.write(evaluate[pardim]())
                    f.write("        self.assertAlmostEqual(np.linalg.norm(du-du2), 0.0)\n")
                    if pardim > 1:
                        f.write("        self.assertAlmostEqual(np.linalg.norm(dv-dv2), 0.0)\n")
                    if pardim > 2:
                        f.write("        self.assertAlmostEqual(np.linalg.norm(dw-dw2), 0.0)\n")

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print("derivative_spline_test.py written")
