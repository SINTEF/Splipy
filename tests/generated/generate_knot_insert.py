from __future__ import annotations

import datetime
import os
import subprocess

import numpy as np
from general import *


def insert_knot(p, knot):
    obj_name = ["", "crv", "surf", "vol"]
    pardim = len(knot)
    result = ""
    for i in range(pardim):
        n = len(knot[i])
        q = p[i]
        dt = knot[i][-q - 1] - knot[i][q]  # knot vector range
        for j in range(q - 1):
            k = np.random.rand() * dt + knot[i][q]  # insert random knot
            n = n + j + 1
            if pardim == 1:
                result += "        %s2.insert_knot([%.3f]*%d) # insert C%d-knot\n" % (
                    obj_name[pardim],
                    k,
                    j + 1,
                    q - j - 2,
                )
            else:
                result += "        %s2.insert_knot([%.3f]*%d, %d) # insert C%d-knot\n" % (
                    obj_name[pardim],
                    k,
                    j + 1,
                    i,
                    q - j - 2,
                )
            result += (
                "        self.assertEqual(len(%s2.knots(direction=%d,with_multiplicities=True)), %d)\n"
                % (obj_name[pardim], i, n)
            )
    return result


# dump large sets of control points
np.set_printoptions(threshold=np.inf, linewidth=100)
np.set_printoptions(threshold=np.inf)
# fetch random state
state = np.random.get_state()  # might consider dumping this for reproducibility
# open file and write headers
f = open("knot_insert_test.py", "w")
f.write("# --- Automatic generated test file  ---\n")
f.write("# Generator    : " + os.path.basename(__file__) + "\n")
f.write("# Date         : " + str(datetime.date.today()) + "\n")
f.write("# Git revision : " + subprocess.check_output(["git", "rev-parse", "HEAD"]))
f.write("""
import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis
from math import sqrt
import unittest


class TestInsertKnot(unittest.TestCase):
""")

evaluate = [None, evaluate_curve, evaluate_surface, evaluate_volume]
for baseP in [2, 5]:
    for dim in [2, 3]:
        for rational in [True, False]:
            for periodic in [-1, 0, 1, 2]:
                for pardim in range(1, 4):
                    p = np.array([baseP] * pardim) + np.random.randint(0, 3, pardim)
                    if periodic >= p[0] - 1:
                        continue
                    if dim < 3 and pardim > 2:
                        continue
                    n = p + 3
                    n += np.random.randint(-2, 3, pardim)
                    n[0] = np.maximum(2 * (p[0] + periodic), n[0])
                    cp = gen_controlpoints(n, dim, rational, periodic)
                    knot = [gen_knot(n[i], p[i], (i == 0) * (periodic + 1) - 1) for i in range(pardim)]
                    f.write("    def test_" + get_name(n, p, dim, rational, periodic) + "(self):\n")
                    if cp.shape[0] > 30:
                        cp_str = repr(cp).replace("\n", "")
                    else:
                        cp_str = repr(cp)
                    f.write("        controlpoints = np." + cp_str + "\n")
                    write_basis(f, p, knot, periodic)
                    write_object_creation(f, rational, pardim)
                    f.write(insert_knot(p, knot))
                    f.write(evaluate[pardim]())
                    f.write("        self.assertAlmostEqual(np.linalg.norm(pt-pt2), 0.0)\n\n")

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print("knot_insert_test.py written")
