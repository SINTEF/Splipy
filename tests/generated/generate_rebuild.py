from __future__ import annotations

import datetime
import os
import subprocess

import numpy as np
from general import *


def rebuild(p, knot):
    obj_name = ["", "crv", "surf", "vol"]
    pardim = len(knot)
    p += np.random.randint(-1, 2, pardim)
    p = [max(q, 2) for q in p]
    n = [len(k) for k in knot]
    if pardim == 1:
        n += np.random.randint(300, 350, pardim)
    elif pardim == 2:
        n += np.random.randint(110, 140, pardim)
    elif pardim == 3:
        n += np.random.randint(70, 90, pardim)
    result = "        %s2 = %s.rebuild(" % (obj_name[pardim], obj_name[pardim])
    if len(p) == 1 or np.random.random() > 0.75:
        result += str(p[0]) + ", " + str(n[0]) + ")"
    else:
        result += repr(list(p)) + ", " + repr(list(n)) + ")"
    return result


# dump large sets of control points
np.set_printoptions(threshold=np.inf, linewidth=100)
np.set_printoptions(threshold=np.inf)
# fetch random state
state = np.random.get_state()  # might consider dumping this for reproducibility
# open file and write headers
f = open("rebuild_test.py", "w")
f.write("# --- Automatic generated test file  ---\n")
f.write("# Generator    : " + os.path.basename(__file__) + "\n")
f.write("# Date         : " + str(datetime.date.today()) + "\n")
f.write("# Git revision : " + subprocess.check_output(["git", "rev-parse", "HEAD"]))
f.write("""
import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis
from math import sqrt
import unittest


class TestRebuild(unittest.TestCase):
""")

evaluate = [None, evaluate_curve, evaluate_surface, evaluate_volume]
for baseP in [2, 4, 5]:
    for dim in [2, 3]:
        for rational in [True, False]:
            for periodic in [-1, 1]:
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
                    write_object_creation(f, rational, pardim, False)
                    f.write(rebuild(p, knot))
                    f.write(evaluate[pardim]())
                    precision = 4 - pardim
                    if periodic > -1:
                        precision -= 2
                    if rational:
                        precision -= 1
                    precision = max(1, precision)
                    f.write(
                        "        self.assertAlmostEqual(np.max(np.abs(pt-pt2)), 0.0, places=%d)\n\n"
                        % (precision)
                    )

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print("rebuild_test.py written")
