from __future__ import annotations

import datetime
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from . import general

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

    from splipy.typing import FloatArray

# Get random generator
rng = np.random.default_rng()


def insert_knot(p: Sequence[int], knot: Sequence[FloatArray]) -> str:
    obj_name = ["crv", "surf", "vol"]
    pardim = len(knot)
    result = ""
    for i in range(pardim):
        n = len(knot[i])
        q = p[i]
        dt = knot[i][-q - 1] - knot[i][q]  # knot vector range
        for j in range(q - 1):
            k = rng.random() * dt + knot[i][q]  # insert random knot
            n = n + j + 1
            if pardim == 1:
                result += (
                    f"        {obj_name[pardim - 1]}2.insert_knot("
                    f"[{k:.3f}]*{j + 1}) # insert C{q - j - 2}-knot\n"
                )
            else:
                result += (
                    f"        {obj_name[pardim - 1]}2.insert_knot([{k:.3f}]*"
                    f"{j + 1}, {i}) # insert C{q - j - 2}-knot\n"
                )
            result += (
                f"        self.assertEqual(len({obj_name[pardim - 1]}2.knots("
                f"direction={i},with_multiplicities=True)), {n})\n"
            )
    return result


# dump large sets of control points
np.set_printoptions(threshold=sys.maxsize, linewidth=100)
# open file and write headers
f = Path("knot_insert_test.py").open("w")
f.write("# --- Automatic generated test file  ---\n")
f.write("# Generator    : " + Path(__file__).name + "\n")
f.write("# Date         : " + str(datetime.date.today()) + "\n")
f.write("# Git revision : " + subprocess.check_output(["git", "rev-parse", "HEAD"]).decode())
f.write("""
import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis, SplineObject
from math import sqrt
import unittest


class TestInsertKnot(unittest.TestCase):
""")

evaluate: list[Callable[[], str]] = [
    general.evaluate_curve,
    general.evaluate_surface,
    general.evaluate_volume,
]
for baseP in [2, 5]:
    for dim in [2, 3]:
        for rational in [True, False]:
            for periodic in [-1, 0, 1, 2]:
                for pardim in range(1, 4):
                    p = tuple(int(x) for x in np.array([baseP] * pardim) + rng.integers(0, 3, pardim))
                    if periodic >= p[0] - 1:
                        continue
                    if dim < 3 and pardim > 2:
                        continue
                    n = [int(x) for x in rng.integers(-2, 3, pardim) + p + 3]
                    n[0] = max(2 * (p[0] + periodic), n[0])
                    cp = general.gen_controlpoints(n, dim, rational, periodic)
                    knot = [
                        general.gen_knot(n[i], p[i], (i == 0) * (periodic + 1) - 1) for i in range(pardim)
                    ]
                    f.write(
                        "    def test_"
                        + general.get_name(n, p, dim, rational, periodic)
                        + "(self) -> None:\n"
                    )
                    cp_str = repr(cp).replace("\n", "") if cp.shape[0] > 30 else repr(cp)
                    f.write("        controlpoints = np." + cp_str + "\n")
                    general.write_basis(f, p, knot, periodic)
                    general.write_object_creation(f, rational, pardim)
                    f.write(insert_knot(p, knot))
                    f.write(evaluate[pardim - 1]())
                    f.write("        self.assertAlmostEqual(float(np.linalg.norm(pt-pt2)), 0.0)\n\n")

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print("knot_insert_test.py written")
