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


def rebuild(p: Sequence[int], knot: Sequence[FloatArray]) -> str:
    obj_name = ["crv", "surf", "vol"]
    pardim = len(knot)
    pp = [int(x) for x in rng.integers(-1, 2, pardim) + p]
    pp = [max(q, 2) for q in pp]
    n = [len(k) for k in knot]
    if pardim == 1:
        n += rng.integers(300, 350, pardim)
    elif pardim == 2:
        n += rng.integers(110, 140, pardim)
    elif pardim == 3:
        n += rng.integers(70, 90, pardim)
    result = f"        {obj_name[pardim - 1]}2 = {obj_name[pardim]}.rebuild("
    if len(pp) == 1 or rng.random() > 0.75:
        result += str(pp[0]) + ", " + str(n[0]) + ")"
    else:
        result += repr(list(pp)) + ", " + repr(list(n)) + ")"
    return result


# dump large sets of control points
np.set_printoptions(threshold=sys.maxsize, linewidth=100)
# open file and write headers
f = Path("rebuild_test.py").open("w")
f.write("# --- Automatic generated test file  ---\n")
f.write("# Generator    : " + Path(__file__).name + "\n")
f.write("# Date         : " + str(datetime.date.today()) + "\n")
f.write("# Git revision : " + subprocess.check_output(["git", "rev-parse", "HEAD"]).decode())
f.write("""
import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis
from math import sqrt
import unittest


class TestRebuild(unittest.TestCase):
""")

evaluate: list[Callable[[], str]] = [
    general.evaluate_curve,
    general.evaluate_surface,
    general.evaluate_volume,
]
for baseP in [2, 4, 5]:
    for dim in [2, 3]:
        for rational in [True, False]:
            for periodic in [-1, 1]:
                for pardim in range(1, 4):
                    p = tuple(int(x) for x in np.array([baseP] * pardim) + rng.integers(0, 3, pardim))
                    if periodic >= p[0] - 1:
                        continue
                    if dim < 3 and pardim > 2:
                        continue
                    n = [int(x) for x in rng.integers(-2, 3, pardim) + p + 3]
                    n[0] = np.maximum(2 * (p[0] + periodic), n[0])
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
                    general.write_object_creation(f, rational, pardim, False)
                    f.write(rebuild(p, knot))
                    f.write(evaluate[pardim - 1]())
                    precision = 4 - pardim
                    if periodic > -1:
                        precision -= 2
                    if rational:
                        precision -= 1
                    precision = max(1, precision)
                    f.write(
                        f"        self.assertAlmostEqual(np.max(np.abs(pt-pt2)), 0.0, places={precision})\n\n"
                    )

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print("rebuild_test.py written")
