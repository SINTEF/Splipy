from __future__ import annotations

import datetime
import subprocess
import sys
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

import numpy as np

from . import general

if TYPE_CHECKING:
    from collections.abc import Callable

    from splipy.typing import FloatArray

# Get random generator
rng = np.random.default_rng()


def do_translate(f: TextIO, dim: int, pardim: int) -> str:
    obj_name = ["crv", "surf", "vol"]
    amount = rng.random(dim) * 10
    amount = np.floor(amount * 10) / 10
    f.write(f"        {obj_name[pardim - 1]}2 += np." + repr(amount) + "\n")
    return expect_translate(amount)


def expect_translate(amount: FloatArray) -> str:
    result = ""
    for i in range(len(amount)):
        result += f"        pt2[...,{i}] -= {amount[i]}\n"
    return result


def do_scale(f: TextIO, dim: int, pardim: int) -> str:
    obj_name = ["crv", "surf", "vol"]
    amount = rng.random(dim) * 10 + 2
    amount = np.floor(amount * 10) / 10
    f.write(f"        {obj_name[pardim - 1]}2 *= np." + repr(amount) + "\n")
    return expect_scale(amount)


def expect_scale(amount: FloatArray) -> str:
    result = ""
    for i in range(len(amount)):
        result += f"        pt2[...,{i}] /= {amount[i]}\n"
    return result


def do_mirror(f: TextIO, dim: int, pardim: int) -> str:
    obj_name = ["", "crv", "surf", "vol"]
    axis = rng.integers(0, 3)
    axes = [0, 0, 0]
    axes[axis] = rng.integers(1, 4)
    f.write(f"        {obj_name[pardim]}2.mirror(" + repr(axes) + ")\n")
    return expect_mirror(axis)


def expect_mirror(axis: int) -> str:
    return f"        pt2[...,{axis}] = -pt2[...,{axis}]\n"


# dump large sets of control points
np.set_printoptions(threshold=sys.maxsize)
# open file and write headers
f = Path("affine_test.py").open("w")
f.write("# --- Automatic generated test file  ---\n")
f.write("# Generator    : " + Path(__file__).name + "\n")
f.write("# Date         : " + str(datetime.date.today()) + "\n")
f.write("# Git revision : " + subprocess.check_output(["git", "rev-parse", "HEAD"]).decode())
f.write("""
import numpy as np
from splipy import Volume, Surface, Curve, BSplineBasis, SplineObject
from math import sqrt
import unittest


class TestAffine(unittest.TestCase):
""")

evaluate: list[Callable[[], str]] = [
    general.evaluate_curve,
    general.evaluate_surface,
    general.evaluate_volume,
]
do: list[Callable[[TextIO, int, int], str]] = [do_translate, do_scale, do_mirror]
name = ["translate", "scale", "mirror"]
for baseP in [2, 5]:
    for dim in [2, 3]:
        for rational in [True, False]:
            for periodic in [-1, 0, 1]:
                for pardim in range(1, 4):
                    for j in range(len(do)):
                        p = tuple(int(x) for x in np.array([baseP] * pardim) + rng.integers(0, 3, pardim))
                        if periodic >= p[0] - 1:
                            continue
                        if dim < 3 and pardim > 2:
                            continue
                        if dim < 3 and name[j] == "mirror":
                            continue
                        n = tuple(int(x) for x in rng.integers(-2, 3, pardim) + p + 3)
                        cp = general.gen_controlpoints(n, dim, rational, periodic)
                        knot = [
                            general.gen_knot(n[i], p[i], (i == 0) * (periodic + 1) - 1) for i in range(pardim)
                        ]
                        f.write(
                            f"    def test_{name[j]}_"
                            + general.get_name(n, p, dim, rational, periodic)
                            + "(self) -> None:\n"
                        )
                        cp_str = repr(cp).replace("\n", "") if cp.shape[0] > 30 else repr(cp)
                        f.write("        controlpoints = np." + cp_str + "\n")
                        general.write_basis(f, p, knot, periodic)
                        general.write_object_creation(f, rational, pardim)
                        new_dim = dim if name[j] == "scale" else min(3, dim + rng.integers(0, 2))
                        expect = do[j](f, new_dim, pardim)
                        f.write(evaluate[pardim - 1]())
                        f.write(expect)
                        if new_dim != dim:
                            f.write("        allZero           = pt2\n")
                            f.write("        allZero[...,:-1] -= pt \n")
                            f.write("        self.assertAlmostEqual(float(np.linalg.norm(allZero)), 0.0)\n\n")
                        else:
                            f.write("        self.assertAlmostEqual(float(np.linalg.norm(pt-pt2)), 0.0)\n\n")

f.write("""
if __name__ == '__main__':
    unittest.main()
""")
print("affine_test.py written")
