from __future__ import annotations

from collections.abc import Sequence
from itertools import chain, repeat
from math import pi
from typing import TextIO

import numpy as np

from splipy.typing import FloatArray

rng = np.random.default_rng()


def gen_knot(n: int, p: int, periodic: int) -> FloatArray:
    m = n + (periodic + 2)
    result = np.fromiter(
        chain(
            repeat(0.0, p),
            range(1, m - p),
            repeat(m - p, p),
        ),
        dtype=np.float64
    )
    result[p : m - 1] += (rng.random(m - p - 1) - 0.5) * 0.8
    result = np.round(result * 10) / 10  # set precision to one digit
    for i in range(periodic + 1):
        result[-periodic + i - 1] += result[p + i]
        result[i] -= m - p - result[-p - periodic + i - 1]
    return result


def gen_cp_curve(n: int, dim: int, periodic: int) -> FloatArray:
    cp = np.zeros((n, dim), dtype=np.float64)
    if periodic > -1:
        t = np.linspace(0, 2 * pi, n + 1)
        cp[:, 0] = 100 * np.cos(t[:-1])
        cp[:, 1] = 100 * np.sin(t[:-1])
    else:
        cp[:, 0] = np.linspace(0, 100, n)
    cp += (rng.random((n, dim)) - 0.5) * 10
    return np.floor(cp)


def gen_cp_surface(n: tuple[int, int], dim: int, periodic: int) -> FloatArray:
    cp = np.zeros((n[0], n[1], dim), dtype=np.float64)
    if periodic > -1:
        t = np.linspace(0, 2 * pi, n[0] + 1)
        r = np.linspace(60, 100, n[1])
        R, T = np.meshgrid(r, t[:-1])
        cp[:, :, 0] = R * np.cos(T)
        cp[:, :, 1] = R * np.sin(T)
    else:
        y, x = np.meshgrid(np.linspace(0, 100, n[1]), np.linspace(0, 100, n[0]))
        cp[:, :, 0] = x
        cp[:, :, 1] = y
    cp += (rng.random((n[0], n[1], dim)) - 0.5) * 10
    return np.floor(cp.transpose(1, 0, 2))


def gen_cp_volume(n: tuple[int, int, int], dim: int, periodic: int) -> FloatArray:
    cp = np.zeros((n[0], n[1], n[2], dim), dtype=np.float64)
    if periodic > -1:
        t = np.linspace(0, 2 * pi, n[0] + 1)
        r = np.linspace(50, 100, n[1])
        z = np.linspace(0, 40, n[2])
        R, T, Z = np.meshgrid(r, t[:-1], z)
        cp[:, :, :, 0] = R * np.cos(T)
        cp[:, :, :, 1] = R * np.sin(T)
        cp[:, :, :, 2] = Z
    else:
        y, x, z = np.meshgrid(np.linspace(0, 100, n[1]), np.linspace(0, 100, n[0]), np.linspace(0, 100, n[2]))
        cp[:, :, :, 0] = x
        cp[:, :, :, 1] = y
        cp[:, :, :, 2] = z
    cp += (rng.random((n[0], n[1], n[2], dim)) - 0.5) * 10
    return np.floor(cp.transpose(2, 1, 0, 3))


def gen_controlpoints(
    n: tuple[int] | tuple[int, int] | tuple[int, int, int],
    dim: int,
    rational: bool,
    periodic: int,
) -> FloatArray:
    if len(n) == 1:  # curve
        cp = gen_cp_curve(n[0], dim, periodic)
        total_n = n[0]
    elif len(n) == 2:  # surface
        cp = gen_cp_surface(n, dim, periodic)
        total_n = n[0] * n[1]
    elif len(n) == 3:  # volume
        cp = gen_cp_volume(n, dim, periodic)
        total_n = n[0] * n[1] * n[2]

    cp = np.reshape(cp, (total_n, dim))

    if rational:
        w = rng.random(total_n) + 0.5
        w = np.round(w * 10) / 10
        cp = np.insert(cp, dim, w, 1)

    return cp


def write_basis(f: TextIO, p: Sequence[int], knot: Sequence[FloatArray], periodic: int) -> None:
    for i in range(len(knot)):
        f.write(f"        basis{i} = BSplineBasis(" + str(p[i]) + ", np." + repr(knot[i]))
        if periodic > -1 and i == 0:  # only consider periodicity in the first direction
            f.write("," + str(periodic))
        f.write(")\n")


def get_name(
    n: tuple[int] | tuple[int, int] | tuple[int, int, int],
    p: Sequence[int],
    dim: int,
    rational: bool,
    periodic: int,
) -> str:
    result = ""
    if len(n) == 1:
        result += "curve"
    elif len(n) == 2:
        result += "surface"
    elif len(n) == 3:
        result += "volume"
    result += "_" + str(dim) + "D"
    result += "_p"
    for q in p:
        result += str(q)
    if rational:
        result += "_rational"
    if periodic > -1:
        result += f"_C{periodic}_periodic"
    return result


def raise_order(p: int) -> str:
    result = "        crv2.raise_order(2)\n"
    result += "        p = crv2.order()\n"
    result += f"        self.assertEqual(p, {p + 2})\n"
    return result


def write_object_creation(f: TextIO, rational: bool, pardim: int, clone: bool = True) -> None:
    if pardim == 1:
        f.write("        crv  = Curve(basis0, controlpoints," + str(rational) + ")\n")
        if clone:
            f.write("        crv2 = crv.clone()\n")
    elif pardim == 2:
        f.write("        surf  = Surface(basis0, basis1, controlpoints," + str(rational) + ")\n")
        if clone:
            f.write("        surf2 = surf.clone()\n")
    elif pardim == 3:
        f.write("        vol  = Volume(basis0, basis1, basis2, controlpoints," + str(rational) + ")\n")
        if clone:
            f.write("        vol2 = vol.clone()\n")


def evaluate_curve() -> str:
    return """
        u    = np.linspace(crv.start(0), crv.end(0), 13)
        pt   = crv(u)
        pt2  = crv2(u)
"""


def evaluate_surface() -> str:
    return """
        u    = np.linspace(surf.start(0), surf.end(0), 9)
        v    = np.linspace(surf.start(1), surf.end(1), 9)
        pt   = surf(u,v)
        pt2  = surf2(u,v)

"""


def evaluate_volume() -> str:
    return """
        u    = np.linspace(vol.start(0), vol.end(0), 7)
        v    = np.linspace(vol.start(1), vol.end(1), 7)
        w    = np.linspace(vol.start(2), vol.end(2), 7)
        pt   = vol(u,v,w)
        pt2  = vol2(u,v,w)
"""
