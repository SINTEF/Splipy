from __future__ import annotations

from math import pi

import numpy as np


def gen_knot(n, p, periodic):
    m = n + (periodic + 2)
    result = [0] * p + range(1, m - p) + [m - p] * p
    result = np.array(result, "float")
    result[p : m - 1] += (np.random.rand(m - p - 1) - 0.5) * 0.8
    result = np.round(result * 10) / 10  # set precision to one digit
    for i in range(periodic + 1):
        result[-periodic + i - 1] += result[p + i]
        result[i] -= m - p - result[-p - periodic + i - 1]
    return result


def gen_cp_curve(n, dim, periodic):
    cp = np.zeros((n, dim))
    if periodic > -1:
        t = np.linspace(0, 2 * pi, n + 1)
        cp[:, 0] = 100 * np.cos(t[:-1])
        cp[:, 1] = 100 * np.sin(t[:-1])
    else:
        cp[:, 0] = np.linspace(0, 100, n)
    cp += (np.random.rand(n, dim) - 0.5) * 10
    return np.floor(cp)


def gen_cp_surface(n, dim, periodic):
    cp = np.zeros((n[0], n[1], dim))
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
    cp += (np.random.rand(n[0], n[1], dim) - 0.5) * 10
    return np.floor(cp.transpose(1, 0, 2))


def gen_cp_volume(n, dim, periodic):
    cp = np.zeros((n[0], n[1], n[2], dim))
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
    cp += (np.random.rand(n[0], n[1], n[2], dim) - 0.5) * 10
    return np.floor(cp.transpose(2, 1, 0, 3))


def gen_controlpoints(n, dim, rational, periodic):
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
        w = np.random.rand(total_n) + 0.5
        w = np.round(w * 10) / 10
        cp = np.insert(cp, dim, w, 1)

    return cp


def write_basis(f, p, knot, periodic):
    for i in range(len(knot)):
        f.write("        basis%d = BSplineBasis(" % (i) + str(p[i]) + ", np." + repr(knot[i]))
        if periodic > -1 and i == 0:  # only consider periodicity in the first direction
            f.write("," + str(periodic))
        f.write(")\n")


def get_name(n, p, dim, rational, periodic):
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
        result += "_C%d_periodic" % periodic
    return result


def raise_order(p):
    result = "        crv2.raise_order(2)\n"
    result += "        p = crv2.order()\n"
    result += "        self.assertEqual(p, %d)\n" % (p + 2)
    return result


def write_object_creation(f, rational, pardim, clone=True):
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


def evaluate_curve():
    result = """
        u    = np.linspace(crv.start(0), crv.end(0), 13)
        pt   = crv(u)
        pt2  = crv2(u)
"""
    return result


def evaluate_surface():
    result = """
        u    = np.linspace(surf.start(0), surf.end(0), 9)
        v    = np.linspace(surf.start(1), surf.end(1), 9)
        pt   = surf(u,v)
        pt2  = surf2(u,v)

"""
    return result


def evaluate_volume():
    result = """
        u    = np.linspace(vol.start(0), vol.end(0), 7)
        v    = np.linspace(vol.start(1), vol.end(1), 7)
        w    = np.linspace(vol.start(2), vol.end(2), 7)
        pt   = vol(u,v,w)
        pt2  = vol2(u,v,w)
"""
    return result
