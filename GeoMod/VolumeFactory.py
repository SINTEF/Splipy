# -*- coding: utf-8 -*-

"""Handy utilities for creating volumes."""

from math import pi, sqrt
import numpy as np
from GeoMod import Surface, Volume, BSplineBasis
import GeoMod.SurfaceFactory as SurfaceFactory

__all__ = ['cube', 'revolve', 'cylinder', 'extrude', 'edge_surfaces']


def cube(size=1):
    """cube([size=1])

    Create a cube with parmetric origin at *(0,0,0)*.

    :param size: Size(s), either a single scalar or a tuple of scalars per axis
    :type size: float or (float)
    :return: A linear parametrized box
    :rtype: Volume
    """
    result = Volume()
    result.scale(size)
    return result


def revolve(surf, theta=2 * pi):
    """revolve(surf, [theta=2pi])

    Revolve a volume by sweeping a surface in a rotational fashion around the
    *z* axis.

    :param Surface surf: Surface to revolve
    :param float theta: Angle to revolve, in radians
    :return: The revolved surface
    :rtype: Volume
    """
    surf = surf.clone()  # clone input surface, throw away old reference
    surf.set_dimension(3)  # add z-components (if not already present)
    surf.force_rational()  # add weight (if not already present)
    n = len(surf)  # number of control points of the surface
    cp = np.zeros((8 * n, 4))
    basis = BSplineBasis(3, [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4], periodic=0)
    basis *= 2 * pi / 4  # set parametric domain to (0,2pi) in w-direction

    # loop around the circle and set control points by the traditional 9-point
    # circle curve with weights 1/sqrt(2), only here C0-periodic, so 8 points
    for i in range(8):
        if i % 2 == 0:
            weight = 1.0
        else:
            weight = 1.0 / sqrt(2)
        cp[i * n:(i + 1) * n, :] = np.reshape(surf.controlpoints.transpose(1, 0, 2), (n, 4))
        cp[i * n:(i + 1) * n, 2] *= weight
        cp[i * n:(i + 1) * n, 3] *= weight
        surf.rotate(pi / 4)
    return Volume(surf.bases[0], surf.bases[1], basis, cp, True)


def cylinder(r=1, h=1):
    """cylinder([r=1], [h=1])

    Create a solid cylinder with the *z* axis as central axis.

    :param float r: Radius
    :param float h: Height
    :return: The cylinder
    :rtype: Volume
    """
    shell = SurfaceFactory.cylinder(r, h)
    cp = []
    for controlpoint in shell:
        cp.append([0, 0, controlpoint[2], controlpoint[3]])  # project to z-axis
    for controlpoint in shell:
        cp.append(list(controlpoint))

    return Volume(shell.bases[0], shell.bases[1], BSplineBasis(), cp, True)


def extrude(surf, h):
    """Extrude a surface by sweeping it in the *z* direction to a given height.

    :param Surface surf: Surface to extrude
    :param float h: Height in the *z* direction
    :return: The extruded surface
    :rtype: Volume
    """
    surf.set_dimension(3)  # add z-components (if not already present)
    cp = []
    for controlpoint in surf:
        cp.append(list(controlpoint))
    surf += (0, 0, h)
    for controlpoint in surf:
        cp.append(list(controlpoint))
    surf -= (0, 0, h)
    return Volume(surf.bases[0], surf.bases[1], BSplineBasis(2), cp, surf.rational)


def edge_surfaces(surfaces):
    """Create the volume defined by the region between the input surfaces.

    In case of six input surfaces, these must be given in the order: bottom,
    top, left, right, back, front. Opposing sides must be parametrized in the
    same directions.

    :param [Surface] surfaces: Two or six edge surfaces
    :return: The enclosed volume
    :rtype: Volume
    :raises ValueError: If the length of *surfaces* is not two or six
    """
    if len(surfaces) == 2:
        surf1 = surfaces[0].clone()
        surf2 = surfaces[1].clone()
        Surface.make_surfaces_identical(surf1, surf2)
        (n1, n2, d) = surf1.controlpoints.shape  # d = dimension + rational

        controlpoints = np.zeros((n1, n2, 2, d))
        controlpoints[:, :, 0, :] = surf1.controlpoints
        controlpoints[:, :, 1, :] = surf2.controlpoints

        # Volume constructor orders control points in a different way, so we
        # create it from scratch here
        result = Volume()
        result.bases = [surf1.bases[0], surf1.bases[1], BSplineBasis(2)]
        result.dimension = surf1.dimension
        result.rational = surf1.rational
        result.controlpoints = controlpoints

        return result
    elif len(surfaces) == 6:
        raise NotImplementedError('Should have been coons patch algorithm here. Come back later')
    else:
        raise ValueError('Requires two or six input surfaces')
