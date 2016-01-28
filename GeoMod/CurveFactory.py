# -*- coding: utf-8 -*-

"""Handy utilities for creating curves."""

from math import pi, cos, sin, sqrt, ceil
from GeoMod import Curve, BSplineBasis
import numpy as np

__all__ = ['line', 'polygon', 'n_gon', 'circle', 'circle_segment', 'interpolate']


def line(a, b):
    """Create a line between two points.

    :param point-like a: Start point
    :param point-like b: End point
    :return: Linear curve from *a* to *b*
    :rtype: Curve
    """
    return Curve(controlpoints=[a, b])


def polygon(points):
    """Create a linear interpolation between input points.

    :param [point-like] points: The points to interpolate
    :return: Linear curve through the input points
    :rtype: Curve
    """
    # establish knot vector based on eucledian length between points
    knot = [0, 0]
    prevPt = points[0]
    dist = 0
    for pt in points[1:]:
        for (x0, x1) in zip(prevPt, pt):  # loop over (x,y) and maybe z-coordinate
            dist += (x1 - x0)**2
        knot.append(sqrt(dist))
        prevPt = pt
    knot.append(knot[-1])

    return Curve(BSplineBasis(2, knot), points)


def n_gon(n=5, r=1):
    """n_gon([n=5], [r=1])

    Create a regular polygon of *n* equal sides centered at the origin.

    :param int n: Number of sides and vertices
    :param float r: Radius
    :return: A linear, periodic, 2D curve
    :rtype: Curve
    :raises ValueError: If radius is not positive
    :raises ValueError: If *n* < 3
    """
    if r <= 0:
        raise ValueError('radius needs to be positive')
    if n < 3:
        raise ValueError('regular polygons need at least 3 sides')

    cp = []
    dt = 2 * pi / n
    knot = [0]
    for i in range(n):
        cp.append([r * cos(i * dt), r * sin(i * dt)])
        knot.append(i)
    knot += [n, n]
    basis = BSplineBasis(2, knot, 0)
    return Curve(basis, cp)


def circle(r=1):
    """circle([r=1])

    Create a circle at the origin.

    :param float r: Radius
    :return: A periodic, quadratic rational curve
    :rtype: Curve
    :raises ValueError: If radius is not positive
    """
    if r <= 0:
        raise ValueError('radius needs to be positive')

    w = 1.0 / sqrt(2)
    controlpoints = [[r, 0, 1],
                     [r*w, r*w, w],
                     [0, r, 1],
                     [-r*w, r*w, w],
                     [-r, 0, 1],
                     [-r*w, -r*w, w],
                     [0, -r, 1],
                     [r*w, -r*w, w]]
    knot = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]) / 4.0 * 2 * pi
    return Curve(BSplineBasis(3, knot, 0), controlpoints, True)


def circle_segment(theta, r=1):
    """circle_segment(theta, [r=1])

    Create a circle segment centered at the origin, starting at *(r,0)*.

    :param float theta: Angle in radians
    :param float r: Radius
    :return: A quadratic rational curve
    :rtype: Curve
    :raises ValueError: If radiusis not positive
    :raises ValueError: If theta is not in the range *[-2pi, 2pi]*
    """
    # error test input
    if abs(theta) > 2 * pi:
        raise ValueError('theta needs to be in range [-2pi,2pi]')
    if r <= 0:
        raise ValueError('radius needs to be positive')

    # build knot vector
    knot_spans = int(ceil(theta / (2 * pi / 3)))
    knot = [0]
    for i in range(knot_spans + 1):
        knot += [i] * 2
    knot += [knot_spans]  # knot vector [0,0,0,1,1,2,2,..,n,n,n]
    knot = np.array(knot) / float(knot[-1]) * theta  # set parametic space to [0,theta]

    n = (knot_spans - 1) * 2 + 3  # number of control points needed
    cp = []
    t = 0  # current angle
    dt = float(theta) / knot_spans / 2  # angle step

    # build control points
    for i in range(n):
        w = 1 - (i % 2) * (1 - cos(dt))  # weights = 1 and cos(dt) every other i
        x = r * cos(t)
        y = r * sin(t)
        cp += [[x, y, w]]
        t += dt

    return Curve(BSplineBasis(3, knot), cp, True)


def interpolate(x_pts, basis):
    """Perform general spline interpolation at the Greville points of a basis.

    :param array-like x_pts: Matrix *X[i,j]* of interpolation points *xj* with
        components *j*
    :param BSplineBasis basis: Basis on which to interpolate
    :return: Interpolated curve
    :rtype: Curve
    """
    # wrap x_pts into a numpy matrix
    x_pts = np.matrix(x_pts)

    # evaluate all basis functions in the interpolation points
    grev_pts = basis.greville()
    N = basis.evaluate(grev_pts)

    # solve interpolation problem
    controlpoints = np.linalg.solve(N, x_pts)

    return Curve(basis, controlpoints)
