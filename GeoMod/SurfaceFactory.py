# -*- coding: utf-8 -*-

"""Handy utilities for creating surfaces."""

from GeoMod import BSplineBasis, Curve, Surface
from math import pi, sqrt
import GeoMod.CurveFactory as CurveFactory
import inspect
import numpy as np

__all__ = ['square', 'disc', 'sphere', 'extrude', 'revolve', 'cylinder', 'torus', 'edge_curves',
           'thicken']


def square(size=1):
    """square([size=1])

    Create a square with parametric origin at *(0,0)*.

    :param size: Size(s), either a single scalar or a tuple of scalars per axis
    :type size: float or (float)
    :return: A linear parametrized square
    :rtype: Surface
    """
    result = Surface()  # unit square
    result.scale(size)
    return result


def disc(r=1, type='radial'):
    """disc([r=1], [type='radial'])

    Create a circular disc with center at (0,0). The *type* parameter
    distinguishes between different parametrizations.

    :param float r: Radius
    :param string type: The type of parametrization ('radial' or 'square')
    :return: The disc
    :rtype: Surface
    """
    if type == 'radial':
        c = CurveFactory.circle(r)
        cp = np.zeros((16, 3))
        cp[:, -1] = 1
        cp[1::2, :] = c.controlpoints
        return Surface(BSplineBasis(), c.bases[0], cp, True)
    elif type == 'square':
        w = 1 / sqrt(2)
        cp = [[-r * w, -r * w, 1],
              [0, -r, w],
              [r * w, -r * w, 1],
              [-r, 0, w],
              [0, 0, 1],
              [r, 0, w],
              [-r * w, r * w, 1],
              [0, r, w],
              [r * w, r * w, 1]]
        basis1 = BSplineBasis(3)
        basis2 = BSplineBasis(3)
        return Surface(basis1, basis2, cp, True)
    else:
        raise ValueError('invalid type argument')


def sphere(r=1):
    """sphere([r=1])

    Create a spherical shell.

    :param float r: Radius
    :return: The spherical shell
    :rtype: Surface
    """
    circle = CurveFactory.circle_segment(pi, r)
    circle.rotate(-pi / 2)
    circle.rotate(pi / 2, (1, 0, 0))  # flip up into xz-plane
    return revolve(circle)


def extrude(curve, h):
    """Extrude a curve by sweeping it in the *z* direction to a given height.

    :param Curve curve: Curve to extrude
    :param float h: Height in the *z* direction
    :return: The extruded curve
    :rtype: Surface
    """
    curve = curve.clone()  # clone input curve, throw away input reference
    curve.set_dimension(3)  # add z-components (if not already present)
    n = len(curve)  # number of control points of the curve
    cp = np.zeros((2 * n, 4))
    cp[:n, :] = curve.controlpoints  # the first control points form the bottom
    curve += (0, 0, h)
    cp[n:, :] = curve.controlpoints  # the last control points form the top
    return Surface(curve.bases[0], BSplineBasis(2), cp, curve.rational)


def revolve(curve, theta=2 * pi):
    """revolve(curve, [theta=2pi])

    Revolve a surface by sweeping a curve in a rotational fashion around the
    *z* axis.

    :param Curve curve: Curve to revolve
    :param float theta: Angle to revolve, in radians
    :return: The revolved curve
    :rtype: Surface
    """
    curve = curve.clone()  # clone input curve, throw away input reference
    curve.set_dimension(3)  # add z-components (if not already present)
    curve.force_rational()  # add weight (if not already present)
    n = len(curve)  # number of control points of the curve
    cp = np.zeros((8 * n, 4))
    basis = BSplineBasis(3, [-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5], periodic=0)
    basis *= 2 * pi / 4  # set parametric domain to (0,2pi) in v-direction

    # loop around the circle and set control points by the traditional 9-point
    # circle curve with weights 1/sqrt(2), only here C0-periodic, so 8 points
    for i in range(8):
        if i % 2 == 0:
            weight = 1.0
        else:
            weight = 1.0 / sqrt(2)
        cp[i * n:(i + 1) * n, :] = curve.controlpoints
        cp[i * n:(i + 1) * n, 2] *= weight
        cp[i * n:(i + 1) * n, 3] *= weight
        curve.rotate(pi / 4)
    return Surface(curve.bases[0], basis, cp, True)


def cylinder(r=1, h=1):
    """cylinder([r=1], [h=1])

    Create a cylinder shell with no top or bottom with the *z* axis as central axis.

    :param float r: Radius
    :param float h: Height
    :return: The cylinder shell
    :rtype: Surface
    """
    return extrude(CurveFactory.circle(r), h)


def torus(minor_r=1, major_r=3):
    """torus([minor_r=1], [major_r=3])

    Create a torus (doughnut) by revolving a circle of size *minor_r* around
    the *z* axis with radius *major_r*.

    :param float minor_r: The thickness of the torus (radius in the *xz* plane)
    :param float major_r: The size of the torus (radius in the *xy* plane)
    :return: A periodic torus
    :rtype: Surface
    """
    circle = CurveFactory.circle(minor_r)
    circle.rotate(pi / 2, (1, 0, 0))  # flip up into xz-plane
    circle.translate((major_r, 0, 0))  # move into position to spin around z-axis
    return revolve(circle)


def edge_curves(*curves):
    """edge_curves(curves...)

    Create the surface defined by the region between the input curves.

    In case of four input curves, these must be given in an ordered directional
    closed loop around the resulting surface.

    :param [Curve] curves: Two or four edge curves
    :return: The enclosed surface
    :rtype: Surface
    :raises ValueError: If the length of *curves* is not two or four
    """
    if len(curves) == 1: # probably gives input as a list-like single variable
        curves = curves[0]
    if len(curves) == 2:
        crv1 = curves[0].clone()
        crv2 = curves[1].clone()
        Curve.make_splines_identical(crv1, crv2)
        (n, d) = crv1.controlpoints.shape  # d = dimension + rational

        controlpoints = np.zeros((2 * n, d))
        controlpoints[:n, :] = crv1.controlpoints
        controlpoints[n:, :] = crv2.controlpoints
        linear = BSplineBasis(2)

        return Surface(crv1.bases[0], linear, controlpoints, crv1.rational)
    elif len(curves) == 4:
        # coons patch (https://en.wikipedia.org/wiki/Coons_patch)
        bottom = curves[0]
        right  = curves[1]
        top    = curves[2].clone()
        left   = curves[3].clone()  # gonna change these two, so make copies
        top.reverse()
        left.reverse()
        # create linear interpolation between opposing sides
        s1 = edge_curves(bottom, top)
        s2 = edge_curves(left, right)
        s2.swap()
        # create (linear,linear) corner parametrization
        linear = BSplineBasis(2)
        rat = s1.rational  # using control-points from top/bottom, so need to know if these are rational
        s3 = Surface(linear, linear, [bottom[0], bottom[-1], top[0], top[-1]], rat)

        # in order to add spline surfaces, they need identical parametrization
        Surface.make_splines_identical(s1, s2)
        Surface.make_splines_identical(s1, s3)
        Surface.make_splines_identical(s2, s3)

        result = s1
        result.controlpoints += s2.controlpoints
        result.controlpoints -= s3.controlpoints
        return result
    else:
        raise ValueError('Requires two or four input curves')


def thicken(curve, amount):
    """Generate a surface by adding thickness to a curve.

    - For 2D curves this will generate a 2D planar surface with the curve
      through the center.

    - For 3D curves this will generate a surface "tube" which is open at both
      ends (that is, for a straight line it will generate a cylinder shell).

    The resulting surface is an approximation generated by interpolating at the
    Greville points. It does not use the same discretization as the input
    curve. It does not check for self-intersecting geometries.

    :param Curve curve: The generating curve
    :param amount: The amount of thickness, either constant or variable (if
        variable, the function must accept parameters named *x*, *y*, *z* and/or *t*)
    :return: Surrounding surface
    :rtype: Surface
    """
    # NOTES: There are several pitfalls with this function
    #  * self intersection:
    #     could be handled more gracefully, but is here ignored
    #  * choice of discretization:
    #     the offset curve is computed from the velocity (tangent) which is of
    #     one continuity less than the original curve. In particular C1
    #     quadratic curves will get very noticable C0-kinks in them. Currently
    #     this is completely ignored and we keep the high continuity of the
    #     original curve.
    #  * width given by function input
    #     could produce wild behaviour. Original discretization might not
    #     produce a satisfactory result
    #  * 3D tube geometries:
    #     unimplemented as of now. Would like to see the three points above
    #     resolved before this is implemented. Rough idea is to compute the
    #     acceleration and binormal vectors to the curve and sketch out a
    #     circle in the plane defined by these two vectors

    curve = curve.clone()  # clone input curve, throw away input reference
    t = curve.basis.greville()
    if curve.dimension == 2:
        # linear parametrization across domain
        n = len(curve)
        left_points = np.matrix(np.zeros((n, 2)))
        right_points = np.matrix(np.zeros((n, 2)))
        linear = BSplineBasis(2)

        x = curve.evaluate(t)  # curve at interpolation points
        v = curve.evaluate_tangent(t)  # velocity at interpolation points
        l = np.sqrt(v[:, 0]**2 + v[:, 1]**2)  # normalizing factor for velocity
        v[:, 0] = v[:, 0] / l
        v[:, 1] = v[:, 1] / l
        v = np.matrix(v)
        if inspect.isfunction(amount):
            arg_names = inspect.getargspec(amount).args
            argc = len(arg_names)
            argv = [0] * argc
            for i in range(n):
                # build up the list of arguments (in case not all of (x,y,t) are specified)
                for j in range(argc):
                    if arg_names[j] == 'x':
                        argv[j] = x[i, 0]
                    elif arg_names[j] == 'y':
                        argv[j] = x[i, 1]
                    elif arg_names[j] == 't':
                        argv[j] = t[i]
                # figure out the distane at this particular point
                dist = amount(*argv)

                # store interpolation points
                right_points[i, 0] = x[i, 0] - v[i, 1] * dist  # x at bottom
                right_points[i, 1] = x[i, 1] + v[i, 0] * dist  # y at bottom
                left_points[i, 0] = x[i, 0] + v[i, 1] * dist  # x at top
                left_points[i, 1] = x[i, 1] - v[i, 0] * dist  # y at top
        else:
            right_points[:, 0] = x[:, 0] - v[:, 1] * amount  # x at bottom
            right_points[:, 1] = x[:, 1] + v[:, 0] * amount  # y at bottom
            left_points[:, 0] = x[:, 0] + v[:, 1] * amount  # x at top
            left_points[:, 1] = x[:, 1] - v[:, 0] * amount  # y at top
        # perform interpolation on each side
        right = CurveFactory.interpolate(right_points, curve.basis)
        left = CurveFactory.interpolate(left_points, curve.basis)

        # draw the linear surface in between
        controlpoints = np.zeros((2 * n, 2))
        controlpoints[:n, :] = right.controlpoints
        controlpoints[n:, :] = left.controlpoints
        return Surface(curve.basis, linear, controlpoints)

    else:  # dimension=3, we will create a surrounding tube
        raise NotImplementedError('Currently only 2D supported. See comments in source code')

# def loft(curves):
    
