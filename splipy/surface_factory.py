# -*- coding: utf-8 -*-

"""Handy utilities for creating surfaces."""

from splipy import BSplineBasis, Curve, Surface
from math import pi, sqrt, atan2
from splipy.utils import flip_and_move_plane_geometry
import splipy.curve_factory as CurveFactory
import inspect
import numpy as np
import os
from os.path import dirname, realpath, join

__all__ = ['square', 'disc', 'sphere', 'extrude', 'revolve', 'cylinder', 'torus', 'edge_curves',
           'thicken', 'sweep', 'loft', 'interpolate', 'least_square_fit', 'teapot']


def square(size=1, lower_left=(0,0)):
    """square([size=1])

    Create a square with parametric origin at *(0,0)*.

    :param size: Size(s), either a single scalar or a tuple of scalars per axis
    :type size: float or (float)
    :return: A linear parametrized square
    :rtype: Surface
    """
    result = Surface()  # unit square
    result.scale(size)
    result += lower_left
    return result


def disc(r=1, center=(0,0,0), normal=(0,0,1), type='radial'):
    """disc([r=1], [type='radial'])

    Create a circular disc. The *type* parameter distinguishes between
    different parametrizations.

    :param float r: Radius
    :param string type: The type of parametrization ('radial' or 'square')
    :return: The disc
    :rtype: Surface
    """
    if type == 'radial':
        c1 = CurveFactory.circle(r)
        c2 = c1*0
        result = edge_curves(c2, c1)
        result.swap()
        result.reparam((0,r), (0,2*pi))
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
        result = Surface(basis1, basis2, cp, True)
    else:
        raise ValueError('invalid type argument')

    return flip_and_move_plane_geometry(result, center, normal)


def sphere(r=1, center=(0,0,0)):
    """sphere([r=1])

    Create a spherical shell.

    :param float r: Radius
    :param point-like center: Local origin of the sphere
    :return: The spherical shell
    :rtype: Surface
    """
    circle = CurveFactory.circle_segment(pi, r)
    circle.rotate(-pi / 2)
    circle.rotate(pi / 2, (1, 0, 0))  # flip up into xz-plane
    return revolve(circle) + center


def extrude(curve, amount):
    """Extrude a curve by sweeping it to a given height.

    :param Curve curve: Curve to extrude
    :param vector-like amount: 3-component vector of sweeping amount and
                               direction
    :return: The extruded curve
    :rtype: Surface
    """
    curve = curve.clone()  # clone input curve, throw away input reference
    curve.set_dimension(3)  # add z-components (if not already present)
    n = len(curve)  # number of control points of the curve
    cp = np.zeros((2 * n, curve.dimension + curve.rational))
    cp[:n, :] = curve.controlpoints  # the first control points form the bottom
    curve += amount
    cp[n:, :] = curve.controlpoints  # the last control points form the top
    return Surface(curve.bases[0], BSplineBasis(2), cp, curve.rational)


def revolve(curve, theta=2 * pi, axis=[0,0,1]):
    """revolve(curve, [theta=2pi], [axis=[0,0,1]])

    Revolve a surface by sweeping a curve in a rotational fashion around the
    *z* axis.

    :param Curve curve: Curve to revolve
    :param float theta: Angle to revolve, in radians
    :param vector-like axis: Axis of rotation
    :return: The revolved surface
    :rtype: Surface
    """
    curve = curve.clone()  # clone input curve, throw away input reference
    curve.set_dimension(3)  # add z-components (if not already present)
    curve.force_rational()  # add weight (if not already present)

    # align axis with the z-axis
    normal_theta = atan2(axis[1], axis[0])
    normal_phi   = atan2(sqrt(axis[0]**2 + axis[1]**2), axis[2])
    curve.rotate(-normal_theta, [0,0,1])
    curve.rotate(-normal_phi,   [0,1,0])

    circle_seg = CurveFactory.circle_segment(theta)

    n = len(curve)      # number of control points of the curve
    m = len(circle_seg) # number of control points of the sweep
    cp = np.zeros((m * n, 4))

    # loop around the circle and set control points by the traditional 9-point
    # circle curve with weights 1/sqrt(2), only here C0-periodic, so 8 points
    dt = 0
    t  = 0
    for i in range(m):
        x,y,w = circle_seg[i]
        dt  = atan2(y,x) - t
        t  += dt
        curve.rotate(dt)
        cp[i * n:(i + 1) * n, :]  = curve[:]
        cp[i * n:(i + 1) * n, 2] *= w
        cp[i * n:(i + 1) * n, 3] *= w
    result = Surface(curve.bases[0], circle_seg.bases[0], cp, True)
    # rotate it back again
    result.rotate(normal_phi,   [0,1,0])
    result.rotate(normal_theta, [0,0,1])
    return result


def cylinder(r=1, h=1, center=(0,0,0), axis=(0,0,1)):
    """cylinder([r=1], [h=1], [center=(0,0,0)], [axis=(0,0,1)])

    Create a cylinder shell with no top or bottom 

    :param float r: Radius
    :param float h: Height
    :param point-like center: The center of the bottom circle
    :param vector-like axis: Cylinder axis
    :return: The cylinder shell
    :rtype: Surface
    """
    return extrude(CurveFactory.circle(r, center, axis), h*axis)


def torus(minor_r=1, major_r=3, center=(0,0,0)):
    """torus([minor_r=1], [major_r=3])

    Create a torus (doughnut) by revolving a circle of size *minor_r* around
    the *z* axis with radius *major_r*.

    :param float minor_r: The thickness of the torus (radius in the *xz* plane)
    :param float major_r: The size of the torus (radius in the *xy* plane)
    :param point-like center: Local origin of the torus
    :return: A periodic torus
    :rtype: Surface
    """
    circle = CurveFactory.circle(minor_r)
    circle.rotate(pi / 2, (1, 0, 0))  # flip up into xz-plane
    circle.translate((major_r, 0, 0))  # move into position to spin around z-axis
    return revolve(circle) + center


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
        if rat:
            bottom = bottom.clone().force_rational() # don't mess with the input curve, make clone
            top.force_rational()                     # this is already a clone
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
    Greville points. It will use the same discretization as the input curve.
    It does not check for self-intersecting geometries.

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
    t = curve.bases[0].greville()
    if curve.dimension == 2:
        # linear parametrization across domain
        n = len(curve)
        left_points = np.matrix(np.zeros((n, 2)))
        right_points = np.matrix(np.zeros((n, 2)))
        linear = BSplineBasis(2)

        x = np.matrix(curve.evaluate(t))      # curve at interpolation points
        v = curve.derivative(t)               # velocity at interpolation points
        l = np.sqrt(v[:, 0]**2 + v[:, 1]**2)  # normalizing factor for velocity
        for i in range(n):
            if l[i] < 1e-13: # in case of zero velocity, use neighbour instead
                if i>0:
                    v[i,:] = v[i-1,:]
                else:
                    v[i,:] = v[i+1,:]
            else:
                v[i,:] /= l[i]
                
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
                    elif arg_names[j] == 'z':
                        argv[j] = 0.0
                    elif arg_names[j] == 't':
                        argv[j] = t[i]
                # figure out the distane at this particular point
                dist = amount(*argv)

                # store interpolation points
                right_points[i, 0] = x[i, 0] - v[i, 1] * dist  # x at bottom
                right_points[i, 1] = x[i, 1] + v[i, 0] * dist  # y at bottom
                left_points[ i, 0] = x[i, 0] + v[i, 1] * dist  # x at top
                left_points[ i, 1] = x[i, 1] - v[i, 0] * dist  # y at top
        else:
            right_points[:, 0] = x[:, 0] - v[:, 1] * amount  # x at bottom
            right_points[:, 1] = x[:, 1] + v[:, 0] * amount  # y at bottom
            left_points[ :, 0] = x[:, 0] + v[:, 1] * amount  # x at top
            left_points[ :, 1] = x[:, 1] - v[:, 0] * amount  # y at top
        # perform interpolation on each side
        right = CurveFactory.interpolate(right_points, curve.bases[0])
        left  = CurveFactory.interpolate(left_points,  curve.bases[0])
        return edge_curves(right, left)

    else:  # dimension=3, we will create a surrounding tube
        raise NotImplementedError('Currently only 2D supported. See comments in source code')

def sweep(path, shape):
    """Generate a surface by sweeping a shape along a path

    The resulting surface is an approximation generated by interpolating at the
    Greville points. It is generated by sweeping a shape curve along a path.

    The *shape* object has to be contained in the 'xy' plane (preferably centered
    around the origin) as its x-coordinate is extruded in the normal direction,
    and its y-coordinate in the binormal direction of the *path* curve.

    :param Curve path:  The path to drag *shape* along
    :param Curve shape: The shape to be dragged out to a surface
    :return: Surrounding surface
    :rtype: Surface
    """
    b1 = path.bases[0]
    b2 = shape.bases[0]
    n1 = b1.num_functions()
    n2 = b2.num_functions()
    # this requires binormals and normals, which only work in 3D, so assume this here
    X  = np.zeros((n1,n2, 3))
    for i in range(n1):
        u = b1.greville(i)
        x = path(u)
        B = path.binormal(u)
        N = path.normal(u)
        for j in range(n2):
            v = b2.greville(j)
            y = shape(v)
            X[i,j,:] = x + N*y[0] + B*y[1]

    return interpolate(X, [b1,b2])




def loft(*curves):
    if len(curves) == 1:
        curves = curves[0]

    # clone input, so we don't change those references
    # make sure everything has the same dimension since we need to compute length
    curves = [c.clone().set_dimension(3) for c in curves]
    if len(curves)==2:
        return edge_curves(curves)
    elif len(curves)==3:
        # can't do cubic spline interpolation, so we'll do quadratic
        basis2 = BSplineBasis(3)
        dist  = basis2.greville()
    else:
        x = [c.center() for c in curves]

        # create knot vector from the euclidian length between the curves
        dist = [0]
        for (x1,x0) in zip(x[1:],x[:-1]):
            dist.append(dist[-1] + np.linalg.norm(x1-x0))

        # using "free" boundary condition by setting N'''(u) continuous at second to last and second knot
        knot = [dist[0]]*4 + dist[2:-2] + [dist[-1]]*4
        basis2 = BSplineBasis(4, knot)

    n = len(curves)
    for i in range(n):
        for j in range(i+1,n):
            Curve.make_splines_identical(curves[i], curves[j])
    
    basis1 = curves[0].bases[0]
    m      = basis1.num_functions()
    u      = basis1.greville() # parametric interpolation points
    v      = dist              # parametric interpolation points
    
    # compute matrices
    Nu     = basis1(u)
    Nv     = basis2(v)
    Nu_inv = np.linalg.inv(Nu)
    Nv_inv = np.linalg.inv(Nv)

    # compute interpolation points in physical space
    x      = np.zeros((m,n, curves[0][0].size))
    for i in range(n):
        x[:,i,:] = Nu * curves[i].controlpoints

    # solve interpolation problem
    cp = np.tensordot(Nv_inv, x,  axes=(1,1))
    cp = np.tensordot(Nu_inv, cp, axes=(1,1))

    # re-order controlpoints so they match up with Surface constructor
    cp = cp.transpose((1, 0, 2))
    cp = cp.reshape(n*m, cp.shape[2])

    return Surface(basis1, basis2, cp, curves[0].rational)


def interpolate(x, bases, u=None):
    """interpolate(x, bases, [u=None])

    Interpolate a surface on a set of regular gridded interpolation points `x`.

    The points can be either a matrix (in which case the first index is
    interpreted as a flat row-first index of the interpolation grid) or a 3D
    tensor. In both cases the last index is the physical coordinates.

    :param x: Grid of interpolation points
    :type x: matrix-like or 3D-tensor-like
    :param [BSplineBasis] bases: The basis to interpolate on
    :param [array-like] u: Parametric interpolation points, defaults to
        Greville points of the basis
    :return: Interpolated surface
    :rtype: Surface
    """
    surf_shape = [b.num_functions() for b in bases]
    dim = x.shape[-1]
    if len(x.shape) == 2:
        x = x.reshape(surf_shape + [dim])
    if u is None:
        u = [b.greville() for b in bases]
    N_all = [b(t) for b,t in zip(bases, u)]
    N_all.reverse()
    cp = x
    for N in N_all:
        cp = np.tensordot(np.linalg.inv(N), cp, axes=(1,1))

    return Surface(bases[0], bases[1], cp.transpose(1,0,2).reshape((np.prod(surf_shape),dim)))


def least_square_fit(x, bases, u):
    """Perform a least-square fit of a point cloud `x` onto a spline basis.

    The points can be either a matrix (in which case the first index is
    interpreted as a flat row-first index of the interpolation grid) or a 3D
    tensor. In both cases the last index is the physical coordinates.

    There must be at least as many points as basis functions.

    :param x: Grid of evaluation points
    :type x: matrix-like or 3D-tensor-like
    :param [BSplineBasis] bases: Basis on which to interpolate
    :param [array-like] u: Parametric values at evaluation points
    :return: Approximated surface
    :rtype: Surface
    """
    surf_shape = [b.num_functions() for b in bases]
    dim = x.shape[-1]
    if len(x.shape) == 2:
        x = x.reshape(surf_shape + [dim])
    N_all = [b(t) for b,t in zip(bases, u)]
    N_all.reverse()
    cp = x
    for N in N_all:
        cp = np.tensordot(N.T, cp, axes=(1,1))
    for N in N_all:
        cp = np.tensordot(np.linalg.inv(N.T*N), cp, axes=(1,1))

    return Surface(bases[0], bases[1], cp.transpose(1,0,2).reshape((np.prod(surf_shape),dim)))


def teapot():
    """Generate the Utah teapot as 32 cubic bezier patches. This teapot has a
    rim, but no bottom. It is also self-intersecting making it unsuitable for
    perfect-match multipatch modeling.

    The data is picked from http://www.holmes3d.net/graphics/teapot/

    :return: The utah teapot
    :rtype: List of Surface
    """
    path = join(dirname(realpath(__file__)), 'templates', 'teapot.bpt')
    with open(path) as f:
        results = []
        numb_patches = int(f.readline())
        for i in range(numb_patches):
            p = np.fromstring(f.readline(), dtype=np.uint8, count=2, sep=' ')
            basis1 = BSplineBasis(p[0]+1)
            basis2 = BSplineBasis(p[1]+1)

            ncp = basis1.num_functions() * basis2.num_functions()
            cp  = [np.fromstring(f.readline(), dtype=np.float, count=3, sep=' ') for j in range(ncp)]
            results.append(Surface(basis1, basis2, cp))

    return results
