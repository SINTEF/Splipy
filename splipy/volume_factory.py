# -*- coding: utf-8 -*-

"""Handy utilities for creating volumes."""

from math import pi, sqrt
import numpy as np
from splipy import Surface, Volume, BSplineBasis
import splipy.surface_factory as SurfaceFactory

__all__ = ['cube', 'revolve', 'cylinder', 'extrude', 'edge_surfaces', 'loft',
           'interpolate', 'least_square_fit']


def cube(size=1, lower_left=(0,0,0)):
    """cube([size=1])

    Create a cube with parmetric origin at *(0,0,0)*.

    :param size: Size(s), either a single scalar or a tuple of scalars per axis
    :type size: float or (float)
    :return: A linear parametrized box
    :rtype: Volume
    """
    result = Volume()
    result.scale(size)
    result += lower_left
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
    basis = BSplineBasis(3, [-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5], periodic=0)
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


def cylinder(r=1, h=1, center=(0,0,0), axis=(0,0,1)):
    """cylinder([r=1], [h=1], [center=(0,0,0)], [axis=(0,0,1)])

    Create a solid cylinder 

    :param float r: Radius
    :param float h: Height
    :param point-like center: The center of the bottom circle
    :param vector-like axis: Cylinder axis
    :return: The cylinder
    :rtype: Volume
    """
    return extrude(SurfaceFactory.disc(r, center, axis), h*axis)


def extrude(surf, amount):
    """Extrude a surface by sweeping it to a given height.

    :param Surface surf: Surface to extrude
    :param vector-like amount: 3-component vector of sweeping amount and
                               direction
    :return: The extruded surface
    :rtype: Volume
    """
    surf.set_dimension(3)  # add z-components (if not already present)
    cp = []
    for controlpoint in surf:
        cp.append(list(controlpoint))
    surf += amount
    for controlpoint in surf:
        cp.append(list(controlpoint))
    surf -= amount
    return Volume(surf.bases[0], surf.bases[1], BSplineBasis(2), cp, surf.rational)


def edge_surfaces(*surfaces):
    """edge_surfaces(surfaces...)

    Create the volume defined by the region between the input surfaces.

    In case of six input surfaces, these must be given in the order: bottom,
    top, left, right, back, front. Opposing sides must be parametrized in the
    same directions.

    :param [Surface] surfaces: Two or six edge surfaces
    :return: The enclosed volume
    :rtype: Volume
    :raises ValueError: If the length of *surfaces* is not two or six
    """
    if len(surfaces) == 1: # probably gives input as a list-like single variable
        surfaces = surfaces[0]
    if len(surfaces) == 2:
        surf1 = surfaces[0].clone()
        surf2 = surfaces[1].clone()
        Surface.make_splines_identical(surf1, surf2)
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
        # coons patch (https://en.wikipedia.org/wiki/Coons_patch)
        umin = surfaces[0]
        umax = surfaces[1]
        vmin = surfaces[2]
        vmax = surfaces[3]
        wmin = surfaces[4]
        wmax = surfaces[5]
        vol1 = edge_surfaces(umin,umax)
        vol2 = edge_surfaces(vmin,vmax)
        vol3 = edge_surfaces(wmin,wmax)
        vol4 = Volume(controlpoints=vol1.corners(order='F'), rational=vol1.rational)
        vol1.swap(0, 2)
        vol1.swap(1, 2)
        vol2.swap(1, 2)
        vol4.swap(1, 2)
        Volume.make_splines_identical(vol1, vol2)
        Volume.make_splines_identical(vol1, vol3)
        Volume.make_splines_identical(vol1, vol4)
        Volume.make_splines_identical(vol2, vol3)
        Volume.make_splines_identical(vol2, vol4)
        Volume.make_splines_identical(vol3, vol4)
        result  = vol1.clone()
        result.controlpoints +=   vol2.controlpoints
        result.controlpoints +=   vol3.controlpoints
        result.controlpoints -= 2*vol4.controlpoints
        return result
    else:
        raise ValueError('Requires two or six input surfaces')

def loft(*surfaces):
    if len(surfaces) == 1:
        surfaces = surfaces[0]

    # clone input, so we don't change those references
    # make sure everything has the same dimension since we need to compute length
    surfaces = [s.clone().set_dimension(3) for s in surfaces]
    if len(surfaces)==2:
        return SurfaceFactory.edge_curves(surfaces)
    elif len(surfaces)==3:
        # can't do cubic spline interpolation, so we'll do quadratic
        basis3 = BSplineBasis(3)
        dist  = basis3.greville()
    else:
        x = [s.center() for s in surfaces]

        # create knot vector from the euclidian length between the surfaces
        dist = [0]
        for (x1,x0) in zip(x[1:],x[:-1]):
            dist.append(dist[-1] + np.linalg.norm(x1-x0))

        # using "free" boundary condition by setting N'''(u) continuous at second to last and second knot
        knot = [dist[0]]*4 + dist[2:-2] + [dist[-1]]*4
        basis3 = BSplineBasis(4, knot)

    n = len(surfaces)
    for i in range(n):
        for j in range(i+1,n):
            Surface.make_splines_identical(surfaces[i], surfaces[j])

    basis1 = surfaces[0].bases[0]
    basis2 = surfaces[0].bases[1]
    m1     = basis1.num_functions()
    m2     = basis2.num_functions()
    dim    = len(surfaces[0][0])
    u      = basis1.greville() # parametric interpolation points
    v      = basis2.greville()
    w      = dist

    # compute matrices
    Nu     = basis1(u)
    Nv     = basis2(v)
    Nw     = basis3(w)
    Nu_inv = np.linalg.inv(Nu)
    Nv_inv = np.linalg.inv(Nv)
    Nw_inv = np.linalg.inv(Nw)

    # compute interpolation points in physical space
    x      = np.zeros((m1,m2,n, dim))
    for i in range(n):
        tmp        = np.tensordot(Nv, surfaces[i].controlpoints, axes=(1,1))
        x[:,:,i,:] = np.tensordot(Nu, tmp                      , axes=(1,1))

    # solve interpolation problem
    cp = np.tensordot(Nw_inv, x,  axes=(1,2))
    cp = np.tensordot(Nv_inv, cp, axes=(1,2))
    cp = np.tensordot(Nu_inv, cp, axes=(1,2))

    # re-order controlpoints so they match up with Surface constructor
    cp = np.reshape(cp.transpose((2, 1, 0, 3)), (m1*m2*n, dim))

    return Volume(basis1, basis2, basis3, cp, surfaces[0].rational)

def interpolate(x, bases, u=None):
    """Interpolate a volume on a set of regular gridded interpolation points

    :param x: Grid of interpolation points. Component *x[i,j,k,l]* denotes
        component *x[l]* at index *(i,j)* in the volume. For matrix input *x[i,j]*
        then index *i* is interpreted as a flat row-first indexing of the 
        interpolation grid with components *x[j]*
    :type x: matrix-like or 4D-tensor-like
    :param [BSplineBasis] bases: The basis to interpolate on
    :param [array-like]   u    : Parametric interpolation points, defaults to
                                 greville points of the basis
    :return: Interpolated volume
    :rtype: Volume
    """
    vol_shape = [b.num_functions() for b in bases]
    dim = x.shape[-1]
    if len(x.shape) == 2:
        x = x.reshape(vol_shape + [dim])
    if u is None:
        u = [b.greville() for b in bases]
    N_all = [b(t) for b,t in zip(bases, u)]
    N_all.reverse()
    cp = x
    for N in N_all:
        cp = np.tensordot(np.linalg.inv(N), cp, axes=(1,2))

    return Volume(bases[0], bases[1], bases[2], cp.transpose(2,1,0,3).reshape((np.prod(vol_shape),dim)))

def least_square_fit(x, bases, u):
    """Perform a least-square fit of a point cloud onto a spline basis

    :param x: Grid of evaluation points. Component *x[i,j,k]* denotes
        component *x[k]* at index *(i,j)* in the volume. For matrix input *x[i,j]*
        then index *i* is interpreted as a flat row-first indexing of the 
        interpolation grid with components *x[j]*. The number of points must be
        equal to or larger than the number of bases functions
    :type x: matrix-like or 3D-tensor-like
    :param [BSplineBasis] bases: Basis on which to interpolate
    :param [array-like]   t    : parametric values at evaluation points
    :return: Approximated volume
    :rtype: Volume
    """
    vol_shape = [b.num_functions() for b in bases]
    dim = x.shape[-1]
    if len(x.shape) == 2:
        x = x.reshape(vol_shape + [dim])
    N_all = [b(t) for b,t in zip(bases, u)]
    N_all.reverse()
    cp = x
    for N in N_all:
        cp = np.tensordot(N.T, cp, axes=(1,2))
    for N in N_all:
        cp = np.tensordot(np.linalg.inv(N.T*N), cp, axes=(1,2))

    return Volume(bases[0], bases[1], bases[2], cp.transpose(2,1,0,3).reshape((np.prod(vol_shape),dim)))
