# -*- coding: utf-8 -*-

"""Handy utilities for creating volumes."""

from math import pi, sqrt
import numpy as np
from splipy import Surface, Volume, BSplineBasis
import splipy.curve_factory as CurveFactory
import splipy.surface_factory as SurfaceFactory

__all__ = ['cube', 'sphere', 'revolve', 'cylinder', 'extrude', 'edge_surfaces',
           'loft', 'interpolate', 'least_square_fit']



def cube(size=1, lower_left=(0,0,0)):
    """  Create a cube with parmetric origin at *(0,0,0)*.

    :param float size: Size(s), either a single scalar or a tuple of scalars per axis
    :param array-like lower_left: local origin, the lower bottom left corner of the cube
    :return: A linear parametrized box
    :rtype: Volume
    """
    result = Volume()
    result.scale(size)
    result += lower_left
    return result

def sphere(r=1, center=(0,0,0), type='radial'):
    """  Create a solid sphere

    :param float r: Radius
    :param array-like center: Local origin of the sphere
    :param string type: The type of parametrization ('radial' or 'square')
    :return: A solid ball
    :rtype: Volume
    """
    if type == 'radial':
        shell    = SurfaceFactory.sphere(r, center)
        midpoint = shell*0 + center
        return edge_surfaces(shell, midpoint)
    elif type == 'square':
        # based on the work of James E.Cobb: "Tiling the Sphere with Rational Bezier Patches"
        # University of Utah, July 11, 1988. UUCS-88-009
        b = BSplineBasis(order=5)
        sr2 = sqrt(2)
        sr3 = sqrt(3)
        sr6 = sqrt(6)
        cp = [[      -4*(sr3-1),       4*(1-sr3),      4*(1-sr3),   4*(3-sr3)  ], # row 0
              [           -sr2 ,     sr2*(sr3-4),    sr2*(sr3-4), sr2*(3*sr3-2)],
              [              0 ,  4./3*(1-2*sr3), 4./3*(1-2*sr3),4./3*(5-sr3)  ],
              [            sr2 ,     sr2*(sr3-4),    sr2*(sr3-4), sr2*(3*sr3-2)],
              [       4*(sr3-1),       4*(1-sr3),      4*(1-sr3),   4*(3-sr3)  ],
              [    -sr2*(4-sr3),            -sr2,    sr2*(sr3-4), sr2*(3*sr3-2)], # row 1
              [    -(3*sr3-2)/2,     (2-3*sr3)/2,     -(sr3+6)/2,     (sr3+6)/2],
              [              0 , sr2*(2*sr3-7)/3,       -5*sr6/3, sr2*(sr3+6)/3],
              [     (3*sr3-2)/2,     (2-3*sr3)/2,     -(sr3+6)/2,     (sr3+6)/2],
              [     sr2*(4-sr3),            -sr2,    sr2*(sr3-4), sr2*(3*sr3-2)],
              [ -4./3*(2*sr3-1),               0, 4./3*(1-2*sr3),   4*(5-sr3)/3], # row 2
              [-sr2/3*(7-2*sr3),               0,       -5*sr6/3, sr2*(sr3+6)/3],
              [              0 ,               0,    4*(sr3-5)/3, 4*(5*sr3-1)/9],
              [ sr2/3*(7-2*sr3),               0,       -5*sr6/3, sr2*(sr3+6)/3],
              [  4./3*(2*sr3-1),               0, 4./3*(1-2*sr3),   4*(5-sr3)/3],
              [    -sr2*(4-sr3),             sr2,    sr2*(sr3-4), sr2*(3*sr3-2)], # row 3
              [    -(3*sr3-2)/2,    -(2-3*sr3)/2,     -(sr3+6)/2,     (sr3+6)/2],
              [              0 ,-sr2*(2*sr3-7)/3,       -5*sr6/3, sr2*(sr3+6)/3],
              [     (3*sr3-2)/2,    -(2-3*sr3)/2,     -(sr3+6)/2,     (sr3+6)/2],
              [     sr2*(4-sr3),             sr2,    sr2*(sr3-4), sr2*(3*sr3-2)],
              [      -4*(sr3-1),      -4*(1-sr3),      4*(1-sr3),   4*(3-sr3)  ], # row 4
              [           -sr2 ,    -sr2*(sr3-4),    sr2*(sr3-4), sr2*(3*sr3-2)],
              [              0 , -4./3*(1-2*sr3), 4./3*(1-2*sr3),4./3*(5-sr3)  ],
              [            sr2 ,    -sr2*(sr3-4),    sr2*(sr3-4), sr2*(3*sr3-2)],
              [       4*(sr3-1),      -4*(1-sr3),      4*(1-sr3),   4*(3-sr3)  ]]
        wmin = Surface(b,b,cp, rational=True)
        wmax = wmin.clone().mirror([0,0,1])
        vmax = wmin.clone().rotate(pi/2, [1,0,0])
        vmin = vmax.clone().mirror([0,1,0])
        umax = vmin.clone().rotate(pi/2, [0,0,1])
        umin = umax.clone().mirror([1,0,0])
        # ideally I would like to call edge_surfaces() now, but that function
        # does not work with rational surfaces, so we'll just manually try
        # and add some inner controlpoints
        cp   = np.zeros((5,5,5,4))
        cp[ :, :, 0,:] = wmin[:,:,:]
        cp[ :, :,-1,:] = wmax[:,:,:]
        cp[ :, 0, :,:] = vmin[:,:,:]
        cp[ :,-1, :,:] = vmax[:,:,:]
        cp[ 0, :, :,:] = umin[:,:,:]
        cp[-1, :, :,:] = umax[:,:,:]
        inner = np.linspace(-.5,.5, 3)
        Y, X, Z = np.meshgrid(inner,inner,inner)
        cp[1:4,1:4,1:4,0] = X
        cp[1:4,1:4,1:4,1] = Y
        cp[1:4,1:4,1:4,2] = Z
        cp[1:4,1:4,1:4,3] = 1
        ball = Volume(b,b,b,cp,rational=True, raw=True)
        return r*ball + center
    else:
        raise ValueError('invalid type argument')


def revolve(surf, theta=2 * pi):
    """  Revolve a volume by sweeping a surface in a rotational fashion around
    the *z* axis.

    :param Surface surf: Surface to revolve
    :param float theta: Angle to revolve, in radians
    :return: The revolved surface
    :rtype: Volume
    """
    surf = surf.clone()  # clone input surface, throw away old reference
    surf.set_dimension(3)  # add z-components (if not already present)
    surf.force_rational()  # add weight (if not already present)
    n = len(surf)  # number of control points of the surface
    path = CurveFactory.circle_segment(theta=theta)
    m = len(path)

    cp = np.zeros((m * n, 4))

    dt = (path.knots(0)[1] - path.knots(0)[0]) / 2.0
    for i in range(m):
        weight = path[i,-1]
        cp[i * n:(i + 1) * n, :] = np.reshape(surf.controlpoints.transpose(1, 0, 2), (n, 4))
        cp[i * n:(i + 1) * n, 2] *= weight
        cp[i * n:(i + 1) * n, 3] *= weight
        surf.rotate(dt)
    return Volume(surf.bases[0], surf.bases[1], path.bases[0], cp, True)


def cylinder(r=1, h=1, center=(0,0,0), axis=(0,0,1)):
    """  Create a solid cylinder

    :param float r: Radius
    :param float h: Height
    :param array-like center: The center of the bottom circle
    :param array-like axis: Cylinder axis
    :return: The cylinder
    :rtype: Volume
    """
    return extrude(SurfaceFactory.disc(r, center, axis), h*np.array(axis))


def extrude(surf, amount):
    """  Extrude a surface by sweeping it to a given height.

    :param Surface surf: Surface to extrude
    :param array-like amount: 3-component vector of sweeping amount and direction
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
    """  Create the volume defined by the region between the input surfaces.

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
        result = Volume(surf1.bases[0], surf1.bases[1], BSplineBasis(2), controlpoints,
                         rational=surf1.rational, raw=True)

        return result
    elif len(surfaces) == 6:
        if any([surf.rational for surf in surfaces]):
            raise RuntimeError('edge_surfaces not supported for rational splines')

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
        result                =   vol1.clone()
        result.controlpoints +=   vol2.controlpoints
        result.controlpoints +=   vol3.controlpoints
        result.controlpoints -= 2*vol4.controlpoints
        return result
    else:
        raise ValueError('Requires two or six input surfaces')

def sweep(path, shape):
    """  Generate a surface by sweeping a shape along a path

    The resulting surface is an approximation generated by interpolating at the
    Greville points. It is generated by sweeping a shape curve along a path.

    The *shape* object has to be contained in the 'xy' plane (preferably centered
    around the origin) as its x-coordinate is extruded in the normal direction,
    and its y-coordinate in the binormal direction of the *path* curve.

    :param Curve path:  The path to drag *shape* along
    :param Surface shape: The shape to be dragged out to a surface
    :return: Surrounding volume
    :rtype: Volume
    """
    b1 = shape.bases[0]
    b2 = shape.bases[1]
    b3 = path.bases[0]
    n1 = b1.num_functions()
    n2 = b2.num_functions()
    n3 = b3.num_functions()
    # this requires binormals and normals, which only work in 3D, so assume this here
    X  = np.zeros((n1,n2,n3, 3))

    # pre-evaluate the surface
    u = b1.greville()
    v = b2.greville()
    y = shape(u,v)

    for k in range(n3):
        w = b3.greville(k)
        x = path(w)
        B = path.binormal(w)
        N = path.normal(w)
        for i in range(n1):
            for j in range(n2):
                X[i,j,k,:] = x + N*y[i,j,0] + B*y[i,j,1]

    return interpolate(X, [b1,b2,b3])


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
    """  Interpolate a volume on a set of regular gridded interpolation points `x`.

    The points can be either a matrix (in which case the first index is
    interpreted as a flat row-first index of the interpolation grid) or a 4D
    tensor. In both cases the last index is the physical coordinates.

    :param numpy.ndarray x: Grid of interpolation points
    :param [BSplineBasis] bases: The basis to interpolate on
    :param [array-like] u: Parametric interpolation points, defaults to
        Greville points of the basis
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
    """  Perform a least-square fit of a point cloud `x` onto a spline basis.

    The points can be either a matrix (in which case the first index is
    interpreted as a flat row-first index of the interpolation grid) or a 4D
    tensor. In both cases the last index is the physical coordinates.

    There must be at least as many points as basis functions.

    :param numpy.ndarray x: Grid of evaluation points
    :param [BSplineBasis] bases: Basis on which to interpolate
    :param [array-like] u: Parametric values at evaluation points
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
