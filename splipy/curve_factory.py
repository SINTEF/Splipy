# -*- coding: utf-8 -*-

"""Handy utilities for creating curves."""

from math import pi, cos, sin, sqrt, ceil, atan2
import copy
import inspect

import numpy as np
from numpy.linalg import norm
import scipy.sparse as sp
import scipy.sparse.linalg as splinalg

from .curve import Curve
from .basis import BSplineBasis
from .utils import flip_and_move_plane_geometry, rotate_local_x_axis
from . import state

__all__ = ['Boundary', 'line', 'polygon', 'n_gon', 'circle', 'ellipse',
           'circle_segment_from_three_points', 'circle_segment', 'interpolate',
           'least_square_fit', 'cubic_curve', 'bezier', 'manipulate', 'fit']

class Boundary:
    """Enumeration representing different boundary conditions used in
    :func:`interpolate`."""

    FREE = 1
    """The curve will be smooth at the second and second-to-last unique knot."""

    NATURAL = 2
    """The curve will have zero second derivatives at the endpoints."""

    HERMITE = 3
    """Specify the derivatives at the knots."""

    PERIODIC = 4
    """The curve will be periodic at the endpoints."""

    TANGENT = 5
    """Specify the tangents at the endpoints."""

    TANGENTNATURAL = 6
    """Use `TANGENT` for the start and `NATURAL` for the end."""


def line(a, b, relative=False):
    """  Create a line between two points.

    :param array-like a: Start point
    :param array-like b: End point
    :param bool relative: Whether *b* is relative to *a*
    :return: Linear curve from *a* to *b*
    :rtype: Curve
    """
    if relative:
        b = tuple(ai + bi for ai, bi in zip(a, b))
    return Curve(controlpoints=[a, b])


def polygon(*points, **keywords):
    """  Create a linear interpolation between input points.

    :param [array-like] points: The points to interpolate
    :param bool relative: If controlpoints are interpreted as relative to the
        previous one
    :param [float] t: specify parametric interpolation points (the knot vector)
    :return: Linear curve through the input points
    :rtype: Curve
    """
    if len(points) == 1:
        points = points[0]

    knot = keywords.get('t', [])
    if len(knot) == 0: # establish knot vector based on eucledian length between points
        knot = [0, 0]
        prevPt = points[0]
        for pt in points[1:]:
            dist = 0
            for (x0, x1) in zip(prevPt, pt):  # loop over (x,y) and maybe z-coordinate
                dist += (x1 - x0)**2
            knot.append(knot[-1] + sqrt(dist))
            prevPt = pt
        knot.append(knot[-1])
    else: # use knot vector given as input argument
        knot = [knot[0]] + list(knot) + [knot[-1]]

    relative = keywords.get('relative', False)
    if relative:
        points = list(points)
        for i in range(1, len(points)):
            points[i] = [x0 + x1 for (x0,x1) in zip(points[i-1], points[i])]

    return Curve(BSplineBasis(2, knot), points)


def n_gon(n=5, r=1, center=(0,0,0), normal=(0,0,1)):
    """  Create a regular polygon of *n* equal sides centered at the origin.

    :param int n: Number of sides and vertices
    :param float r: Radius
    :param array-like center: local origin
    :param array-like normal: local normal
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
    knot = [-1]
    for i in range(n):
        cp.append([r * cos(i * dt), r * sin(i * dt)])
        knot.append(i)
    knot += [n, n+1]
    basis = BSplineBasis(2, knot, 0)

    result =  Curve(basis, cp)
    return flip_and_move_plane_geometry(result, center, normal)

def circle(r=1, center=(0,0,0), normal=(0,0,1), type='p2C0', xaxis=(1,0,0)):
    """  Create a circle.

    :param float r: Radius
    :param array-like center: local origin
    :param array-like normal: local normal
    :param string type: The type of parametrization ('p2C0' or 'p4C1')
    :param array-like xaxis: direction of sem, i.e. parametric start point t=0
    :return: A periodic, quadratic rational curve
    :rtype: Curve
    :raises ValueError: If radius is not positive
    """
    if r <= 0:
        raise ValueError('radius needs to be positive')

    if type == 'p2C0' or type == 'C0p2':
        w = 1.0 / sqrt(2)
        controlpoints = [[1, 0, 1],
                         [w, w, w],
                         [0, 1, 1],
                         [-w, w, w],
                         [-1, 0, 1],
                         [-w, -w, w],
                         [0, -1, 1],
                         [w, -w, w]]
        knot = np.array([-1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5]) / 4.0 * 2 * pi

        result = Curve(BSplineBasis(3, knot, 0), controlpoints, True)
    elif type.lower() == 'p4c1' or type.lower() == 'c1p4':
        w = 2*sqrt(2)/3
        a = 1.0/2/sqrt(2)
        b = 1.0/6 * (4*sqrt(2)-1)
        controlpoints = [[ 1,-a, 1],
                         [ 1, a, 1],
                         [ b, b, w],
                         [ a, 1, 1],
                         [-a, 1, 1],
                         [-b, b, w],
                         [-1, a, 1],
                         [-1,-a, 1],
                         [-b,-b, w],
                         [-a,-1, 1],
                         [ a,-1, 1],
                         [ b,-b, w]]
        knot = np.array([ -1, -1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5]) / 4.0 * 2 * pi
        result = Curve(BSplineBasis(5, knot, 1), controlpoints, True)
    else:
        raise ValueError('Unkown type: %s' %(type))

    result *= r
    result.rotate(rotate_local_x_axis(xaxis, normal))
    return flip_and_move_plane_geometry(result, center, normal)

def ellipse(r1=1, r2=1, center=(0,0,0), normal=(0,0,1), type='p2C0', xaxis=(1,0,0)):
    """  Create an ellipse

    :param float r1: Radius along xaxis
    :param float r2: Radius orthogonal to xaxis
    :param array-like center: local origin
    :param array-like normal: local normal
    :param string type: The type of parametrization ('p2C0' or 'p4C1')
    :param array-like xaxis: direction of sem, i.e. parametric start point t=0
    :return: A periodic, quadratic rational curve
    :rtype: Curve
    :raises ValueError: If radius is not positive
    """
    result = circle(type=type)
    result *= [r1,r2,1]
    result.rotate(rotate_local_x_axis(xaxis, normal))
    return flip_and_move_plane_geometry(result, center, normal)

def circle_segment_from_three_points(x0, x1, x2):
    """circle_segment_from_three_points(x0, x1, x2)

    Create a circle segment going from the point x0 to x2 through x1

    :param array-like x0: The start point (2D or 3D point)
    :param array-like x1: An intermediate point (2D or 3D)
    :param array-like x2: The end point (2D or 3D)
    :rtype: Curve
    """
    # wrap input into 3d numpy arrays
    pt0 = np.array([0,0,0], dtype='float')
    pt1 = np.array([0,0,0], dtype='float')
    pt2 = np.array([0,0,0], dtype='float')
    pt0[:len(x0)] = x0
    pt1[:len(x1)] = x1
    pt2[:len(x2)] = x2

    # figure out normal, center and radius
    normal = np.cross(pt1-pt0, pt2-pt0)
    A = np.vstack((2*(pt1-pt0), 2*(pt2-pt0), normal))
    b = np.array([ np.dot(pt1,pt1) - np.dot(pt0,pt0),
                   np.dot(pt2,pt2) - np.dot(pt0,pt0),
                   np.dot(normal,pt0)])
    center = np.linalg.solve(A,b)
    radius = norm(pt2-center)
    v2 = pt2-center
    v1 = pt1-center
    v0 = pt0-center
    w0 = pt0-pt2
    w1 = pt1-pt2
    w2 = np.cross(w0, w1)
    normal = np.cross(v0,v2)
    len_v2 = norm(v2)
    len_v0 = norm(v0)
    theta  = np.arccos(np.dot(v2,v0) / len_v2 / len_v0)
    if not np.all([np.sign(i)==np.sign(j) or abs(i-j) < state.controlpoint_absolute_tolerance for (i,j) in zip(w2,normal)]):
        theta = 2*pi - theta
        normal = -normal

    result = circle_segment(theta, radius, center, np.cross(v0,v1), v0)

    # spit out 2D curve if all input points were 2D, otherwise return 3D
    result.set_dimension(np.max([len(x0), len(x1), len(x2)]))
    return result


def circle_segment(theta, r=1, center=(0,0,0), normal=(0,0,1), xaxis=(1,0,0)):
    """  Create a circle segment starting parallel to the rotated x-axis.

    :param float theta: Angle in radians
    :param float r: Radius
    :param array-like center: circle segment center
    :param array-like normal: normal vector to the plane that contains circle
    :param array-like xaxis: direction of the parametric start point t=0
    :return: A quadratic rational curve
    :rtype: Curve
    :raises ValueError: If radius is not positive
    :raises ValueError: If theta is not in the range *[-2pi, 2pi]*
    """
    # error test input
    if abs(theta) > 2 * pi:
        raise ValueError('theta needs to be in range [-2pi,2pi]')
    if r <= 0:
        raise ValueError('radius needs to be positive')
    if theta == 2*pi:
        return circle(r, center, normal)

    # build knot vector
    knot_spans = int(ceil(abs(theta) / (2 * pi / 3)))
    knot = [0]
    for i in range(knot_spans + 1):
        knot += [i] * 2
    knot += [knot_spans]  # knot vector [0,0,0,1,1,2,2,..,n,n,n]
    knot = np.array(knot) / float(knot[-1]) * theta  # set parametric space to [0,theta]

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

    if theta < 0:
        cp.reverse()
        result = Curve(BSplineBasis(3, np.flip(knot,0)), cp, True)
    else:
        result = Curve(BSplineBasis(3, knot), cp, True)
    result.rotate(rotate_local_x_axis(xaxis, normal))
    return flip_and_move_plane_geometry(result, center, normal)

def interpolate(x, basis, t=None):
    """  Perform general spline interpolation on a provided basis.

    :param matrix-like x: Matrix *X[i,j]* of interpolation points *xi* with
        components *j*
    :param BSplineBasis basis: Basis on which to interpolate
    :param array-like t: parametric values at interpolation points; defaults to
        Greville points if not provided
    :return: Interpolated curve
    :rtype: Curve
    """
    # wrap input into an array
    x = np.array(x)

    # evaluate all basis functions in the interpolation points
    if t is None:
        t = basis.greville()
    N = basis.evaluate(t, sparse=True)

    # solve interpolation problem
    cp = splinalg.spsolve(N, x)
    cp = cp.reshape(x.shape)

    return Curve(basis, cp)

def least_square_fit(x, basis, t):
    """  Perform a least-square fit of a point cloud onto a spline basis

    :param matrix-like x: Matrix *X[i,j]* of interpolation points *xi* with
        components *j*. The number of points must be equal to or larger than
        the number of basis functions in *basis*
    :param BSplineBasis basis: Basis on which to interpolate
    :param array-like t: parametric values at evaluation points
    :return: Approximated curve
    :rtype: Curve
    """

    # evaluate all basis functions at evaluation points
    N = basis.evaluate(t)

    # solve interpolation problem
    controlpoints,_,_,_ = np.linalg.lstsq(N, x, rcond=None)

    return Curve(basis, controlpoints)


def cubic_curve(x, boundary=Boundary.FREE, t=None, tangents=None):
    """  Perform cubic spline interpolation on a provided basis.

    The valid boundary conditions are enumerated in :class:`Boundary`. The
    meaning of the `tangents` parameter depends on the specified boundary
    condition:

    - `TANGENT`: two points,
    - `TANGENTNATURAL`: one point,
    - `HERMITE`: *n* points

    :param matrix-like x: Matrix *X[i,j]* of interpolation points *x_i* with
        components *j*
    :param int boundary: Any value from :class:`Boundary`.
    :param array-like t: parametric values at interpolation points, defaults
        to Euclidean distance between evaluation points
    :param matrix-like tangents: Tangent information according to the boundary
        conditions.
    :return: Interpolated curve
    :rtype: Curve
    """

    # if periodic input is not closed, make sure we do it now
    if boundary == Boundary.PERIODIC and not (
       np.allclose(x[0,:], x[-1,:],
                   rtol = state.controlpoint_relative_tolerance,
                   atol = state.controlpoint_absolute_tolerance)):
        x = np.append(x, [x[0,:]], axis=0)
        if t is not None:
            # augment interpolation knot by euclidian distance to end
            t = list(t) + [t[-1] + norm(np.array(x[0,:])- np.array(x[-2,:]))]

    n = len(x)
    if t is None:
        t = [0.0]
        for (x0,x1) in zip(x[:-1,:], x[1:,:]):
            # eucledian distance between two consecutive points
            dist = norm(np.array(x1)-np.array(x0))
            t.append(t[-1]+dist)

    # modify knot vector for chosen boundary conditions
    knot = [t[0]]*3 + list(t) + [t[-1]]*3
    if boundary == Boundary.FREE:
        del knot[-5]
        del knot[4]
    elif boundary == Boundary.HERMITE:
        knot = sorted(list(knot) + list(t[1:-1]))

    # create the interpolation basis and interpolation matrix on this
    if boundary == Boundary.PERIODIC:
        # C2-periodic knots
        knot[0]  = t[0]  + t[-4] - t[-1]
        knot[1]  = t[0]  + t[-3] - t[-1]
        knot[2]  = t[0]  + t[-2] - t[-1]
        knot[-3] = t[-1] + t[1]  - t[0]
        knot[-2] = t[-1] + t[2]  - t[0]
        knot[-1] = t[-1] + t[3]  - t[0]

        basis = BSplineBasis(4, knot, 2)

        # do not duplicate the interpolation at the sem (start=end is the same point)
        # identical points equal singular interpolation matrix
        t = t[:-1]
        x = x[:-1,:]
    else:
        basis = BSplineBasis(4, knot)
    N = basis(t, sparse=True)  # left-hand-side matrix

    # add derivative boundary conditions if applicable
    if boundary in [Boundary.TANGENT, Boundary.HERMITE, Boundary.TANGENTNATURAL]:
        if boundary == Boundary.TANGENT:
            dn = basis([t[0], t[-1]], d=1)
        elif boundary == Boundary.TANGENTNATURAL:
            dn = basis(t[0], d=1)
        elif boundary == Boundary.HERMITE:
            dn = basis(t, d=1)
        N  = sp.vstack([N, sp.csr_matrix(dn)])
        x  = np.vstack([x, tangents])

    # add double derivative boundary conditions if applicable
    if boundary in [Boundary.NATURAL, Boundary.TANGENTNATURAL]:
        if boundary == Boundary.NATURAL:
            ddn  = basis([t[0], t[-1]], d=2)
            new  = 2
        elif boundary == Boundary.TANGENTNATURAL:
            ddn  = basis(t[-1], d=2)
            new  = 1
        N  = sp.vstack([N, sp.csr_matrix(ddn)])
        x  = np.vstack([x, np.zeros((new, x.shape[1]))])

    # solve system to get controlpoints
    cp = splinalg.spsolve(N,x)
    cp = cp.reshape(x.shape)

    # wrap it all into a curve and return
    return Curve(basis, cp)

def bezier(pts, quadratic=False, relative=False):
    """  Generate a cubic or quadratic bezier curve from a set of control points

    :param [array-like] pts: list of control-points. In addition to a starting
        point we need three points per bezier interval for cubic splines and
        two points for quadratic splines
    :param bool quadratic: True if a quadratic spline is to be returned, False
        if a cubic spline is to be returned
    :param bool relative: If controlpoints are interpreted as relative to the
        previous one
    :return: Bezier curve
    :rtype: Curve

    """
    if quadratic:
        p = 3
    else:
        p = 4
    # compute number of intervals
    n = int((len(pts)-1)/(p-1))
    # generate uniform knot vector of repeated integers
    knot = list(range(n+1)) * (p-1) + [0, n]
    knot.sort()
    if relative:
        pts = copy.deepcopy(pts)
        for i in range(1, len(pts)):
            pts[i] = [x0 + x1 for (x0,x1) in zip(pts[i-1], pts[i])]
    return Curve(BSplineBasis(p, knot), pts)

def manipulate(crv, f, normalized=False, vectorized=False):
    """  Create a new curve based on an expression-evaluation of an existing one
    :param Curve crv: original curve on which f is to be applied
    :param function f: expression of the physical point *x*, the velocity
        (tangent) *v*, parametric point *t* and/or acceleration *a*.
    :param normalized: If velocity and acceleration terms should be normalized
        to have length 1
    :param vectorized: True if *f* is expressed in terms of vectorized
        operations. The method is sped up whenever this is used.

    Examples:

    .. code:: python

        def scale_by_two(x):
            return 2*x

        new_curve = manipulate(old_curve, scale_by_two)
        new_curve = old_curve * 2 # will give the same result

        def move_3_to_right(x, v):
            result = x
            # *v* is velocity. Rotate this 90 to the right and it points (+v[1], -v[0])
            result[0] += 3*v[1]
            result[1] -= 3*v[0]
            return result

        # note that the velocity vector will not have length one unless normalized is passed
        new_curve = manipulate(old_curve, move_3_to_right, normalized=True)

        def move_3_to_right_fast(x, v):
            result = x
            result[:,0] += 3*v[:,1]
            result[:,1] -= 3*v[:,0]
            return result

        new_curve = manipulate(old_curve, move_3_to_right_fast, normalized=True, vectorized=True)
    """

    b = crv.bases[0]
    t = np.array(b.greville())
    n = len(crv)

    if vectorized:
        x = crv(t)
        arg_names = inspect.getargspec(f).args
        argc = len(arg_names)
        argv = [0] * argc
        for j in range(argc):
            if arg_names[j] == 'x':
                argv[j] = x
            elif arg_names[j] == 't':
                argv[j] = t
            elif arg_names[j] == 'v':
                c0 = np.array([i for i in range(n) if b.continuity(t[i]) == 0])
                v = crv.derivative(t, 1)
                if len(c0)>0:
                    v[c0,:] = (v[c0,:] + crv.derivative(t[c0], 1, above=False)) / 2.0
                if normalized:
                    v[:] = [vel / norm(vel) for vel in v]
                argv[j] = v
            elif arg_names[j] == 'a':
                c1 = np.array([i for i in range(n) if b.continuity(t[i]) < 2])
                a = crv.derivative(t, 2)
                if len(c1)>0:
                    a[c1,:] = (a[c1,:] + crv.derivative(t[c1], 2, above=False)) / 2.0
                if normalized:
                    a[:] = [acc / norm(acc) for acc in a]
                argv[j] = a
        destination   = f(*argv)
    else:
        destination = np.zeros((len(crv), crv.dimension))
        for (t1, i) in zip(t, range(len(t))):
            x = crv(t1)
            arg_names = inspect.getargspec(f).args
            argc = len(arg_names)
            argv = [0] * argc
            for j in range(argc):
                if arg_names[j] == 'x':
                    argv[j] = x
                elif arg_names[j] == 't':
                    argv[j] = t1
                elif arg_names[j] == 'v':
                    v = crv.derivative(t1, 1)
                    if b.continuity(t1) < 1:
                        v += crv.derivative(t1, 1, above=False)
                        v /= 2.0
                    if normalized:
                        v /= norm(v)
                    argv[j] = v
                elif arg_names[j] == 'a':
                    a = crv.derivative(t1, 2)
                    if b.continuity(t1) < 2:
                        a += crv.derivative(t1, 1, above=False)
                        a /= 2.0
                    if normalized:
                        a /= norm(a)
                    argv[j] = a
            destination[i] = f(*argv)

    N = b(t, sparse=True)
    controlpoints = splinalg.spsolve(N, destination)
    return Curve(b, controlpoints)

def fit(x, t0, t1, rtol=1e-4, atol=0.0):
    """  Computes an interpolation for a parametric curve up to a specified tolerance.
    The method will iteratively refine parts where needed resulting in a non-uniform
    knot vector with as optimized knot locations as possible.

    :param function x: callable function which takes as input a vector of evaluation
        points t and gives as output a matrix x where x[i,j] is component j evaluated
        at point t[i]
    :param float t0: start of parametric domain
    :param float t1: end of parametric domain
    :param float rtol: relative tolerance for stopping criterium. It is defined
        to be ||e||_L2 / D, where D is the length of the curve and ||e||_L2 is the
        L2-error (see Curve.error)
    :param float atol: absolute tolerance for stopping criterium. It is defined to
        be the maximal distance between the curve approximation and the exact curve
    :return: Curve Non-uniform cubic B-spline curve

    Examples:

    .. code:: python

        import numpy as np
        import splipy.curve_factory as curve_factory

        # gives a B-spline approximation to the circle with arclength parametrization
        # unlike curve_factory.circle which is exact, but not arclength
        def arclength_circle(t):
            return np.array( [np.cos(t), np.sin(t)] ).T
        crv = curve_factory.fit(arclength_circle, 0, 2*np.pi)
        print(crv)

        # approximates a difficult function with wild behaviour around t=0, but
        # this is overcome by a higher knot density around this point
        def one_over_t(t):
            eps = 1e-8 # to avoid 1/0 we add a small epsilon
            return np.array( [t, 1.0/(t+eps)] ).T
        crv = curve_factory.fit(one_over_t, 0, 1, rtol=1e-6)
        print(crv) # first knot span is ~1e-9, last knot is ~1e-1

        # one can specify the target curve in terms of existing Curve objects
        crv = curve_factory.circle(r=1)     # Curve-object, quadratic NURBS
        def move_along_tangent(t):
            return crv(t) + crv.tangent(t)  # can evaluate curve or its derivatives
        # fit() will create a B-spline approximation using non-uniform refinement
        crv2 = curve_factory.fit(move_along_tangent, crv.start(0), crv.end(0))
    """

    b = BSplineBasis(4, [t0,t0,t0,t0, t1,t1,t1,t1])
    t = np.array(b.greville())
    crv = interpolate(x(t), b, t)
    (err2, maxerr) = crv.error(x)
    # polynomial input (which can be exactly represented) only use one knot span
    if maxerr < 1e-13:
        return crv

    # for all other curves, start with 4 knot spans
    knot_vector = [t0,t0,t0,t0] + [i/5.0*(t1-t0)+t0 for i in range(1,5)] + [t1,t1,t1,t1]
    b = BSplineBasis(4, knot_vector)
    t = np.array(b.greville())
    crv = interpolate(x(t), b, t)
    (err2, maxerr) = crv.error(x)
    # this is technically false since we need the length of the target function *x*
    # and not our approximation *crv*, but we don't have the derivative of *x*, so
    # we can't compute it. This seems like a healthy compromise
    length = crv.length()
    while np.sqrt(np.sum(err2))/length > rtol and maxerr > atol:
        knot_span    = crv.knots(0) # knot vector without multiplicities
        target_error = (rtol*length)**2 / len(err2) # equidistribute error among all knot spans
        refinements  = []
        for i in range(len(err2)):
            # figure out how many new knots we require in this knot interval:
            # if we converge with *scale* and want an error of *target_error*
            # |e|^2 * (1/n)^scale = target_error^2

            conv_order = 4                   # cubic interpolateion is order=4
            square_conv_order = 2*conv_order # we are computing with square of error
            scale = square_conv_order + 4    # don't want to converge too quickly in case of highly non-uniform mesh refinement is required
            n = int(np.ceil(np.exp((np.log(err2[i]) - np.log(target_error))/scale)))

            # add *n* new interior knots to this knot span
            new_knots = np.linspace(knot_span[i], knot_span[i+1], n+1)
            knot_vector = knot_vector + list(new_knots[1:-1])

        # build new refined knot vector
        knot_vector.sort()
        b = BSplineBasis(4, knot_vector)
        # do interpolation and return result
        t = np.array(b.greville())
        crv = interpolate(x(t), b, t)
        (err2, maxerr) = crv.error(x)
        length = crv.length()

    return crv

def fit_points(x, t=[], rtol=1e-4, atol=0.0):
    """ Computes an approximation for a list of points up to a specified tolerance.
    The method will iteratively refine parts where needed resulting in a non-uniform
    knot vector with as optimized knot locations as possible. The target curve is the
    linear interpolation of the input points

    :param [point-like] x: The points to approximate
    :param [float] t: parametric values for curve points (same length as x)
    :param float rtol: relative tolerance for stopping criterium. It is defined
        to be ||e||_L2 / D, where D is the length of the curve and ||e||_L2 is the
        L2-error (see Curve.error)
    :param float atol: absolute tolerance for stopping criterium. It is defined to
        be the maximal distance between the curve approximation and the exact curve
    :return: Curve Non-uniform cubic B-spline curve
    """

    if len(t)>0:
        linear = polygon(x, t=t)
    else:
        linear = polygon(x)
    return fit(linear, linear.start(0), linear.end(0), rtol=rtol, atol=atol)
