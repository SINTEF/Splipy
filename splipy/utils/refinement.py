from __future__ import division

__doc__ = 'Implementation of various refinement schemes.'

from splipy.utils import ensure_listlike, check_direction
from math import atan, pi
import numpy as np

# TODO: put control over these tolerances somewhere. Modstate in splipy seems
#       to be the place for it, but we can't let splipy.utils influence the
#       structure of splipy.
def knot_exists(existing_knots, new_knot):
    return np.any(np.isclose(existing_knots, new_knot, atol=1e-7, rtol=1e-10))


def geometric_refine(obj, alpha, n, direction=0, reverse=False):
    """geometric_refine(obj, alpha, n, [direction=0], [reverse=False])

    Refine a spline object by making a geometric distribution of element sizes.

    :param obj: The object to refine
    :type obj: :class:`splipy.SplineObject`
    :param float alpha: The length ratio between two sequential knot segments
    :param int n: The number of knots to insert
    :param direction: The direction to refine in
    :param bool reverse: Set to `True` to refine towards the other end
    """
    # some error tests on input
    if n<=0:
        raise ValueError('n should be greater than 0')

    direction = check_direction(direction, obj.pardim)
    if reverse:
        obj.reverse(direction)

    # fetch knots
    knots = obj.knots()
    knot_start = knots[direction][0]
    knot_end   = knots[direction][-1]
    dk = knot_end - knot_start

    # evaluate the factors
    n = n+1 # redefine n to be knot spans instead of new internal knots
    totProd = 1.0
    totSum  = 0.0
    for i in range(n):
        totSum  += totProd
        totProd *= alpha
    d1 = 1.0 / totSum
    knot = d1

    # compute knot locations
    new_knots = []
    for i in range(n-1):
        k = knot_start + knot*dk
        if not knot_exists(knots[direction], k):
            new_knots.append(k)
        knot += alpha*d1
        d1   *= alpha

    # do the actual knot insertion
    obj.insert_knot(new_knots, direction)

    if reverse:
        obj.reverse(direction)
    return obj


def edge_refine(obj, S, n, direction=0):
    """edge_refine(obj, S, n, [direction=0])

    Refine an object towards both edges in a direction, by sampling an
    arctan function.

    :param obj: The object to refine
    :type obj: :class:`splipy.SplineObject`
    :param float S: The slope of the arctan function.
    :param int n: The number of knots to insert
    :param direction: The direction to refine in
    """
    # some error tests on input
    if n<=0:
        raise ValueError('n should be greater than 0')

    direction = check_direction(direction, obj.pardim)

    # fetch knots
    knots = obj.knots()
    knot_start = knots[direction][0]
    knot_end   = knots[direction][-1]
    dk = knot_end - knot_start

    # compute knot locations
    new_knots = []
    max_atan  = atan(S)
    for i in range(1,n+1):
        xi  = -1.0 + 2.0*i/(n+1)
        xi *= S
        k   = knot_start + (atan(xi)+max_atan)/2/max_atan*dk
        if not knot_exists(knots[direction], k):
            new_knots.append(k)

    # do the actual knot insertion
    obj.insert_knot(new_knots, direction)
    return obj


def _splitvector(len, parts):
    delta = len // parts
    sizes = [delta for i in range(parts)]
    remainder = len-parts*delta
    for i in range(parts-remainder+1, parts):
        sizes[i] = sizes[i]+1
    result = [0]
    for i in range(1,parts):
        result.append(sizes[i]+result[i-1])
    return result


def subdivide(objs, n):
    """Subdivide a list of objects by splitting them up along existing knot
    lines. The resulting partition will roughly the same number of elements on
    all pieces. By splitting along *n* lines, we generate *n* + 1 new blocks.

    The number of subdivisions can be a list of integers: one for each
    direction, or a single integer for uniform sudivision.

    :param objs: Objects to split
    :param n: Number of subdivisions to perform
    :type n: int or [int]
    :return: New objects
    :rtype: [:class:`splipy.SplineObject`]
    """
    pardim = objs[0].pardim # 1 for curves, 2 for surfaces, 3 for volumes
    n = ensure_listlike(n, pardim)

    result = objs
    for d in range(pardim):
        # split all objects so far along direction d
        new_results = []
        for obj in result:
            splitting_points = [obj.knots(d)[i] for i in _splitvector(len(obj.knots(d)), n[d]+1)]
            new_results += obj.split(splitting_points[1:], d)

        # only keep the smallest pieces in our result list
        result = new_results

    return result

