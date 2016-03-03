from __future__ import division

__doc__ = 'Implementation of various refinement schemes.'

from GeoMod.Utils import ensure_listlike
from math import atan, pi
import numpy as np

# TODO: put control over these tolerances somewhere. Modstate in GeoMod seems
#       to be the place for it, but we can't let GeoUtils influence the
#       structure of GeoMod.
def knot_exists(existing_knots, new_knot):
    return np.any(np.isclose(existing_knots, new_knot, atol=1e-7, rtol=1e-10))

# Geometric distribution of knots
def geometric_refine(obj, alpha, n, direction=1):
    """Refine a SplineObject by making a geometric distribution of element sizes.
    To refine close to the other side, give direction as a negative integer
    :param SplineObject obj: The object to refine
    :param float  alpha    : The ratio between two sequential knot segments
    :param int    n        : The number of knots to insert
    :param int    direction: The direction to refine in (u=1, v=2 or w=3) 
    :return                : None
    """
    
    # some error tests on input
    if n<=0:
        raise ValueError('n should be greater than 0')

    flipBack = False
    if direction < 0:
        obj.reverse(-direction-1)
        direction = -direction
        flipBack = True

    # fetch knots
    knots = obj.knots()
    knot_start = knots[direction-1][0]
    knot_end   = knots[direction-1][-1]
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
        if not knot_exists(knots[direction-1], k):
            new_knots.append(k)
        knot += alpha*d1
        d1   *= alpha

    # do the actual knot insertion
    obj.insert_knot(new_knots, direction-1)

    if flipBack:
        obj.reverse(direction-1)


# Edge refinement
def edge_refine(obj, S, n, direction=1):
    """Refine an object by both edges, by sampling a arctan-function
    :param SplineObject obj      : The object to refine
    :param float        S        : The slope of the atan-function
    :param int          n        : The number of knots to insert
    :param int          direction: The direction to refine in (u=1, v=2, w=3) 
    :return                 : None
    """
    
    # some error tests on input
    if n<=0:
        raise ValueError('n should be greater than 0')
    
    # fetch knots
    knots = obj.knots()
    knot_start = knots[direction-1][0]
    knot_end   = knots[direction-1][-1]
    dk = knot_end - knot_start
    
    # compute knot locations
    new_knots = []
    max_atan  = atan(S)
    for i in range(1,n+1):
        xi  = -1.0 + 2.0*i/(n+1)
        xi *= S
        k   = knot_start + (atan(xi)+max_atan)/2/max_atan*dk
        if not knot_exists(knots[direction-1], k):
            new_knots.append(k)
    
    # do the actual knot insertion
    obj.insert_knot(new_knots, direction-1)

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
       lines. The resulting partition will roughly the same number of elements
       on all pieces. By splitting along *n* lines, we generate *n+1* new blocks
       :param [SplineObject] srfs : SplineObject to split
       :param int or [int]   n    : Number of subdivisions to perform, uniformly
                                    in all directions, or a list of each direction
                                    separately.
       :return                    : List of new SplineObjects
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

