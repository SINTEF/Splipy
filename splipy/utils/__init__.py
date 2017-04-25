# -*- coding: utf-8 -*-

from __future__ import division
from itertools import combinations, product
from math import atan2, sqrt
import numpy as np

try:
    from collections.abc import Sized
except ImportError:
    from collections import Sized

def sections(src_dim, tgt_dim):
    """Generate all boundary sections from a source dimension to a target
    dimension. For example, `sections(3,1)` generates all edges on a volume.

    The return values are lists of length `src_dim` with each element either 0,
    --1 or `None`, which are suitable for passing to
    :func:`splipy.SplineObject.section`.
    """
    # Enumerate all combinations of fixed directions
    nfixed = src_dim - tgt_dim
    for fixed in combinations(range(src_dim), r=nfixed):
        # Enumerate all {0,-1}^n over the fixed directions
        for indices in product([0, -1], repeat=nfixed):
            args = [None] * src_dim
            for f, i in zip(fixed, indices[::-1]):
                args[f] = i
            yield args

def section_from_index(src_dim, tgt_dim, i):
    """Return the i'th section from a source dimension to a target dimension.

    See :func:`splipy.Utils.sections` for more information.
    """
    for j, s in sections(src_dim, tgt_dim):
        if i == j:
            return s

def section_to_index(section):
    """Return the index corresponding to a section."""
    src_dim = len(section)
    tgt_dim = sum(1 for s in section if s is None)
    for i, t in enumerate(sections(src_dim, tgt_dim)):
        if tuple(section) == tuple(t):
            return i

def check_section(*args, **kwargs):
    """check_section(u, v, ...)

    Parse arguments and return a section spec.

    The keyword argument `pardim` *must* be provided. The return value is a
    section as described in :func:`splipy.Utils.sections`.
    """
    pardim = kwargs['pardim']
    args = list(args)
    while len(args) < pardim:
        args.append(None)
    for k in set(kwargs.keys()) & set('uvw'):
        index = 'uvw'.index(k)
        args[index] = kwargs[k]
    return args

def check_direction(direction, pardim):
    if   direction in {0, 'u', 'U'} and 0 < pardim:
        return 0
    elif direction in {1, 'v', 'V'} and 1 < pardim:
        return 1
    elif direction in {2, 'w', 'W'} and 2 < pardim:
        return 2
    raise ValueError('Invalid direction')

def ensure_flatlist(x):
    """Flattens a multi-list x to a single index list."""
    if isinstance(x[0], Sized):
        return x[0]
    return x

def is_singleton(x):
    """Checks if x is list-like."""
    return not isinstance(x, Sized)

def ensure_listlike(x, dups=1):
    """Wraps x in a list if it's not list-like."""
    try:
        while len(x) < dups:
            x = list(x)
            x.append(x[-1])
        return x
    except TypeError:
        return [x] * dups
    except IndexError:
        return []

def flip_and_move_plane_geometry(obj, center=(0,0,0), normal=(0,0,1)):
    """re-orients a planar geometry by moving it to a different location and
    tilting it"""
    # don't touch it if not needed. translate or scale operations may force
    # object into 3D space
    if normal != (0,0,1):
        theta = atan2(normal[1], normal[0])
        phi   = atan2(sqrt(normal[0]**2+normal[1]**2), normal[2])
        obj.rotate(phi,   (0,1,0))
        obj.rotate(theta, (0,0,1))
    if center != (0,0,0):
        obj.translate(center)
    return obj

def reshape(cps, newshape, order='C', ncomps=None):
    """Like numpy's reshape, but preserves control points of several dimensions
    that are stored contiguously.

    The return value has shape (*newshape, ncomps), where ncomps is the number
    of components per control point, as inferred by the size of `cps` and the
    desired shape.

    The `order` argument ('C' or 'F') determines the order in which control
    points are read, but does *not* affect the order in which each component of
    a control point is read.
    """
    npts = np.prod(newshape)
    if ncomps is None:
        try:
            ncomps = cps.size // npts
        except AttributeError:
            ncomps = len(cps) // npts

    if order == 'C':
        shape = list(newshape) + [ncomps]
    elif order == 'F':
        shape = list(newshape[::-1]) + [ncomps]
    cps = np.reshape(cps, shape)
    if order == 'F':
        spec = list(range(len(newshape)))[::-1] + [len(newshape)]
        cps = cps.transpose(spec)
    return cps
