# -*- coding: utf-8 -*-

from itertools import combinations, product

def sections(src_dim, tgt_dim):
    """Generate all boundary sections from a source dimension to a target
    dimension. For example, `sections(3,1)` generates all edges on a volume.

    The return values are lists of length `src_dim` with each element either 0,
    --1 or `None`, which are suitable for passing to
    :func:`GeoMod.SplineObject.section`.
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

    See :func:`GeoMod.Utils.sections` for more information.
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
    section as described in :func:`GeoMod.Utils.sections`.
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
    try:
        len(x[0])
        return x[0]
    except TypeError:
        return x

def is_singleton(x):
    """Checks if x is list-like."""
    try:
        len(x)
        return False
    except TypeError:
        return True

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
