# -*- coding: utf-8 -*-

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

def ensure_listlike(x, dups=1):
    """Wraps x in a list if it's not list-like."""
    try:
        while len(x) < dups:
            x = list(x)
            x.append(x[-1])
        return x
    except TypeError:
        return [x] * dups
