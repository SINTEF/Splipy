# -*- coding: utf-8 -*-

def check_direction(direction, pardim):
    if   (direction == 0 or direction == 'u' or direction == 'U') and 0<pardim:
        return 0
    elif (direction == 1 or direction == 'v' or direction == 'V') and 1<pardim:
        return 1
    elif (direction == 2 or direction == 'w' or direction == 'W') and 2<pardim:
        return 2
    else:
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
