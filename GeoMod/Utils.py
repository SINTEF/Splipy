# -*- coding: utf-8 -*-

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
