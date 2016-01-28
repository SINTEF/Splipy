def ensure_listlike(x, dups=1):
    """Wraps x in a list if it's not list-like."""
    try:
        len(x)
        return x
    except TypeError:
        return [x] * dups
