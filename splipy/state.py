"""This module handles the global Splipy state."""

from contextlib import contextmanager
import sys

states = ['controlpoint_relative_tolerance',
          'controlpoint_absolute_tolerance',
          'knot_tolerance']
__all__ = states + ['state']


controlpoint_relative_tolerance = 0.0
"""Relative tolerance used for matching control points."""

controlpoint_absolute_tolerance = 1e-8
"""Absolute tolerance used for matching control points."""

knot_tolerance = 1e-10
"""Absolute tolerance used for matching parametric values (knot vectors)."""


@contextmanager
def state(**kwargs):
    """A context manager for running code in a modified state.

    This takes an arbitrary number of keyword arguments, which correspond to
    state variables. The managed block runs with these values, which are then
    restored to their previous values afterwards.

    For example:

    .. code:: python

       from splipy.state import state

       with state(controlpoint_relative_tolerance=1e-2):
           # Code that runs with relative tolerance 1e-2
       # Relative tolerance restored at this point
    """
    module = sys.modules[__name__]

    # Store the current values
    before = {k: getattr(module, k) for k in states}

    # Set values according to the input
    for k, v in kwargs.items():
        setattr(module, k, v)

    yield

    # Return settings to their previous values
    for k, v in before.items():
        setattr(module, k, v)
