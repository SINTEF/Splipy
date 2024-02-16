"Implementation of various smoothing operations on a per-controlpoint level."

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np
from scipy import ndimage

from . import check_direction

if TYPE_CHECKING:
    from splipy.splineobject import SplineObject
    from splipy.types import Direction


def smooth(obj: SplineObject, comp: Optional[Direction] = None) -> None:
    """Smooth an object by setting the interior control points to the average of
    itself and all neighbours (e.g. 9 for surfaces, 27 for volumes). The edges
    are kept unchanged, and any rational weights are kept unchanged.

    :param obj: The object to smooth
    :type obj: :class:`splipy.SplineObject`
    :param comp: which component to smooth 0=x-component, 1=y-component, 2=z-component, None=all components
    :type comp: int or None
    """
    n = obj.shape

    averaging_mask = np.ones([3] * len(n) + [1])
    averaging_mask /= averaging_mask.size

    new_controlpoints = ndimage.convolve(obj.controlpoints, averaging_mask, mode="wrap")

    # build up the indexing for the domain 'interior'. This would be
    # controlpoints[1:-1, 1:-1 ,:]        for non-rational surface
    # controlpoints[ :  , 1:-1 ,:]        for surfaces which is periodic in 'u'
    # controlpoints[1:-1, 1:-1 ,:-1]      for rational surfaces
    # controlpoints[1:-1, 1:-1 , 1:-1, :] for non-rational volumes
    # controlpoints[ :  ,  :   ,  :  , :] for non-rational volumes
    #                                     which are periodic in all three parametric directions
    interior = []
    for pardim in range(len(n)):
        if obj.periodic(pardim):
            interior.append(slice(None, None, None))
        else:
            interior.append(slice(1, -1, None))
    if obj.rational:
        interior.append(slice(0, -1, None))
    elif comp is not None:
        cix = check_direction(comp, obj.dimension)
        interior.append(slice(cix, cix + 1, None))
    else:
        interior.append(slice(None, None, None))

    ix = tuple(interior)
    obj.controlpoints[ix] = new_controlpoints[ix]
