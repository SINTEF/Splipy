__doc__ = 'Implementation of various smoothing operations on a per-controlpoint level.'

from scipy import ndimage
import numpy as np
import splipy.curve_factory as cf

def smooth(obj):
    """Smooth an object by setting the interior control points to the average of
    itself and all neighbours (e.g. 9 for surfaces, 27 for volumes). The edges
    are kept unchanged, and any rational weights are kept unchanged.

    :param obj: The object to smooth
    :type obj: :class:`splipy.SplineObject`
    """
    n = obj.shape
    averaging_mask  = np.ones([3]*len(n)+[1])
    averaging_mask /= averaging_mask.size

    new_controlpoints = ndimage.convolve(obj.controlpoints, averaging_mask)

    if obj.rational:
        interior = tuple([slice(1,-1,None)]*len(n) + [slice(0,-1,None)])
    else:
        interior = tuple([slice(1,-1,None)]*len(n) + [slice(None,None,None)])

    obj.controlpoints[interior] =  new_controlpoints[interior]
