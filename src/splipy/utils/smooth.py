__doc__ = 'Implementation of various smoothing operations on a per-controlpoint level.'

from scipy import ndimage
import numpy as np

from . import check_direction


def smooth(obj, comp=None):
    """Smooth an object by setting the interior control points to the average of
    itself and all neighbours (e.g. 9 for surfaces, 27 for volumes). The edges
    are kept unchanged, and any rational weights are kept unchanged.

    :param obj: The object to smooth
    :type obj: :class:`splipy.SplineObject`
    :param comp: which component to smooth 0=x-component, 1=y-component, 2=z-component, None=all components
    :type comp: int or None
    """
    n = obj.shape
    if comp is not None:
        comp = check_direction(comp, len(obj))

    averaging_mask  = np.ones([3]*len(n)+[1])
    averaging_mask /= averaging_mask.size

    new_controlpoints = ndimage.convolve(obj.controlpoints, averaging_mask, mode='wrap')

    # build up the indexing for the domain 'interior'. This would be
    # controlpoints[1:-1, 1:-1 ,:]        for non-rational surface
    # controlpoints[ :  , 1:-1 ,:]        for surfaces which is periodic in 'u'
    # controlpoints[1:-1, 1:-1 ,:-1]      for rational surfaces
    # controlpoints[1:-1, 1:-1 , 1:-1, :] for non-rational volumes
    # controlpoints[ :  ,  :   ,  :  , :] for non-rational volumes which is periodic in all three parametric directions
    interior = []
    for pardim in range(len(n)):
        if obj.periodic(pardim):
            interior.append(slice(None,None,None))
        else:
            interior.append(slice(1,-1,None))
    if obj.rational:
        interior.append(slice(0,-1,None))
    elif comp is not None:
        interior.append(slice(comp,comp+1,None))
    else:
        interior.append(slice(None,None,None))

    interior = tuple(interior)
    obj.controlpoints[interior] =  new_controlpoints[interior]
