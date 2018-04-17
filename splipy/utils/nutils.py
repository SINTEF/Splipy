from __future__ import division

__doc__ = 'Implementation of convenience methods with respect to nutils integration.'

from splipy import Curve, Surface, Volume
import numpy as np


def splipy_to_nutils(spline):
    """ Return controlpoints according to nutils ordering """
    n = len(spline)
    dim = spline.dimension
    if(type(spline) == Curve):
        return np.reshape(spline[:,:]                  , (n, dim), order='F')
    elif(type(spline) == Surface):
        return np.reshape(spline[:,:,:].swapaxes(0,1)  , (n, dim), order='F')
    elif(type(spline) == Volume):
        return np.reshape(spline[:,:,:,:].swapaxes(0,2), (n, dim), order='F')

def multiplicities(spline):
    """  Returns the multiplicity of the knots at all knot values as a 2D array for
    all parametric directions, for all knots """
    return [[spline.order(d) - spline.bases[d].continuity(k) - 1 for k in spline.knots(d)] for d in range(spline.pardim)]

def degree(spline):
    """ Returns polynomial degree (splipy order - 1) for all parametric directions """
    return [p-1 for p in spline.order()]

