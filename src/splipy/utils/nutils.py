__doc__ = 'Implementation of convenience methods with respect to nutils integration.'

import numpy as np

from ..curve import Curve
from ..surface import Surface
from ..volume import Volume


def controlpoints(spline):
    """ Return controlpoints according to nutils ordering """
    n = len(spline)
    dim = spline.dimension
    if isinstance(spline, Curve):
        return np.reshape(spline[:,:]                  , (n, dim), order='F')
    elif isinstance(spline, Surface):
        return np.reshape(spline[:,:,:].swapaxes(0,1)  , (n, dim), order='F')
    elif isinstance(spline, Volume):
        return np.reshape(spline[:,:,:,:].swapaxes(0,2), (n, dim), order='F')
    raise RuntimeError('Non-spline argument detected')

def multiplicities(spline):
    """  Returns the multiplicity of the knots at all knot values as a 2D array for
    all parametric directions, for all knots """
    return [[spline.order(d) - spline.bases[d].continuity(k) - 1 for k in spline.knots(d)] for d in range(spline.pardim)]

def degree(spline):
    """ Returns polynomial degree (splipy order - 1) for all parametric directions """
    return [p-1 for p in spline.order()]

def splipy_to_nutils(spline):
    """ Returns nutils domain and geometry object for spline mapping given by the argument """
    from nutils import mesh, function
    domain, geom = mesh.rectilinear(spline.knots())
    cp    = controlpoints(spline)
    basis = domain.basis('spline', degree=degree(spline), knotmultiplicities=multiplicities(spline))
    geom  = function.matmat(basis, cp)
    #TODO: add correct behaviour for rational and/or periodic geometries
    return domain, geom

