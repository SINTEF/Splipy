from __future__ import annotations

__doc__ = 'Implementation of convenience methods with respect to nutils integration.'

from typing import TYPE_CHECKING

import numpy as np

from ..splineobject import SplineObject
from ..curve import Curve
from ..surface import Surface
from ..volume import Volume
from ..types import FArray

if TYPE_CHECKING:
    from nutils import mesh, function


def controlpoints(spline: SplineObject) -> FArray:
    """Return controlpoints according to nutils ordering."""
    n = len(spline)
    dim = spline.dimension
    if isinstance(spline, Curve):
        return np.reshape(spline[:,:]                  , (n, dim), order='F')
    elif isinstance(spline, Surface):
        return np.reshape(spline[:,:,:].swapaxes(0,1)  , (n, dim), order='F')
    elif isinstance(spline, Volume):
        return np.reshape(spline[:,:,:,:].swapaxes(0,2), (n, dim), order='F')
    raise RuntimeError('Non-spline argument detected')


def multiplicities(spline: SplineObject) -> list[list[int]]:
    """Return the multiplicity of the knots at all knot values as a 2D array for
    all parametric directions, for all knots.
    """
    return [[spline.order(d) - spline.bases[d].continuity(k) - 1 for k in spline.knots(d)] for d in range(spline.pardim)]


def degree(spline: SplineObject) -> list[int]:
    """ Returns polynomial degree (splipy order - 1) for all parametric directions """
    return [p-1 for p in spline.order()]


def splipy_to_nutils(spline: SplineObject) -> tuple[mesh.Domain, function.Function]:
    """ Returns nutils domain and geometry object for spline mapping given by the argument """
    from nutils import mesh, function
    domain, geom = mesh.rectilinear(spline.knots())
    cp    = controlpoints(spline)
    basis = domain.basis('spline', degree=degree(spline), knotmultiplicities=multiplicities(spline))
    geom  = function.matmat(basis, cp)
    #TODO: add correct behaviour for rational and/or periodic geometries
    return domain, geom
