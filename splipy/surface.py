# -*- coding: utf-8 -*-

from bisect import bisect_left
from itertools import chain

import numpy as np

from .basis import BSplineBasis
from .curve import Curve
from .splineobject import SplineObject, evaluate
from .utils import is_singleton, ensure_listlike, check_direction, sections

__all__ = ['Surface']


class Surface(SplineObject):
    """Surface()

    Represents a surface: an object with a two-dimensional parameter space."""

    _intended_pardim = 2

    def __init__(self, basis1=None, basis2=None, controlpoints=None, rational=False, **kwargs):
        """  Construct a surface with the given basis and control points.

        The default is to create a linear one-element mapping from and to the
        unit square.

        :param BSplineBasis basis1: The basis of the first parameter direction
        :param BSplineBasis basis2: The basis of the second parameter direction
        :param array-like controlpoints: An *n1* × *n2* × *d* matrix of control points
        :param bool rational: Whether the surface is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
            which is the last coordinate)
        """
        super(Surface, self).__init__([basis1, basis2], controlpoints, rational, **kwargs)

    def normal(self, u, v, above=(True,True), tensor=True):
        """  Evaluate the normal of the surface at given parametric values.

        This is equal to the cross-product between tangents. The return value
        is normalized.

        If *tensor* is true, evaluation will take place on a tensor product
        grid, i.e. it will return an *n1* × *n2* × ... × *dim* array, where
        *ni* is the number of evaluation points in direction *i*, and *dim* is
        the physical dimension of the object.

        If *tensor* is false, there must be an equal number *n* of evaluation
        points in all directions, and the return value will be an *n* × *dim*
        array.

        :param u: Parametric coordinate(s) in the first direction
        :type u: float or [float]
        :param v: Parametric coordinate(s) in the second direction
        :type v: float or [float]
        :param (bool) above: Evaluation in the limit from above
        :param tensor: Whether to evaluate on a tensor product grid
        :type tensor: bool
        :return: Normal array *X[i,j,k]* of component *xk* evaluated at *(u[i], v[j])*
        :rtype: numpy.array
        :raises RuntimeError: If the physical dimension is not 2 or 3
        """
        if not tensor and len(u) != len(v):
            raise ValueError('Parametes must have same length')

        if self.dimension == 2:
            try:
                shape = (len(u), len(v), 3) if tensor else (len(u), 3)
                result = np.zeros(shape)
                result[..., 2] = 1
                return result
            except TypeError:  # single valued input u, fails on len(u)
                return np.array([0, 0, 1])
        elif self.dimension == 3:
            # fetch the tangent vectors
            (du, dv) = self.tangent(u, v, above=above, tensor=tensor)

            # compute normals
            normals = np.cross(du,dv)

            # normalize output
            if len(du.shape) == 1:
                return normals / np.linalg.norm(normals)
            magnitude = np.linalg.norm( normals, axis=-1)
            magnitude = magnitude.reshape(magnitude.shape + (1,))
            return normals / magnitude
        else:
            raise RuntimeError('Normal evaluation only defined for 2D and 3D geometries')


    def derivative(self, u, v, d=(1,1), above=True, tensor=True):
        """  Evaluate the derivative of the surface at the given parametric values.

        This function returns an *n* × *m* x *dim* array, where *n* is the number of
        evaluation points in u, *m* is the number of evaluation points in v, and
        *dim* is the physical dimension of the curve.

        If there is only one evaluation point, a vector of length *dim* is
        returned instead.

        :param u: Parametric coordinate(s) in the first direction
        :type u: float or [float]
        :param v: Parametric coordinate(s) in the second direction
        :type v: float or [float]
        :param int d: Number of derivatives to compute in each direction
        :type d: [int]
        :param (bool) above: Evaluation in the limit from above
        :param tensor: Whether to evaluate on a tensor product grid
        :type tensor: bool
        :return: Derivative array *X[i,j,k]* of component *xk* evaluated at *(u[i], v[j])*
        :rtype: numpy.array
        """

        squeeze = all(is_singleton(t) for t in [u,v])
        derivs = ensure_listlike(d, self.pardim)
        if not self.rational or np.sum(derivs) < 2 or np.sum(derivs) > 3:
            return super(Surface, self).derivative(u,v, d=derivs, above=above, tensor=tensor)

        u = ensure_listlike(u)
        v = ensure_listlike(v)
        result = np.zeros((len(u), len(v), self.dimension))
        # dNus = [self.bases[0].evaluate(u, d, above) for d in range(derivs[0]+1)]
        # dNvs = [self.bases[1].evaluate(v, d, above) for d in range(derivs[1]+1)]
        dNus = [self.bases[0].evaluate(u, d, above) for d in range(np.sum(derivs)+1)]
        dNvs = [self.bases[1].evaluate(v, d, above) for d in range(np.sum(derivs)+1)]

        d0ud0v = evaluate([dNus[0], dNvs[0]], self.controlpoints, tensor)
        d1ud0v = evaluate([dNus[1], dNvs[0]], self.controlpoints, tensor)
        d0ud1v = evaluate([dNus[0], dNvs[1]], self.controlpoints, tensor)
        d1ud1v = evaluate([dNus[1], dNvs[1]], self.controlpoints, tensor)
        d2ud0v = evaluate([dNus[2], dNvs[0]], self.controlpoints, tensor)
        d0ud2v = evaluate([dNus[0], dNvs[2]], self.controlpoints, tensor)
        W     = d0ud0v[:,:,-1]
        dWdu  = d1ud0v[:,:,-1]
        dWdv  = d0ud1v[:,:,-1]
        d2Wduv= d1ud1v[:,:,-1]
        d2Wdu = d2ud0v[:,:,-1]
        d2Wdv = d0ud2v[:,:,-1]

        for i in range(self.dimension):
            H1   = d1ud0v[:,:,i] * W - d0ud0v[:,:,i] * dWdu
            H2   = d0ud1v[:,:,i] * W - d0ud0v[:,:,i] * dWdv
            dH1du= d2ud0v[:,:,i] * W - d0ud0v[:,:,i] * d2Wdu
            dH1dv= d1ud1v[:,:,i] * W + d1ud0v[:,:,i] * dWdv - d0ud1v[:,:,i] * dWdu - d0ud0v[:,:,i] * d2Wduv
            dH2du= d1ud1v[:,:,i] * W + d0ud1v[:,:,i] * dWdu - d1ud0v[:,:,i] * dWdv - d0ud0v[:,:,i] * d2Wduv
            dH2dv= d0ud2v[:,:,i] * W - d0ud0v[:,:,i] * d2Wdv
            G1   = dH1du*W - 2*H1*dWdu
            G2   = dH2dv*W - 2*H2*dWdv
            if derivs == (1,0):
                result[:,:,i] = H1 / W/W
            elif derivs == (0,1):
                result[:,:,i] = H2 / W/W
            elif derivs == (1,1):
                result[:,:,i] = (dH1dv*W - 2*H1*dWdv) /W/W/W
            elif derivs == (2,0):
                result[:,:,i] = G1 /W/W/W
            elif derivs == (0,2):
                result[:,:,i] = G2 /W/W/W
            if np.sum(derivs) > 2:
                d2ud1v = evaluate([dNus[2], dNvs[1]], self.controlpoints, tensor)
                d1ud2v = evaluate([dNus[1], dNvs[2]], self.controlpoints, tensor)
                d3ud0v = evaluate([dNus[3], dNvs[0]], self.controlpoints, tensor)
                d0ud3v = evaluate([dNus[0], dNvs[3]], self.controlpoints, tensor)
                d3Wdu   = d3ud0v[:,:,-1]
                d3Wdv   = d0ud3v[:,:,-1]
                d3Wduuv = d2ud1v[:,:,-1]
                d3Wduvv = d1ud2v[:,:,-1]
                d2H1du  = d3ud0v[:,:,i]*W + d2ud0v[:,:,i]*dWdu - d1ud0v[:,:,i]*d2Wdu - d0ud0v[:,:,i]*d3Wdu
                d2H1duv = d2ud1v[:,:,i]*W + d2ud0v[:,:,i]*dWdv - d0ud1v[:,:,i]*d2Wdu - d0ud0v[:,:,i]*d3Wduuv
                d2H2dv  = d0ud3v[:,:,i]*W + d0ud2v[:,:,i]*dWdv - d0ud1v[:,:,i]*d2Wdv - d0ud0v[:,:,i]*d3Wdv
                d2H2duv = d1ud2v[:,:,i]*W + d0ud2v[:,:,i]*dWdu - d1ud0v[:,:,i]*d2Wdv - d0ud0v[:,:,i]*d3Wduvv
                dG1du   = d2H1du *W + dH1du*dWdu - 2*dH1du*dWdu - 2*H1*d2Wdu
                dG1dv   = d2H1duv*W + dH1du*dWdv - 2*dH1dv*dWdu - 2*H1*d2Wduv
                dG2du   = d2H2duv*W + dH2dv*dWdu - 2*dH2du*dWdv - 2*H2*d2Wduv
                dG2dv   = d2H2dv *W + dH2dv*dWdv - 2*dH2dv*dWdv - 2*H2*d2Wdv

                if derivs == (3,0):
                    result[:,:,i] = (dG1du*W -3*G1*dWdu) /W/W/W/W
                elif derivs == (0,3):
                    result[:,:,i] = (dG2dv*W -3*G2*dWdv) /W/W/W/W
                elif derivs == (2,1):
                    result[:,:,i] = (dG1dv*W -3*G1*dWdv) /W/W/W/W
                elif derivs == (1,2):
                    result[:,:,i] = (dG2du*W -3*G2*dWdu) /W/W/W/W

        # Squeeze the singleton dimensions if we only have one point
        if squeeze:
            result = result.reshape(self.dimension)

        return result


    def area(self):
        """ Computes the area of the surface in geometric space """
        # fetch integration points
        (x1,w1) = np.polynomial.legendre.leggauss(self.order(0)+1)
        (x2,w2) = np.polynomial.legendre.leggauss(self.order(1)+1)
        # map points to parametric coordinates (and update the weights)
        (knots1,knots2) = self.knots()
        u  = np.array([ (x1+1)/2*(t1-t0)+t0 for t0,t1 in zip(knots1[:-1], knots1[1:]) ])
        w1 = np.array([     w1/2*(t1-t0)    for t0,t1 in zip(knots1[:-1], knots1[1:]) ])
        v  = np.array([ (x2+1)/2*(t1-t0)+t0 for t0,t1 in zip(knots2[:-1], knots2[1:]) ])
        w2 = np.array([     w2/2*(t1-t0)    for t0,t1 in zip(knots2[:-1], knots2[1:]) ])

        # wrap everything to vectors
        u = np.ndarray.flatten(u)
        v = np.ndarray.flatten(v)
        w1 = np.ndarray.flatten(w1)
        w2 = np.ndarray.flatten(w2)

        # compute all quantities of interest (i.e. the jacobian)
        du = self.derivative(u,v, d=(1,0))
        dv = self.derivative(u,v, d=(0,1))
        J  = np.cross(du,dv)

        if self.dimension == 3:
            J = np.sqrt(np.sum(J**2, axis=2))
        else:
            J = np.abs(J)
        return w1.dot(J).dot(w2)

    def edges(self):
        """Return the four edge curves in (parametric) order: umin, umax, vmin, vmax

        :return: Edge curves
        :rtype: (Curve)
        """
        return tuple(self.section(*args) for args in sections(2, 1))

    def const_par_curve(self, knot, direction):
        """  Get a Curve representation of the parametric line of some constant
        knot value.
        :param float knot: The constant knot value to sample the surface
        :param int direction: The parametric direction for the constant value
        :return: curve on this surface
        :rtype: Curve
        """
        direction = check_direction(direction, 2)

        # clone basis since we need to augment this by knot insertion
        b    = self.bases[direction].clone()

        # compute mapping matrix C which is the knotinsertion operator
        mult = min(b.continuity(knot), b.order-1)
        C    = np.identity(self.shape[direction])
        for i in range(mult):
            C = b.insert_knot(knot) @ C

        # at this point we have a C0 basis, find the right interpolating index
        i  = max(bisect_left(b.knots, knot) - 1,0)

        # compute the controlpoints and return Curve
        cp = np.tensordot(C[i,:], self.controlpoints, axes=(0, direction))
        return Curve(self.bases[1-direction], cp, self.rational)

    def rebuild(self, p, n):
        """  Creates an approximation to this surface by resampling it using
        uniform knot vectors of order *p* with *n* control points.

        :param (int) p: Tuple of polynomial discretization order in each direction
        :param (int) n: Tuple of number of control points in each direction
        :return: A new approximate surface
        :rtype: Surface
        """
        p = ensure_listlike(p, dups=2)
        n = ensure_listlike(n, dups=2)

        old_basis = self.bases
        basis = []
        u = []
        N = []
        # establish uniform open knot vectors
        for i in range(2):
            knot = [0] * p[i] + list(range(1, n[i] - p[i] + 1)) + [n[i] - p[i] + 1] * p[i]
            basis.append(BSplineBasis(p[i], knot))

            # make these span the same parametric domain as the old ones
            basis[i].normalize()
            t0 = old_basis[i].start()
            t1 = old_basis[i].end()
            basis[i] *= (t1 - t0)
            basis[i] += t0

            # fetch evaluation points and evaluate basis functions
            u.append(basis[i].greville())
            N.append(basis[i].evaluate(u[i]))

        # find interpolation points as evaluation of existing surface
        x = self.evaluate(u[0], u[1])

        # solve interpolation problem
        cp = np.tensordot(np.linalg.inv(N[1]), x, axes=(1, 1))
        cp = np.tensordot(np.linalg.inv(N[0]), cp, axes=(1, 1))

        # re-order controlpoints so they match up with Surface constructor
        cp = cp.transpose((1, 0, 2))
        cp = cp.reshape(n[0] * n[1], cp.shape[2])

        # return new resampled curve
        return Surface(basis[0], basis[1], cp)

    def __repr__(self):
        result = str(self.bases[0]) + '\n' + str(self.bases[1]) + '\n'
        # print legacy controlpoint enumeration
        n1, n2, n3 = self.controlpoints.shape
        for j in range(n2):
            for i in range(n1):
                result += str(self.controlpoints[i, j, :]) + '\n'
        return result

    get_derivative_surface = SplineObject.get_derivative_spline
