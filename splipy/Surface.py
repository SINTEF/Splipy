# -*- coding: utf-8 -*-

from splipy import BSplineBasis, Curve, SplineObject
from splipy.utils import ensure_listlike, check_direction, sections
from bisect import bisect_left
from itertools import chain
import numpy as np

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

    def normal(self, u, v, above=(True,True)):
        """  Evaluate the normal of the surface at given parametric values.

        This is equal to the cross-product between tangents. The return value
        is normalized.

        :param u: Parametric coordinate(s) in the first direction
        :type u: float or [float]
        :param v: Parametric coordinate(s) in the second direction
        :type v: float or [float]
        :param (bool) above: Evaluation in the limit from above
        :return: Normal array *X[i,j,k]* of component *xj* evaluated at *(u[i], v[j])*
        :rtype: numpy.array
        :raises RuntimeError: If the physical dimension is not 2 or 3
        """
        if self.dimension == 2:
            try:
                result = np.zeros((len(u), len(v), 3))
                result[:, :, 2] = 1
                return result
            except TypeError:  # single valued input u, fails on len(u)
                return np.array([0, 0, 1])
        elif self.dimension == 3:
            # fetch the tangent vectors
            (du, dv) = self.tangent(u, v, above=above)

            # compute normals
            normals = np.cross(du,dv)

            # normalize output
            if len(du.shape) == 1:
                return normals / np.linalg.norm(normals)
            magnitude = np.apply_along_axis(np.linalg.norm, -1, normals)
            magnitude = magnitude.reshape(magnitude.shape + (1,))
            return normals / magnitude
        else:
            raise RuntimeError('Normal evaluation only defined for 2D and 3D geometries')

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
        C    = np.matrix(np.identity(self.shape[direction]))
        for i in range(mult):
            C = b.insert_knot(knot) * C

        # at this point we have a C0 basis, find the right interpolating index
        i  = max(bisect_left(b.knots, knot) - 1,0)

        # compute the controlpoints and return Curve
        cp = np.tensordot(C[i,:], self.controlpoints, axes=(1, direction))
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
