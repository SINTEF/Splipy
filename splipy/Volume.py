# -*- coding: utf-8 -*-

from splipy import BSplineBasis, SplineObject
from splipy.utils import ensure_listlike, check_direction, sections
from itertools import chain
import numpy as np

__all__ = ['Volume']


class Volume(SplineObject):
    """Volume()

    Represents a volume: an object with a three-dimensional parameter space."""

    _intended_pardim = 3

    def __init__(self, basis1=None, basis2=None, basis3=None, controlpoints=None, rational=False, **kwargs):
        """  Construct a volume with the given basis and control points.

        The default is to create a linear one-element mapping from and to the
        unit cube.

        :param BSplineBasis basis1: The basis of the first parameter direction
        :param BSplineBasis basis2: The basis of the second parameter direction
        :param BSplineBasis basis3: The basis of the third parameter direction
        :param array-like controlpoints: An *n1* × *n2* × *n3* × *d* matrix of
            control points
        :param bool rational: Whether the volume is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
            which is the last coordinate)
        """
        super(Volume, self).__init__([basis1, basis2, basis3], controlpoints, rational, **kwargs)

    def edges(self):
        """Return the twelve edges of this volume in order:

        - umin, vmin
        - umax, vmin
        - umin, vmax
        - umax, vmax
        - umin, wmin
        - umax, wmin
        - umin, wmax
        - umax, wmax
        - vmin, wmin
        - vmax, wmin
        - vmin, wmax
        - vmax, wmax

        :return: Edges
        :rtype: (Curve)
        """
        return tuple(self.section(*args) for args in sections(3, 1))

    def faces(self):
        """Return the six faces of this volume in order: umin, umax, vmin, vmax, wmin, wmax.

        :return: Boundary faces
        :rtype: (Surface)
        """
        return tuple(self.section(*args) for args in sections(3, 2))

    def volume(self):
        """ Computes the volume of the object in geometric space """
        # fetch integration points
        (x1,w1) = np.polynomial.legendre.leggauss(self.order(0)+1)
        (x2,w2) = np.polynomial.legendre.leggauss(self.order(1)+1)
        (x3,w3) = np.polynomial.legendre.leggauss(self.order(2)+1)
        # map points to parametric coordinates (and update the weights)
        (knots1,knots2,knots3) = self.knots()
        u  = np.array([ (x1+1)/2*(t1-t0)+t0 for t0,t1 in zip(knots1[:-1], knots1[1:]) ])
        w1 = np.array([     w1/2*(t1-t0)    for t0,t1 in zip(knots1[:-1], knots1[1:]) ])
        v  = np.array([ (x2+1)/2*(t1-t0)+t0 for t0,t1 in zip(knots2[:-1], knots2[1:]) ])
        w2 = np.array([     w2/2*(t1-t0)    for t0,t1 in zip(knots2[:-1], knots2[1:]) ])
        w  = np.array([ (x3+1)/2*(t1-t0)+t0 for t0,t1 in zip(knots3[:-1], knots3[1:]) ])
        w3 = np.array([     w3/2*(t1-t0)    for t0,t1 in zip(knots3[:-1], knots3[1:]) ])

        # wrap everything to vectors
        u = np.ndarray.flatten(u)
        v = np.ndarray.flatten(v)
        w = np.ndarray.flatten(w)
        w1 = np.ndarray.flatten(w1)
        w2 = np.ndarray.flatten(w2)
        w3 = np.ndarray.flatten(w3)

        # compute all quantities of interest (i.e. the jacobian)
        du = self.derivative(u,v,w, d=(1,0,0))
        dv = self.derivative(u,v,w, d=(0,1,0))
        dw = self.derivative(u,v,w, d=(0,0,1))

        J  = du[:,:,:,0] * np.cross(dv[:,:,:,1:],   dw[:,:,:,1:]  ) -  \
             du[:,:,:,1] * np.cross(dv[:,:,:,0::2], dw[:,:,:,0::2]) +  \
             du[:,:,:,2] * np.cross(dv[:,:,:,:-1],  dw[:,:,:,:-1] )

        return np.abs(J).dot(w3).dot(w2).dot(w1)

    def rebuild(self, p, n):
        """  Creates an approximation to this volume by resampling it using
        uniform knot vectors of order *p* with *n* control points.

        :param (int) p: Tuple of polynomial discretization order in each direction
        :param (int) n: Tuple of number of control points in each direction
        :return: A new approximate volume
        :rtype: Volume
        """
        p = ensure_listlike(p, dups=3)
        n = ensure_listlike(n, dups=3)

        old_basis = [self.bases[0], self.bases[1], self.bases[2]]
        basis = []
        u = []
        N = []
        # establish uniform open knot vectors
        for i in range(3):
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

        # find interpolation points as evaluation of existing volume
        x = self.evaluate(u[0], u[1], u[2])

        # solve interpolation problem
        cp = np.tensordot(np.linalg.inv(N[2]), x, axes=(1, 2))
        cp = np.tensordot(np.linalg.inv(N[1]), cp, axes=(1, 2))
        cp = np.tensordot(np.linalg.inv(N[0]), cp, axes=(1, 2))

        # re-order controlpoints so they match up with Volume constructor
        cp = cp.transpose((2, 1, 0, 3))
        cp = cp.reshape(n[0] * n[1] * n[2], cp.shape[3])

        # return new resampled curve
        return Volume(basis[0], basis[1], basis[2], cp)

    def __repr__(self):
        result = str(self.bases[0]) + '\n'
        result += str(self.bases[1]) + '\n'
        result += str(self.bases[2]) + '\n'
        # print legacy controlpoint enumeration
        n1, n2, n3, dim = self.controlpoints.shape
        for k in range(n3):
            for j in range(n2):
                for i in range(n1):
                    result += str(self.controlpoints[i, j, k, :]) + '\n'
        return result

    get_derivative_volume = SplineObject.get_derivative_spline
