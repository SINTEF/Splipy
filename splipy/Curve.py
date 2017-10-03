# -*- coding: utf-8 -*-

from splipy import BSplineBasis
from splipy.SplineObject import SplineObject
from splipy.utils import ensure_listlike, is_singleton
from itertools import chain
from bisect import bisect_left, bisect_right
import numpy as np
import scipy.sparse.linalg as splinalg

__all__ = ['Curve']


class Curve(SplineObject):
    """Curve()

    Represents a curve: an object with a one-dimensional parameter space."""

    _intended_pardim = 1

    def __init__(self, basis=None, controlpoints=None, rational=False, **kwargs):
        """  Construct a curve with the given basis and control points.

        The default is to create a linear one-element mapping from (0,1) to the
        unit interval.

        :param BSplineBasis basis: The underlying B-Spline basis
        :param array-like controlpoints: An *n* × *d* matrix of control points
        :param bool rational: Whether the curve is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
            which is the last coordinate)
        """
        super(Curve, self).__init__([basis], controlpoints, rational, **kwargs)

    def evaluate(self, *params):
        """  Evaluate the object at given parametric values.

        This function returns an *n1* × *n2* × ... × *dim* array, where *ni* is
        the number of evaluation points in direction *i*, and *dim* is the
        physical dimension of the object.

        If there is only one evaluation point, a vector of length *dim* is
        returned instead.

        :param u,v,...: Parametric coordinates in which to evaluate
        :type u,v,...: float or [float]
        :return: Geometry coordinates
        :rtype: numpy.array
        """
        squeeze = is_singleton(params[0])
        params = [ensure_listlike(p) for p in params]

        self._validate_domain(*params)

        # Evaluate the derivatives of the corresponding bases at the corresponding points
        # and build the result array
        N = self.bases[0].evaluate(params[0], sparse=True)
        result = N*self.controlpoints

        # For rational objects, we divide out the weights, which are stored in the
        # last coordinate
        if self.rational:
            for i in range(self.dimension):
                result[..., i] /= result[..., -1]
            result = np.delete(result, self.dimension, -1)

        # Squeeze the singleton dimensions if we only have one point
        if squeeze:
            result = result.reshape(self.dimension)

        return result

    def derivative(self, t, d=1, above=True, tensor=True):
        """  Evaluate the derivative of the curve at the given parametric values.

        This function returns an *n* × *dim* array, where *n* is the number of
        evaluation points, and *dim* is the physical dimension of the curve.

        If there is only one evaluation point, a vector of length *dim* is
        returned instead.

        :param t: Parametric coordinates in which to evaluate
        :type t: float or [float]
        :param int d: Number of derivatives to compute
        :param bool above: Evaluation in the limit from above
        :param bool tensor: Not used in this method
        :return: Derivative array
        :rtype: numpy.array
        """
        if not self.rational or d != 2:
            return super(Curve, self).derivative(t, d=d, above=above, tensor=tensor)

        t = ensure_listlike(t)
        dN = self.bases[0].evaluate(t, d, above)
        result = np.array(dN * self.controlpoints)

        d2 = result
        d1 = np.array(self.bases[0].evaluate(t, 1, above) * self.controlpoints)
        d0 = np.array(self.bases[0].evaluate(t) * self.controlpoints)
        W = d0[:, -1]   # W(t)
        W1 = d1[:, -1]  # W'(t)
        W2 = d2[:, -1]  # W''(t)
        for i in range(self.dimension):
            result[:, i] = (d2[:, i] * W * W - 2 * W1 *
                            (d1[:, i] * W - d0[:, i] * W1) - d0[:, i] * W2 * W) / W / W / W

        result = np.delete(result, self.dimension, 1)  # remove the weight column

        if result.shape[0] == 1:  # in case of single value input t, return vector instead of matrix
            result = np.array(result[0, :]).reshape(self.dimension)

        return result

    def binormal(self, t, above=True):
        """  Evaluate the normalized binormal of the curve at the given parametric value(s).

        This function returns an *n* × 3 array, where *n* is the number of
        evaluation points.

        The binormal is computed as the normalized cross product between the
        velocity and acceleration of the curve.

        :param t: Parametric coordinates in which to evaluate
        :type t: float or [float]
        :param bool above: Evaluation in the limit from above
        :return: Derivative array
        :rtype: numpy.array
        """
        # error test input
        if self.dimension != 3:
            raise ValueError('Binormals require dimension = 3')

        # compute derivative
        dx    = self.derivative(t, d=1, above=above)
        ddx   = self.derivative(t, d=2, above=above)

        result = np.cross(dx, ddx)

        # in case of single value input t, return vector instead of matrix
        if len(dx.shape) == 1:
            return result / np.linalg.norm(result)

        # normalize
        magnitude = np.apply_along_axis(np.linalg.norm, 1, result)
        magnitude = magnitude.reshape(-1,1)

        return result / magnitude

    def normal(self, t, above=True):
        """  Evaluate the normal of the curve at the given parametric value(s).

        This function returns an *n* × 3 array, where *n* is the number of
        evaluation points.

        The normal is computed as the cross product between the binormal and
        the tangent of the curve.

        :param t: Parametric coordinates in which to evaluate
        :type t: float or [float]
        :param bool above: Evaluation in the limit from above
        :return: Derivative array
        :rtype: numpy.array
        """
        # error test input
        if self.dimension != 3:
            raise RuntimeError('Normals require dimension = 3')

        # compute derivative
        T = self.tangent(t,  above=above)
        B = self.binormal(t, above=above)

        return np.cross(B,T)

    def curvature(self, t, above=True):
        """  Evaluate the curvaure at specified point(s). The curvature is defined as

        .. math:: \\frac{|\\boldsymbol{v}\\times \\boldsymbol{a}|}{|\\boldsymbol{v}|^3}

        :param t: Parametric coordinates in which to evaluate
        :type t: float or [float]
        :param bool above: Evaluation in the limit from above
        :return: Derivative array
        :rtype: numpy.array
        """
        # compute derivative
        v = self.derivative(t, d=1, above=above)
        a = self.derivative(t, d=2, above=above)
        w = np.cross(v,a)

        if len(v.shape) == 1: # single evaluation point
            magnitude = np.linalg.norm(w)
            speed     = np.linalg.norm(v)
        else:                 # multiple evaluation points
            if self.dimension == 2:
                magnitude = w # for 2D-cases np.cross() outputs scalars
            else:             # for 3D, it is vectors
                magnitude = np.apply_along_axis(np.linalg.norm, -1, w)
            speed     = np.apply_along_axis(np.linalg.norm, -1, v)

        return magnitude / np.power(speed,3)

    def torsion(self, t, above=True):
        """  Evaluate the torsion for a 3D curve at specified point(s). The torsion is defined as

        .. math:: \\frac{(\\boldsymbol{v}\\times \\boldsymbol{a})\\cdot (d\\boldsymbol{a}/dt)}{|\\boldsymbol{v}\\times \\boldsymbol{a}|^2}

        :param t: Parametric coordinates in which to evaluate
        :type t: float or [float]
        :param bool above: Evaluation in the limit from above
        :return: Derivative array
        :rtype: numpy.array
        """
        if self.dimension == 2:
            # no torsion for 2D curves
            t = ensure_listlike(t)
            return np.zeros(len(t))
        elif self.dimension == 3:
            # only allow 3D curves
            pass
        else:
            raise ValueError('dimension must be 2 or 3')

        # compute derivative
        v  = self.derivative(t, d=1, above=above)
        a  = self.derivative(t, d=2, above=above)
        da = self.derivative(t, d=3, above=above)
        w = np.cross(v,a)

        if len(v.shape) == 1: # single evaluation point
            magnitude = np.linalg.norm(w)
            nominator = np.dot(w, a)
        else:                 # multiple evaluation points
            magnitude = np.apply_along_axis(np.linalg.norm, -1, w)
            nominator = np.array([np.dot(w1,da1) for (w1,da1) in zip(w, da)])

        return nominator / np.power(magnitude, 2)


    def raise_order(self, amount):
        """  Raise the polynomial order of the curve.

        :param int amount: Number of times to raise the order
        :return: self
        """
        if amount < 0:
            raise ValueError('Raise order requires a non-negative parameter')
        elif amount == 0:
            return

        # create the new basis
        newBasis = self.bases[0].raise_order(amount)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_t = newBasis.greville()  # parametric interpolation points (t)
        N_old = self.bases[0].evaluate(interpolation_pts_t)
        N_new = newBasis.evaluate(interpolation_pts_t, sparse=True)
        interpolation_pts_x = N_old * self.controlpoints  # projective interpolation points (x,y,z,w)

        # solve the interpolation problem
        self.controlpoints = np.array(splinalg.spsolve(N_new, interpolation_pts_x))
        self.bases = [newBasis]

        return self

    def append(self, curve):
        """  Extend the curve by merging another curve to the end of it.

        The curves are glued together in a C0 fashion with enough repeated
        knots. The function assumes that the end of this curve perfectly
        matches the start of the input curve.

        :param Curve curve: Another curve
        :raises RuntimeError: If either curve is periodic
        :return: self
        """
        # ASSUMPTION: open knot vectors

        # error test input
        if self.bases[0].periodic > -1 or curve.bases[0].periodic > -1:
            raise RuntimeError('Cannot append with periodic curves')

        # copy input curve so we don't change that one directly
        extending_curve = curve.clone()

        # make sure both are in the same space, and (if needed) have rational weights
        Curve.make_splines_compatible(self, extending_curve)

        # make sure both have the same discretization order
        p1 = self.order(0)
        p2 = extending_curve.order(0)
        if p1 < p2:
            self.raise_order(p2 - p1)
        else:
            extending_curve.raise_order(p1 - p2)
        p = max(p1, p2)

        # build new knot vector by merging the two existing ones
        old_knot = self.knots(direction=0, with_multiplicities=True)
        add_knot = extending_curve.knots(direction=0, with_multiplicities=True)
        # make sure that the new one starts where the old one stops
        add_knot -= add_knot[0]
        add_knot += old_knot[-1]
        new_knot = np.zeros(len(add_knot) + len(old_knot) - p - 1)
        new_knot[:len(old_knot) - 1] = old_knot[:-1]
        new_knot[len(old_knot) - 1:] = add_knot[p:]

        # build new control points by merging the two existing matrices
        n1 = len(self)
        n2 = len(extending_curve)
        new_controlpoints = np.zeros((n1 + n2 - 1, self.dimension + self.rational))
        new_controlpoints[:n1, :] = self.controlpoints
        new_controlpoints[n1:, :] = extending_curve.controlpoints[1:, :]

        # update basis and controlpoints
        self.bases = [BSplineBasis(p, new_knot)]
        self.controlpoints = new_controlpoints

        return self

    def continuity(self, knot):
        """  Get the parametric continuity of the curve at a given point. Will
        return p-1-m, where m is the knot multiplicity and inf between knots"""
        return self.bases[0].continuity(knot)

    def get_kinks(self):
        """  Get the parametric coordinates at all points which have C0-
        continuity"""
        return [k for k in self.knots(0) if self.continuity(k)<1]

    def length(self, t0=None, t1=None):
        """ Computes the euclidian length of the curve in geometric space

        .. math:: \\int_{t_0}^{t_1}\\sqrt{x(t)^2 + y(t)^2 + z(t)^2} dt

        """
        (x,w) = np.polynomial.legendre.leggauss(self.order(0)+1)
        knots = self.knots(0)
        # keep only integration boundaries within given start (t0) and stop (t1) interval
        if t0 is not None:
            i = bisect_left(knots, t0)
            knots = np.insert(knots, i, t0)
            knots = knots[i:]
        if t1 is not None:
            i = bisect_right(knots, t1)
            knots = knots[:i]
            knots = np.insert(knots, i, t1)

        t = np.array([ (x+1)/2*(t1-t0)+t0 for t0,t1 in zip(knots[:-1], knots[1:]) ])
        w = np.array([     w/2*(t1-t0)    for t0,t1 in zip(knots[:-1], knots[1:]) ])
        t = np.ndarray.flatten(t)
        w = np.ndarray.flatten(w)
        dx = self.derivative(t)
        detJ = np.sqrt(np.sum(dx**2, axis=1))
        return np.dot(detJ, w)

    def rebuild(self, p, n):
        """  Creates an approximation to this curve by resampling it using a
        uniform knot vector of order *p* with *n* control points.

        :param int p: Polynomial discretization order
        :param int n: Number of control points
        :return: A new approximate curve
        :rtype: Curve
        """
        # establish uniform open knot vector
        knot = [0] * p + list(range(1, n - p + 1)) + [n - p + 1] * p
        basis = BSplineBasis(p, knot)
        # set parametric range of the new basis to be the same as the old one
        basis.normalize()
        t0 = self.bases[0].start()
        t1 = self.bases[0].end()
        basis *= (t1 - t0)
        basis += t0
        # fetch evaluation points and solve interpolation problem
        t = basis.greville()
        N = basis.evaluate(t, sparse=True)
        controlpoints = splinalg.spsolve(N, self.evaluate(t))

        # return new resampled curve
        return Curve(basis, controlpoints)

    def error(self, target):
        """  Computes the L2 (squared and per knot span) and max error between
        this curve and a target curve

        .. math:: ||\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)||_{L^2(t_1,t_2)}^2 = \\int_{t_1}^{t_2}
            |\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)|^2 dt, \\quad \\forall \;\\text{knots}\;t_1 < t_2

        .. math:: ||\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)||_{L^\infty} = \\max_t |\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)|

        :param function target: callable function which takes as input a vector
            of evaluation points t and gives as output a matrix x where
            x[i,j] is component j evaluated at point t[i]
        :return: L2 error per knot-span and the maximum error
        :rtype:  tuple(list(float), float)

        Examples:

        .. code:: python

            import numpy as np

            def arclength_circle(t):
                return np.array( [np.cos(t), np.sin(t)] ).T

            crv = curve_factory.circle(r=1)
            (err2, maxerr) = crv.error(arclength_circle)
            print '|| e ||_L2  = ', np.sqrt(np.sum(err2))
            print '|| e ||_max = ', maxerr
        """
        knots = self.knots(0)
        (x,w) = np.polynomial.legendre.leggauss(self.order(0)+1)
        err2    = []
        err_inf = 0.0
        for t0,t1 in zip(knots[:-1], knots[1:]): # for all knot spans
            tg  = (x+1)/2*(t1-t0)+t0 # evaluation points
            wg  =   w  /2*(t1-t0)    # integration weights
            error = self(tg) - target(tg)    # [x-xh, y-yh, z-zh]
            error = np.sum(error**2, axis=1) # |x-xh|^2
            err2.append(np.dot(error, wg))   # integrate over domain
            err_inf = max(np.max(np.sqrt(error)), err_inf)
        return (err2, err_inf)

    def __repr__(self):
        return str(self.bases[0]) + '\n' + str(self.controlpoints)

    def get_derivative_curve(self):
        return super(Curve, self).get_derivative_spline(0)
