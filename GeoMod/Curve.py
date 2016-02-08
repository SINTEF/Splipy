# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis
from GeoMod.SplineObject import SplineObject
from GeoMod.Utils import ensure_listlike
from bisect import bisect_left
import numpy as np

__all__ = ['Curve']


class Curve(SplineObject):
    """Curve()

    Represents a curve: an object with a one-dimensional parameter space."""

    def __init__(self, basis=None, controlpoints=None, rational=False):
        """__init__([basis=None], [controlpoints=None], [rational=False])

        Construct a curve with the given basis and control points.

        The default is to create a linear one-element mapping from (0,1) to the
        unit interval.

        :param BSplineBasis basis: The underlying B-Spline basis
        :param array-like controlpoints: An *n* × *d* matrix of control points
        :param bool rational: Whether the curve is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
            which is the last coordinate)
        """
        super(Curve, self).__init__([basis], controlpoints, rational)

    def evaluate_tangent(self, t):
        """Evaluate the tangent of the curve at given parametric values.

        This is equivalent to :func:`GeoMod.Curve.evaluate_derivative` with the
        default value of *d* = 1.

        :param t: Parametric coordinate(s) in which to evaluate
        :type t: float or [float]
        :return: Tangent matrix *X[i,j]* of component *xj'(t)* evaluated at *t(i)*
        :rtype: numpy.array
        """
        return self.evaluate_derivative(t, d=1)

    def evaluate_derivative(self, t, d=1):
        """evaluate_derivative(u, [d=1])

        Evaluate the derivative of the curve at the given parametric values.

        This function returns an *n* × *dim* array, where *n* is the number of
        evaluation points, and *dim* is the physical dimension of the curve.

        If there is only one evaluation point, a vector of length *dim* is
        returned instead.

        :param u: Parametric coordinates in which to evaluate
        :type u: float or [float]
        :param int d: Number of derivatives to compute
        :return: Derivative array
        :rtype: numpy.array
        """
        if not self.rational or d != 2:
            return super(Curve, self).evaluate_derivative(t, d=d)

        t = ensure_listlike(t)
        dN = self.bases[0].evaluate(t, d)
        result = np.array(dN * self.controlpoints)

        d2 = result
        d1 = np.array(self.bases[0].evaluate(t, 1) * self.controlpoints)
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

    def raise_order(self, amount):
        """Raise the order of the curve.

        :param int amount: Number of degrees to increase
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
        N_new = newBasis.evaluate(interpolation_pts_t)
        interpolation_pts_x = N_old * self.controlpoints  # projective interpolation points (x,y,z,w)

        # solve the interpolation problem
        self.controlpoints = np.linalg.solve(N_new, interpolation_pts_x)
        self.bases = [newBasis]

    def append(self, curve):
        """Extend the curve by merging another curve to the end of it.

        The curves are glued together in a C0 fashion with enough repeated
        knots. The function assumes that the end of this curve perfectly
        matches the start of the input curve.

        :param Curve curve: Another curve
        :raises RuntimeError: If either curve is periodic
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
        old_knot = self.knots(with_multiplicities=True)
        add_knot = extending_curve.knots(with_multiplicities=True)
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
        self.basis = [BSplineBasis(p, new_knot)]
        self.controlpoints = new_controlpoints

        return self

    def continuity(self, knot):
        """Get the parametric continuity of the curve at a given point. Will
        return p-1-m, where m is the knot multiplicity and inf between knots"""
        return self.bases[0].continuity(knot)

    def split(self, knots):
        """Split a curve into two or more separate representations with C0
        continuity between them.

        :param knots: The splitting points
        :type knots: float or [float]
        :return: The new curves
        :rtype: [Curve]
        """
        # for single-value input, wrap it into a list
        knots = ensure_listlike(knots)

        p = self.order(0)
        results = []
        splitting_curve = self.clone()
        # insert knots to produce C{-1} at all splitting points
        for k in knots:
            continuity = splitting_curve.continuity(k)
            if continuity == np.inf:
                continuity = p - 1
            splitting_curve.insert_knot([k] * (continuity + 1))

        b = splitting_curve.bases[0]
        if b.periodic > -1:
            mu = bisect_left(b.knots, knots[0])
            b.roll(mu)
            splitting_curve.controlpoints = np.roll(splitting_curve.controlpoints, -mu, 0)
            b.knots = b.knots[:-b.periodic-1]
            b.periodic = -1
            if len(knots) > 1:
                return splitting_curve.split(knots[1:])
            else:
                return splitting_curve


        # everything is available now, just have to find the right index range
        # in the knot vector and controlpoints to store in each separate curve
        # piece
        last_cp_i = 0
        last_knot_i = 0
        for k in knots:
            mu = bisect_left(splitting_curve.bases[0].knots, k)
            n_cp = mu - last_knot_i
            basis = BSplineBasis(p, splitting_curve.bases[0].knots[last_knot_i:mu + p])
            controlpoints = splitting_curve.controlpoints[last_cp_i:last_cp_i + n_cp, :]

            results.append(Curve(basis, controlpoints, self.rational))
            last_knot_i = mu
            last_cp_i += n_cp
        # with n splitting points, we're getting n+1 pieces. Add the final one:
        basis = BSplineBasis(p, splitting_curve.bases[0].knots[last_knot_i:])
        controlpoints = splitting_curve.controlpoints[last_cp_i:, :]
        results.append(Curve(basis, controlpoints, self.rational))

        return results

    def rebuild(self, p, n):
        """Creates an approximation to this curve by resampling it using a
        uniform knot vector of order *p* with *n* control points.

        :param int p: Polynomial discretization order
        :param int n: Number of control points
        :return: A new approximate curve
        :rtype: Curve
        """
        # establish uniform open knot vector
        knot = [0] * p + range(1, n - p + 1) + [n - p + 1] * p
        basis = BSplineBasis(p, knot)
        # set parametric range of the new basis to be the same as the old one
        basis.normalize()
        t0 = self.bases[0].start()
        t1 = self.bases[0].end()
        basis *= (t1 - t0)
        basis += t0
        # fetch evaluation points and solve interpolation problem
        t = basis.greville()
        N = basis.evaluate(t)
        controlpoints = np.linalg.solve(N, self.evaluate(t))

        # return new resampled curve
        return Curve(basis, controlpoints)

    def write_g2(self, outfile):
        """Write the curve in GoTools format.

        :param file-like outfile: The file to write to
        """
        outfile.write('100 1 0 0\n')  # surface header, gotools version 1.0.0
        outfile.write('%i %i\n' % (self.dimension, int(self.rational)))
        self.bases[0].write_g2(outfile)

        (n1, n2) = self.controlpoints.shape
        for i in range(n1) + range(self.bases[0].periodic + 1):
            for j in range(n2):
                outfile.write('%f ' % self.controlpoints[i, j])
            outfile.write('\n')

    def __repr__(self):
        return str(self.bases[0]) + '\n' + str(self.controlpoints)
