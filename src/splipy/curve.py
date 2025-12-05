from __future__ import annotations

from bisect import bisect_left, bisect_right
from typing import TYPE_CHECKING, ClassVar, Self, cast

import numpy as np
import scipy.sparse.linalg as splinalg

from . import state
from .basis import BSplineBasis
from .splineobject import SplineObject
from .utils import ensure_listlike, is_singleton

if TYPE_CHECKING:
    from collections.abc import Sequence

    from .typing import ArrayLike, Direction, FloatArray, Scalar

__all__ = ["Curve"]


class Curve(SplineObject):
    """Curve()

    Represents a curve: an object with a one-dimensional parameter space."""

    _intended_pardim: ClassVar[int] = 1

    def __init__(
        self,
        basis: BSplineBasis | None = None,
        controlpoints: ArrayLike | None = None,
        rational: bool = False,
        raw: bool = False,
    ) -> None:
        """Construct a curve with the given basis and control points.

        The default is to create a linear one-element mapping from (0,1) to the
        unit interval.

        :param BSplineBasis basis: The underlying B-Spline basis
        :param array-like controlpoints: An *n* × *d* matrix of control points
        :param bool rational: Whether the curve is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
            which is the last coordinate)
        """
        super().__init__([basis], controlpoints, rational, raw=raw)

    def evaluate(self, *params: ArrayLike | Scalar, tensor: bool = True) -> FloatArray:
        """Evaluate the object at given parametric values.

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
        params_np = [np.atleast_1d(np.asarray(t, dtype=np.float64)) for t in params]

        self._validate_domain(*params_np)

        # Evaluate the derivatives of the corresponding bases at the corresponding points
        # and build the result array
        N = self.bases[0].evaluate(params_np[0], sparse=True)
        result: FloatArray = N @ self.controlpoints

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

    def derivative(
        self,
        *params: ArrayLike | Scalar,
        d: int | Sequence[int] = 1,
        above: bool | Sequence[bool] = True,
        tensor: bool = True,
    ) -> FloatArray:
        """Evaluate the derivative of the curve at the given parametric values.

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
        d = ensure_listlike(d)[0]
        if not self.rational or d < 2 or d > 3:
            return super().derivative(params[0], d=d, above=above, tensor=tensor)

        above = ensure_listlike(above)[0]
        t = np.atleast_1d(np.asarray(params[0], dtype=np.float64))
        result: FloatArray = np.zeros((len(t), self.dimension), dtype=np.float64)

        d2 = np.array(self.bases[0].evaluate(t, 2, above) @ self.controlpoints)
        d1 = np.array(self.bases[0].evaluate(t, 1, above) @ self.controlpoints)
        d0 = np.array(self.bases[0].evaluate(t) @ self.controlpoints)
        W = d0[:, -1]  # W(t)
        W1 = d1[:, -1]  # W'(t)
        W2 = d2[:, -1]  # W''(t)
        if d == 2:
            for i in range(self.dimension):
                result[:, i] = (
                    (d2[:, i] * W * W - 2 * W1 * (d1[:, i] * W - d0[:, i] * W1) - d0[:, i] * W2 * W)
                    / W
                    / W
                    / W
                )
        elif d == 3:
            d3 = np.array(self.bases[0].evaluate(t, 3, above) @ self.controlpoints)
            W3 = d3[:, -1]  # W'''(t)
            for i in range(self.dimension):
                H = d1[:, i] * W - d0[:, i] * W1
                H1 = d2[:, i] * W - d0[:, i] * W2
                H2 = d3[:, i] * W + d2[:, i] * W1 - d1[:, i] * W2 - d0[:, i] * W3
                G = H1 * W - 2 * H * W1
                G1 = H2 * W - 2 * H * W2 - H1 * W1
                result[:, i] = (G1 * W - 3 * G * W1) / W / W / W / W

        if result.shape[0] == 1:  # in case of single value input t, return vector instead of matrix
            result = np.array(result[0, :]).reshape(self.dimension)

        return result

    def binormal(self, t: ArrayLike | Scalar, above: bool = True) -> FloatArray:
        """Evaluate the normalized binormal of the curve at the given parametric value(s).

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
            raise ValueError("Binormals require dimension = 3")

        # compute derivative
        dx = self.derivative(t, d=1, above=above)
        ddx = self.derivative(t, d=2, above=above)

        # in case of vanishing acceleration, colinear velocity and acceleration,
        # such as linear curves we guess an appropriate binbormal (multiple choice available)
        if len(dx.shape) == 1:
            if np.allclose(ddx, 0):
                ddx = np.array([1, 0, 0]) if np.allclose(dx[:2], 0) else np.array([0, 0, 1])
        else:
            for i in range(ddx.shape[0]):
                if np.allclose(ddx[i, :], 0):
                    if np.allclose(dx[i, :2], 0):  # dx = [0,0,1]
                        ddx[i, :] = np.array([1, 0, 0])
                    else:
                        ddx[i, :] = np.array([0, 0, 1])

        result = np.cross(dx, ddx)

        # in case of single value input t, return vector instead of matrix
        if len(dx.shape) == 1:
            return result / np.linalg.norm(result)

        # normalize
        magnitude: FloatArray = np.linalg.norm(result, axis=1)
        magnitude = magnitude.reshape(-1, 1)

        return result / magnitude

    def normal(self, t: ArrayLike | Scalar, above: bool = True) -> FloatArray:
        """Evaluate the normal of the curve at the given parametric value(s).

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
            raise RuntimeError("Normals require dimension = 3")

        # compute derivative
        T = self.tangent(t, above=above)
        B = self.binormal(t, above=above)

        return np.cross(B, T)

    def curvature(self, t: ArrayLike | Scalar, above: bool = True) -> FloatArray | float:
        """Evaluate the curvaure at specified point(s). The curvature is defined as

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

        w = v[..., 0] * a[..., 1] - v[..., 1] * a[..., 0] if self.dimension == 2 else np.cross(v, a)

        if len(v.shape) == 1:  # single evaluation point
            return float(np.linalg.norm(w) / np.linalg.norm(v))

        magnitude: FloatArray = np.abs(w) if self.dimension == 2 else np.linalg.norm(w, axis=-1)
        speed: FloatArray = np.linalg.norm(v, axis=-1)
        return magnitude / speed**3

    def torsion(self, t: ArrayLike | Scalar, above: bool = True) -> FloatArray | float:
        """Evaluate the torsion for a 3D curve at specified point(s). The torsion is defined as

        .. math:: \\frac{(\\boldsymbol{v}\\times \\boldsymbol{a})\\cdot
            (d\\boldsymbol{a}/dt)}{|\\boldsymbol{v}\\times \\boldsymbol{a}|^2}

        :param t: Parametric coordinates in which to evaluate
        :type t: float or [float]
        :param bool above: Evaluation in the limit from above
        :return: Derivative array
        :rtype: numpy.array
        """
        if self.dimension == 2:
            # no torsion for 2D curves
            t = np.atleast_1d(np.asarray(t, dtype=np.float64))
            return np.zeros(len(t), dtype=np.float64)
        if self.dimension != 3:
            raise ValueError("dimension must be 2 or 3")

        # compute derivative
        v = self.derivative(t, d=1, above=above)
        a = self.derivative(t, d=2, above=above)
        da = self.derivative(t, d=3, above=above)
        w = np.cross(v, a)

        # magnitude: float = float(np.linalg.norm(w))
        # magnitude: FloatArray = np.linalg.norm(w) if v.ndim == 1 else np.linalg.norm(w, axis=-1)

        if v.ndim == 1:  # single evaluation point
            # magnitude = np.linalg.norm(w)
            dot: FloatArray = np.dot(w, a)
            return dot / np.linalg.norm(w) ** 2

        magnitude: FloatArray = np.linalg.norm(w, axis=-1)
        nominator: FloatArray = np.einsum("ik,ik->i", w, da)
        return nominator / magnitude**2

    def raise_order(self, diff: int, *args: int, direction: Direction | None = None) -> Self:
        """Raise the polynomial order of the curve.

        :param int amount: Number of times to raise the order
        :return: self
        """
        if diff < 0:
            raise ValueError("Raise order requires a non-negative parameter")
        if diff == 0:
            return self

        # create the new basis
        newBasis = self.bases[0].raise_order(diff)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_t = newBasis.greville()  # parametric interpolation points (t)
        N_old = self.bases[0].evaluate(interpolation_pts_t)
        N_new = newBasis.evaluate(interpolation_pts_t, sparse=True)
        interpolation_pts_x = N_old @ self.controlpoints  # projective interpolation points (x,y,z,w)

        # solve the interpolation problem
        self.controlpoints = np.array(splinalg.spsolve(N_new, interpolation_pts_x), dtype=np.float64)
        self.bases = [newBasis]

        return self

    def append(self, curve: Curve) -> Self:
        """Extend the curve by merging another curve to the end of it.

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
            raise RuntimeError("Cannot append with periodic curves")

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
        new_knot[: len(old_knot) - 1] = old_knot[:-1]
        new_knot[len(old_knot) - 1 :] = add_knot[p:]

        # build new control points by merging the two existing matrices
        n1 = len(self)
        n2 = len(extending_curve)
        new_controlpoints = np.zeros((n1 + n2 - 1, self.dimension + self.rational), dtype=np.float64)
        new_controlpoints[:n1, :] = self.controlpoints
        new_controlpoints[n1:, :] = extending_curve.controlpoints[1:, :]

        # update basis and controlpoints
        self.bases = [BSplineBasis(p, new_knot)]
        self.controlpoints = new_controlpoints

        return self

    def continuity(self, knot: Scalar) -> int | float:
        """Get the parametric continuity of the curve at a given point. Will
        return p-1-m, where m is the knot multiplicity and inf between knots"""
        return self.bases[0].continuity(knot)

    def get_kinks(self) -> FloatArray:
        """Get the parametric coordinates at all points which have C0-
        continuity"""
        return np.array([k for k in self.knots(0) if self.continuity(k) < 1], dtype=np.float64)

    def length(self, t0: Scalar | None = None, t1: Scalar | None = None) -> float:
        """Computes the euclidian length of the curve in geometric space

        .. math:: \\int_{t_0}^{t_1}\\sqrt{x(t)^2 + y(t)^2 + z(t)^2} dt

        """
        knots = self.knots(0)
        quadrature_points = self.order(0) + 1
        if len(knots) == 2:
            quadrature_points *= 2
        (x, w) = np.polynomial.legendre.leggauss(quadrature_points)
        # keep only integration boundaries within given start (t0) and stop (t1) interval
        if t0 is not None:
            t0 = float(t0)
            i = bisect_left(knots, t0)
            knots = np.insert(knots, i, t0)
            knots = knots[i:]
        if t1 is not None:
            t1 = float(t1)
            i = bisect_right(knots, t1)
            knots = knots[:i]
            knots = np.insert(knots, i, t1)

        t = np.array(
            [(x + 1) / 2 * (t1 - t0) + t0 for t0, t1 in zip(knots[:-1], knots[1:])], dtype=np.float64
        ).flatten()
        w = np.array([w / 2 * (t1 - t0) for t0, t1 in zip(knots[:-1], knots[1:])], dtype=np.float64).flatten()
        dx = self.derivative(t)
        detJ = np.sqrt(np.sum(dx**2, axis=1))
        return float(np.dot(detJ, w))

    def rebuild(self, p: int, n: int) -> Curve:
        """Creates an approximation to this curve by resampling it using a
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
        basis *= t1 - t0
        basis += t0
        # fetch evaluation points and solve interpolation problem
        t = basis.greville()
        N = basis.evaluate(t, sparse=True)
        controlpoints = splinalg.spsolve(N, self.evaluate(t))

        # return new resampled curve
        return Curve(basis, controlpoints)

    def _closest_point_linear_curve(self, pt: ArrayLike) -> tuple[FloatArray, float]:
        """Computes the closest point on a linear curve to a given point.
        :param array-like pt: point to which the closest point on the curve is sought
        :return: the closest point on the curve and its parametric location
        :rtype: tuple(numpy.array, float)

        """
        knots = self.knots(0)
        mindist_squared = np.linalg.norm(pt - self.controlpoints[0]) ** 2
        t = knots[0]
        for p1, p0, t1, t0 in zip(self.controlpoints[1:], self.controlpoints[:-1], knots[2:-1], knots[1:-2]):
            b = p1 - p0
            a = pt - p0
            if 0 <= np.dot(a, b) <= np.dot(b, b):
                dist_squared = np.dot(a, a) - np.dot(a, b) ** 2 / np.dot(b, b)
                if dist_squared < mindist_squared:
                    mindist_squared = dist_squared
                    t = t0 + (np.dot(a, b) / np.dot(b, b)) * (t1 - t0)
                if np.dot(p1 - pt, p1 - p0) < mindist_squared:
                    mindist_squared = np.dot(p1 - pt, p1 - pt)
                    t = t1
        return self(t), t

    def closest_point(self, pt: ArrayLike, t0: Scalar = None) -> tuple[FloatArray, float]:
        """Computes the closest point on this curve to a given point. This is done by newton iteration
        and is using the state variables `controlpoint_absolute_tolerance`
        to determine convergence; but limited to 15 iterations.
        :param array-like pt: point to which the closest point on the curve is sought
        :param float t0: optional starting guess for the parametric location of the closest point
        :return: the closest point on the curve and its parametric location
        :rtype: tuple(numpy.array, float)

        """
        if self.order(0) == 1:
            return self._closest_point_linear_curve(pt)

        if t0 is None:
            dist = [np.linalg.norm(cp - pt) for cp in self.controlpoints]
        i = np.argmin(dist)
        t0 = self.bases[0].greville(i)
        t = t0
        iter = 0
        atol = state.controlpoint_absolute_tolerance
        F = np.dot(self(t) - pt, self.derivative(t))
        while np.abs(F) > atol:
            x = self(t)
            dx = self.derivative(t)
            ddx = self.derivative(t, d=2)
            e = x - pt
            dF = np.dot(dx, dx) + np.dot(e, ddx)
            dt = -F / dF
            # closest point outside curve definition. Return the closest endpoint
            if (t == self.bases[0].start() and dt < 0) or (t == self.bases[0].end() and dt > 0):
                break
            t += dt
            t = np.clip(t, self.bases[0].start(), self.bases[0].end())
            F = np.dot(self(t) - pt, self.derivative(t))
            iter += 1
            if iter > 15:
                # print(f'Warning: did not converge in 15 iterations, returning last {t=}, {F=}')
                break
        return self(t), t

    def error(self, target: Curve) -> tuple[FloatArray, float]:
        """Computes the L2 (squared and per knot span) and max error between
        this curve and a target curve

        .. math:: ||\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)||_{L^2(t_1,t_2)}^2 = \\int_{t_1}^{t_2}
            |\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)|^2 dt, \\quad \\forall \\;\\text{knots}\\;t_1 < t_2

        .. math:: ||\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)||_{L^\\infty} =
            \\max_t |\\boldsymbol{x_h}(t)-\\boldsymbol{x}(t)|

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
            print('|| e ||_L2  = ', np.sqrt(np.sum(err2)))
            print('|| e ||_max = ', maxerr)
        """
        knots = self.knots(0)
        (x, w) = np.polynomial.legendre.leggauss(self.order(0) + 1)
        err2: list[float] = []
        err_inf = 0.0
        for t0, t1 in zip(knots[:-1], knots[1:]):  # for all knot spans
            tg = (x + 1) / 2 * (t1 - t0) + t0  # evaluation points
            wg = w / 2 * (t1 - t0)  # integration weights
            error = self(tg) - target(tg)  # [x-xh, y-yh, z-zh]
            error = np.sum(error**2, axis=1)  # |x-xh|^2
            err2.append(np.dot(error, wg))  # integrate over domain
            err_inf = max(np.max(np.sqrt(error)), err_inf)
        return (np.array(err2, dtype=np.float64), err_inf)

    def __repr__(self) -> str:
        return str(self.bases[0]) + "\n" + str(self.controlpoints)

    def get_derivative_curve(self) -> Curve:
        return cast("Curve", super().get_derivative_spline(0))
