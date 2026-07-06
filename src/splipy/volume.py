from __future__ import annotations

from bisect import bisect_left
from typing import TYPE_CHECKING, ClassVar, cast

import numpy as np

from .basis import BSplineBasis
from .splineobject import SplineObject
from .surface import Surface
from .utils import check_direction, ensure_listlike, sections

if TYPE_CHECKING:
    from collections.abc import Sequence

    from splipy.typing import Point, Scalar

    from .curve import Curve
    from .typing import ArrayLike, Direction

__all__ = ["Volume"]


class Volume(SplineObject):
    """Volume()

    Represents a volume: an object with a three-dimensional parameter space."""

    _intended_pardim: ClassVar[int] = 3

    def __init__(
        self,
        basis1: BSplineBasis | None = None,
        basis2: BSplineBasis | None = None,
        basis3: BSplineBasis | None = None,
        controlpoints: ArrayLike | None = None,
        rational: bool = False,
        raw: bool = False,
    ):
        """Construct a volume with the given basis and control points.

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
        super().__init__([basis1, basis2, basis3], controlpoints, rational, raw=raw)

    def edges(
        self,
    ) -> tuple[Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve]:
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
        return cast(
            "tuple[Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve, Curve]",
            tuple(self.section(*args) for args in sections(3, 1)),
        )

    def faces(
        self,
    ) -> tuple[
        Surface | None, Surface | None, Surface | None, Surface | None, Surface | None, Surface | None
    ]:
        """Return the six faces of this volume in order: umin, umax, vmin, vmax, wmin, wmax.

        :return: Boundary faces
        :rtype: (Surface)
        """
        boundary_faces: list[Surface | None] = cast(
            "list[Surface | None]",
            [self.section(*args) for args in sections(3, 2)],
        )
        for i, b in enumerate(self.bases):
            if b.periodic > -1:
                boundary_faces[2 * i] = None
                boundary_faces[2 * i + 1] = None
        return cast(
            "tuple[Surface | None, Surface | None, Surface | None"
            ", Surface | None, Surface | None, Surface | None]",
            tuple(boundary_faces),
        )

    def const_par_surface(self, knot: Scalar, direction: Direction) -> Surface:
        """Get a Surface representation of the parametric volume at some constant
        knot value.
        :param float knot: The constant knot value to sample the volume
        :param int direction: The parametric direction for the constant value
        :return: surface in this volume
        :rtype: Surface
        """
        direction = check_direction(direction, 3)

        # clone basis since we need to augment this by knot insertion
        b = self.bases[direction].clone()

        # snap to existing knot if close enough
        knot = b.snap_point(knot)

        # compute mapping matrix C which is the knot insertion operator
        mult = b.min_continuity(knot, b.order - 1)
        C = np.identity(self.shape[direction])
        for i in range(mult):
            C = b.insert_knot(knot) @ C

        # at this point we have a C0 basis, find the right interpolating index
        i = max(bisect_left(b.knots, knot) - 1, 0)

        # compute the controlpoints and return Surface
        cp = np.tensordot(C[i, :], self.controlpoints, axes=(0, direction))
        cp = cp.transpose(1, 0, 2)

        # return surface with the two remaining directions
        remaining_dirs = [d for d in range(3) if d != direction]
        return Surface(self.bases[remaining_dirs[0]], self.bases[remaining_dirs[1]], cp, self.rational)

    def volume(self) -> Scalar:
        """Computes the volume of the object in geometric space"""
        # fetch integration points
        (x1, w1) = np.polynomial.legendre.leggauss(self.order(0) + 1)
        (x2, w2) = np.polynomial.legendre.leggauss(self.order(1) + 1)
        (x3, w3) = np.polynomial.legendre.leggauss(self.order(2) + 1)

        # map points to parametric coordinates (and update the weights)
        (knots1, knots2, knots3) = self.knots()
        u = np.array([(x1 + 1) / 2 * (t1 - t0) + t0 for t0, t1 in zip(knots1[:-1], knots1[1:])])
        w1 = np.array([w1 / 2 * (t1 - t0) for t0, t1 in zip(knots1[:-1], knots1[1:])])
        v = np.array([(x2 + 1) / 2 * (t1 - t0) + t0 for t0, t1 in zip(knots2[:-1], knots2[1:])])
        w2 = np.array([w2 / 2 * (t1 - t0) for t0, t1 in zip(knots2[:-1], knots2[1:])])
        w = np.array([(x3 + 1) / 2 * (t1 - t0) + t0 for t0, t1 in zip(knots3[:-1], knots3[1:])])
        w3 = np.array([w3 / 2 * (t1 - t0) for t0, t1 in zip(knots3[:-1], knots3[1:])])

        # wrap everything to vectors
        u = np.ndarray.flatten(u)
        v = np.ndarray.flatten(v)
        w = np.ndarray.flatten(w)
        w1 = np.ndarray.flatten(w1)
        w2 = np.ndarray.flatten(w2)
        w3 = np.ndarray.flatten(w3)

        # compute all quantities of interest (i.e. the jacobian)
        du = self.derivative(u, v, w, d=(1, 0, 0))
        dv = self.derivative(u, v, w, d=(0, 1, 0))
        dw = self.derivative(u, v, w, d=(0, 0, 1))

        c1 = dv[..., 1] * dw[..., 2] - dv[..., 2] * dw[..., 1]
        c2 = dv[..., 0] * dw[..., 2] - dv[..., 2] * dw[..., 0]
        c3 = dv[..., 0] * dw[..., 1] - dv[..., 1] * dw[..., 0]

        J = du[:, :, :, 0] * c1 - du[:, :, :, 1] * c2 + du[:, :, :, 2] * c3

        return np.abs(J).dot(w3).dot(w2).dot(w1)  # type: ignore[no-any-return]

    def rebuild(self, p: int | Sequence[int], n: int | Sequence[int]) -> Volume:
        """Creates an approximation to this volume by resampling it using
        uniform knot vectors of order *p* with *n* control points.

        :param (int) p: Tuple of polynomial discretization order in each direction
        :param (int) n: Tuple of number of control points in each direction
        :return: A new approximate volume
        :rtype: Volume
        """
        ps = ensure_listlike(p, dups=3)
        ns = ensure_listlike(n, dups=3)

        old_basis = [self.bases[0], self.bases[1], self.bases[2]]
        basis = []
        u = []
        N = []
        # establish uniform open knot vectors
        for i in range(3):
            knot = [0] * ps[i] + list(range(1, ns[i] - ps[i] + 1)) + [ns[i] - ps[i] + 1] * ps[i]
            basis.append(BSplineBasis(ps[i], knot))

            # make these span the same parametric domain as the old ones
            basis[i].normalize()
            t0 = old_basis[i].start()
            t1 = old_basis[i].end()
            basis[i] *= t1 - t0
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
        cp = cp.reshape(ns[0] * ns[1] * ns[2], cp.shape[3])

        # return new resampled curve
        return Volume(basis[0], basis[1], basis[2], cp)

    def __repr__(self) -> str:
        result = str(self.bases[0]) + "\n"
        result += str(self.bases[1]) + "\n"
        result += str(self.bases[2]) + "\n"
        # print legacy controlpoint enumeration
        n1, n2, n3, dim = self.controlpoints.shape
        for k in range(n3):
            for j in range(n2):
                for i in range(n1):
                    result += str(self.controlpoints[i, j, k, :]) + "\n"
        return result

    def get_antiderivative_volume(self, direction: Direction, constant: Point | None = None) -> Volume:
        """Compute the antiderivative (integral) of the volume in a given parametric direction.

        The antiderivative is computed by inverting the derivative operator on
        the spline space in the specified parametric direction. The result is a
        new volume of order p+1 in that direction (where p is the current order)
        whose derivative in that direction equals this volume.

        The antiderivative is only unique up to an additive constant volume. By
        default, the constant is chosen such that the antiderivative evaluates to
        zero at the start of the parametric domain in the given direction. You can
        specify a different constant to shift the result.

        :param direction: The parametric direction to integrate in (0, 1, or 2; or 'u', 'v', 'w')
        :type direction: int or str
        :param array-like constant: Optional constant vector to add to the result.
            If not provided, defaults to zero (antiderivative is zero at parameter start).
            Must have the same dimension as the volume's physical space.
        :type constant: array-like or None
        :return: A new volume whose derivative in the given direction equals self
        :rtype: Volume
        :raises RuntimeError: If the volume is rational (not supported)

        Examples:

        .. code:: python

            import splipy as sp
            import numpy as np

            # Create a simple trilinear volume
            vol = sp.volume_factory.cube()

            # Compute antiderivative in u-direction
            integral_u = vol.get_antiderivative_volume('u')

            # The derivative of integral should equal the original volume
            u = np.linspace(0, 1, 11)
            v = np.linspace(0, 1, 11)
            w = np.linspace(0, 1, 11)
            diff = np.linalg.norm(integral_u.derivative(u, v, w, d=(1,0,0)) - vol(u, v, w))
            print(f"Error: {diff}")  # Should be near machine precision

        """
        return cast("Volume", super().get_antiderivative_spline(direction, constant))

    get_derivative_volume = SplineObject.get_derivative_spline
