from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar, cast

import numpy as np

from .basis import BSplineBasis
from .splineobject import SplineObject
from .utils import ensure_listlike, sections

if TYPE_CHECKING:
    from collections.abc import Sequence

    from .curve import Curve
    from .surface import Surface
    from .typing import ArrayLike, FloatArray

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

    def volume(self) -> float:
        """Computes the volume of the object in geometric space"""

        w1: FloatArray
        w2: FloatArray
        w3: FloatArray

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

        return float(np.abs(J).dot(w3).dot(w2).dot(w1))

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

    get_derivative_volume = SplineObject.get_derivative_spline
