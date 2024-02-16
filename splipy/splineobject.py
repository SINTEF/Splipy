from __future__ import annotations

import copy
from bisect import bisect_left
from itertools import product
from operator import attrgetter, methodcaller
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    Generic,
    Literal,
    MutableSequence,
    Optional,
    Protocol,
    Sequence,
    SupportsFloat,
    SupportsIndex,
    TypeVar,
    Union,
    cast,
    overload,
)

import numpy as np
from numpy.typing import NDArray
from typing_extensions import Self, Unpack

from .basis import BSplineBasis
from .utils import (
    check_direction,
    check_section,
    ensure_listlike,
    ensure_scalars,
    is_singleton,
    raise_order_1D,
    reshape,
    rotation_matrix,
    sections,
)

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from .types import Direction, FArray, Scalar, ScalarOrScalars, Scalars, SectionElt, SectionKwargs

__all__ = ["SplineObject"]


IPArray = NDArray[np.intp]


def _transpose_fix(pardim: int, direction: int) -> tuple[int, ...]:
    ret = list(range(1, pardim + 1))
    ret.insert(direction, 0)
    return tuple(ret)


def _evaluate(bases: Sequence[FArray], cps: FArray, tensor: bool = True) -> FArray:
    if tensor:
        idx = len(bases) - 1
        for N in bases[::-1]:
            cps = np.tensordot(N, cps, axes=(1, idx))
    else:
        cps = np.einsum("ij,j...->i...", bases[0], cps)
        for N in bases[1:]:
            cps = np.einsum("ij,ij...->i...", N, cps)
    return cps


# Standard constructor for SplineObjects
T = TypeVar("T", bound="SplineObject", covariant=True)


class Constructor(Protocol, Generic[T]):
    def __call__(
        self,
        bases: Sequence[BSplineBasis],
        controlpoints: Any = None,
        rational: bool = False,
        raw: bool = False,
    ) -> T:
        ...


class SplineObject:
    """Master class for spline objects with arbitrary dimensions.

    This class should be subclassed instead of used directly.

    All SplineObjects support basic arithmetic operators, which are interpreted
    as translation and scaling. In-place operators (e.g. ``+=``) mutate the
    object, while infix operators (e.g. ``+``) create new objects.
    """

    _intended_pardim: ClassVar[Optional[int]] = None

    dimension: int
    bases: list[BSplineBasis]
    controlpoints: FArray
    rational: bool

    @staticmethod
    def constructor(pardim: int) -> Constructor[SplineObject]:
        constructor = next(iter([c for c in SplineObject.__subclasses__() if c._intended_pardim == pardim]))

        def wrapped_constructor(
            bases: Sequence[BSplineBasis],
            controlpoints: Any = None,
            rational: bool = False,
            raw: bool = False,
        ) -> SplineObject:
            return constructor(*bases, controlpoints, rational=rational, raw=raw)  # type: ignore[arg-type, misc, call-arg]

        return wrapped_constructor

    @property
    def self_constructor(self) -> Constructor[Self]:
        return cast(Constructor[Self], SplineObject.constructor(self.pardim))

    def __init__(
        self,
        bases: Sequence[Optional[BSplineBasis]],
        controlpoints: Any = None,
        rational: bool = False,
        raw: bool = False,
    ) -> None:
        """Construct a spline object with the given bases and control points.

        The default is to create a linear one-element mapping from and to the
        unit (hyper)cube.

        :param [BSplineBasis] bases: The basis of each parameter direction
        :param array-like controlpoints: An *n1* × *n2* × ... × *d* matrix of
            control points
        :param bool rational: Whether the object is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
            which is the last coordinate)
        :param bool raw: If True, skip any control point reordering.
            (For internal use.)
        """

        self.bases = [(b.clone() if b else BSplineBasis()) for b in bases]

        if controlpoints is None:
            # `product' produces tuples in row-major format (the last input varies quickest)
            # We want them in column-major format, so we reverse the basis orders, and then
            # also reverse the output tuples
            cps = np.array(
                [c[::-1] for c in product(*(b.greville() for b in self.bases[::-1]))],
                dtype=float,
            )

            # Minimum two dimensions
            m, n = cps.shape
            if n == 1:
                zeros = np.zeros_like(cps, shape=(m, 1))
                cps = np.concatenate((cps, zeros), axis=1)

            # Add weight = 1 for identiy-mapping rational splines
            if rational:
                ones = np.ones_like(cps, shape=(m, 1))
                controlpoints = np.concatenate((cps, ones), axis=1)

            self.controlpoints = cps

        else:
            self.controlpoints = np.array(controlpoints, dtype=float)

        self.dimension = self.controlpoints.shape[-1] - rational
        self.rational = rational

        if not raw:
            shape = tuple(b.num_functions() for b in self.bases)
            ncomps = self.dimension + rational
            self.controlpoints = reshape(self.controlpoints, shape, order="F", ncomps=ncomps)

    def _validate_domain(self, *params: MutableSequence[float]) -> None:
        """Check whether the given evaluation parameters are valid.

        :raises ValueError: If the parameters are outside the domain
        """
        for b, p in zip(self.bases, params):
            b.snap(p)
            if b.periodic < 0 and (min(p) < b.start() or b.end() < max(p)):
                raise ValueError("Evaluation outside parametric domain")

    def evaluate(
        self,
        *params: ScalarOrScalars,
        tensor: bool = True,
    ) -> FArray:
        """Evaluate the object at given parametric values.

        If *tensor* is true, evaluation will take place on a tensor product
        grid, i.e. it will return an *n1* × *n2* × ... × *dim* array, where
        *ni* is the number of evaluation points in direction *i*, and *dim* is
        the physical dimension of the object.

        If *tensor* is false, there must be an equal number *n* of evaluation
        points in all directions, and the return value will be an *n* × *dim*
        array.

        If there is only one evaluation point, a vector of length *dim* is
        returned instead.

        :param u,v,...: Parametric coordinates in which to evaluate
        :type u,v,...: float or [float]
        :param tensor: Whether to evaluate on a tensor product grid
        :type tensor: bool
        :return: Geometry coordinates
        :rtype: numpy.array
        """
        squeeze = all(is_singleton(p) for p in params)
        params_list = [ensure_scalars(p) for p in params]

        if not tensor and len({len(p) for p in params_list}) != 1:
            raise ValueError("Parameters must have same length")

        self._validate_domain(*params_list)

        # Evaluate the corresponding bases at the corresponding points
        # and build the result array
        Ns = [b.evaluate(p) for b, p in zip(self.bases, params_list)]
        result = _evaluate(Ns, self.controlpoints, tensor)

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
        *params: ScalarOrScalars,
        d: Union[int, Sequence[int]] = 1,
        above: Union[bool, Sequence[bool]] = True,
        tensor: bool = True,
    ) -> FArray:
        """Evaluate the derivative of the object at the given parametric values.

        If *tensor* is true, evaluation will take place on a tensor product
        grid, i.e. it will return an *n1* × *n2* × ... × *dim* array, where
        *ni* is the number of evaluation points in direction *i*, and *dim* is
        the physical dimension of the object.

        If *tensor* is false, there must be an equal number *n* of evaluation
        points in all directions, and the return value will be an *n* × *dim*
        array.

        If there is only one evaluation point, a vector of length *dim* is
        returned instead.

        Examples:

        .. code:: python

           # Tangent of curve at single point
           curve.derivative(1.0)

           # Double derivative of curve at single point:
           curve.derivative(1.0, d=2)

           # Third derivative of curve at several points:
           curve.derivative([0.0, 1.0, 2.0], d=3)

           # Tangents of surface:
           surface.derivative(0.5, 0.7, d=(1,0))
           surface.derivative(0.5, 0.7, d=(0,1))

           # Cross-derivative of surface:
           surface.derivative(0.5, 0.7, d=(1,1))

        :param u,v,...: Parametric coordinates in which to evaluate
        :type u,v,...: float or [float]
        :param (int) d: Order of derivative to compute
        :param (bool) above: Evaluation in the limit from above
        :param tensor: Whether to evaluate on a tensor product grid
        :type tensor: bool
        :return: Derivatives
        :rtype: numpy.array
        """
        squeeze = all(is_singleton(p) for p in params)
        params_list = [ensure_scalars(p) for p in params]
        derivs = ensure_listlike(d, self.pardim)
        above = ensure_listlike(above, self.pardim)

        if not tensor and len({len(p) for p in params_list}) != 1:
            raise ValueError("Parameters must have same length")

        self._validate_domain(*params_list)

        # Evaluate the derivatives of the corresponding bases at the corresponding points
        # and build the result array
        dNs = [
            b.evaluate_dense(p, d=d, from_right=from_right)
            for b, p, d, from_right in zip(self.bases, params_list, derivs, above)
        ]
        result = _evaluate(dNs, self.controlpoints, tensor)

        # For rational curves, we need to use the quotient rule
        # (n/W)' = (n' W - n W') / W^2 = n'/W - nW'/W^2
        # * n'(i) = result[..., i]
        # * W'(i) = result[..., -1]
        # We evaluate in the regular way to compute n and W.
        if self.rational:
            if sum(derivs) > 1:
                raise RuntimeError("Rational derivative not implemented for order %i" % sum(derivs))
            Ns = [b.evaluate(p) for b, p in zip(self.bases, params_list)]
            non_derivative = _evaluate(Ns, self.controlpoints, tensor)
            W = non_derivative[..., -1]  # W
            Wd = result[..., -1]  # W'
            for i in range(self.dimension):
                result[..., i] = result[..., i] / W - non_derivative[..., i] * Wd / W / W
            result = np.delete(result, self.dimension, -1)

        # Squeeze the singleton dimensions if we only have one point
        if squeeze:
            result = result.reshape(self.dimension)

        return result

    @overload
    def get_derivative_spline(self) -> tuple[Self, ...]:
        ...

    @overload
    def get_derivative_spline(self, direction: Direction) -> Self:
        ...

    def get_derivative_spline(self, direction=None):  # type: ignore[no-untyped-def]
        """Compute the controlpoints associated with the derivative spline object

        If `direction` is given, only the derivatives in that direction are
        returned.

        If `direction` is not given, this function returns a tuple of all
        partial derivatives

        .. code:: python

           # Create a 4x4 element cubic spline surface
           surf = Surface()
           surf.raise_order(2,2)
           surf.refine(3,3)
           surf[1:4,1:4,:] += 0.1 # make the surface non-trivial by moving controlpoints

           # Create the derivative surface
           du = surf.get_derivative_spline(direction='u')

           # evaluation is identical
           print(du.evaluate(0.3, 0.4))
           print(surf.derivative(0.3, 0.4, d=(1,0)))

           print(surf.order()) # prints (3,3)
           print(du.order())   # prints (2,3)

        :param int direction: The tangential direction
        :return: Derivative spline
        :rtype: SplineObject
        """

        if self.rational:
            raise RuntimeError("Not working for rational splines")

        # if no direction is specified, return a tuple with all derivatives
        if direction is None:
            return tuple([self.get_derivative_spline(dim) for dim in range(self.pardim)])

        d = check_direction(direction, self.pardim)
        k = self.knots(d, with_multiplicities=True)
        p = self.order(d) - 1
        n = self.shape[d]
        if self.bases[d].periodic < 0:
            C = np.zeros((n - 1, n))
            for i in range(n - 1):
                C[i, i] = -float(p) / (k[i + p + 1] - k[i + 1])
                C[i, i + 1] = float(p) / (k[i + p + 1] - k[i + 1])
        else:
            C = np.zeros((n, n))
            for i in range(n):
                ip1 = np.mod(i + 1, n)
                C[i, i] = -float(p) / (k[i + p + 1] - k[i + 1])
                C[i, ip1] = float(p) / (k[i + p + 1] - k[i + 1])

        derivative_cps = np.tensordot(C, self.controlpoints, axes=(1, d))
        derivative_cps = derivative_cps.transpose(_transpose_fix(self.pardim, d))
        bases = list(self.bases)
        bases[d] = BSplineBasis(p, k[1:-1], bases[d].periodic - 1)

        return self.self_constructor(bases, derivative_cps, rational=self.rational, raw=True)

    @overload
    def tangent(
        self,
        *params: ScalarOrScalars,
        direction: Direction,
        above: Union[bool, Sequence[bool]] = True,
        tensor: bool = True,
    ) -> FArray:
        ...

    @overload
    def tangent(
        self,
        *params: ScalarOrScalars,
        direction: None = None,
        above: Union[bool, Sequence[bool]] = True,
        tensor: bool = True,
    ) -> tuple[FArray, ...]:
        ...

    def tangent(self, *params, direction=None, above=True, tensor=True):  # type: ignore[no-untyped-def]
        """Evaluate the tangents of the object at the given parametric values.

        If `direction` is given, only the derivatives in that direction are
        evaluated. This is equivalent to calling
        :func:`splipy.SplineObject.derivative` with
        `d=(0,...,0,1,0,...,0)`, the unit vector corresponding to the given
        direction.

        If `direction` is not given, this function returns a tuple of all
        tangential derivatives at the given points.

        :param u,v,...: Parametric coordinates in which to evaluate
        :type u,v,...: float or [float]
        :param int direction: The tangential direction
        :param (bool) above: Evaluation in the limit from above
        :param tensor: Whether to evaluate on a tensor product grid
        :type tensor: bool
        :return: Tangents
        :rtype: tuple<numpy.array>
        """

        derivative = [0] * self.pardim
        above = ensure_listlike(above, self.pardim)

        if direction is None:
            result: tuple[FArray, ...] = ()
            for i in range(self.pardim):
                derivative[i] = 1

                # Compute velocity in this direction
                v = self.derivative(*params, d=derivative, above=above, tensor=tensor)

                # Normalize
                if v.ndim == 1:
                    v /= np.linalg.norm(v)
                else:
                    speed = np.linalg.norm(v, axis=-1)
                    v /= np.reshape(speed, speed.shape + (1,))

                # Store in result tuple
                result += (v,)
                derivative[i] = 0
            return result

        if self.pardim == 1:  # curves
            direction = 0

        i = check_direction(direction, self.pardim)
        derivative[i] = 1

        # Compute velocity in this direction
        v = self.derivative(*params, d=derivative, above=above, tensor=tensor)

        # Normalize
        if v.ndim == 1:
            v /= np.linalg.norm(v)
        else:
            speed = np.linalg.norm(v, axis=-1)
            v /= np.reshape(speed, speed.shape + (1,))

        return v

    @overload
    def section(
        self, *args: SectionElt, unwrap_points: Literal[True] = True, **kwargs: Unpack[SectionKwargs]
    ) -> Union[SplineObject, FArray]:
        ...

    @overload
    def section(
        self, *args: SectionElt, unwrap_points: Literal[False], **kwargs: Unpack[SectionKwargs]
    ) -> SplineObject:
        ...

    def section(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        """Return a section from the object. A section can be any sub-object of
        parametric dimension not exceeding that of the object. E.g. for a
        volume, sections include vertices, edges, faces, etc.

        The arguments are control point indices for each direction. `None`
        means that direction should be variable in the returned object.

        .. code:: python

           # Get the face with u=max
           vol.section(-1, None, None)

           # Keyword arguments are supported for u, v and w
           # This is the same as the above
           vol.section(u=-1)

           # Get the edge with u=min, v=max
           vol.section(0, -1, None)

           # This is equivalent to vol.clone()
           vol.section()

        If a specific subclass of `SplineObject` is found that handles the
        requested number of variable directions (parametric dimension), then
        the return value is of that type. If not, it will be a generic `SplineObject`.

        If the section has no variable directions (it is a point), then the
        return value will be an array, unless the keyword argument
        `unwrap_points` is true, in which case it will return a
        zero-dimensional `SplineObject`.

        :param u,v,...: Control point indices
        :type u,v,...: int or None
        :return: Section
        :rtype: SplineObject or np.array
        """
        section = check_section(*args, pardim=self.pardim, **kwargs)
        unwrap_points = kwargs.get("unwrap_points", True)

        slices = tuple(slice(None) if p is None else p for p in section)
        bases = [b for b, p in zip(self.bases, section) if p is None]
        if bases or not unwrap_points:
            if 1 <= len(bases) <= 3:
                return SplineObject.constructor(len(bases))(
                    bases, self.controlpoints[slices], rational=self.rational, raw=True
                )
            return SplineObject(bases, self.controlpoints[slices], rational=self.rational, raw=True)
        return self.controlpoints[slices]

    def set_order(self, *order: int) -> Self:
        """Set the polynomial order of the object. If only one argument is
        given, the order is set uniformly over all directions.

        :param int u,v,...: The new order in a given direction.
        :raises ValueError: If the order is reduced in any direction.
        :return: self
        """
        if len(order) == 1:
            order = (order[0],) * self.pardim
        if not all(new >= old for new, old in zip(order, self.order())):
            raise ValueError("Cannot lower order using set_order")

        diff = [new - old for new, old in zip(order, self.order())]
        return self.raise_order(*diff)

    def raise_order(self, *raises: int, direction: Optional[Direction] = None) -> Self:
        """Raise the polynomial order of the object. If only one
        argument is given, the order is raised equally over all
        directions, unless the `direction` argument is also given. The
        explicit version is only implemented on open knot vectors. The
        method raise_order_implicit is used otherwise.

        :param int u,v,...: Number of times to raise the order in a given
            direction.
        :param int direction: The direction in which to raise the order.
        :return: self
        """
        if len(raises) == 1 and direction is None:
            raises_list = [raises[0]] * self.pardim
        elif len(raises) == 1 and direction is not None:
            raises_list = [0] * self.pardim
            raises_list[check_direction(direction, self.pardim)] = raises[0]
        else:
            raises_list = list(raises)

        if not all(r >= 0 for r in raises_list):
            raise ValueError("Cannot lower order using raise_order")
        if all(r == 0 for r in raises_list):
            return self

        if any(b.continuity(b.knots[0]) < b.order or b.periodic > -1 for b in self.bases):
            self.raise_order_implicit(*raises_list)
            return self

        new_bases = [b.raise_order(r) for b, r in zip(self.bases, raises_list)]

        d_p = self.pardim

        controlpoints = self.controlpoints
        for i in range(0, d_p):
            dimensions = np.array(controlpoints.shape)
            indices = np.array(range(0, d_p + 1))
            indices[i], indices[d_p] = d_p, i
            controlpoints = np.transpose(controlpoints, indices)
            controlpoints = np.reshape(controlpoints, (np.prod(dimensions[indices[:-1]]), dimensions[i]))
            controlpoints = raise_order_1D(
                controlpoints.shape[1] - 1,
                self.order(i),
                self.bases[i].knots,
                controlpoints,
                raises_list[i],
                self.bases[i].periodic,
            )
            controlpoints = np.reshape(
                controlpoints, np.append(dimensions[indices[:-1]], controlpoints.shape[1])
            )
            controlpoints = np.transpose(controlpoints, indices)

        self.controlpoints = controlpoints
        self.bases = new_bases
        return self

    def raise_order_implicit(self, *raises: int) -> Self:
        """Raise the polynomial order of the object. If only one argument is
        given, the order is raised equally over all directions.

        :param int u,v,...: Number of times to raise the order in a given
            direction.
        :return: self
        """

        new_bases = [b.raise_order(r) for b, r in zip(self.bases, raises)]

        # Set up an interpolation problem
        # This works in projective space, so no special handling for rational objects
        interpolation_pts = [b.greville() for b in new_bases]
        N_old = [b(pts) for b, pts in zip(self.bases, interpolation_pts)]
        N_new = [b(pts) for b, pts in zip(new_bases, interpolation_pts)]

        # Calculate the projective interpolation points
        result = self.controlpoints
        for n in N_old[::-1]:
            result = np.tensordot(n, result, axes=(1, self.pardim - 1))

        # Solve the interpolation problem
        for n in N_new[::-1]:
            result = np.tensordot(np.linalg.inv(n), result, axes=(1, self.pardim - 1))

        self.controlpoints = result
        self.bases = new_bases

        return self

    def lower_order(self, *lowers: int) -> Self:
        """Lower the polynomial order of the object. If only one argument is
        given, the order is lowered equally over all directions.

        :param int u,v,...: Number of times to lower the order in a given
            direction.
        :return SplineObject: Approximation of the current object on a lower
            order basis
        """
        if len(lowers) == 1:
            lowers = (lowers[0],) * self.pardim
        if all(l == 0 for l in lowers):
            return self.clone()

        new_bases = [b.lower_order(l) for b, l in zip(self.bases, lowers)]

        # Set up an interpolation problem
        # This works in projective space, so no special handling for rational objects
        interpolation_pts = [b.greville() for b in new_bases]
        N_old = [b(pts) for b, pts in zip(self.bases, interpolation_pts)]
        N_new = [b(pts) for b, pts in zip(new_bases, interpolation_pts)]

        # Calculate the projective interpolation points
        new_controlpts = self.controlpoints
        for n in N_old[::-1]:
            new_controlpts = np.tensordot(n, new_controlpts, axes=(1, self.pardim - 1))

        # Solve the interpolation problem
        for n in N_new[::-1]:
            new_controlpts = np.tensordot(np.linalg.inv(n), new_controlpts, axes=(1, self.pardim - 1))

        # search for the right subclass constructor, i.e. Volume, Surface or Curve
        return self.self_constructor(new_bases, new_controlpts, rational=self.rational, raw=True)

    @overload
    def start(self) -> tuple[int, ...]:
        ...

    @overload
    def start(self, direction: Direction, /) -> int:
        ...

    def start(self, direction=None):  # type: ignore[no-untyped-def]
        """Return the start of the parametric domain.

        If `direction` is given, returns the start of that direction, as a
        float. If it is not given, returns the start of all directions, as a
        tuple.

        :param int direction: Direction in which to get the start.
        :raises ValueError: For invalid direction
        """
        if direction is None:
            return tuple(b.start() for b in self.bases)
        direction = check_direction(direction, self.pardim)
        return self.bases[direction].start()

    @overload
    def end(self) -> tuple[int, ...]:
        ...

    @overload
    def end(self, direction: Direction, /) -> int:
        ...

    def end(self, direction=None):  # type: ignore[no-untyped-def]
        """Return the end of the parametric domain.

        If `direction` is given, returns the end of that direction, as a float.
        If it is not given, returns the end of all directions, as a tuple.

        :param int direction: Direction in which to get the end.
        :raises ValueError: For invalid direction
        """
        if direction is None:
            return tuple(b.end() for b in self.bases)
        direction_index = check_direction(direction, self.pardim)
        return self.bases[direction_index].end()

    @overload
    def order(self) -> tuple[int, ...]:
        ...

    @overload
    def order(self, direction: Direction, /) -> int:
        ...

    def order(self, direction=None):  # type: ignore[no-untyped-def]
        """Return polynomial order (degree + 1).

        If `direction` is given, returns the order of that direction, as an
        int. If it is not given, returns the order of all directions, as a
        tuple.

        :param int direction: Direction in which to get the order.
        :raises ValueError: For invalid direction
        """
        if direction is None:
            return tuple(b.order for b in self.bases)
        direction = check_direction(direction, self.pardim)
        return self.bases[direction].order

    @overload
    def knots(self, /, with_multiplicities: bool = False) -> tuple[FArray, ...]:
        ...

    @overload
    def knots(self, direction: Direction, /, with_multiplicities: bool = False) -> FArray:
        ...

    def knots(self, direction=None, with_multiplicities=False):  # type: ignore[no-untyped-def]
        """Return knots vector

        If `direction` is given, returns the knots in that direction, as a
        list. If it is not given, returns the knots of all directions, as a
        tuple.

        :param int direction: Direction in which to get the knots.
        :param bool with_multiplicities: If true, return knots with
            multiplicities (i.e. repeated).
        :raises ValueError: For invalid direction
        """
        getter: Callable[[BSplineBasis], FArray]
        getter = attrgetter("knots") if with_multiplicities else methodcaller("knot_spans")  # type: ignore[assignment]

        if direction is None:
            return tuple(getter(b) for b in self.bases)
        direction_index = check_direction(direction, self.pardim)
        return getter(self.bases[direction_index])

    def reverse(self, direction: Direction = 0) -> Self:
        """Swap the direction of a parameter by making it go in the reverse
        direction. The parametric domain remains unchanged.

        :param int direction: The direction to flip.
        :return: self
        """
        direction_index = check_direction(direction, self.pardim)
        self.bases[direction_index].reverse()

        # This creates the following slice programmatically
        # array[:, :, :, ..., ::-1,]
        # index=direction -----^
        # :    => slice(None, None, None)
        # ::-1 => slice(None, None, -1)
        slices = [slice(None, None, None) for _ in range(direction_index)] + [slice(None, None, -1)]
        self.controlpoints = self.controlpoints[tuple(slices)]

        return self

    def swap(self, dir1: Direction = 0, dir2: Direction = 1, /) -> Self:
        """Swap two parameter directions.

        This function silently passes for curves.

        :param direction dir1: The first direction (default u)
        :param direction dir2: The second direction (default v)
        :return: self
        """
        if self.pardim == 1:
            return self

        dir1_index = check_direction(dir1, self.pardim)
        dir2_index = check_direction(dir2, self.pardim)

        # Reorder control points
        new_directions = list(range(self.pardim + 1))
        new_directions[dir1_index] = dir2_index
        new_directions[dir2_index] = dir1_index
        self.controlpoints = self.controlpoints.transpose(new_directions)

        # Swap knot vectors
        self.bases[dir1_index], self.bases[dir2_index] = self.bases[dir2_index], self.bases[dir1_index]

        return self

    def insert_knot(self, knot: ScalarOrScalars, direction: Direction = 0) -> Self:
        """Insert a new knot into the spline.

        :param int direction: The direction to insert in
        :param knot: The new knot(s) to insert
        :type knot: float or [float]
        :raises ValueError: For invalid direction
        :return: self
        """
        shape = self.controlpoints.shape

        # for single-value input, wrap it into a list
        knot_list = ensure_scalars(knot)

        direction_index = check_direction(direction, self.pardim)

        C = np.identity(shape[direction_index])
        for k in knot_list:
            C = self.bases[direction_index].insert_knot(k) @ C
        self.controlpoints = np.tensordot(C, self.controlpoints, axes=(1, direction_index))
        self.controlpoints = self.controlpoints.transpose(_transpose_fix(self.pardim, direction_index))

        return self

    @overload
    def refine(self, n: int, /, direction: Direction) -> Self:
        ...

    @overload
    def refine(self, *args: int) -> Self:
        ...

    def refine(self, *ns, **kwargs):  # type: ignore[no-untyped-def]
        """Enrich the spline space by inserting knots into each existing knot
        span.

        This method supports three different usage patterns:

        .. code:: python

           # Refine each direction by a factor n
           obj.refine(n)

           # Refine a single direction by a factor n
           obj.refine(n, direction='v')

           # Refine all directions by given factors
           obj.refine(nu, nv, ...)

        :param int nu,nv,...: Number of new knots to insert into each span
        :param int direction: Direction to refine in
        :return: self
        """
        direction = kwargs.get("direction", None)

        if len(ns) == 1 and direction is not None:
            directions = iter([check_direction(direction, self.pardim)])
        else:
            directions = iter(range(self.pardim))

        factors: Sequence[int] = ns if len(ns) > 1 else [ns[0]] * self.pardim

        for n, d in zip(factors, directions):
            knots = self.knots(d)  # excluding multiple knots
            new_knots = []
            for k0, k1 in zip(knots[:-1], knots[1:]):
                new_knots.extend(np.linspace(k0, k1, n + 2)[1:-1])
            self.insert_knot(new_knots, d)

        return self

    @overload
    def reparam(self, *args: Union[FArray, tuple[Scalar, Scalar]]) -> Self:
        ...

    @overload
    def reparam(self, arg: Union[FArray, tuple[Scalar, Scalar]], /, direction: Direction) -> Self:
        ...

    @overload
    def reparam(self, /, direction: Direction) -> Self:
        ...

    def reparam(self, *args, **kwargs):  # type: ignore[no-untyped-def]
        """Redefine the parametric domain. This function accepts two calling
        conventions:

        `reparametrize(u, v, ...)` reparametrizes each direction to the domains
        given by the tuples *u*, *v*, etc. It is equivalent to calling
        `reparametrize(u[0], u[1])` on each basis. The default domain for
        directions not given is (0,1). In particular, if no arguments are
        given, the new parametric domain will be the unit (hyper)cube.

        `reparametrize(u, direction=d)` reparametrizes just the direction given
        by *d* and leaves the others untouched.

        :param tuple u, v, ...: New parametric domains, default to (0,1)
        :param int direction: The direction to reparametrize
        :return: self
        """
        if "direction" not in kwargs:
            # Pad the args with (0, 1) for the extra directions
            intervals: list[tuple[Scalar, Scalar]] = list(args) + [(0, 1)] * (len(self.bases) - len(args))
            for b, (start, end) in zip(self.bases, intervals):
                b.reparam(start, end)
        else:
            direction = check_direction(kwargs["direction"], self.pardim)
            if len(args) == 0:
                self.bases[direction].reparam(0, 1)
            else:
                start, end = args[0]
                self.bases[direction].reparam(start, end)

        return self

    def translate(self, x: Scalars) -> Self:
        """Translate (i.e. move) the object by a given distance.

        :param array-like x: The vector to translate by.
        :return: self
        """
        # 3D rational example: create a 4x4 translation matrix
        #
        #  |xw|      |  1   0   0  x1 |   |xw|
        #  |yw|   =  |  0   1   0  x2 | * |yw|
        #  |zw|      |  0   0   1  x3 |   |zw|
        #  | w|_new  |  0   0   0   1 |   | w|_old
        #
        #  PS: we even need a rational representation for non-rational splines
        #      in order to formulate translation as a matrix-matrix product
        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points
        if len(x) > dim:  # typical case: requesting movement in z-direction for 2D geometries
            self.set_dimension(len(x))
            dim = self.dimension

        # set up the translation matrix
        translation_matrix = np.identity(dim + 1)
        for i in range(dim):
            translation_matrix[i, -1] = x[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        if not self.rational:
            cp = np.ones((n, dim + 1))  # pad with weights=1
            cp[:, :-1] = np.reshape(self.controlpoints, (n, dim))
        else:
            cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual scaling by matrix-matrix multiplication
        cp = cp @ translation_matrix.T  # right-mult, so we need transpose

        # store results
        if self.rational:
            self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)
        else:
            self.controlpoints = np.reshape(np.array(cp[:, :-1]), self.controlpoints.shape)

        return self

    def scale(self, *args: Scalar) -> Self:
        """Scale, or magnify the object by a given amount.

        In case of one input argument, the scaling is uniform.

        :param args: Scaling factors, possibly different in each direction.
        :type args: array-like or float
        :return: self
        """
        # 3D rational example: create a 4x4 scaling matrix
        #
        #  |xw|      |  sx    0    0   0 |   |xw|
        #  |yw|   =  |   0   sy    0   0 | * |yw|
        #  |zw|      |   0    0   sz   0 |   |zw|
        #  | w|_new  |   0    0    0   1 |   | w|_old
        #
        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points
        s = ensure_scalars(args, dups=3)

        # set up the scaling matrix
        scale_matrix = np.identity(dim + rat)
        for i in range(dim):
            scale_matrix[i, i] = s[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual scaling by matrix-matrix multiplication
        cp = cp @ scale_matrix

        # store results
        self.controlpoints = np.reshape(cp, self.controlpoints.shape)

        return self

    def rotate(self, theta: Scalar, normal: Scalars = (0, 0, 1)) -> Self:
        """Rotate the object around an axis.

        :param float theta: Angle to rotate about, measured in radians
        :param array-like normal: The normal axis (if 3D) to rotate about
        :raises RuntimeError: If the physical dimension is not 2 or 3
        :return: self
        """
        # 2D rational example: create a 3x3 rotation matrix
        #
        #  |xw|     | cos(t)  -sin(t)   0 |   |xw|
        #  |yw|   = | sin(t)   cos(t)   0 | * |yw|
        #  | w|_new |   0        0      1 |   | w|_old
        #

        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points
        if not (normal[0] == 0 and normal[1] == 0):  # rotating a 2D geometry out of the xy-plane
            self.set_dimension(3)
            dim = self.dimension

        # set up the rotation matrix
        if dim == 2:
            R = np.array(
                [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]
            ).T  # we do right-multiplication, so we need a transpose
        elif dim == 3:
            normal = np.array(normal)
            R = rotation_matrix(theta, normal)
        else:
            raise RuntimeError("rotation undefined for geometries other than 2D and 3D")

        rot_matrix = np.identity(dim + rat)
        rot_matrix[0:dim, 0:dim] = R

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual rotation by matrix-matrix multiplication
        cp = cp @ rot_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def mirror(self, normal: Scalars) -> Self:
        """Mirror the object around a plane through the origin.

        :param array-like normal: The plane normal to mirror about.
        :raises RuntimeError: If the physical dimension is not 2 or 3
        :return: self
        """
        # 3D rational example: create a 4x4 reflection matrix
        #
        #  normal = [a,b,c]
        #
        #  |xw|     |  1-2a^2   -2ab    -2ac   0 |   |xw|
        #  |yw|   = |   -2ab   1-2b^2   -2bc   0 | * |yw|
        #  |zw|   = |   -2ac    -2bc   1-2c^2  0 | * |zw|
        #  | w|_new |    0       0         0   1 |   | w|_old
        #
        #  PS: A reflection about a line is not a linear transformation; it is
        #      an affine transformation.

        dim = self.dimension
        rat = self.rational
        n = len(self)  # number of control points

        if dim != 3:
            raise RuntimeError("reflection undefined for geometries other than 3D")

        # fixup the input normal to right form
        normal_vec = np.asarray(normal, dtype=float)
        normal_vec /= np.sqrt(np.dot(normal_vec, normal_vec))  # normalize it

        # set up the reflection matrix
        reflection_matrix = np.identity(dim + rat)
        reflection_matrix[0:dim, 0:dim] -= 2 * np.outer(normal_vec, normal_vec)

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual rotation by matrix-matrix multiplication
        cp = cp @ reflection_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def project(self, plane: Literal["", "x", "y", "z", "xy", "yz", "xz", "xyz"]) -> Self:
        """Project the geometry onto a plane or axis.

        - `project('xy')` will project the object onto the *xy* plane, setting
          all *z* components to zero.
        - `project('y')` will project the object onto the *y* axis, setting all
          *x* and *z* components to zero.

        :param string plane: Any combination of 'x', 'y' and 'z'
        :return: self
        """
        keep = [c in plane.lower() for c in "xyz"]

        dim = self.dimension
        for i in range(dim):
            if not keep[i]:
                self.controlpoints[..., i] = 0

        return self

    def bounding_box(self) -> list[tuple[float, float]]:
        """Get the bounding box of a spline object, computed from the
        control-point values. Could be inaccurate for rational splines.

        Returns the minima and maxima for each direction:
        [(xmin, xmax), (ymin, ymax), ...]

        :return: Bounding box
        :rtype: [(float)]
        """
        dim = self.dimension

        result = []
        for i in range(dim):
            result.append(
                (
                    np.min(self.controlpoints[..., i]),
                    np.max(self.controlpoints[..., i]),
                )
            )
        return result

    def center(self) -> FArray:
        """Get the center of the domain

        For curves this will return :math:`(\\tilde{x}, \\tilde{y},...)`, where

        .. math:: \\tilde{x} = \\frac{1}{L} \\int_{t_0}^{t_1} x(t) \\; dt

        and :math:`L=t_1-t_0` is the length of the parametric domain :math:`[t_0,t_1]`.

        For surfaces this will return :math:`(\\tilde{x}, \\tilde{y},...)`, where

        .. math:: \\tilde{x} = \\frac{1}{A} \\int_{v_0}^{v_1} \\int_{u_0}^{u_1} x(u,v) \\; du \\; dv

        and :math:`A=(u_1-u_0)(v_1-v_0)` is the area of the parametric domain
        :math:`[u_0,u_1]\\times[v_0,v_1]`.

        .. warning:: For rational splines, this will integrate in projective
            coordinates, then project the centerpoint. This is as opposed to
            integrating the rational functions :math:`\\frac{N_i(t)w_i}{\\sum_j
            N_j(t)w_j}`.
        """

        # compute integration of basis functions
        Ns = [b.integrate(b.start(), b.end()) for b in self.bases]

        # compute parametric size
        par_size = np.prod([t1 - t0 for (t0, t1) in zip(self.start(), self.end())])

        # multiply basis functions with control points
        idx = self.pardim - 1
        result = self.controlpoints
        for N in Ns[::-1]:
            result = np.tensordot(N, result, axes=(0, idx))
            idx -= 1

        result /= par_size

        # project to physical space
        if self.rational:
            result[:-1] /= result[-1]
            result = np.delete(result, self.dimension)

        return result

    def corners(self, order: Literal["C", "F"] = "C") -> FArray:
        """Return the corner control points.

        The `order` parameter determines which order to use, either ``'F'`` or
        ``'C'``, for row-major or column-major ordering. E.g. for a volume, in
        parametric coordinates,

        - ``'C'`` gives (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), etc.
        - ``'F'`` gives (0,0,0), (0,0,1), (0,1,0), (0,1,1), (1,0,0), etc.

        :param str order: The ordering to use
        :return: Corners
        :rtype: np.array

        .. warning:: For rational splines, this will return the corners in
            projective coordinates, including weights.
        """
        result = np.zeros((2**self.pardim, self.dimension + int(self.rational)), dtype=float)
        for i, args in enumerate(sections(self.pardim, 0)):
            result[i, :] = self.section(*(args[::-1] if order == "F" else args))
        return result

    def lower_periodic(self, periodic: int, direction: Direction = 0) -> Self:
        """Set the periodicity of the spline object in the given direction,
        keeping the geometry unchanged.

        :param int periodic: new periodicity, i.e. the basis is C^k over the start/end
        :param int direction: the parametric direction of the basis to modify
        :return: self
        """
        direction_index = check_direction(direction, self.pardim)

        b = self.bases[direction_index]
        while periodic < b.periodic:
            self.insert_knot(self.start(direction_index), direction_index)
            self.controlpoints = np.roll(self.controlpoints, -1, direction_index)
            b.roll(1)
            b.periodic -= 1
            b.knots = b.knots[:-1]
        if periodic > b.periodic:
            raise ValueError("Cannot raise periodicity")

        return self

    def set_dimension(self, new_dim: int) -> Self:
        """Set the physical dimension of the object. If increased, the new
        components are set to zero.

        :param int new_dim: New dimension.
        :return: self
        """
        dim = self.dimension
        shape = self.controlpoints.shape
        while new_dim > dim:
            self.controlpoints = np.insert(self.controlpoints, dim, np.zeros(shape[:-1]), self.pardim)
            dim += 1
        while new_dim < dim:
            self.controlpoints = np.delete(self.controlpoints, -2 if self.rational else -1, -1)
            dim -= 1
        self.dimension = new_dim

        return self

    def periodic(self, direction: Direction = 0) -> bool:
        """Return true if the spline object is periodic in the given parametric direction."""
        direction = check_direction(direction, self.pardim)
        return self.bases[direction].periodic > -1

    def force_rational(self) -> Self:
        """Force a rational representation of the object.

        The weights of a non-rational object will be set to 1.

        :return: self
        """
        if not self.rational:
            dim = self.dimension
            shape = self.controlpoints.shape
            self.controlpoints = np.insert(self.controlpoints, dim, np.ones(shape[:-1]), self.pardim)
            self.rational = True

        return self

    def split_periodic(self, knot: Scalar, direction: Direction = 0) -> Self:
        """Split a periodic object along one of its periodic axes.

        :param knot: The splitting knot
        :type knot: float
        :param direction: Parametric direction
        :type direction: int
        :return: The new object
        :rtype: SplineObject
        """

        direction_index = check_direction(direction, self.pardim)
        assert self.periodic(direction_index)

        splitting_obj = self.clone()

        continuity = splitting_obj.bases[direction_index].continuity(knot)
        continuity = min(continuity, self.order(direction_index) - 1)
        splitting_obj.insert_knot([float(knot)] * (continuity + 1), direction_index)

        basis = splitting_obj.bases[direction_index]
        mu = bisect_left(basis.knots, knot)
        basis.roll(mu)
        splitting_obj.controlpoints = np.roll(splitting_obj.controlpoints, -mu, direction_index)
        basis.knots = basis.knots[: -basis.periodic - 1]
        basis.periodic = -1

        return splitting_obj

    def split_nonperiodic(self, knot: ScalarOrScalars, direction: Direction = 0) -> list[Self]:
        """Split an object into two or more separate representations with C0
        continuity between them.

        :param knots: The splitting points
        :type knots: float or [float]
        :param direction: Parametric direction
        :type direction: int
        :return: The new objects
        :rtype: [SplineObject]
        """

        # for single-value input, wrap it into a list
        knots_list = ensure_scalars(knot)

        # error test input
        direction_index = check_direction(direction, self.pardim)
        assert not self.periodic(direction_index)

        p = self.order(direction_index)
        results: list[Self] = []
        splitting_obj = self.clone()
        bases = self.bases
        # insert knots to produce C{-1} at all splitting points
        for k in knots_list:
            continuity = bases[direction_index].continuity(k)
            if continuity == np.inf:
                continuity = p - 1
            splitting_obj.insert_knot([k] * (continuity + 1), direction_index)

        # everything is available now, just have to find the right index range
        # in the knot vector and controlpoints to store in each separate curve
        # piece
        last_cp_i = 0
        last_knot_i = 0

        bases = splitting_obj.bases
        b = bases[direction_index]
        cp_slice = [slice(None, None, None)] * len(self.controlpoints.shape)
        for k in knots_list:
            if self.start(direction_index) < k < self.end(direction_index):  # skip start/end points
                mu = bisect_left(b.knots, k)
                n_cp = mu - last_knot_i
                knot_slice = slice(last_knot_i, mu + p, None)
                cp_slice[direction_index] = slice(last_cp_i, last_cp_i + n_cp, None)

                cp = splitting_obj.controlpoints[tuple(cp_slice)]
                bases[direction_index] = BSplineBasis(p, b.knots[knot_slice])

                results.append(self.self_constructor(bases, cp, rational=splitting_obj.rational, raw=True))

                last_knot_i = mu
                last_cp_i += n_cp

        # with n splitting points, we're getting n+1 pieces. Add the final one:
        knot_slice = slice(last_knot_i, None, None)
        cp_slice[direction_index] = slice(last_cp_i, None, None)
        bases[direction_index] = BSplineBasis(p, b.knots[knot_slice])
        cp = splitting_obj.controlpoints[tuple(cp_slice)]
        results.append(self.self_constructor(bases, cp, rational=splitting_obj.rational, raw=True))

        return results

    def split(self, knots: ScalarOrScalars, direction: Direction = 0) -> Union[Self, list[Self]]:
        """Split an object into two or more separate representations with C0
        continuity between them.

        :param knots: The splitting points
        :type knots: float or [float]
        :param direction: Parametric direction
        :type direction: int
        :return: The new objects
        :rtype: [SplineObject]
        """
        direction_index = check_direction(direction, self.pardim)
        if self.periodic(direction_index):
            knots_list = ensure_scalars(knots)
            split_obj = self.split_periodic(knots_list[0], direction)
            if len(knots_list) > 1:
                return split_obj.split_nonperiodic(knots_list[1:], direction)
            return split_obj
        return self.split_nonperiodic(knots, direction)

    def split_many(self, knots: ScalarOrScalars, direction: Direction = 0) -> list[Self]:
        """Like split, but always returns a list."""
        split_obj = self.split(knots, direction)
        if isinstance(split_obj, list):
            return split_obj
        return [split_obj]

    def make_periodic(self, continuity: Optional[int] = None, direction: Direction = 0) -> Self:
        """Make the spline object periodic in a given parametric direction.

        :param int continuity: The continuity along the boundary (default max).
        :param int direction: The direction to ensure continuity in.
        """

        direction = check_direction(direction, self.pardim)
        basis = self.bases[direction]
        if continuity is None:
            continuity = basis.order - 2
        if not -1 <= continuity <= basis.order - 2:
            raise ValueError(f"Illegal continuity for basis of order {continuity}: {basis.order}")
        if continuity == -1:
            raise ValueError(
                "Operation not supported. For discontinuous spline spaces, consider SplineObject.split()."
            )
        if basis.periodic >= 0:
            raise ValueError("Basis is already periodic")

        basis = basis.make_periodic(continuity)

        # Merge control points
        index_beg: list[Union[slice, int]] = [slice(None, None, None)] * (self.pardim + 1)
        index_end: list[Union[slice, int]] = [slice(None, None, None)] * (self.pardim + 1)
        cps = np.array(self.controlpoints)
        weights = iter(np.linspace(0, 1, continuity + 1)) if continuity > 0 else iter([0.5])
        for i, j, t in zip(range(continuity + 1), range(-continuity - 1, 0), weights):
            # Weighted average between cps[..., i, ..., :] and cps[..., -c-1+i, ..., :]
            # The weights are chosen so that, for periodic c, the round trip
            # c.split().make_periodic() with suitable arguments produces an
            # object identical to c. (Mostly black magic.)
            index_beg[direction] = i
            index_end[direction] = j
            cps[tuple(index_beg)] = t * cps[tuple(index_beg)] + (1 - t) * cps[tuple(index_end)]

        # cps[..., :-(continuity+1), ..., :]
        index_beg[direction] = slice(None, -(continuity + 1), None)
        cps = cps[tuple(index_beg)]

        bases = list(self.bases)
        bases[direction] = basis

        return self.self_constructor(bases, cps, rational=self.rational, raw=True)

    @property
    def pardim(self) -> int:
        """The number of parametric dimensions: 1 for curves, 2 for surfaces, 3
        for volumes, etc.
        """
        return self.controlpoints.ndim - 1

    def clone(self) -> Self:
        """Clone the object."""
        return copy.deepcopy(self)

    __call__ = evaluate

    def __len__(self) -> int:
        """Return the number of control points (basis functions) for the object."""
        n = 1
        for b in self.bases:
            n *= b.num_functions()
        return n

    def _unravel_flat_index(self, i: Union[slice, SupportsIndex]) -> tuple[IPArray, ...]:
        """Unravels a flat index i to multi-indexes.

        :param i: Flat index
        :type i: int or slice
        :rtype: Tuple of np.array
        :raises IndexError: If the index is out of bounds
        """
        # i is int => make sure we deal with negative i properly
        # i is slice => use i.indices to compute the actual indices
        total = len(self)
        if isinstance(i, SupportsIndex):
            j = i.__index__()
            indexes = [j] if j >= 0 else [total + j]
        else:
            indexes = list(range(*i.indices(total)))

        # Convert to multi-indexes
        try:
            unraveled = np.unravel_index(indexes, self.controlpoints.shape[:-1], order="F")
        except ValueError:
            raise IndexError

        return unraveled

    def __getitem__(self, i: Union[slice, SupportsIndex, tuple[Union[slice, SupportsIndex], ...]]) -> FArray:
        """Get the control point at a given index.

        Indexing is in column-major order. Examples of supported indexing
        modes:

        .. code:: python

           # Flat indexing with an int
           obj[4]

           # Flat indexing from the end
           obj[-1]

           # Flat indexing with a slice
           obj[2:5]

           # Multi-indexing with ints, negative ints and slices
           obj[0,-1,:]

        :rtype: numpy.array
        """
        if isinstance(i, tuple):
            return self.controlpoints[i]

        unraveled = self._unravel_flat_index(i)

        # Singleton dimensions should be squeezed if the input was an int
        if isinstance(i, int):
            return self.controlpoints[unraveled][0]  # type: ignore[no-any-return]
        return self.controlpoints[unraveled]  # type: ignore[no-any-return]

    def __setitem__(
        self,
        i: Union[int, slice, SupportsIndex, tuple[Union[int, slice, SupportsIndex], ...]],
        cp: ArrayLike,
    ) -> None:
        """Set the control points at given indices.

        This function supports the same indexing modes as
        :func:`SplineObject.__getitem__`

        :param int i: Index or indices
        :param numpy.array cp: New control point(s)
        """
        if isinstance(i, tuple):
            self.controlpoints[i] = cp
            return

        unraveled = self._unravel_flat_index(i)
        self.controlpoints[unraveled] = cp

    @property
    def shape(self) -> tuple[int, ...]:
        """The dimensions of the control point array."""
        return self.controlpoints.shape[:-1]

    # TODO(Eivind): Py310 from types import NotImplementedType
    # Then change all the return types to Self | NotImplementedType

    def __iadd__(self, x: Any) -> Self:
        if isinstance(x, (Sequence, np.ndarray)):
            self.translate(x)
            return self
        return NotImplemented

    def __isub__(self, x: Any) -> Self:
        if isinstance(x, (Sequence, np.ndarray)):
            self.translate(-np.array(x, dtype=float))  # can't do -x if x is a list, so we rewrap it here
            return self
        return NotImplemented

    def __imul__(self, x: Any) -> Self:
        if isinstance(x, (Sequence, np.ndarray)):
            self.scale(*x)
            return self
        if isinstance(x, float):
            self.scale(x)
            return self
        if isinstance(x, SupportsFloat):
            self.scale(float(x))
            return self
        return NotImplemented

    def __itruediv__(self, x: Any) -> Self:
        if isinstance(x, (Sequence, np.ndarray)):
            y = 1 / np.ndarray(x, dtype=float)
            self.scale(*y)
            return self
        if isinstance(x, float):
            self.scale(1 / x)
            return self
        if isinstance(x, SupportsFloat):
            self.scale(1 / float(x))
            return self
        return NotImplemented

    __ifloordiv__ = __itruediv__  # integer division (should not distinguish)

    def __add__(self, x: Any) -> Self:
        return self.clone().__iadd__(x)

    def __radd__(self, x: Any) -> Self:
        return self.__add__(x)

    def __sub__(self, x: Any) -> Self:
        return self.clone().__isub__(x)

    def __mul__(self, x: Any) -> Self:
        return self.clone().__imul__(x)

    def __rmul__(self, x: Any) -> Self:
        return self.__mul__(x)

    def __div__(self, x: Any) -> Self:
        return self.clone().__itruediv__(x)

    def flip_and_move_plane_geometry(self, center: Scalars = (0, 0, 0), normal: Scalars = (0, 0, 1)) -> Self:
        """Re-orient a planar geometry by moving it to a different location and
        tilting it.

        Don't call unless necessary. Translate or scale operations may force
        an object into 3D space.
        """
        if not np.allclose(np.asarray(normal), np.array([0, 0, 1])):
            theta = np.arctan2(normal[1], normal[0])
            phi = np.arctan2(np.sqrt(normal[0] ** 2 + normal[1] ** 2), normal[2])
            self.rotate(phi, (0, 1, 0))
            self.rotate(theta, (0, 0, 1))
        if not np.allclose(np.asarray(center), 0):
            self.translate(center)
        return self

    @classmethod
    def make_splines_compatible(cls, spline1: SplineObject, spline2: SplineObject) -> None:
        """Ensure that two splines are compatible.

        This will manipulate one or both to ensure that they are both rational
        or nonrational, and that they lie in the same physical space.

        :param SplineObject spline1: The first spline
        :param SplineObject spline2: The second spline
        """
        # make both rational (if needed)
        if spline1.rational:
            spline2.force_rational()
        elif spline2.rational:
            spline1.force_rational()

        # make both in the same geometric space
        if spline1.dimension > spline2.dimension:
            spline2.set_dimension(spline1.dimension)
        else:
            spline1.set_dimension(spline2.dimension)

    @classmethod
    def make_splines_identical(
        cls, spline1: SplineObject, spline2: SplineObject, direction: Optional[Direction] = None
    ) -> None:
        """Ensure that two splines have identical discretization.

        This will first make them compatible (see
        :func:`splipy.SplineObject.make_curves_compatible`), reparametrize them, and
        possibly raise the order and insert knots as required.

        :param SplineObject spline1: The first spline
        :param SplineObject spline2: The second spline
        :param int direction: The direction to make identical. If
            None, make all directions identical.
        """

        # make sure that rational/dimension is the same
        SplineObject.make_splines_compatible(spline1, spline2)

        # If all directions, just call the same method several times
        if direction is None:
            for i in range(spline1.pardim):
                cls.make_splines_identical(spline1, spline2, direction=i)
            return

        # From this point, assume we're running on a single direction
        i = check_direction(direction, spline1.pardim)

        # make both have knot vectors in domain (0,1)
        spline1.reparam(direction=i)
        spline2.reparam(direction=i)

        # settle on the lowest periodicity if different appear
        if spline1.bases[i].periodic < spline2.bases[i].periodic:
            spline2.lower_periodic(spline1.bases[i].periodic, i)
        elif spline2.bases[i].periodic < spline1.bases[i].periodic:
            spline1.lower_periodic(spline2.bases[i].periodic, i)

        # make sure both have the same order
        p1 = spline1.order(i)
        p2 = spline2.order(i)
        p = max(p1, p2)
        spline1.raise_order(p - p1, direction=i)
        spline2.raise_order(p - p2, direction=i)

        # make sure both have the same knot vectors
        knot1 = spline1.knots(i)
        knot2 = spline2.knots(i)
        b1 = spline1.bases[i]
        b2 = spline2.bases[i]

        inserts = []
        for k in knot1:
            c1 = b1.continuity(k)
            c2 = b2.continuity(k)
            if c2 > c1:
                m = min(c2 - c1, p - 1 - c1)  # c2=np.inf if knot does not exist
                inserts.extend([k] * m)
        spline2.insert_knot(inserts, direction=i)

        inserts = []
        for k in knot2:
            c1 = b1.continuity(k)
            c2 = b2.continuity(k)
            if c1 > c2:
                m = min(c1 - c2, p - 1 - c2)  # c1=np.inf if knot does not exist
                inserts.extend([k] * m)
        spline1.insert_knot(inserts, direction=i)
