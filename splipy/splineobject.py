# -*- coding: utf-8 -*-

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray, ArrayLike
from numpy import floating
import copy
from operator import attrgetter, methodcaller
from itertools import chain, product
from bisect import bisect_left
from typing import (
    Sequence, Optional, List, overload, Union, Tuple,
    TypedDict, ClassVar, Protocol, cast, TypeVar, Literal,
    MutableSequence,
)
from typing_extensions import Unpack, Self

from . import state
from .basis import BSplineBasis
from .types import Direction, Seq, Float, OneOrMore
from .utils import (
    reshape, rotation_matrix, is_singleton, ensure_listlike,
    check_direction, ensure_flatlist, check_section, sections,
    raise_order_1D
)

__all__ = ['SplineObject']


def transpose_fix(pardim: int, direction: int) -> Tuple[int, ...]:
    ret = list(range(1, pardim+1))
    ret.insert(direction, 0)
    return tuple(ret)


def evaluate(bases: Sequence[NDArray[floating]], cps: NDArray[floating], tensor: bool = True) -> NDArray[floating]:
    if tensor:
        idx = len(bases) - 1
        for N in bases[::-1]:
            cps = np.tensordot(N, cps, axes=(1, idx))
    else:
        cps = np.einsum('ij,j...->i...', bases[0], cps)
        for N in bases[1:]:
            cps = np.einsum('ij,ij...->i...', N, cps)
    return cps


class SectionKwargs(TypedDict, total=False):
    u: Optional[int]
    v: Optional[int]
    w: Optional[int]


T = TypeVar("T", covariant=True)

class Constructor(Protocol[T]):
    def __call__(
        self,
        *bases: Optional[BSplineBasis],
        controlpoints: Optional[ArrayLike] = None,
        rational: bool = False,
        raw: bool = False,
    ) -> T: ...


class SplineObject:
    """Master class for spline objects with arbitrary dimensions.

    This class should be subclassed instead of used directly.

    All SplineObjects support basic arithmetic operators, which are interpreted
    as translation and scaling. In-place operators (e.g. ``+=``) mutate the
    object, while infix operators (e.g. ``+``) create new objects.
    """

    bases: List[BSplineBasis]
    controlpoints: NDArray[floating]
    dimension: int
    rational: bool

    _intended_pardim: ClassVar[int]

    def __init__(
        self,
        bases: Sequence[Optional[BSplineBasis]],
        controlpoints: Optional[ArrayLike] = None,
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
            cps = [c[::-1] for c in product(*(b.greville() for b in self.bases[::-1]))]

            # Minimum two dimensions
            if len(cps[0]) == 1:
                cps = [tuple(list(c) + [0.0]) for c in cps]

            # Add weight = 1 for identiy-mapping rational splines
            if rational:
                cps = [tuple(list(c) + [1.0]) for c in cps]

            self.controlpoints = np.array(cps)

        else:
            self.controlpoints = np.array(controlpoints)

        self.dimension = self.controlpoints.shape[-1] - rational
        self.rational = rational

        if not raw:
            shape = tuple(b.num_functions() for b in self.bases)
            ncomps = self.dimension + rational
            self.controlpoints = reshape(self.controlpoints, shape, order='F', ncomps=ncomps)

    @staticmethod
    def _subclass_constructor(pardim: int) -> Optional[Constructor[SplineObject]]:
        classes = [c for c in SplineObject.__subclasses__() if c._intended_pardim == pardim]
        if classes:
            return cast(Constructor[SplineObject], classes[0])
        return None

    @property
    def _construct(self) -> Constructor[Self]:
        return cast(Constructor[Self], self.__class__)

    def _validate_domain(self, *params: Seq[Float]) -> None:
        """Check whether the given evaluation parameters are valid.

        :raises ValueError: If the parameters are outside the domain
        """
        for b, p in zip(self.bases, params):
            if b.periodic < 0:
                if min(p) < b.start() - state.knot_tolerance or b.end() + state.knot_tolerance < max(p):
                    raise ValueError('Evaluation outside parametric domain')

    def evaluate(self, *params: OneOrMore[Float], tensor: bool = True) -> NDArray[floating]:
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
        eval_params = [ensure_listlike(p) for p in params]

        if not tensor and len({len(p) for p in eval_params}) != 1:
            raise ValueError('Parameters must have same length')

        self._validate_domain(*eval_params)

        # Evaluate the corresponding bases at the corresponding points
        # and build the result array
        Ns = [b.evaluate(p) for b, p in zip(self.bases, eval_params)]
        result = evaluate(Ns, self.controlpoints, tensor)

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
        *params: OneOrMore[Float],
        d: OneOrMore[int] = 1,
        above: OneOrMore[bool] = True,
        tensor: bool = True,
    ) -> NDArray[floating]:
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
        eval_params = [ensure_listlike(p) for p in params]

        derivs = ensure_listlike(d, self.pardim)
        eval_above = ensure_listlike(above, self.pardim)

        if not tensor and len({len(p) for p in eval_params}) != 1:
            raise ValueError('Parameters must have same length')

        self._validate_domain(*eval_params)

        # Evaluate the derivatives of the corresponding bases at the corresponding points
        # and build the result array
        dNs = [b.evaluate(p, d, from_right) for b, p, d, from_right in zip(self.bases, eval_params, derivs, eval_above)]
        result = evaluate(dNs, self.controlpoints, tensor)

        # For rational curves, we need to use the quotient rule
        # (n/W)' = (n' W - n W') / W^2 = n'/W - nW'/W^2
        # * n'(i) = result[..., i]
        # * W'(i) = result[..., -1]
        # We evaluate in the regular way to compute n and W.
        if self.rational:
            if sum(derivs) > 1:
                raise RuntimeError('Rational derivative not implemented for order %i' % sum(derivs))
            Ns = [b.evaluate(p) for b, p in zip(self.bases, eval_params)]
            non_derivative = evaluate(Ns, self.controlpoints, tensor)
            W = non_derivative[..., -1]  # W
            Wd = result[..., -1]         # W'
            for i in range(self.dimension):
                result[..., i] = result[..., i] / W - non_derivative[..., i] * Wd / W / W
            result = np.delete(result, self.dimension, -1)

        # Squeeze the singleton dimensions if we only have one point
        if squeeze:
            result = result.reshape(self.dimension)

        return result

    @overload
    def get_derivative_spline(self, direction: Direction) -> Self: ...

    @overload
    def get_derivative_spline(self) -> Tuple[Self]: ...

    def get_derivative_spline(self, direction = None):   # type: ignore[no-untyped-def]
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
            raise RuntimeError('Not working for rational splines')

        # if no direction is specified, return a tuple with all derivatives
        if direction is None:
            return tuple([self.get_derivative_spline(dim) for dim in range(self.pardim)])
        else:
            d = check_direction(direction, self.pardim)
            k = self.knots(d, with_multiplicities=True)
            p = self.order(d)-1
            n = self.shape[d]
            if self.bases[d].periodic < 0:
                C = np.zeros((n-1, n))
                for i in range(n-1):
                    C[i,i]   = -float(p) / (k[i+p+1] - k[i+1])
                    C[i,i+1] =  float(p) / (k[i+p+1] - k[i+1])
            else:
                C = np.zeros((n, n))
                for i in range(n):
                    ip1 = np.mod(i+1,n)
                    C[i,i]   = -float(p) / (k[i+p+1] - k[i+1])
                    C[i,ip1] =  float(p) / (k[i+p+1] - k[i+1])

            derivative_cps = np.tensordot(C, self.controlpoints, axes=(1, d))
            derivative_cps = derivative_cps.transpose(transpose_fix(self.pardim, d))
            bases    = [b for b in self.bases]
            bases[d] = BSplineBasis(p, k[1:-1], bases[d].periodic-1)

            return self._construct(*bases, controlpoints=derivative_cps, rational=self.rational, raw=True)

    @overload
    def tangent(
        self,
        *params: OneOrMore[Float],
        direction: Direction,
        above: OneOrMore[bool] = True,
        tensor: bool = True,
    ) -> NDArray[floating]: ...

    @overload
    def tangent(
        self,
        *params: OneOrMore[Float],
        above: OneOrMore[bool] = True,
        tensor: bool = True,
    ) -> NDArray[floating]: ...

    def tangent(self, *params, direction = None, above = True, tensor = True):  # type: ignore[no-untyped-def]
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
        eval_above = ensure_listlike(above, self.pardim)

        if self.pardim == 1: # curves
            direction = 0

        if direction is None:
            result = ()
            for i in range(self.pardim):
                derivative[i] = 1
                # compute velocity in this direction
                v = self.derivative(*params, d=derivative, above=eval_above, tensor=tensor)
                # normalize
                if len(v.shape)==1:
                    speed = np.linalg.norm(v)
                else:
                    speed = np.linalg.norm( v, axis=-1)
                    speed = np.reshape(speed, speed.shape +(1,))
                # store in result tuple
                result += (v/speed,)
                derivative[i] = 0
            return result

        i = check_direction(direction, self.pardim)
        derivative[i] = 1
        # compute velocity in this direction
        v = self.derivative(*params, d=derivative, above=eval_above, tensor=tensor)
        # normalize
        if len(v.shape)==1:
            speed = np.linalg.norm(v)
        else:
            speed = np.linalg.norm( v, axis=-1)
            speed = np.reshape(speed, speed.shape +(1,))

        return v / speed

    @overload
    def section(
        self,
        *args: Optional[int],
        unwrap_points: Literal[False] = False,
        **kwargs: Unpack[SectionKwargs],
    ) -> SplineObject: ...

    @overload
    def section(
        self,
        *args: Optional[int],
        unwrap_points: Literal[True] = True,
        **kwargs: Unpack[SectionKwargs],
    ) -> Union[SplineObject, NDArray[floating]]: ...

    def section(self, *args, unwrap_points=True, **kwargs):  # type: ignore[no-untyped-def]
        """Returns a section from the object. A section can be any sub-object of
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

        slices = tuple(slice(None) if p is None else p for p in section)
        bases = [b for b, p in zip(self.bases, section) if p is None]
        if bases or not unwrap_points:
            constructor = SplineObject._subclass_constructor(len(bases))
            if constructor:
                return constructor(*bases, controlpoints=self.controlpoints[slices], rational=self.rational, raw=True)
            return SplineObject(bases, self.controlpoints[slices], self.rational, raw=True)
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
        function raise_order_implicit is used otherwise.

        :param int u,v,...: Number of times to raise the order in a given
            direction.
        :param int direction: The direction in which to raise the order.
        :return: self
        """
        if len(raises) == 1 and direction is None:
            raises = (raises[0],) * self.pardim
        elif len(raises) == 1:
            newraises = [0] * self.pardim
            newraises[check_direction(direction, self.pardim)] = raises[0]
            raises = tuple(newraises)
        if not all(r >= 0 for r in raises):
            raise ValueError("Cannot lower order using raise_order")
        if all(r == 0 for r in raises):
            return self

        if any(b.continuity(b.knots[0]) < b.order or b.periodic > -1 for b in self.bases):
            self.raise_order_implicit(*raises)
            return self

        new_bases = [b.raise_order(r) for b, r in zip(self.bases, raises)]

        d_p = self.pardim

        controlpoints = self.controlpoints
        for i in range(0,d_p):
            dimensions = np.array(controlpoints.shape)
            indices = np.array(range(0,d_p+1))
            indices[i],indices[d_p] = d_p,i
            controlpoints = np.transpose(controlpoints,indices)
            controlpoints = np.reshape(controlpoints,(np.prod(dimensions[indices[:-1]]),dimensions[i]))
            controlpoints = raise_order_1D(controlpoints.shape[1]-1,self.order(i),
                                           self.bases[i].knots,controlpoints,raises[i],self.bases[i].periodic)
            controlpoints = np.reshape(controlpoints,np.append(dimensions[indices[:-1]],controlpoints.shape[1]))
            controlpoints = np.transpose(controlpoints,indices)

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
            result = np.tensordot(n, result, axes=(1, self.pardim-1))

        # Solve the interpolation problem
        for n in N_new[::-1]:
            result = np.tensordot(np.linalg.inv(n), result, axes=(1, self.pardim-1))

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
            new_controlpts = np.tensordot(n, new_controlpts, axes=(1, self.pardim-1))

        # Solve the interpolation problem
        for n in N_new[::-1]:
            new_controlpts = np.tensordot(np.linalg.inv(n), new_controlpts, axes=(1, self.pardim-1))

        return self._construct(*new_bases, controlpoints=new_controlpts, rational=self.rational, raw=True)

    @overload
    def start(self) -> Tuple[float]: ...

    @overload
    def start(self, direction: Direction) -> float: ...

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
    def end(self) -> Tuple[float]: ...

    @overload
    def end(self, direction: Direction) -> float: ...

    def end(self, direction=None):  # type: ignore[no-untyped-def]
        """  Return the end of the parametric domain.

        If `direction` is given, returns the end of that direction, as a float.
        If it is not given, returns the end of all directions, as a tuple.

        :param int direction: Direction in which to get the end.
        :raises ValueError: For invalid direction
        """
        if direction is None:
            return tuple(b.end() for b in self.bases)
        direction = check_direction(direction, self.pardim)
        return self.bases[direction].end()

    @overload
    def order(self) -> Tuple[int]: ...

    @overload
    def order(self, direction: Direction) -> int: ...

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
    def knots(self, *, with_multiplicities: bool = False) -> Tuple[NDArray[floating], ...]: ...

    @overload
    def knots(self, direction: Direction, *, with_multiplicities: bool = False) -> NDArray[floating]: ...

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
        getter = attrgetter('knots') if with_multiplicities else methodcaller('knot_spans')
        if direction is None:
            return tuple(getter(b) for b in self.bases)
        direction = check_direction(direction, self.pardim)
        return getter(self.bases[direction])

    def reverse(self, direction: Direction = 0) -> Self:
        """Swap the direction of a parameter by making it go in the reverse
        direction. The parametric domain remains unchanged.

        :param int direction: The direction to flip.
        :return: self
        """
        index = check_direction(direction, self.pardim)
        self.bases[index].reverse()

        # This creates the following slice programmatically
        # array[:, :, :, ..., ::-1,]
        # index=direction -----^
        # :    => slice(None, None, None)
        # ::-1 => slice(None, None, -1)
        slices = [slice(None, None, None) for _ in range(index)] + [slice(None, None, -1)]
        self.controlpoints = self.controlpoints[tuple(slices)]

        return self

    def swap(self, dir1: Direction = 0, dir2: Direction = 1) -> Self:
        """Swaps two parameter directions.

        This function silently passes for curves.

        :param direction dir1: The first direction (default u)
        :param direction dir2: The second direction (default v)
        :return: self
        """
        if self.pardim == 1:
            return self

        idx1 = check_direction(dir1, self.pardim)
        idx2 = check_direction(dir2, self.pardim)

        # Reorder control points
        new_directions = list(range(self.pardim + 1))
        new_directions[idx1] = idx2
        new_directions[idx2] = idx1
        self.controlpoints = self.controlpoints.transpose(new_directions)

        # Swap knot vectors
        self.bases[idx1], self.bases[idx2] = self.bases[idx2], self.bases[idx1]

        return self

    def insert_knot(self, knot: OneOrMore[Float], direction: Direction = 0) -> Self:
        """Insert a new knot into the spline.

        :param int direction: The direction to insert in
        :param knot: The new knot(s) to insert
        :type knot: float or [float]
        :raises ValueError: For invalid direction
        :return: self
        """
        shape  = self.controlpoints.shape

        # for single-value input, wrap it into a list
        knots = ensure_listlike(knot)

        index = check_direction(direction, self.pardim)

        C = np.identity(shape[index])
        for k in knots:
            C = self.bases[index].insert_knot(k) @ C
        self.controlpoints = np.tensordot(C, self.controlpoints, axes=(1, index))
        self.controlpoints = self.controlpoints.transpose(transpose_fix(self.pardim, index))

        return self

    @overload
    def refine(self, n: int, direction: Direction) -> Self: ...

    @overload
    def refine(self, *n: int) -> Self: ...

    def refine(self, *ns, direction = None):  # type: ignore[no-untyped-def]
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
        if len(ns) == 1 and direction is not None:
            directions = [check_direction(direction, self.pardim)]
        else:
            directions = range(self.pardim)

        if len(ns) == 1:
            ns = (ns[0],) * self.pardim

        for n, d in zip(ns, directions):
            knots = self.knots(direction=d)  # excluding multiple knots
            new_knots = []
            for (k0, k1) in zip(knots[:-1], knots[1:]):
                new_knots.extend(np.linspace(k0, k1, n+2)[1:-1])
            self.insert_knot(new_knots, d)

        return self

    @overload
    def reparam(self, *args: Tuple[float, float]) -> Self: ...

    @overload
    def reparam(self, args: Tuple[float, float], direction: Direction) -> Self: ...

    @overload
    def reparam(self, direction: Direction) -> Self: ...

    def reparam(self, *args, direction=None):  # type: ignore[no-untyped-def]
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
        if direction is None:
            # Pad the args with (0, 1) for the extra directions
            args = list(args) + [(0, 1)] * (len(self.bases) - len(args))
            for b, (start, end) in zip(self.bases, args):
                b.reparam(start, end)
        else:
            index = check_direction(direction, self.pardim)
            if len(args) == 0:
                self.bases[index].reparam(0,1)
            else:
                start, end = args[0]
                self.bases[index].reparam(start, end)

        return self

    def translate(self, x: ArrayLike) -> Self:
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

        xx = np.array(x)
        if len(xx) > dim:  # typical case: requesting movement in z-direction for 2D geometries
            self.set_dimension(len(xx))
            dim = self.dimension

        # set up the translation matrix
        translation_matrix = np.identity(dim + 1)
        for i in range(dim):
            translation_matrix[i, -1] = xx[i]

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

    def scale(self, *args: float) -> Self:
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
        s = ensure_flatlist(args)
        s = ensure_listlike(s, dups=3)

        # set up the scaling matrix
        scale_matrix = np.identity(dim + rat)
        for i in range(dim):
            scale_matrix[i, i] = s[i]

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual scaling by matrix-matrix multiplication
        cp = cp @ scale_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def rotate(self, theta: float, normal: ArrayLike = (0, 0, 1)) -> Self:
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

        norm = np.array(normal)
        if not (norm[0] == 0 and norm[1] == 0):  # rotating a 2D geometry out of the xy-plane
            self.set_dimension(3)
            dim = self.dimension

        # set up the rotation matrix
        if dim == 2:
            R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]
                 ]).T  # we do right-multiplication, so we need a transpose
        elif dim == 3:
            R = rotation_matrix(theta, norm)
        else:
            raise RuntimeError('rotation undefined for geometries other than 2D and 3D')

        rot_matrix = np.identity(dim + rat)
        rot_matrix[0:dim, 0:dim] = R

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual rotation by matrix-matrix multiplication
        cp = cp @ rot_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def mirror(self, normal: ArrayLike) -> Self:
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
            raise RuntimeError('reflection undefined for geometries other than 3D')

        # fixup the input normal to right form
        norm = np.array(normal)
        norm = norm / np.sqrt(np.dot(norm, norm))  # normalize it

        # set up the reflection matrix
        reflection_matrix = np.identity(dim + rat)
        reflection_matrix[0:dim, 0:dim] -= 2 * np.outer(norm, norm)

        # wrap out the controlpoints to a matrix (down from n-D tensor)
        cp = np.reshape(self.controlpoints, (n, dim + rat))

        # do the actual rotation by matrix-matrix multiplication
        cp = cp @ reflection_matrix

        # store results
        self.controlpoints = np.reshape(np.array(cp), self.controlpoints.shape)

        return self

    def project(self, plane: str) -> Self:
        """  Projects the geometry onto a plane or axis.

        - `project('xy')` will project the object onto the *xy* plane, setting
          all *z* components to zero.
        - `project('y')` will project the object onto the *y* axis, setting all
          *x* and *z* components to zero.

        :param string plane: Any combination of 'x', 'y' and 'z'
        :return: self
        """
        keep = [c in plane.lower() for c in 'xyz']

        dim = self.dimension
        for i in range(dim):
            if not keep[i]:
                self.controlpoints[..., i] = 0

        return self

    def bounding_box(self) -> List[Tuple[float, float]]:
        """Gets the bounding box of a spline object, computed from the
        control-point values. Could be inaccurate for rational splines.

        Returns the minima and maxima for each direction:
        [(xmin, xmax), (ymin, ymax), ...]

        :return: Bounding box
        :rtype: [(float)]
        """
        dim = self.dimension

        result = []
        for i in range(dim):
            result.append((np.min(self.controlpoints[..., i]),
                           np.max(self.controlpoints[..., i])))
        return result

    def center(self) -> NDArray[floating]:
        """Gets the center of the domain

        For curves this will return :math:`(\\tilde{x}, \\tilde{y},...)`, where

        .. math:: \\tilde{x} = \\frac{1}{L} \\int_{t_0}^{t_1} x(t) \\; dt

        and :math:`L=t_1-t_0` is the length of the parametric domain :math:`[t_0,t_1]`.

        For surfaces this will return :math:`(\\tilde{x}, \\tilde{y},...)`, where

        .. math:: \\tilde{x} = \\frac{1}{A} \\int_{v_0}^{v_1} \\int_{u_0}^{u_1} x(u,v) \\; du \\; dv

        and :math:`A=(u_1-u_0)(v_1-v_0)` is the area of the parametric domain :math:`[u_0,u_1]\\times[v_0,v_1]`.

        .. warning:: For rational splines, this will integrate in projective
            coordinates, then project the centerpoint. This is as opposed to
            integrating the rational functions :math:`\\frac{N_i(t)w_i}{\\sum_j
            N_j(t)w_j}`.
        """

        # compute integration of basis functions
        Ns = [b.integrate(b.start(), b.end()) for b in self.bases]

        # compute parametric size
        par_size = np.prod([t1-t0 for (t0,t1) in zip(self.start(), self.end())])

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

    def corners(self, order: Literal['F', 'C'] = 'C') -> NDArray[floating]:
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
        result = np.zeros((2**self.pardim, self.dimension + int(self.rational)))
        for i, args in enumerate(sections(self.pardim, 0)):
            result[i,:] = self.section(*(args[::-1] if order == 'F' else args))
        return result

    def lower_periodic(self, periodic: int, direction: Direction = 0) -> Self:
        """Sets the periodicity of the spline object in the given direction,
        keeping the geometry unchanged.

        :param int periodic: new periodicity, i.e. the basis is C^k over the start/end
        :param int direction: the parametric direction of the basis to modify
        :return: self
        """
        index = check_direction(direction, self.pardim)

        b  = self.bases[index]
        while periodic < b.periodic:
            self.insert_knot(self.start(index), index)
            self.controlpoints = np.roll(self.controlpoints, -1, index)
            b.roll(1)
            b.periodic -= 1
            b.knots = b.knots[:-1]
        if periodic > b.periodic:
            raise ValueError('Cannot raise periodicity')

        return self

    def set_dimension(self, new_dim: int) -> Self:
        """Sets the physical dimension of the object. If increased, the new
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
        """Returns true if the spline object is periodic in the given parametric direction"""
        index = check_direction(direction, self.pardim)
        return self.bases[index].periodic > -1

    def force_rational(self) -> Self:
        """Force a rational representation of the object.

        The weights of a non-rational object will be set to 1.

        :return: self
        """
        if self.rational:
            return self

        dim = self.dimension
        shape = self.controlpoints.shape
        self.controlpoints = np.insert(self.controlpoints, dim, np.ones(shape[:-1]), self.pardim)
        self.rational = True
        return self

    def split(self, knots: OneOrMore[Float], direction: Direction = 0) -> Union[Self, List[Self]]:
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
        vec_knots = ensure_listlike(knots)
        index = check_direction(direction, self.pardim)

        p = self.order(index)
        results = []
        splitting_obj = self.clone()
        bases = self.bases
        # insert knots to produce C{-1} at all splitting points
        for k in vec_knots:
            continuity = bases[index].continuity(k)
            if continuity == np.inf:
                continuity = p - 1
            splitting_obj.insert_knot([k] * (continuity + 1), index)

        b = splitting_obj.bases[index]
        if b.periodic > -1:
            mu = bisect_left(b.knots, vec_knots[0])
            b.roll(mu)
            splitting_obj.controlpoints = np.roll(splitting_obj.controlpoints, -mu, index)
            b.knots = b.knots[:-b.periodic-1]
            b.periodic = -1
            if len(vec_knots) > 1:
                return splitting_obj.split(vec_knots[1:], index)
            else:
                return splitting_obj

        # everything is available now, just have to find the right index range
        # in the knot vector and controlpoints to store in each separate curve
        # piece
        last_cp_i = 0
        last_knot_i = 0

        bases = splitting_obj.bases
        b     = bases[index]
        cp_slice = [slice(None, None, None)] * len(self.controlpoints.shape)
        for k in vec_knots:
            if self.start(index) < k < self.end(index): # skip start/end points
                mu = bisect_left(b.knots, k)
                n_cp = mu - last_knot_i
                knot_slice      = slice(last_knot_i, mu+p, None)
                cp_slice[index] = slice(last_cp_i,   last_cp_i+n_cp,  None)

                cp = splitting_obj.controlpoints[ tuple(cp_slice) ]
                bases[index] = BSplineBasis(p, b.knots[knot_slice])

                results.append(self._construct(*bases, controlpoints=cp, rational=splitting_obj.rational, raw=True))

                last_knot_i = mu
                last_cp_i += n_cp

        # with n splitting points, we're getting n+1 pieces. Add the final one:
        knot_slice      = slice(last_knot_i, None, None)
        cp_slice[index] = slice(last_cp_i,   None, None)
        bases[index] = BSplineBasis(p, b.knots[knot_slice])
        cp = splitting_obj.controlpoints[ tuple(cp_slice) ]
        results.append(self._construct(*bases, controlpoints=cp, rational=splitting_obj.rational, raw=True))

        return results

    def make_periodic(self, continuity: Optional[int] = None, direction: Direction = 0) -> Self:
        """Make the spline object periodic in a given parametric direction.

        :param int continuity: The continuity along the boundary (default max).
        :param int direction: The direction to ensure continuity in.
        """

        index = check_direction(direction, self.pardim)
        basis = self.bases[index]
        cont = continuity if continuity is not None else basis.order - 2
        if not -1 <= cont <= basis.order - 2:
            raise ValueError('Illegal continuity for basis of order {}: {}'.format(
                cont, basis.order
            ))
        if cont == -1:
            raise ValueError(
                'Operation not supported. '
                'For discontinuous spline spaces, consider SplineObject.split().'
            )
        if basis.periodic >= 0:
            raise ValueError('Basis is already periodic')

        basis = basis.make_periodic(cont)

        # Merge control points
        index_beg: List[Union[slice, int]] = [slice(None,None,None)] * (self.pardim + 1)
        index_end: List[Union[slice, int]] = [slice(None,None,None)] * (self.pardim + 1)
        cps = np.array(self.controlpoints)
        weights = np.linspace(0, 1, cont + 1) if cont > 0 else np.array([0.5])
        for i, j, t in zip(range(cont + 1), range(-cont-1, 0), weights):
            # Weighted average between cps[..., i, ..., :] and cps[..., -c-1+i, ..., :]
            # The weights are chosen so that, for periodic c, the round trip
            # c.split().make_periodic() with suitable arguments produces an
            # object identical to c. (Mostly black magic.)
            index_beg[index] = i
            index_end[index] = j
            cps[tuple(index_beg)] = t * cps[tuple(index_beg)] + (1 - t) * cps[tuple(index_end)]

        # cps[..., :-(continuity+1), ..., :]
        index_beg[index] = slice(None, -(cont + 1), None)
        cps = cps[tuple(index_beg)]

        bases = list(self.bases)
        bases[index] = basis
        return self._construct(*bases, controlpoints=cps, rational=self.rational, raw=True)

    @property
    def pardim(self) -> int:
        """The number of parametric dimensions: 1 for curves, 2 for surfaces, 3
        for volumes, etc.
        """
        return len(self.controlpoints.shape)-1

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

    def _unravel_flat_index(self, i: Union[int, slice]) -> Tuple[NDArray[np.intp], ...]:
        """Unravels a flat index i to multi-indexes.

        :param i: Flat index
        :type i: int or slice
        :rtype: Tuple of np.array
        :raises IndexError: If the index is out of bounds
        """
        # i is int => make sure we deal with negative i properly
        # i is slice => use i.indices to compute the actual indices
        total = len(self)
        if isinstance(i, int):
            indexes = [i] if i >= 0 else [total + i]
        else:
            indexes = list(range(*i.indices(total)))

        # Convert to multi-indexes
        try:
            unraveled = np.unravel_index(indexes, self.controlpoints.shape[:-1], order='F')
        except ValueError:
            raise IndexError

        return unraveled

    def __getitem__(self, i: Union[int, slice, Tuple[Union[int, slice], ...]]) -> NDArray[floating]:
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
            return self.controlpoints[unraveled][0]
        return self.controlpoints[unraveled]

    def __setitem__(self, i: Union[int, slice, Tuple[Union[int, slice], ...]], cp: ArrayLike) -> None:
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
    def shape(self) -> Tuple[int, ...]:
        """The dimensions of the control point array."""
        return self.controlpoints.shape[:-1]

    def __iadd__(self, x: ArrayLike) -> Self:
        self.translate(x)
        return self

    def __isub__(self, x: ArrayLike) -> Self:
        self.translate(-np.array(x))  # can't do -x if x is a list, so we rewrap it here
        return self

    def __imul__(self, x: float) -> Self:
        self.scale(x)
        return self

    def __itruediv__(self, x: float) -> Self:
        self.scale(1.0 / x)
        return self

    __ifloordiv__ = __itruediv__  # integer division (should not distinguish)
    __idiv__ = __itruediv__  # python2 compatibility

    def __add__(self, x: ArrayLike) -> Self:
        new_obj = copy.deepcopy(self)
        new_obj += x
        return new_obj

    def __radd__(self, x: ArrayLike) -> Self:
        return self + x

    def __sub__(self, x: ArrayLike) -> Self:
        new_obj = copy.deepcopy(self)
        new_obj -= x
        return new_obj

    def __mul__(self, x: float) -> Self:
        new_obj = copy.deepcopy(self)
        new_obj *= x
        return new_obj

    def __rmul__(self, x: float) -> Self:
        return self * x

    def __truediv__(self, x: float) -> Self:
        new_obj = copy.deepcopy(self)
        new_obj /= x
        return new_obj

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
        cls,
        spline1: SplineObject,
        spline2: SplineObject,
        direction: Optional[Direction] = None,
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
        knot1 = spline1.knots(direction=i)
        knot2 = spline2.knots(direction=i)
        b1    = spline1.bases[i]
        b2    = spline2.bases[i]

        inserts = []
        for k in knot1:
            c1 = b1.continuity(k)
            c2 = b2.continuity(k)
            if c2 > c1:
                m = cast(int, min(c2-c1, p-1-c1)) # c2=np.inf if knot does not exist
                inserts.extend([k] * m)
        spline2.insert_knot(inserts, direction=i)

        inserts = []
        for k in knot2:
            c1 = b1.continuity(k)
            c2 = b2.continuity(k)
            if c1 > c2:
                m = cast(int, min(c1-c2, p-1-c2)) # c1=np.inf if knot does not exist
                inserts.extend([k]*m)
        spline1.insert_knot(inserts, direction=i)
