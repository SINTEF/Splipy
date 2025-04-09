from __future__ import annotations

import copy
from bisect import bisect_left, bisect_right
from itertools import chain, repeat
from typing import TYPE_CHECKING, Literal, Self, cast, overload

import numpy as np
import numpy.typing as npt
import splipy_core
from deprecated import deprecated
from scipy.sparse import csr_matrix

from . import state
from .utils import ensure_listlike_old

if TYPE_CHECKING:
    from .typing import ArrayLike, FloatArray, ScalarLike

__all__ = ["BSplineBasis"]


class NotAKnotError(ValueError):
    pass


class BSplineBasis:
    """BSplineBasis()

    Represents a one-dimensional B-Spline basis.

    BSplineBasis objects support basic arithmetic operators, which are
    interpreted as acting on the parametric domain.
    """

    __slots__ = ["knots", "order", "periodic"]
    knots: FloatArray
    order: int
    periodic: int

    def __init__(
        self,
        order: int = 2,
        knots: ArrayLike | None = None,
        periodic: int = -1,
    ) -> None:
        """Construct a B-Spline basis with a given order and knot vector.

        :param int order: Spline order, i.e. one greater than the polynomial degree.
        :param [float] knots: Knot vector of non-decreasing components.
            Defaults to open knot vector on domain [0,1].
        :param int periodic: Number of continuous derivatives at start and end.
            --1 is not periodic, 0 is continuous, etc.
        :raises ValueError: for inapproriate knot vectors
        """
        if order < 1:
            raise ValueError("invalid spline order")

        periodic = max(periodic, -1)
        if knots is None:
            knots = [0] * order + [1] * order
            for i in range(periodic + 1):
                knots[i] = -1
                knots[-i - 1] = 2

        self.knots = np.asarray(knots, dtype=np.float64)
        self.order = order
        self.periodic = periodic

        # error test input
        if len(self.knots) < 2 * order:
            raise ValueError("knot vector has too few elements")
        for i in range(len(self.knots) - 1):
            if self.knots[i + 1] - self.knots[i] < -state.knot_tolerance:
                raise ValueError("knot vector needs to be non-decreasing")
        if periodic >= 0:
            kts = self.knots
            for i in range(order + periodic - 1):
                diff = abs(
                    (kts[i + 1] - kts[i]) - (kts[-order - periodic + i] - kts[-order - periodic - 1 + i])
                )
                if diff > state.knot_tolerance:
                    raise ValueError("periodic knot vector is mis-matching at the start/end")

    def num_functions(self) -> int:
        """Returns the number of basis functions in the basis.

        .. warning:: This is different from :func:`splipy.BSplineBasis.__len__`."""
        return len(self.knots) - self.order - (self.periodic + 1)

    def start(self) -> float:
        """Start point of parametric domain. For open knot vectors, this is the
        first knot.

        :return: Knot number *p*, where *p* is the spline order
        :rtype: float
        """
        return float(self.knots.flat[self.order - 1])

    def end(self) -> float:
        """End point of parametric domain. For open knot vectors, this is the
        last knot.

        :return: Knot number *n*--*p*, where *p* is the spline order and *n* is
            the number of knots
        :rtype: Float
        """
        return float(self.knots.flat[-self.order])

    def greville_all(self) -> FloatArray:
        """Fetch all greville points, also known as knot averages:

        .. math:: \\sum_{j=i+1}^{i+p-1} \\frac{t_j}{p-1}

        :return: Greville points
        :rtype: [float]
        """
        operator = np.ones((self.order - 1,), dtype=np.float64) / (self.order - 1)
        return np.convolve(self.knots[1 : -1 - (self.periodic + 1)], operator, mode="valid")

    def greville_single(self, index: int) -> float:
        """Fetch a greville point, also known as a knot averages:

        .. math:: \\sum_{j=i+1}^{i+p-1} \\frac{t_j}{p-1}

        :param index: The index of the Greville point to calculate
        :type index: int
        :return: A Greville point
        :rtype: float
        """
        return float(np.sum(self.knots[index + 1 : index + self.order]) / (self.order - 1))

    @overload
    def greville(self, index: int) -> float: ...

    @overload
    def greville(self) -> FloatArray: ...

    def greville(self, index: int | None = None) -> float | FloatArray:
        """Fetch greville points, also known as knot averages:

        .. math:: \\sum_{j=i+1}^{i+p-1} \\frac{t_j}{p-1}

        :return: One, or all of the Greville points
        :rtype: [float] (if *index* is ``None``) or float
        """
        if index is not None:
            return self.greville_single(index)
        return self.greville_all()

    @overload
    def evaluate(
        self,
        t: ArrayLike | ScalarLike,
        d: int = 0,
        from_right: bool = ...,
    ) -> npt.NDArray[np.double]: ...

    @overload
    def evaluate(
        self,
        t: ArrayLike | ScalarLike,
        d: int = 0,
        from_right: bool = ...,
        sparse: Literal[False] = ...,
    ) -> npt.NDArray[np.double]: ...

    @overload
    def evaluate(
        self,
        t: ArrayLike | ScalarLike,
        d: int = 0,
        from_right: bool = ...,
        sparse: Literal[True] = ...,
    ) -> csr_matrix[np.float64]: ...

    def evaluate(
        self,
        t: ArrayLike | ScalarLike,
        d: int = 0,
        from_right: bool = True,
        sparse: bool = False,
    ) -> npt.NDArray[np.double] | csr_matrix[np.float64]:
        """Evaluate all basis functions in a given set of points.

        :param t: The parametric coordinate(s) in which to evaluate
        :type t: float or [float]
        :param int d: Number of derivatives to compute
        :param bool from_right: True if evaluation should be done in the limit
            from above
        :param bool sparse: True if computed matrix should be returned as sparse
        :return: A matrix *N[i,j]* of all basis functions *j* evaluated in all
            points *i*
        :rtype: numpy.array
        """
        t_arr = np.atleast_1d(np.asarray(t, dtype=np.float64))
        splipy_core.snap(self.knots, t_arr, state.knot_tolerance)

        if self.order <= d:  # requesting more derivatives than polynomial degree: return all zeros
            return np.zeros((len(t_arr), self.num_functions()))

        (data, size) = splipy_core.evaluate(
            self.knots,
            self.order,
            t_arr,
            self.periodic,
            state.knot_tolerance,
            d,
            from_right,
        )

        N = csr_matrix(data, size)
        return N.toarray() if not sparse else N

    @deprecated(version="1.11", reason="Use .evaluate()")
    def evaluate_old(self, t, d=0, from_right=True, sparse=False):  # type: ignore[no-untyped-def]
        """Evaluate all basis functions in a given set of points.
        :param t: The parametric coordinate(s) in which to evaluate
        :type t: float or [float]
        :param int d: Number of derivatives to compute
        :param bool from_right: True if evaluation should be done in the limit
            from above
        :param bool sparse: True if computed matrix should be returned as sparse
        :return: A matrix *N[i,j]* of all basis functions *j* evaluated in all
            points *i*
        :rtype: numpy.array
        """
        # for single-value input, wrap it into a list so it don't crash on the loop below
        t = ensure_listlike_old(t)
        self.snap(t)

        p = self.order  # knot vector order
        n_all = len(self.knots) - p  # number of basis functions (without periodicity)
        n = len(self.knots) - p - (self.periodic + 1)  # number of basis functions (with periodicity)
        m = len(t)
        data = np.zeros(m * p)
        indices = np.zeros(m * p, dtype="int32")
        indptr = np.array(range(0, m * p + 1, p), dtype="int32")
        if p <= d:  # requesting more derivatives than polymoial degree: return all zeros
            return np.zeros((m, n))
        if self.periodic >= 0:
            t = copy.deepcopy(t)
            # Wrap periodic evaluation into domain
            for i in range(len(t)):
                if t[i] < self.start() or t[i] > self.end():
                    t[i] = (t[i] - self.start()) % (self.end() - self.start()) + self.start()
        for i in range(len(t)):
            right = from_right
            evalT = t[i]
            # Special-case the endpoint, so the user doesn't need to
            if abs(t[i] - self.end()) < state.knot_tolerance:
                right = False
            # Skip non-periodic evaluation points outside the domain
            if t[i] < self.start() or t[i] > self.end():
                continue

            # mu = index of last non-zero basis function
            mu = bisect_right(self.knots, evalT) if right else bisect_left(self.knots, evalT)
            mu = min(mu, n_all)

            M = np.zeros(p)  # temp storage to keep all the function evaluations
            M[-1] = 1  # the last entry is a dummy-zero which is never used
            for q in range(1, p - d):
                for j in range(p - q - 1, p):
                    k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                    if j != p - q - 1:
                        M[j] = M[j] * float(evalT - self.knots[k]) / (self.knots[k + q] - self.knots[k])

                    if j != p - 1:
                        M[j] = M[j] + M[j + 1] * float(self.knots[k + q + 1] - evalT) / (
                            self.knots[k + q + 1] - self.knots[k + 1]
                        )

            for q in range(p - d, p):
                for j in range(p - q - 1, p):
                    k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                    if j != p - q - 1:
                        M[j] = M[j] * float(q) / (self.knots[k + q] - self.knots[k])
                    if j != p - 1:
                        M[j] = M[j] - M[j + 1] * float(q) / (self.knots[k + q + 1] - self.knots[k + 1])

            data[i * p : (i + 1) * p] = M
            indices[i * p : (i + 1) * p] = np.arange(mu - p, mu) % n

        N = csr_matrix((data, indices, indptr), (m, n))
        return N.toarray() if sparse else N

    def integrate(self, t0: ScalarLike, t1: ScalarLike) -> FloatArray:
        """Integrate all basis functions over a given domain

        :param float t0: The parametric starting point
        :param float t1: The parametric end point
        :return: The integration of all functions over the input domain
        :rtype: list
        """
        start = np.double(t0)
        end = np.double(t1)

        if self.periodic > -1 and (start < self.start() or end > self.end()):
            raise NotImplementedError("Periodic functions integrated across seam")

        start = np.fmax(start, self.start())
        end = np.fmin(end, self.end())
        p = self.order

        knots = np.empty_like(self.knots, shape=(2 + len(self.knots),))
        knots[1:-1] = self.knots
        knots[0] = self.knots[0]
        knots[-1] = self.knots[-1]

        integration_basis = BSplineBasis(p + 1, knots)
        ib_eval = integration_basis.evaluate([start, end])
        sumdiff: npt.NDArray[np.double] = ib_eval[1] - ib_eval[0]
        sumdiff = np.cumsum(sumdiff[::-1])[::-1]

        N = (knots[p + 1 : -1] - knots[1 : -p - 1]) / p * sumdiff[1:]

        # collapse periodic functions onto themselves
        if self.periodic > -1:
            N[: self.periodic + 1] += N[-self.periodic - 1 :]
            N = N[: -self.periodic - 1]

        return N

    def normalize(self) -> None:
        """Set the parametric domain to be (0,1)."""
        self -= self.start()  # set start-point to 0
        self /= self.end()  # set end-point to 1

    def reparam(self, start: ScalarLike = 0, end: ScalarLike = 1) -> None:
        """Set the parametric domain to be (start, end)

        :raises ValueError: If *end* ≤ *start*
        """
        start = float(start)
        end = float(end)

        if end <= start:
            raise ValueError("end must be larger than start")
        self.normalize()
        self *= end - start
        self += start

    def reverse(self) -> None:
        """Reverse parametric domain, keeping start/end values unchanged."""
        a = self.start()
        b = self.end()
        self.knots = (self.knots[::-1] - a) / (b - a) * (a - b) + b

    def knot_continuity(self, knot: ScalarLike) -> int:
        """Get the continuity of the basis functions at a given point.

        This is similar to continuity(), but throws an error if the input is not
        a knot value. Consequently, the return type is guaranteed to be an integer.

        :return: *p*--*m*--1 at a knot with multiplicity *m*, or ``inf`` between
            knots.
        :rtype: int or float
        """
        knot = float(knot)

        if self.periodic >= 0:
            if knot < self.start() or knot > self.end():
                knot = (knot - self.start()) % (self.end() - self.start()) + self.start()
        elif knot < self.start() or self.end() < knot:
            raise ValueError("out of range")

        # First knot that is larger than the right tolerance point
        hi = bisect_left(self.knots, knot + state.knot_tolerance)

        # Last knot that is smaller than the left tolerance point
        lo = bisect_left(self.knots, knot - state.knot_tolerance)

        if hi == lo:
            raise NotAKnotError
        return self.order - (hi - lo) - 1

    def continuity(self, knot: ScalarLike) -> int | float:
        """Get the continuity of the basis functions at a given point.

        :return: *p*--*m*--1 at a knot with multiplicity *m*, or ``inf``
            between knots.
        :rtype: int or float
        """
        try:
            return self.knot_continuity(knot)
        except NotAKnotError:
            return np.inf

    def make_periodic(self, continuity: int) -> BSplineBasis:
        """Create a periodic basis with a given continuity."""
        deg = self.order - 1
        new_knots = self.knots[deg:-deg]

        diff = self.end() - self.start()
        n_reps = deg - continuity - 1
        n_copy = deg - n_reps

        head = new_knots[-n_copy - 1 : -1] - diff
        tail = new_knots[1 : n_copy + 1] + diff

        new_knots = np.hstack((head, [self.start()] * n_reps, new_knots, [self.end()] * n_reps, tail))
        return BSplineBasis(self.order, new_knots, continuity)

    def knot_spans(self, include_ghost_knots: bool = False) -> FloatArray:
        """Return the set of unique knots in the knot vector.

        :param bool include_ghost_knots: if knots outside start/end are to be
            included. These knots are used by periodic bases.
        :return: List of unique knots
        :rtype: [float]"""
        p = self.order
        haystack = self.knots if include_ghost_knots else self.knots[p - 1 : -p + 1]

        result: list[np.double] = [haystack[0]]
        for k in haystack[1:]:
            if abs(k - result[-1]) > state.knot_tolerance:
                result.append(k)

        return np.array(result, dtype=np.double)

    def raise_order(self, amount: int) -> BSplineBasis:
        """Create a knot vector with higher order.

        The continuity at the knots are kept unchanged by increasing their
        multiplicities.

        :return: New knot vector
        :rtype: [float]
        :raises ValueError: If `amount` is negative
        """
        if amount < 0:
            raise ValueError("amount needs to be a non-negative integer")
        if amount == 0:
            return self.clone()

        # Accumulate the knots in order, by picking lowest of the next uniqe
        # know (repeated) or the next non-unique knot (not repeated)
        knot_spans: list[float] = list(self.knot_spans(include_ghost_knots=True))
        knot_spans_queue = knot_spans[::-1]
        my_knots_queue: list[float] = cast("list[float]", list(self.knots)[::-1])
        new_knots = []

        while knot_spans_queue or my_knots_queue:
            if my_knots_queue and not knot_spans_queue:
                new_knots.extend(my_knots_queue[::-1])
                break
            if (knot_spans_queue and not my_knots_queue) or (knot_spans_queue[-1] < my_knots_queue[-1]):
                new_knots.extend(repeat(knot_spans_queue.pop(), amount))
            else:
                new_knots.append(my_knots_queue.pop())

        # Remove excessive ghost knots which appear at both ends of the knot vector
        if self.periodic > -1:
            n0 = bisect_left(knot_spans, self.start())
            n1 = len(knot_spans) - bisect_left(knot_spans, self.end()) - 1
            new_knots = new_knots[n0 * amount : -n1 * amount]

        return BSplineBasis(self.order + amount, new_knots, self.periodic)

    def lower_order(self, amount: int) -> BSplineBasis:
        """Create a knot vector with lower order.

        The continuity at the knots are kept unchanged by decreasing their
        multiplicities.

        :return: New knot vector
        :rtype: [float]
        :raises TypeError: If `amount` is not an int
        :raises ValueError: If `amount` is negative
        """
        if amount < 0:
            raise ValueError("amount needs to be a non-negative integer")
        if self.order - amount < 2:
            raise ValueError("cannot lower order to less than linears")

        p = self.order - amount
        knots = list(
            chain.from_iterable(
                repeat(k, max(1, p - 1 - self.knot_continuity(k)))
                for k in self.knot_spans(include_ghost_knots=True)
            )
        )

        # Remove excessive ghost knots which appear at both ends of the knot vector
        if self.periodic > -1:
            n0 = bisect_left(knots, self.start())
            n1 = len(knots) - bisect_left(knots, self.end()) - 1
            knots = knots[n0 * amount : -n1 * amount]

        return BSplineBasis(p, knots, self.periodic)

    def insert_knot(self, new_knot: ScalarLike) -> FloatArray:
        """Inserts a knot in the knot vector.

        The return value is a sparse matrix *C* (actually, a dense matrix with
        lots of zeros), such that *N_new* = *N_old* x *C*, where *N* are row
        vectors of basis functions.

        :param float new_knot: The parametric coordinate of the point to insert
        :return: Transformation matrix *C*
        :rtype: numpy.array
        :raises ValueError: If the new knot is outside the domain
        """
        new_knot = float(new_knot)

        if self.periodic >= 0:
            if new_knot < self.start() or new_knot > self.end():
                new_knot = (new_knot - self.start()) % (self.end() - self.start()) + self.start()
        elif new_knot < self.start() or self.end() < new_knot:
            raise ValueError("new_knot out of range")

        # mu is the index of last non-zero (old) basis function
        mu = bisect_right(self.knots, new_knot)
        n = self.num_functions()
        p = self.order
        C = np.zeros((n + 1, n))
        # the modulus operator i%n in the C-matrix is needed for periodic basis functions
        for i in range(mu - p):
            C[i % (n + 1), i % n] = 1
        for i in range(mu - p, mu):
            if self.knots[i + p - 1] <= new_knot and new_knot <= self.knots[i + p]:
                C[i % (n + 1), i % n] = 1
            else:
                C[i % (n + 1), i % n] = (new_knot - self.knots[i]) / (self.knots[i + p - 1] - self.knots[i])
            if self.knots[i] <= new_knot and new_knot <= self.knots[i + 1]:
                C[(i + 1) % (n + 1), i % n] = 1
            else:
                C[(i + 1) % (n + 1), i % n] = (self.knots[i + p] - new_knot) / (
                    self.knots[i + p] - self.knots[i + 1]
                )
        for i in range(mu, n + 1):
            C[i % (n + 1), (i - 1) % n] = 1

        self.knots = np.insert(self.knots, mu, new_knot)

        # make sure that it is correct periodic after knot insertion
        if self.periodic > -1:
            m = len(self.knots)
            r = self.periodic
            if mu <= p + r:  # need to fix ghost knots on right side
                k0 = self.knots[0]
                k1 = self.knots[-p - r - 1]
                for i in range(p + r + 1):
                    self.knots[m - p - r - 1 + i] = k1 + (self.knots[i] - k0)
            elif mu >= m - p - r - 1:  # need to fix ghost knots on left side
                k0 = self.knots[p + r]
                k1 = self.knots[-1]
                for i in range(p + r + 1):
                    self.knots[i] = k0 - (k1 - self.knots[m - p - r - 1 + i])
        return C

    def roll(self, new_start: int) -> None:
        """Rotate a periodic knot vector by setting a new starting index.

        :param int new_start: The index of to the new first knot
        """
        if self.periodic < 0:
            raise RuntimeError("roll only applicable for periodic knot vectors")

        p = self.order
        k = self.periodic
        n = len(self.knots)
        t1 = self.knots[0] - self.knots[-p - k - 1]
        left = slice(new_start, n - p - k - 1, None)
        len_left = left.stop - left.start
        right = slice(0, n - len_left, None)
        (self.knots[:len_left], self.knots[len_left:]) = (self.knots[left], self.knots[right] - t1)

    def matches(self, bspline: BSplineBasis, reverse: bool = False) -> bool:
        """Checks if this basis equals another basis.

        This test disregards scaling and translation of the knots vector. I.e.
        will this basis and *bspline* yield the same spline object if paired
        with identical controlpoints.
        """
        if self.order != bspline.order or self.periodic != bspline.periodic:
            return False
        dt = self.knots[-1] - self.knots[0]
        dt2 = bspline.knots[-1] - bspline.knots[0]
        if reverse:
            return np.allclose(
                (self.knots[-1] - self.knots[::-1]) / dt,
                (bspline.knots - bspline.knots[0]) / dt2,
                atol=state.knot_tolerance,
            )
        return np.allclose(
            (self.knots - self.knots[0]) / dt,
            (bspline.knots - bspline.knots[0]) / dt2,
            atol=state.knot_tolerance,
        )

    def snap_points(self, t: FloatArray) -> None:
        """Snap evaluation points to knots if they are sufficiently close
        as given in by state.state.knot_tolerance. This will modify the input
        vector t

        :param t: evaluation points
        :type t: [float]
        :return: none
        """
        splipy_core.snap(self.knots, t, state.knot_tolerance)

    # TODO(Eivind): Deprecate this later.
    def snap(self, t):  # type: ignore[no-untyped-def]
        """Snap evaluation points to knots if they are sufficiently close
        as given in by state.state.knot_tolerance. This will modify the input
        vector t

        :param t: evaluation points
        :type t: [float]
        :return: none
        """
        n = len(self.knots)
        for j in range(len(t)):
            i = bisect_left(self.knots, t[j])
            if i < n and abs(self.knots[i] - t[j]) < state.knot_tolerance:
                t[j] = self.knots[i]
            elif i > 0 and abs(self.knots[i - 1] - t[j]) < state.knot_tolerance:
                t[j] = self.knots[i - 1]

    def clone(self) -> BSplineBasis:
        """Clone the object."""
        return copy.deepcopy(self)

    __call__ = evaluate

    def __len__(self) -> int:
        """Returns the number of knots in this basis."""
        return len(self.knots)

    def __getitem__(self, i: int) -> float:
        """Returns the knot at a given index."""
        return float(self.knots[i])

    def __iadd__(self, a: ScalarLike) -> Self:
        self.knots += float(a)
        return self

    def __isub__(self, a: ScalarLike) -> Self:
        self.knots -= float(a)
        return self

    def __imul__(self, a: ScalarLike) -> Self:
        self.knots *= float(a)
        return self

    def __itruediv__(self, a: ScalarLike) -> Self:
        self.knots /= float(a)
        return self

    __ifloordiv__ = __itruediv__  # integer division (should not distinguish)

    def __repr__(self) -> str:
        result = "p=" + str(self.order) + ", " + str(self.knots)
        if self.periodic > -1:
            result += ", C" + str(self.periodic) + "-periodic"
        return result
