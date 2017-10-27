# -*- coding: utf-8 -*-

from splipy.utils import ensure_listlike
import splipy.state as state
from bisect import bisect_right, bisect_left
import numpy as np
import copy
from scipy.sparse import csr_matrix

__all__ = ['BSplineBasis']


class BSplineBasis:
    """BSplineBasis()

    Represents a one-dimensional B-Spline basis.

    BSplineBasis objects support basic arithmetic operators, which are
    interpreted as acting on the parametric domain.
    """
    knots = [0, 0, 1, 1]
    order = 2
    periodic = -1

    def __init__(self, order=2, knots=None, periodic=-1):
        """  Construct a B-Spline basis with a given order and knot vector.

        :param int order: Spline order, i.e. one greater than the polynomial degree.
        :param [float] knots: Knot vector of non-decreasing components.
            Defaults to open knot vector on domain [0,1].
        :param int periodic: Number of continuous derivatives at start and end.
            --1 is not periodic, 0 is continuous, etc.
        :raises ValueError: for inapproriate knot vectors
        """

        periodic = max(periodic, -1)
        if knots is None:
            knots = [0] * order + [1] * order
            for i in range(periodic+1):
                knots[   i] = -1
                knots[-i-1] =  2

        self.knots = np.array(knots)
        self.knots = self.knots.astype(float)
        self.order = order
        self.periodic = periodic

        # error test input
        p          = order
        k          = periodic
        n          = len(knots)
        if p < 1:
            raise ValueError('invalid spline order')
        if n < 2*p:
            raise ValueError('knot vector has too few elements')
        if periodic >= 0:
            for i in range(p + k - 1):
                if abs((knots[i + 1] - knots[i]) - (knots[-p - k + i ] - knots[-p - k - 1 + i])) > state.knot_tolerance:
                    raise ValueError('periodic knot vector is mis-matching at the start/end')
        for i in range(len(knots) - 1):
            if knots[i + 1] - knots[i] < -state.knot_tolerance:
                raise ValueError('knot vector needs to be non-decreasing')

    def num_functions(self):
        """  Returns the number of basis functions in the basis.

        .. warning:: This is different from :func:`splipy.BSplineBasis.__len__`."""
        return len(self.knots) - self.order - (self.periodic + 1)

    def start(self):
        """Start point of parametric domain. For open knot vectors, this is the
        first knot.

        :return: Knot number *p*, where *p* is the spline order
        :rtype: float
        """
        return self.knots[self.order - 1]

    def end(self):
        """End point of parametric domain. For open knot vectors, this is the
        last knot.

        :return: Knot number *n*--*p*, where *p* is the spline order and *n* is
            the number of knots
        :rtype: Float
        """
        return self.knots[-self.order]

    def greville(self, index=None):
        """  Fetch greville points, also known as knot averages:

        .. math:: \sum_{j=i+1}^{i+p-1} \\frac{t_j}{p-1}

        :return: One, or all of the Greville points
        :rtype: [float] (if *index* is ``None``) or float
        """
        result = []
        p = self.order
        n = self.num_functions()
        if index is None:
            for i in range(n):
                result.append(float(np.sum(self.knots[i + 1:i + p])) / (p - 1))
        else:
            result = float(np.sum(self.knots[index + 1:index + p])) / (p - 1)
        return result

    def evaluate(self, t, d=0, from_right=True, sparse=False):
        """  Evaluate all basis functions in a given set of points.

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
        t = ensure_listlike(t)
        self.snap(t)

        p = self.order  # knot vector order
        n_all = len(self.knots) - p  # number of basis functions (without periodicity)
        n = len(self.knots) - p - (self.periodic+1)  # number of basis functions (with periodicity)
        m = len(t)
        data    = np.zeros(m*p)
        indices = np.zeros(m*p, dtype='int32')
        indptr  = np.array(range(0,m*p+1,p), dtype='int32')
        if p <= d: # requesting more derivatives than polymoial degree: return all zeros
            return np.matrix(np.zeros((m,n)))
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
            if right:
                mu = bisect_right(self.knots, evalT)
            else:
                mu = bisect_left(self.knots, evalT)
            mu = min(mu, n_all)

            M = np.zeros(p)  # temp storage to keep all the function evaluations
            M[-1] = 1  # the last entry is a dummy-zero which is never used
            for q in range(1, p-d):
                for j in range(p - q - 1, p):
                    k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                    if j != p-q-1:
                        M[j] = M[j] * float(evalT - self.knots[k]) / (self.knots[k + q] - self.knots[k])

                    if j != p-1:
                        M[j] = M[j] + M[j + 1] * float(self.knots[k + q + 1] - evalT) / (self.knots[k + q + 1] - self.knots[k + 1])


            for q in range(p-d, p):
                for j in range(p - q - 1, p):
                    k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                    if j != p-q-1:
                        M[j] = M[j] * float(q) / (self.knots[k + q] - self.knots[k])
                    if j != p-1:
                        M[j] = M[j] - M[j + 1] * float(q) / (self.knots[k + q + 1] - self.knots[k + 1])


            data[i*p:(i+1)*p]    = M
            indices[i*p:(i+1)*p] = np.arange(mu-p, mu) % n

        N = csr_matrix((data, indices, indptr), (m,n))
        if not sparse:
            N = N.todense()
        return N

    def integrate(self, t0, t1):
        """  Integrate all basis functions over a given domain

        :param float t0: The parametric starting point
        :param float t1: The parametric end point
        :return: The integration of all functions over the input domain
        :rtype: list
        """
        if self.periodic > -1 and (t0<self.start() or t1>self.end()):
            raise NotImplemented('Periodic functions integrated across sem')

        t0 = max(t0, self.start())
        t1 = min(t1, self.end()  )
        p  = self.order
        knot = [self.knots[0]] + list(self.knots) + [self.knots[-1]]
        integration_basis = BSplineBasis(p + 1, knot)
        N0 = np.array(integration_basis.evaluate(t0)).flatten()
        N1 = np.array(integration_basis.evaluate(t1)).flatten()
        N  = [(knot[i+p]-knot[i])*1.0/p * np.sum(N1[i:]-N0[i:]) for i in range(N0.size)]
        N  = N[1:]

        # collapse periodic functions onto themselves
        if self.periodic > -1:
            for j in range(self.periodic + 1):
                N[j] += N[-self.periodic - 1 + j]
            N = N[:-self.periodic-1]

        return N

    def normalize(self):
        """Set the parametric domain to be (0,1)."""
        self -= self.start()  # set start-point to 0
        self /= self.end()  # set end-point to 1

    def reparam(self, start=0, end=1):
        """  Set the parametric domain to be (start, end)

        :raises ValueError: If *end* ≤ *start*"""
        if end <= start:
            raise ValueError('end must be larger than start')
        self.normalize()
        self *= (end - start)
        self += start

    def reverse(self):
        """Reverse parametric domain, keeping start/end values unchanged."""
        a = float(self.start())
        b = float(self.end())
        self.knots = (self.knots[::-1] - a) / (b - a) * (a - b) + b

    def continuity(self, knot):
        """Get the continuity of the basis functions at a given point.

        :return: *p*--*m*--1 at a knot with multiplicity *m*, or ``inf``
            between knots.
        :rtype: int or float
        """
        if self.periodic >= 0:
            if knot < self.start() or knot > self.end():
                knot = (knot - self.start()) % (self.end() - self.start()) + self.start()
        elif knot < self.start() or self.end() < knot:
            raise ValueError('out of range')

        p = self.order
        mu = bisect_left(self.knots, knot)

        # Pick the knot to the left if it exists and is closer
        if mu > 0 and abs(self.knots[mu-1] - knot) < abs(self.knots[mu] - knot):
            mu -= 1

        if abs(self.knots[mu] - knot) > state.knot_tolerance:
            return np.inf
        continuity = p - 1
        while mu < len(self.knots) and abs(self.knots[mu] - knot) < state.knot_tolerance:
            continuity -= 1
            mu += 1
        return continuity

    def make_periodic(self, continuity):
        """Create a periodic basis with a given continuity."""
        deg = self.order - 1
        new_knots = self.knots[deg:-deg]

        diff = self.end() - self.start()
        n_reps = deg - continuity - 1
        n_copy = deg - n_reps

        head = new_knots[-n_copy-1:-1] - diff
        tail = new_knots[1:n_copy+1] + diff

        new_knots = np.hstack((head, [self.start()] * n_reps, new_knots, [self.end()] * n_reps, tail))
        return BSplineBasis(self.order, new_knots, continuity)

    def knot_spans(self, include_ghost_knots=False):
        """Return the set of unique knots in the knot vector.

        :param bool include_ghost_knots: if knots outside start/end are to be
            included. These knots are used by periodic basis.
        :return: List of unique knots
        :rtype: [float]"""
        p = self.order
        if include_ghost_knots:
            result = [self.knots[0]]
            for k in self.knots:
                if abs(k - result[-1]) > state.knot_tolerance:
                    result.append(k)
        else:
            result = [self.knots[p-1]]
            for k in self.knots[p-1:-p+1]:
                if abs(k - result[-1]) > state.knot_tolerance:
                    result.append(k)
        return result

    def raise_order(self, amount):
        """Create a knot vector with higher order.

        The continuity at the knots are kept unchanged by increasing their
        multiplicities.

        :return: New knot vector
        :rtype: [float]
        :raises TypeError: If `amount` is not an int
        :raises ValueError: If `amount` is negative
        """
        if type(amount) is not int:
            raise TypeError('amount needs to be a non-negative integer')
        if amount < 0:
            raise ValueError('amount needs to be a non-negative integer')
        if amount == 0:
            return self.clone()
        knot_spans = list(self.knot_spans(True))  # list of unique knots
        # For every degree we raise, we need to increase the multiplicity by one
        knots = list(self.knots) + knot_spans * amount
        # make it a proper knot vector by ensuring that it is non-decreasing
        knots.sort()
        if self.periodic > -1:
            # remove excessive ghost knots which appear at both ends of the knot vector
            n0 =                   bisect_left(knot_spans, self.start())
            n1 = len(knot_spans) - bisect_left(knot_spans, self.end())   - 1
            knots = knots[n0*amount : -n1*amount]

        return BSplineBasis(self.order + amount, knots, self.periodic)

    def lower_order(self, amount):
        """Create a knot vector with lower order.

        The continuity at the knots are kept unchanged by decreasing their
        multiplicities.

        :return: New knot vector
        :rtype: [float]
        :raises TypeError: If `amount` is not an int
        :raises ValueError: If `amount` is negative
        """
        if type(amount) is not int:
            raise TypeError('amount needs to be a non-negative integer')
        if amount < 0:
            raise ValueError('amount needs to be a non-negative integer')
        if self.order - amount < 2:
            raise ValueError('cannot lower order to less than linears')

        p = self.order - amount
        knots = [ [k] * max(p-1-self.continuity(k), 1) for k in self.knot_spans(True)]
        knots = [ k for sublist in knots for k in sublist]

        if self.periodic > -1:
            # remove excessive ghost knots which appear at both ends of the knot vector
            n0 =                   bisect_left(knot_spans, self.start())
            n1 = len(knot_spans) - bisect_left(knot_spans, self.end())   - 1
            knots = knots[n0*amount : -n1*amount]

        return BSplineBasis(p, knots, self.periodic)

    def insert_knot(self, new_knot):
        """Inserts a knot in the knot vector.

        The return value is a sparse matrix *C* (actually, a dense matrix with
        lots of zeros), such that *N_new* = *N_old* x *C*, where *N* are row
        vectors of basis functions.

        :param float new_knot: The parametric coordinate of the point to insert
        :return: Transformation matrix *C*
        :rtype: numpy.array
        :raises ValueError: If the new knot is outside the domain
        """
        if self.periodic >= 0:
            if new_knot < self.start() or new_knot > self.end():
                new_knot = (new_knot - self.start()) % (self.end() - self.start()) + self.start()
        elif new_knot < self.start() or self.end() < new_knot:
            raise ValueError('new_knot out of range')
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
                C[i % (n + 1), i % n] = (new_knot - self.knots[i]) / (
                    self.knots[i + p - 1] - self.knots[i])
            if self.knots[i] <= new_knot and new_knot <= self.knots[i + 1]:
                C[(i + 1) % (n + 1), i % n] = 1
            else:
                C[(i + 1) % (n + 1), i % n] = (self.knots[i + p] - new_knot) / (
                    self.knots[i + p] - self.knots[i + 1])
        for i in range(mu, n + 1):
            C[i % (n + 1), (i - 1) % n] = 1

        self.knots = np.insert(self.knots, mu, new_knot)

        # make sure that it is correct periodic after knot insertion
        if self.periodic > -1:
            m  = len(self.knots)
            r  = self.periodic
            if mu <= p+r: # need to fix ghost knots on right side
                k0 = self.knots[0]
                k1 = self.knots[-p-r-1]
                for i in range(p+r+1):
                    self.knots[m-p-r-1+i] = k1 + (self.knots[i]-k0)
            elif mu >= m-p-r-1: # need to fix ghost knots on left side
                k0 = self.knots[p+r]
                k1 = self.knots[-1]
                for i in range(p+r+1):
                    self.knots[i] = k0 - (k1-self.knots[m-p-r-1+i])
        return C

    def roll(self, new_start):
        """rotate a periodic knot vector by setting a new starting index.

        :param int new_start: The index of to the new first knot
        """
        if self.periodic < 0:
            raise  RuntimeError("roll only applicable for periodic knot vectors")

        p = self.order
        k = self.periodic
        n = len(self.knots)
        t1 = self.knots[0] - self.knots[-p - k - 1]
        left  = slice(new_start, n-p-k-1, None)
        len_left = left.stop - left.start
        right = slice(0, n-len_left, None)
        (self.knots[:len_left], self.knots[len_left:]) = (self.knots[left], self.knots[right] - t1)

    def matches(self, bspline, reverse=False):
        """ Checks if this basis equals another basis, when disregarding
        scaling and translation of the knots vector. I.e. will this basis and
        *bspline* yield the same spline object if paired with identical
        controlpoints """
        if self.order != bspline.order or self.periodic != bspline.periodic:
            return False
        dt  = self.knots[-1]    - self.knots[0]
        dt2 = bspline.knots[-1] - bspline.knots[0]
        if reverse:
            return np.allclose( (self.knots[-1]-self.knots[::-1]) / dt,
                                (bspline.knots-bspline.knots[0]) / dt2,
                                atol=state.knot_tolerance)
        else:
            return np.allclose( (self.knots-self.knots[0]) / dt,
                                (bspline.knots-bspline.knots[0]) / dt2,
                                atol=state.knot_tolerance)

    def snap(self, t):
        """  Snap evaluation points to knots if they are sufficiently close
        as given in by state.state.knot_tolerance. This will modify the input vector t

        :param t: evaluation points
        :type  t: [float]
        :return: none
        """
        n = len(self.knots)
        for j in range(len(t)):
            i = bisect_left(self.knots, t[j])
            if i < n and abs(self.knots[i]-t[j]) < state.knot_tolerance:
                t[j] = self.knots[i]
            elif i > 0 and abs(self.knots[i-1]-t[j]) < state.knot_tolerance:
                t[j] = self.knots[i-1]


    def clone(self):
        """Clone the object."""
        return copy.deepcopy(self)

    __call__ = evaluate

    def __len__(self):
        """Returns the number of knots in this basis."""
        return len(self.knots)

    def __getitem__(self, i):
        """Returns the knot at a given index."""
        return self.knots[i]

    def __iadd__(self, a):
        self.knots += a
        return self

    def __isub__(self, a):
        self.knots -= a
        return self

    def __imul__(self, a):
        self.knots *= a
        return self

    def __itruediv__(self, a):
        self.knots /= a
        return self

    __ifloordiv__ = __itruediv__  # integer division (should not distinguish)
    __idiv__ = __itruediv__  # python2 compatibility

    def __repr__(self):
        result = 'p=' + str(self.order) + ', ' + str(self.knots)
        if self.periodic > -1:
            result += ', C' + str(self.periodic) + '-periodic'
        return result
