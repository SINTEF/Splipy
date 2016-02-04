# -*- coding: utf-8 -*-

from GeoMod.Utils import ensure_listlike
from bisect import bisect_right, bisect_left
import numpy as np

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
    tol = 1e-10

    def __init__(self, order=2, knots=None, periodic=-1):
        """__init__([order=2], [knots=None], [periodic=-1])

        Construct a B-Spline basis with a given order and knot vector.

        :param int order: Spline order, i.e. one greater than the polynomial degree.
        :param [float] knots: Knot vector of non-decreasing components.
            Defaults to open knot vector on domain [0,1].
        :param int periodic: Number of continuous derivatives at start and end.
            --1 is not periodic, 0 is continuous, etc.
        :raises ValueError: for inapproriate knot vectors
        """

        if knots is None:
            knots = [0] * order + [1] * order
            for i in range(periodic+1):
                knots[   i] = -1
                knots[-i-1] =  2

        self.knots = np.array(knots)
        self.knots = self.knots.astype(float)
        self.order = order
        p          = order
        self.periodic = periodic

        # error test input
        if len(knots) < order:
            raise ValueError('knot vector has too few elements')
        if len(knots) < 2*(order+periodic):
            raise ValueError('knot vector has too few elements')
        if periodic >= 0:
            for i in range(periodic + 1):
                if abs((knots[i + 1] - knots[i]) - (knots[-p - periodic + i ] - knots[-p - periodic - 1 + i])) > self.tol:
                    raise ValueError('periodic knot vector is mis-matching at the start/end')
        for i in range(len(knots) - 1):
            if knots[i + 1] - knots[i] < -self.tol:
                raise ValueError('knot vector needs to be non-decreasing')

    def num_functions(self):
        """Returns the number of basis functions in the basis.

        .. warning:: This is different from :func:`GeoMod.BSplineBasis.__len__`."""
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
        """greville([index=None])

        Fetch greville points, also known as knot averages:

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

    def evaluate(self, t, d=0, from_right=True):
        """evaluate(t, [d=0], [from_right=True])

        Evaluate all basis functions in a given set of points.

        :param t: The parametric coordinate(s) in which to evaluate
        :type t: float or [float]
        :param int d: Number of derivatives to compute
        :param bool from_right: True if evaluation should be done in the limit
            from above
        :return: A matrix *N[i,j]* of all basis functions *j* evaluated in all
            points *i*
        :rtype: numpy.array
        """

        # for single-value input, wrap it into a list so it don't crash on the loop below
        t = ensure_listlike(t)

        p = self.order  # knot vector order
        n = len(self.knots) - p  # number of basis functions (without periodicity)
        N = np.matrix(np.zeros((len(t), n)))
        for i in range(len(t)):
            right = from_right
            if p <= d:
                continue  # requesting more derivatives than polymoial degree: return all zeros
            evalT = t[i]
            # Wrap periodic evaluation into domain
            if self.periodic >= 0:
                if t[i] < self.start() or t[i] > self.end():
                    evalT = (t[i] - self.start()) % (self.end() - self.start()) + self.start()
            else:
                # Special-case the endpoint, so the user doesn't need to
                if abs(t[i] - self.end()) < self.tol:
                    right = False
                # Skip non-periodic evaluation points outside the domain
                if t[i] < self.start() or t[i] > self.end():
                    continue

            # mu = index of last non-zero basis function
            if right:
                mu = bisect_right(self.knots, evalT)
            else:
                mu = bisect_left(self.knots, evalT)
            mu = min(mu, n)

            M = np.zeros(p + 1)  # temp storage to keep all the function evaluations
            M[-2] = 1  # the last entry is a dummy-zero which is never used
            for q in range(1, p):
                for j in range(p - q - 1, p):
                    k = mu - p + j  # 'i'-index in global knot vector (ref Hughes book pg.21)
                    if abs(self.knots[k + q] - self.knots[k]) > self.tol:
                        if q < p - d:  # in case of normal evaluation
                            M[j] = M[j] * float(evalT - self.knots[k]) / (
                                self.knots[k + q] - self.knots[k])
                        else:  # in case of derivative evaluation
                            M[j] = M[j] * float(q) / (self.knots[k + q] - self.knots[k])
                    else:  # in case of multiplicative knot
                        M[j] = 0
                    if abs(self.knots[k + q + 1] - self.knots[
                            k + 1]) > self.tol:  # and the same for the second term in the sum
                        if q < p - d:
                            M[j] = M[j] + M[j + 1] * float(self.knots[k + q + 1] - evalT) / (
                                self.knots[k + q + 1] - self.knots[k + 1])
                        else:
                            M[j] = M[j] - M[j + 1] * float(q) / (self.knots[k + q + 1] -
                                                                 self.knots[k + 1])
            N[i, (mu - p):mu] = M[0:-1]

        # collapse periodic functions onto themselves
        for j in range(self.periodic + 1):
            N[:, j] += N[:, -self.periodic - 1 + j]
        N = np.delete(N, range(n - self.periodic - 1, n), 1)

        return N

    def integrate(self, t0, t1):
        """integrate(t0, t1)

        Integrate all basis functions over a given domain

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
        """reparam([start=0], [end=1])

        Set the parametric domain to be (start, end)

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
        p = self.order
        mu = bisect_left(self.knots, knot)
        if abs(self.knots[mu] - knot) > self.tol:
            return np.inf
        continuity = p - 1
        while mu < len(self.knots) and abs(self.knots[mu] - knot) < self.tol:
            continuity -= 1
            mu += 1
        return continuity

    def knot_spans(self):
        """Return the set of unique knots in the knot vector.

        :return: List of unique knots
        :rtype: [float]"""
        result = [self.knots[0]]
        for k in self.knots:
            if abs(k - result[-1]) > self.tol:
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
        knot_spans = list(self.knot_spans())  # list of unique knots
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
        if new_knot < self.start() or self.end() < new_knot:
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

    def write_g2(self, outfile):
        """Write the basis in GoTools format.

        :param file-like outfile: The file to write to.
        """
        outfile.write('%i %i\n' % (len(self.knots) - self.order, self.order))
        for k in self.knots:
            outfile.write('%f ' % k)
        outfile.write('\n')

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
