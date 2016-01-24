from bisect import bisect_right, bisect_left
import numpy as np

__all__ = ['BSplineBasis']


class BSplineBasis:

    knots = [0, 0, 1, 1]
    order = 2
    periodic = -1
    tol = 1e-10

    def __init__(self, order=2, knots=None, periodic=-1):
        """Constructor
        @param order:    Spline order, i.e. one greater than the polynomial degere.
        @type  order:    Int
        @param knots:    Knot vector of non-decreasing components. Defaults to open knot vector on domain [0,1]
        @type  knots:    List of floats
        @param periodic: Number of continuous derivatives at start t0 and end tn. -1 is not periodic, 0 is N(t0)=N(tn), 1 is N'(t0)=N'(tn) etc.
        @type  periodic: Int
        """

        if knots is None:
            knots = [0] * order + [1] * order

        self.knots = np.array(knots)
        self.knots = self.knots.astype(float)
        self.order = order
        self.periodic = periodic

        # error test input
        if len(knots) < order:
            raise ValueError('knot vector has too few elements')
        if periodic >= 0:
            for i in range(periodic + 1):
                if abs((knots[i + 1] - knots[i]) - (knots[-i - 1] - knots[-i - 2])) > self.tol:
                    raise ValueError('periodic knot vector is mis-matching at the start/end')
        for i in range(len(knots) - 1):
            if knots[i + 1] - knots[i] < -self.tol:
                raise ValueError('knot vector needs to be non-decreasing')

    def start(self):
        """Start point of parametric domain. For open knot vectors, this is the first knot.
        @return: Knot number p, where p is the spline order
        @rtype : Float
        """
        return self.knots[self.order - 1]

    def end(self):
        """End point of parametric domain. For open knot vectors, this is the last knot.
        @return: Knot number n-p, where p is the spline order and n is the number of knots
        @rtype : Float
        """
        return self.knots[-self.order]

    def greville(self, index=None):
        """Fetch greville point(s), also known as knot averages, corresponding to this basis: sum_{j=index+1}^{index+p-1} t_j / (p-1)
        @return: One, or all of the Greville points
        @rtype : Float or List of Floats
        """
        result = []
        p = self.order
        n = len(self)
        if index is None:
            for i in range(n):
                result.append(float(np.sum(self.knots[i + 1:i + p])) / (p - 1))
        else:
            result = float(np.sum(self.knots[index + 1:index + p])) / (p - 1)
        return result

    def evaluate(self, t, d=0, from_right=True):
        """Get all basis functions evaluated in a given set of points
        @param t:          The parametric coordinate(s) to perform evaluation
        @type  t:          Float or List of Floats
        @param d:          Number of derivatives
        @type  d:          Int
        @param from_right: True if evaluation should be done in the limit from above
        @type  from_right: Boolean
        @return:           Matrix N[i,j] of all basis functions j in all points i
        @rtype:            numpy.array
        """

        # for single-value input, wrap it into a list so it don't crash on the loop below
        try:
            len(t)
        except TypeError:
            t = [t]

        p = self.order  # knot vector order
        n = len(self.knots) - p  # number of basis functions (without periodicity)
        N = np.matrix(np.zeros((len(t), n)))
        for i in range(len(t)):
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
                    from_right = False
                # Skip non-periodic evaluation points outside the domain
                if t[i] < self.start() or t[i] > self.end():
                    continue

            # mu = index of last non-zero basis function
            if from_right:
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
            N[:, j] += N[:, -j - 1]
        N = np.delete(N, range(n - self.periodic - 1, n), 1)

        return N

    def normalize(self):
        """Set the parametric domain to be (0,1)"""
        self -= self.start()  # set start-point to 0
        self /= self.end()  # set end-point to 1

    def reverse(self):
        """Reverse parametric domain, keeping start/end values unchanged"""
        a = float(self.start())
        b = float(self.end())
        self.knots = (self.knots[::-1] - a) / (b - a) * (a - b) + b

    def get_continuity(self, knot):
        """Get the continuity of the basis functions at a given point. Will
        return p-1-m, where m is the knot multiplicity and inf between knots"""
        p = self.order
        mu = bisect_left(self.knots, knot)
        if abs(self.knots[mu] - knot) > self.tol:
            return np.inf
        continuity = p - 1
        while mu < len(self.knots) and abs(self.knots[mu] - knot) < self.tol:
            continuity -= 1
            mu += 1
        return continuity

    def get_knot_spans(self):
        """Return the set of unique knots in the knot vector"""
        result = [self.knots[0]]
        for k in self.knots:
            if abs(k - result[-1]) > self.tol:
                result.append(k)
        return result

    def raise_order(self, amount):
        """Return a basis corresponding to a raise_order operation,
        keeping the continuity at the knots unchanged by increasing their
        multiplicity"""
        if type(amount) is not int:
            raise TypeError('amount needs to be a positive integer')
        if amount < 0:
            raise ValueError('amount needs to be a positive integer')
        knot_spans = list(self.get_knot_spans())  # list of unique knots
        # For every degree we raise, we need to increase the multiplicity by one
        knots = list(self.knots) + knot_spans * amount
        knots.sort()  # make it a proper knot vector by ensuring that it is non-decreasing

        return BSplineBasis(self.order + amount, knots, self.periodic)

    def insert_knot(self, new_knot):
        """ Returns an extremely sparse matrix C such that N_new = N_old*C,
        where N is a (row)vector of the basis functions. Also inserts new_knot
        in the knot vector.
        @param new_knot: The parametric coordinate point to insert
        @type  new_knot: Float
        @return:         Matrix C
        @rtype:          numpy.array
        """
        if new_knot < self.knots[0] or self.knots[-1] < new_knot:
            raise ValueError('new_knot out of range')
        # mu is the index of last non-zero (old) basis function
        mu = bisect_right(self.knots, new_knot)
        n = len(self)
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

        return C

    def write_g2(self, outfile):
        """write GoTools formatted BSplineBasis to file"""
        outfile.write('%i %i\n' % (len(self.knots) - self.order, self.order))
        for k in self.knots:
            outfile.write('%f ' % k)
        outfile.write('\n')

    def __call__(self, t):
        """see evaluate(t)"""
        return self.evaluate(t)

    def __len__(self):
        """returns the number of functions in this basis"""
        return len(self.knots) - self.order - (self.periodic + 1)

    def __getitem__(self, i):
        """returns knot i, including multiplicities"""
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
