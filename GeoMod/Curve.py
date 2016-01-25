from GeoMod import BSplineBasis
from GeoMod.ControlPointOperations import ControlPointOperations
from GeoMod.Utils import ensure_listlike
from bisect import bisect_left
import numpy as np

__all__ = ['Curve']


class Curve(ControlPointOperations):
    def __init__(self, basis=None, controlpoints=None, rational=False):
        super(Curve, self).__init__([basis], controlpoints, rational)
        self.basis = self.bases[0]

    def evaluate_tangent(self, t):
        """Evaluate the tangent of the curve at given parametric values
        @param t: Parametric coordinate point(s)
        @type  t: Float or list of Floats
        @return : Tangent evaluation. Matrix DX(i,j) of component x(j) evaluated at t(i)
        @rtype  : numpy.array
        """
        return self.evaluate_derivative(t, d=1)

    def evaluate_derivative(self, t, d=1):
        """Evaluate the derivative of the curve at given parametric values
        @param t: Parametric coordinate point(s)
        @type  t: Float or list of Floats
        @param d: Number of derivatives
        @type  d: Int
        @return : Matrix D^n X(i,j) of component x(j) differentiated d times at point t(i)
        @rtype  : numpy.array
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

    def flip_parametrization(self):
        """Swap direction of the curve by making it go in the reverse direction. Parametric domain remain unchanged"""
        self.bases[0].reverse()
        self.controlpoints = self.controlpoints[::-1]

    def get_knots(self, with_multiplicities=False):
        """Get the knots of a spline curve
        @param with_multiplicities: Set to true to obtain the knot vector with multiplicities
        @type with_multiplicities : Boolean
        @return:                    List with the knot values
        @rtype :                    List of float
        """
        if with_multiplicities:
            return self.bases[0].knots
        else:
            return self.bases[0].get_knot_spans()

    def reparametrize(self, start=0, end=1):
        """Redefine the parametric domain to be (start,end)"""
        if end <= start:
            raise ValueError('end must be larger than start')
        self.bases[0].normalize()  # set domain to (0,1)
        self.bases[0] *= (end - start)
        self.bases[0] += start

    def raise_order(self, amount):
        """Raise the order of a spline curve
        @param amount: Number of polynomial degrees to increase
        @type  amount: Int
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

    def refine(self, n):
        """Enrich spline space by inserting n knots into each existing knot
        span
        @param n: The number of new knots to insert into each span
        @type  n: Int
        """
        new_knots = []
        knots = self.get_knots()  # excluding multiple knots
        for (k0, k1) in zip(knots[:-1], knots[1:]):
            element_knots = np.linspace(k0, k1, n + 2)
            new_knots += list(element_knots[1:-1])
        self.insert_knot(new_knots)

    def insert_knot(self, knot):
        """Insert a knot into this spline curve
        @param knot: The knot(s) to insert
        @type  knot: Float or list of Floats
        """
        # for single-value input, wrap it into a list
        knot = ensure_listlike(knot)

        C = np.matrix(np.identity(len(self)))
        for k in knot:
            C = self.bases[0].insert_knot(k) * C

        self.controlpoints = C * self.controlpoints

    def append(self, curve):
        """ Extend this curve by merging another curve to the end of it. The
        curves are glued together in a C0 fashion with enough repeated knots.
        The function assumes that the end of this curve perfectly matches the
        start of the input curve.
        @param curve: Another curve
        @type  curve: Curve
        """
        # ASSUMPTION: open knot vectors

        # error test input
        if self.bases[0].periodic > -1 or curve.bases[0].periodic > -1:
            raise RuntimeError('Cannot append with periodic curves')

        # copy input curve so we don't change that one directly
        extending_curve = curve.clone()

        # make sure both are in the same space, and (if needed) have rational weights
        Curve.make_curves_compatible(self, extending_curve)

        # make sure both have the same discretization order
        p1 = self.order(0)
        p2 = extending_curve.order(0)
        if p1 < p2:
            self.raise_order(p2 - p1)
        else:
            extending_curve.raise_order(p1 - p2)
        p = max(p1, p2)

        # build new knot vector by merging the two existing ones
        old_knot = self.get_knots(True)
        add_knot = extending_curve.get_knots(True)
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

    def get_continuity(self, knot):
        """Get the parametric continuity of the curve at a given point. Will
        return p-1-m, where m is the knot multiplicity and inf between knots"""
        return self.bases[0].get_continuity(knot)

    def split(self, knots):
        """ Split a curve into two or more separate representations with C0
        continuity between them.
        @param knots: splitting point(s)
        @type  knots: Float or list of Floats
        @return     : The curve split into multiple pieces
        @rtype      : List of Curves
        """
        # for single-value input, wrap it into a list
        knots = ensure_listlike(knots)

        p = self.order(0)
        results = []
        splitting_curve = self.clone()
        # insert knots to produce C{-1} at all splitting points
        for k in knots:
            continuity = splitting_curve.get_continuity(k)
            if continuity == np.inf:
                continuity = p - 1
            splitting_curve.insert_knot([k] * (continuity + 1))

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
        """ Creates an approximation to this curve by resampling it using a
        uniform knot vector of order p and with n control points.
        @param p: Discretization order
        @type  p: Int
        @param n: Number of control points
        @type  n: Int
        @return : Approximation of this curve on a different basis
        @rtype  : Curve
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
        """write GoTools formatted SplineCurve to file"""
        outfile.write('100 1 0 0\n')  # surface header, gotools version 1.0.0
        outfile.write('%i %i\n' % (self.dimension, int(self.rational)))
        self.bases[0].write_g2(outfile)

        (n1, n2) = self.controlpoints.shape
        for i in range(n1) + range(self.bases[0].periodic + 1):
            for j in range(n2):
                outfile.write('%f ' % self.controlpoints[i, j])
            outfile.write('\n')

    def __len__(self):
        """return the number of control points (basis functions) for this curve"""
        return self.bases[0].num_functions()

    def __getitem__(self, i):
        return self.controlpoints[i, :]

    def __setitem__(self, i, newCP):
        self.controlpoints[i, :] = newCP
        return self

    def __repr__(self):
        return str(self.bases[0]) + '\n' + str(self.controlpoints)

    @classmethod
    def make_curves_compatible(cls, crv1, crv2):
        """ Make sure that curves are compatible, i.e. for merging. This
        will manipulate one or both to make sure that they are both rational
        and in the same geometric space (2D/3D).
        @param crv1: first curve
        @type  crv1: Curve
        @param crv2: second curve
        @type  crv2: Curve
        """
        # make both rational
        if crv1.rational:
            crv2.force_rational()
        if crv2.rational:
            crv1.force_rational()

        # make both in the same geometric space
        if crv1.dimension > crv2.dimension:
            crv2.set_dimension(crv1.dimension)
        else:
            crv1.set_dimension(crv2.dimension)

    @classmethod
    def make_curves_identical(cls, crv1, crv2):
        """ Make sure that curves have identical discretization, i.e. same
        knot vector and order. May be used to draw a linear surface
        interpolation between them, or to add curves together.
        @param crv1: first curve
        @type  crv1: Curve
        @param crv2: second curve
        @type  crv2: Curve
        """
        # make sure that rational/dimension is the same
        Curve.make_curves_compatible(crv1, crv2)

        # make both have knot vectors in domain (0,1)
        crv1.reparametrize()
        crv2.reparametrize()

        # make sure both have the same order
        p1 = crv1.order(0)
        p2 = crv2.order(0)
        if p1 < p2:
            crv1.raise_order(p2 - p1)
        else:
            crv2.raise_order(p1 - p2)

        # make sure both have the same knot vector
        knot1 = crv1.get_knots(True)
        knot2 = crv2.get_knots(True)
        i1 = 0
        i2 = 0
        while i1 < len(knot1) and i2 < len(knot2):
            if abs(knot1[i1] - knot2[i2]) < crv1.bases[0].tol:
                i1 += 1
                i2 += 1
            elif knot1[i1] < knot2[i2]:
                crv2.insert_knot(knot1[i1])
                i1 += 1
            else:
                crv1.insert_knot(knot2[i2])
                i2 += 1
