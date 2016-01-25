from GeoMod import Curve, BSplineBasis
from GeoMod.ControlPointOperations import ControlPointOperations
from GeoMod.Utils import ensure_listlike
from bisect import bisect_left
import numpy as np

__all__ = ['Surface']


class Surface(ControlPointOperations):
    def __init__(self, basis1=None, basis2=None, controlpoints=None, rational=False):
        super(Surface, self).__init__([basis1, basis2], controlpoints, rational)

    def evaluate_normal(self, u, v):
        """Evaluate the normal vector of the surface at given parametric value(s). The returned values are not normalized
        @param t : Parametric coordinate point(s)
        @type  t : Float or list of Floats
        @param u : Parametric coordinate point(s) in first direction
        @type  u : Float or list of Floats
        @param v : Parametric coordinate point(s) in second direction
        @type  v : Float or list of Floats
        @return  : 3D-array X(i,j,k) of the normal component x(k) evaluated at (u[i],v[j])
        @rtype   : numpy.ndarray
        """
        if self.dimension == 2:
            try:
                result = np.zeros((len(u), len(v), 3))
                result[:, :, 2] = 1
                return result
            except TypeError:  # single valued input u, fails on len(u)
                return np.array([0, 0, 1])
        elif self.dimension == 3:
            (du, dv) = self.evaluate_tangent(u, v)
            result = np.zeros(du.shape)
            # the cross product of the tangent is the normal
            if len(du.shape) == 1:
                result[0] = du[1] * dv[2] - du[2] * dv[1]
                result[1] = du[2] * dv[0] - du[0] * dv[2]
                result[2] = du[0] * dv[1] - du[1] * dv[0]
            else:  # grid evaluation
                result[:, :, 0] = du[:, :, 1] * dv[:, :, 2] - du[:, :, 2] * dv[:, :, 1]
                result[:, :, 1] = du[:, :, 2] * dv[:, :, 0] - du[:, :, 0] * dv[:, :, 2]
                result[:, :, 2] = du[:, :, 0] * dv[:, :, 1] - du[:, :, 1] * dv[:, :, 0]
            return result
        else:
            raise RuntimeError('Normal evaluation only defined for 2D and 3D geometries')

    def evaluate_tangent(self, u, v):
        """Evaluate the two tangent vectors of the surface at given parametric values
        @param t : Parametric coordinate point(s)
        @type  t : Float or list of Floats
        @param u : Parametric coordinate point(s) in first direction
        @type  u : Float or list of Floats
        @param v : Parametric coordinate point(s) in second direction
        @type  v : Float or list of Floats
        @return  : two 3D-arrays (dX/du, dX/dv) of the tangential component x(k) evaluated at (u[i],v[j])
        @rtype   : tuple of two numpy.ndarray
        """
        return (self.evaluate_derivative(u, v, d=(1, 0)),
                self.evaluate_derivative(u, v, d=(0, 1)))

    def swap_parametrization(self):
        """Swaps the two surface parameter directions"""
        self.controlpoints = self.controlpoints.transpose((1, 0, 2))
        self.bases = self.bases[::-1]

    def reparametrize(self, umin=0, umax=1, vmin=0, vmax=1):
        """Redefine the parametric domain to be (umin,umax) x (vmin,vmax)"""
        if umax <= umin or vmax <= vmin:
            raise ValueError('end must be larger than start')
        self.bases[0].normalize()  # set domain to (0,1)
        self.bases[0] *= (umax - umin)
        self.bases[0] += umin
        self.bases[1].normalize()
        self.bases[1] *= (vmax - vmin)
        self.bases[1] += vmin

    def get_edges(self):
        """Return the four edge curves in (parametric) order: bottom, right, top, left"""
        # ASSUMPTION: open knot vectors
        (p1, p2) = self.order()
        (n1, n2, dim) = self.controlpoints.shape
        rat = self.rational
        umin = Curve(self.bases[1], np.reshape(self.controlpoints[0, :, :], (n2, dim)), rat)
        umax = Curve(self.bases[1], np.reshape(self.controlpoints[-1, :, :], (n2, dim)), rat)
        vmin = Curve(self.bases[0], np.reshape(self.controlpoints[:, 0, :], (n1, dim)), rat)
        vmax = Curve(self.bases[0], np.reshape(self.controlpoints[:, -1, :], (n1, dim)), rat)
        # make the curves form a clockwise oriented closed loop around surface
        umax.flip_parametrization()
        vmax.flip_parametrization()
        return (vmin, umax, vmax, umin)

    def raise_order(self, raise_u, raise_v):
        """Raise the order of a spline surface
        @param raise_u: Number of polynomial degrees to increase in u
        @type  raise_u: Int
        @param raise_v: Number of polynomial degrees to increase in v
        @type  raise_v: Int
        """
        # create the new basis
        newBasis1 = self.bases[0].raise_order(raise_u)
        newBasis2 = self.bases[1].raise_order(raise_v)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_u = newBasis1.greville()  # parametric interpolation points u
        interpolation_pts_v = newBasis2.greville()  # parametric interpolation points v
        N_u_old = self.bases[0].evaluate(interpolation_pts_u)
        N_u_new = newBasis1.evaluate(interpolation_pts_u)
        N_v_old = self.bases[1].evaluate(interpolation_pts_v)
        N_v_new = newBasis2.evaluate(interpolation_pts_v)
        tmp = np.tensordot(N_v_old, self.controlpoints, axes=(1, 1))
        tmp = np.tensordot(N_u_old, tmp, axes=(1, 1))  # projective interpolation points (x,y,z,w)
        interpolation_pts_x = tmp  # 3D-tensor with elements (i,j,k) of component x[k] evaluated at u[i] v[j]

        # solve the interpolation problem
        N_u_inv = np.linalg.inv(N_u_new)
        N_v_inv = np.linalg.inv(
            N_v_new
        )  # these are inverses of the 1D problems, and small compared to the total number of unknowns
        tmp = np.tensordot(N_v_inv, interpolation_pts_x, axes=(1, 1))
        tmp = np.tensordot(N_u_inv, tmp, axes=(1, 1))

        # update the basis and controlpoints of the surface
        self.controlpoints = tmp
        self.bases = [newBasis1, newBasis2]

    def refine(self, n):
        """Enrich spline space by inserting n knots into each existing knot
        span
        @param n: The number of new knots to insert into each span
        @type  n: Int
        """
        (knots1, knots2) = self.knots()  # excluding multiple knots

        # insert new knots in the u-direction
        new_knots = []
        for (k0, k1) in zip(knots1[:-1], knots1[1:]):
            element_knots = np.linspace(k0, k1, n + 2)
            new_knots += list(element_knots[1:-1])
        self.insert_knot(0, new_knots)

        # insert new knots in the v-direction
        new_knots = []
        for (k0, k1) in zip(knots2[:-1], knots2[1:]):
            element_knots = np.linspace(k0, k1, n + 2)
            new_knots += list(element_knots[1:-1])
        self.insert_knot(1, new_knots)

    def insert_knot(self, direction, knot):
        """Insert a knot into this spline surface
        @param direction: The parametric direction (u=0, v=1)
        @type  direction: Int
        @param knot:      The knot(s) to insert
        @type  knot:      Float or list of Floats
        """
        # for single-value input, wrap it into a list
        knot = ensure_listlike(knot)
        if direction != 0 and direction != 1:
            raise ValueError('direction must be 0 or 1')

        (n1, n2, dim) = self.controlpoints.shape
        if direction == 0:
            C = np.matrix(np.identity(n1))
            for k in knot:
                C = self.bases[0].insert_knot(k) * C
            self.controlpoints = np.tensordot(C, self.controlpoints, axes=(1, 0))
        else:
            C = np.matrix(np.identity(n2))
            for k in knot:
                C = self.bases[1].insert_knot(k) * C
            self.controlpoints = np.tensordot(C,
                                              self.controlpoints,
                                              axes=(1, 1)).transpose((1, 0, 2))

    def split(self, direction, knots):
        """ Split a surface into two or more separate representations with C0
        continuity between them.
        @param direction: The parametric direction (u=0, v=1)
        @type  direction: Int
        @param knots    : splitting point(s)
        @type  knots    : Float or list of Floats
        @return         : The surface split into multiple pieces
        @rtype          : List of Surface
        """
        # for single-value input, wrap it into a list
        knots = ensure_listlike(knots)
        # error test input
        if direction != 0 and direction != 1:
            raise ValueError('direction must be 0 or 1')

        p = self.order()
        results = []
        splitting_surf = self.clone()
        basis = self.bases
        # insert knots to produce C{-1} at all splitting points
        for k in knots:
            continuity = basis[direction].get_continuity(k)
            if continuity == np.inf:
                continuity = p[direction] - 1
            splitting_surf.insert_knot(direction, [k] * (continuity + 1))

        # everything is available now, just have to find the right index range
        # in the knot vector and controlpoints to store in each separate curve
        # piece
        last_cp_i = 0
        last_knot_i = 0
        (n1, n2, dim) = splitting_surf.controlpoints.shape
        if direction == 0:
            for k in knots:
                mu = bisect_left(splitting_surf.bases[0].knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[0], splitting_surf.bases[0].knots[last_knot_i:mu + p[0]])
                cp = splitting_surf.controlpoints[last_cp_i:last_cp_i + n_cp, :, :]
                cp = np.reshape(cp.transpose((1, 0, 2)), (n_cp * n2, dim))

                results.append(Surface(basis, self.bases[1], cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp

            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[0], splitting_surf.bases[0].knots[last_knot_i:])
            n_cp = basis.num_functions()
            cp = splitting_surf.controlpoints[last_cp_i:, :, :]
            cp = np.reshape(cp.transpose((1, 0, 2)), (n_cp * n2, dim))
            results.append(Surface(basis, self.bases[1], cp, self.rational))
        else:
            for k in knots:
                mu = bisect_left(splitting_surf.bases[1].knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[1], splitting_surf.bases[1].knots[last_knot_i:mu + p[1]])
                cp = splitting_surf.controlpoints[:, last_cp_i:last_cp_i + n_cp, :]
                cp = np.reshape(cp.transpose((1, 0, 2)), (n_cp * n1, dim))

                results.append(Surface(self.bases[0], basis, cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp
            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[1], splitting_surf.bases[1].knots[last_knot_i:])
            n_cp = basis.num_functions()
            cp = splitting_surf.controlpoints[:, last_cp_i:, :]
            cp = np.reshape(cp.transpose((1, 0, 2)), (n_cp * n1, dim))
            results.append(Surface(self.bases[0], basis, cp, self.rational))

        return results

    def rebuild(self, p, n):
        """ Creates an approximation to this surface by resampling it using a
        uniform knot vector of order p and with n control points.
        @param p: Discretization order
        @type  p: Int or list of two int
        @param n: Number of control points
        @type  n: Int or list of two int
        @return : Approximation of this surface on a different basis
        @rtype  : Surface
        """
        p = ensure_listlike(p, dups=2)
        n = ensure_listlike(n, dups=2)

        old_basis = self.bases
        basis = []
        u = []
        N = []
        # establish uniform open knot vectors
        for i in range(2):
            knot = [0] * p[i] + range(1, n[i] - p[i] + 1) + [n[i] - p[i] + 1] * p[i]
            basis.append(BSplineBasis(p[i], knot))

            # make these span the same parametric domain as the old ones
            basis[i].normalize()
            t0 = old_basis[i].start()
            t1 = old_basis[i].end()
            basis[i] *= (t1 - t0)
            basis[i] += t0

            # fetch evaluation points and evaluate basis functions
            u.append(basis[i].greville())
            N.append(basis[i].evaluate(u[i]))

        # find interpolation points as evaluation of existing surface
        x = self.evaluate(u[0], u[1])

        # solve interpolation problem
        cp = np.tensordot(np.linalg.inv(N[1]), x, axes=(1, 1))
        cp = np.tensordot(np.linalg.inv(N[0]), cp, axes=(1, 1))

        # re-order controlpoints so they match up with Surface constructor
        cp = cp.transpose((1, 0, 2))
        cp = cp.reshape(n[0] * n[1], cp.shape[2])

        # return new resampled curve
        return Surface(basis[0], basis[1], cp)

    def write_g2(self, outfile):
        """write GoTools formatted SplineSurface to file"""
        outfile.write('200 1 0 0\n')  # surface header, gotools version 1.0.0
        outfile.write('%i %i\n' % (self.dimension, int(self.rational)))
        self.bases[0].write_g2(outfile)
        self.bases[1].write_g2(outfile)

        (n1, n2, n3) = self.controlpoints.shape
        for j in range(n2) + range(self.bases[1].periodic + 1):
            for i in range(n1) + range(self.bases[0].periodic + 1):
                for k in range(n3):
                    outfile.write('%f ' % self.controlpoints[i, j, k])
                outfile.write('\n')

    def __len__(self):
        """return the number of control points (basis functions) for this surface"""
        return self.bases[0].num_functions() * self.bases[1].num_functions()

    def __getitem__(self, i):
        (n1, n2, dim) = self.controlpoints.shape
        i1 = i % n1
        i2 = int(i / n1)
        return self.controlpoints[i1, i2, :]

    def __setitem__(self, i, newCP):
        (n1, n2, dim) = self.controlpoints.shape
        i1 = i % n1
        i2 = int(i / n1)
        self.controlpoints[i1, i2, :] = newCP
        return self

    def __repr__(self):
        result = str(self.bases[0]) + '\n' + str(self.bases[1]) + '\n'
        # print legacy controlpoint enumeration
        n1, n2, n3 = self.controlpoints.shape
        for j in range(n2):
            for i in range(n1):
                result += str(self.controlpoints[i, j, :]) + '\n'
        return result

    @classmethod
    def make_surfaces_compatible(cls, surf1, surf2):
        """ Make sure that surfaces are compatible, i.e. for merging. This
        will manipulate one or both to make sure that they are both rational
        and in the same geometric space (2D/3D).
        @param surf1: first surface
        @type  surf1: Surface
        @param surf2: second surface
        @type  surf2: Surface
        """
        # make both rational (if needed)
        if surf1.rational:
            surf2.force_rational()
        if surf2.rational:
            surf1.force_rational()

        # make both in the same geometric space
        if surf1.dimension > surf2.dimension:
            surf2.set_dimension(surf1.dimension)
        else:
            surf1.set_dimension(surf2.dimension)

    @classmethod
    def make_surfaces_identical(cls, surf1, surf2):
        """ Make sure that surfaces have identical discretization, i.e. same
        knot vector and order. May be used to draw a linear surface
        interpolation between them, or to add surfaces together.
        @param surf1: first surface
        @type  surf1: Surface
        @param surf2: second surface
        @type  surf2: Surface
        """
        # make sure that rational/dimension is the same
        Surface.make_surfaces_compatible(surf1, surf2)

        # make both have knot vectors in domain (0,1)
        surf1.reparametrize()
        surf2.reparametrize()

        # make sure both have the same order
        p1 = surf1.order()
        p2 = surf2.order()
        p = (max(p1[0], p2[0]), max(p1[1], p2[1]))
        surf1.raise_order(p[0] - p1[0], p[1] - p1[1])
        surf2.raise_order(p[0] - p2[0], p[1] - p2[1])

        # make sure both have the same knot vector in u-direction
        knot1 = surf1.knots(with_multiplicities=True)
        knot2 = surf2.knots(with_multiplicities=True)
        i1 = 0
        i2 = 0
        while i1 < len(knot1[0]) and i2 < len(knot2[0]):
            if abs(knot1[0][i1] - knot2[0][i2]) < surf1.bases[0].tol:
                i1 += 1
                i2 += 1
            elif knot1[0][i1] < knot2[0][i2]:
                surf2.insert_knot(0, knot1[0][i1])
                i1 += 1
            else:
                surf1.insert_knot(0, knot2[0][i2])
                i2 += 1

        # make sure both have the same knot vector in v-direction
        i1 = 0
        i2 = 0
        while i1 < len(knot1[1]) and i2 < len(knot2[1]):
            if abs(knot1[1][i1] - knot2[1][i2]) < surf1.bases[1].tol:
                i1 += 1
                i2 += 1
            elif knot1[1][i1] < knot2[1][i2]:
                surf2.insert_knot(1, knot1[1][i1])
                i1 += 1
            else:
                surf1.insert_knot(1, knot2[1][i2])
                i2 += 1
