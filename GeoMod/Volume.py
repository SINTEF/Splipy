from GeoMod import BSplineBasis, Surface
from GeoMod.ControlPointOperations import ControlPointOperations
from bisect import bisect_left
import numpy as np

__all__ = ['Volume']


class Volume(ControlPointOperations):
    def __init__(self, basis1=None, basis2=None, basis3=None, controlpoints=None, rational=False):

        if basis1 is None:
            basis1 = BSplineBasis()
        if basis2 is None:
            basis2 = BSplineBasis()
        if basis3 is None:
            basis3 = BSplineBasis()

        self.basis1 = basis1
        self.basis2 = basis2
        self.basis3 = basis3

        # if none provided, create the default geometry which is the linear mapping onto the unit cube (0,1)^3
        if controlpoints is None:
            controlpoints = []
            greville_points1 = self.basis1.greville()
            greville_points2 = self.basis2.greville()
            greville_points3 = self.basis3.greville()
            for p3 in greville_points3:
                for p2 in greville_points2:
                    for p1 in greville_points1:
                        controlpoints.append([p1, p2, p3])

        self.dimension = len(controlpoints[0]) - rational
        self.rational = rational

        # controlpoints are given in as 2-index (kji,l) for u[i], v[j], w[k], x[l]
        # reshape these into 4-index (k,j,i,l)
        self.controlpoints = np.reshape(controlpoints, (
            len(self.basis3), len(self.basis2), len(self.basis1), self.dimension + self.rational))
        # swap axis 0 and 2, to make it (i,j,k,l)
        self.controlpoints = self.controlpoints.transpose((2, 1, 0, 3))

    def evaluate(self, u, v, w):
        """Evaluate the volume at given parametric values
        @param u: Parametric coordinate point(s) in first direction
        @type  u: Float or list of Floats
        @param v: Parametric coordinate point(s) in second direction
        @type  v: Float or list of Floats
        @param w: Parametric coordinate point(s) in second direction
        @type  w: Float or list of Floats
        @return : Geometry coordinates. 4D-array X(i,j,k,l) of component x(l) evaluated at (u[i],v[j],w[k])
        @rtype  : numpy.array
        """
        # for single-value input, wrap it into a list
        try:
            len(u)
        except TypeError:
            u = [u]
        try:
            len(v)
        except TypeError:
            v = [v]
        try:
            len(w)
        except TypeError:
            w = [w]

        # error test input
        if self.basis1.periodic < 0:  # periodic functions can evaluate everywhere
            if min(u) < self.basis1.start() or self.basis1.end() < max(u):
                raise ValueError('evaluation outside parametric domain')
        if self.basis2.periodic < 0:
            if min(v) < self.basis2.start() or self.basis2.end() < max(v):
                raise ValueError('evaluation outside parametric domain')
        if self.basis3.periodic < 0:
            if min(w) < self.basis3.start() or self.basis3.end() < max(w):
                raise ValueError('evaluation outside parametric domain')

        # compute basis functions for all points t. Nu(i,j) is a matrix of all functions j for all points u[i]
        Nu = self.basis1.evaluate(u)
        Nv = self.basis2.evaluate(v)
        Nw = self.basis3.evaluate(w)

        # compute physical points [x,y,z] for all points (u[i],v[j],w[k]). For rational volumes, compute [X,Y,Z,W] (in projective space)
        result = np.tensordot(Nw, self.controlpoints, axes=(1, 2))
        result = np.tensordot(Nv, result, axes=(1, 2))
        result = np.tensordot(Nu, result, axes=(1, 2))

        # Project rational volumes down to geometry space: x = X/W, y=Y/W, z=Z/W
        if self.rational:
            for i in range(self.dimension):
                result[:, :, :, i] /= result[:, :, :, -1]
            result = np.delete(result, self.dimension, 3)  # remove all weight evaluations

        # in case of single value input (u,v,w), return vector instead of 4D-matrix
        if result.shape[0] == 1 and result.shape[1] == 1 and result.shape[2] == 1:
            result = np.array(result[0, 0, 0, :]).reshape(self.dimension)
        return result

    def evaluate_derivative(self, u, v, w, du=1, dv=1, dw=1):
        """Evaluate the derivative of the volume at given parametric values
        @param u : Parametric coordinate point(s) in first direction
        @type  u : Float or list of Floats
        @param v : Parametric coordinate point(s) in second direction
        @type  v : Float or list of Floats
        @param w : Parametric coordinate point(s) in third direction
        @type  w : Float or list of Floats
        @param du: Number of derivatives in u
        @type  du: Int
        @param dv: Number of derivatives in v
        @type  dv: Int
        @param dw: Number of derivatives in w
        @type  dw: Int
        @return  : 4D-array D^(du,dv,dw) X(i,j) of component x(l) differentiated (du,dv,dw) times at (u[i],v[j],w[k])
        @rtype   : numpy.ndarray
        """
        # for single-value input, wrap it into a list
        try:
            len(u)
        except TypeError:
            u = [u]
        try:
            len(v)
        except TypeError:
            v = [v]
        try:
            len(w)
        except TypeError:
            w = [w]

        # error test input
        if self.basis1.periodic < 0:  # periodic functions can evaluate everywhere
            if min(u) < self.basis1.start() or self.basis1.end() < max(u):
                raise ValueError('evaluation outside parametric domain')
        if self.basis2.periodic < 0:
            if min(v) < self.basis2.start() or self.basis2.end() < max(v):
                raise ValueError('evaluation outside parametric domain')
        if self.basis3.periodic < 0:
            if min(w) < self.basis3.start() or self.basis3.end() < max(w):
                raise ValueError('evaluation outside parametric domain')

        # compute basis functions for all points t. dNu(i,j) is a matrix of the derivative of all functions j for all points u[i]
        dNu = self.basis1.evaluate(u, du)
        dNv = self.basis2.evaluate(v, dv)
        dNw = self.basis3.evaluate(w, dw)

        # compute physical points [dx/dt,dy/dt,dz/dt] for all points (u[i],v[j],w[k])
        result = np.tensordot(dNw, self.controlpoints, axes=(1, 2))
        result = np.tensordot(dNv, result, axes=(1, 2))
        result = np.tensordot(dNu, result, axes=(1, 2))
        result = np.array(result)

        # Rational surfaces need the quotient rule to compute derivatives (controlpoints are stored as x_i*w_i)
        # x(u,v) = sum_ijk N_i(u) * N_j(v) * N_k(w) * (w_ijk*x_ijk) / W(u,v,w)
        # W(u,v) = sum_ijk N_i(u) * N_j(v) * N_k(w) * w_ijk
        # dx/du =  sum_ijk N_i'(u)*N_j(v)*N_k(w) * (w_ijk*x_ijk) / W(u,v,w) - sum_ijk N_i(u)N_j(v)_N_k(w)*w_ijk*x_ijk* W'(u,v,w)/W(u,v,w)^2
        if self.rational:
            if du + dv + dw > 1:
                raise RuntimeError(
                    'Rational derivatives not implemented for derivatives larger than 1')
            Nu = self.basis.evaluate(u)
            Nv = self.basis.evaluate(v)
            Nw = self.basis.evaluate(w)
            non_derivative = np.tensordot(Nw, self.controlpoints, axes=(1, 2))
            non_derivative = np.tensordot(Nv, non_derivative, axes=(1, 2))
            non_derivative = np.tensordot(Nu, non_derivative, axes=(1, 2))
            non_derivative = np.array(non_derivative)
            W = non_derivative[:, :, :, -1]  # W(u,v,w)
            Wder = result[:, :, :, -1]  # dW/du or dW/dv or dW/dw
            for i in range(self.dimension):
                result[:, :, :, i] = result[:, :, :, i] / W - non_derivative[:, :, :,
                                                                             i] * Wder / W / W

            result = np.delete(result, self.dimension, 1)  # remove the weight column

        return result

    def flip_parametrization(self, direction):
        """Swap direction of the volume by making it go in the reverse direction. Parametric domain remain unchanged
           @param direction: The parametric direction to flip (0=u, 1=v, 2=w)
           @type  direction: Int
        """
        if direction == 0:
            self.basis1.reverse()
            self.controlpoints = self.controlpoints[::-1, :, :, :]
        elif direction == 1:
            self.basis2.reverse()
            self.controlpoints = self.controlpoints[:, ::-1, :, :]
        elif direction == 2:
            self.basis3.reverse()
            self.controlpoints = self.controlpoints[:, :, ::-1, :]
        else:
            raise ValueError('direction must be 0,1 or 2')

    def get_order(self):
        """Return spline volume order (polynomial degree + 1) in all parametric directions"""
        return (self.basis1.order, self.basis2.order, self.basis3.order)

    def get_knots(self, with_multiplicities=False):
        """Get the knots of the spline volume
        @param with_multiplicities: Set to true to obtain the knot vector with multiplicities
        @type with_multiplicities : Boolean
        @return:                    List with the knot values
        @rtype :                    Tuple with List of float
        """
        if with_multiplicities:
            return (self.basis1.knots, self.basis2.knots, self.basis3.knots)
        else:
            return (self.basis1.get_knot_spans(), self.basis2.get_knot_spans(),
                    self.basis3.get_knot_spans())

    def start(self):
        """Return the start of the parametric domain"""
        return (self.basis1.start(), self.basis2.start(), self.basis3.start())

    def end(self):
        """Return the end of the parametric domain"""
        return (self.basis1.end(), self.basis2.end(), self.basis3.end())

    def swap_parametrization(self, pardir1, pardir2):
        """Swaps the two volume parameter directions"""
        if (pardir1 == 0 and pardir2 == 1) or (pardir1 == 1 and pardir2 == 0):
            self.controlpoints = self.controlpoints.transpose(
                (1, 0, 2, 3))  # re-order controlpoints
            self.basis1, self.basis2 = self.basis2, self.basis1  # swap knot vectors
        elif (pardir1 == 0 and pardir2 == 2) or (pardir1 == 2 and pardir2 == 0):
            self.controlpoints = self.controlpoints.transpose(
                (2, 1, 0, 3))  # re-order controlpoints
            self.basis1, self.basis3 = self.basis3, self.basis1  # swap knot vectors
        elif (pardir1 == 1 and pardir2 == 2) or (pardir1 == 2 and pardir2 == 1):
            self.controlpoints = self.controlpoints.transpose(
                (0, 2, 1, 3))  # re-order controlpoints
            self.basis2, self.basis3 = self.basis3, self.basis2  # swap knot vectors
        else:
            raise ValueError(
                'pardir1 and pardir2 must be different from each other and either 0,1 or 2')

    def reparametrize(self, umin=0, umax=1, vmin=0, vmax=1, wmin=0, wmax=1):
        """Redefine the parametric domain to be (umin,umax) x (vmin,vmax) x (wmin,wmax)"""
        if umax <= umin or vmax <= vmin or wmax <= wmin:
            raise ValueError('end must be larger than start')
        self.basis1.normalize()  # set domain to (0,1)
        self.basis1 *= (umax - umin)
        self.basis1 += umin
        self.basis2.normalize()
        self.basis2 *= (vmax - vmin)
        self.basis2 += vmin
        self.basis3.normalize()
        self.basis3 *= (wmax - wmin)
        self.basis3 += wmin

    def get_faces(self):
        """Return a list of the 6 boundary faces of this volume (with outward normal vector). They are ordered as (umin,umax,vmin,vmax,wmin,wmax)"""
        # ASSUMPTION: open knot vectors
        (p1, p2, p3) = self.get_order()
        (n1, n2, n3, dim) = self.controlpoints.shape
        rat = self.rational
        umin = Surface(p3, p2, self.basis3, self.basis2, np.reshape(self.controlpoints[0, :, :, :],
                                                                    (n2 * n3, dim), rat))
        umax = Surface(p3, p2, self.basis3, self.basis2, np.reshape(self.controlpoints[-1, :, :, :],
                                                                    (n2 * n3, dim), rat))
        vmin = Surface(p3, p1, self.basis3, self.basis1, np.reshape(self.controlpoints[:, 0, :, :],
                                                                    (n1 * n3, dim), rat))
        vmax = Surface(p3, p1, self.basis3, self.basis1, np.reshape(self.controlpoints[:, -1, :, :],
                                                                    (n1 * n3, dim), rat))
        wmin = Surface(p2, p1, self.basis2, self.basis1, np.reshape(self.controlpoints[:, :, 0, :],
                                                                    (n1 * n2, dim), rat))
        wmax = Surface(p2, p1, self.basis2, self.basis1, np.reshape(self.controlpoints[:, :, -1, :],
                                                                    (n1 * n2, dim), rat))
        umax.swap_parametrization()
        vmax.swap_parametrization()
        wmax.swap_parametrization()
        return [umin, umax, vmin, vmax, wmin, wmax]

    def raise_order(self, raise_u, raise_v, raise_w):
        """Raise the order of a spline volume
        @param raise_u: Number of polynomial degrees to increase in u
        @type  raise_u: Int
        @param raise_v: Number of polynomial degrees to increase in v
        @type  raise_v: Int
        @param raise_w: Number of polynomial degrees to increase in w
        @type  raise_w: Int
        """
        # create the new basis
        newBasis1 = self.basis1.raise_order(raise_u)
        newBasis2 = self.basis2.raise_order(raise_v)
        newBasis3 = self.basis3.raise_order(raise_w)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_u = newBasis1.greville()  # parametric interpolation points u
        interpolation_pts_v = newBasis2.greville()  # parametric interpolation points v
        interpolation_pts_w = newBasis3.greville()  # parametric interpolation points w
        N_u_old = self.basis1.evaluate(interpolation_pts_u)
        N_u_new = newBasis1.evaluate(interpolation_pts_u)
        N_v_old = self.basis2.evaluate(interpolation_pts_v)
        N_v_new = newBasis2.evaluate(interpolation_pts_v)
        N_w_old = self.basis3.evaluate(interpolation_pts_w)
        N_w_new = newBasis3.evaluate(interpolation_pts_w)
        tmp = np.tensordot(N_w_old, self.controlpoints, axes=(1, 2))
        tmp = np.tensordot(N_v_old, tmp, axes=(1, 2))
        tmp = np.tensordot(N_u_old, tmp, axes=(1, 2))  # projective interpolation points (x,y,z,w)
        interpolation_pts_x = tmp

        # solve the interpolation problem
        N_u_inv = np.linalg.inv(N_u_new)
        N_v_inv = np.linalg.inv(N_v_new)
        N_w_inv = np.linalg.inv(
            N_w_new
        )  # these are inverses of the 1D problems, and small compared to the total number of unknowns
        tmp = np.tensordot(N_w_inv, interpolation_pts_x, axes=(1, 2))
        tmp = np.tensordot(N_v_inv, tmp, axes=(1, 2))
        tmp = np.tensordot(N_u_inv, tmp, axes=(1, 2))

        # update the basis and controlpoints of the volume
        self.controlpoints = tmp
        self.basis1 = newBasis1
        self.basis2 = newBasis2
        self.basis3 = newBasis3

    def refine(self, n):
        """Enrich spline space by inserting n knots into each existing knot
        span
        @param n: The number of new knots to insert into each span
        @type  n: Int
        """
        (knots1, knots2, knots3) = self.get_knots()  # excluding multiple knots

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

        # insert new knots in the w-direction
        new_knots = []
        for (k0, k1) in zip(knots3[:-1], knots3[1:]):
            element_knots = np.linspace(k0, k1, n + 2)
            new_knots += list(element_knots[1:-1])
        self.insert_knot(2, new_knots)

    def insert_knot(self, direction, knot):
        """Insert a knot into this spline volume
        @param direction: The parametric direction (u=0, v=1, w=2)
        @type  direction: Int
        @param knot:      The knot(s) to insert
        @type  knot:      Float or list of Floats
        """
        # for single-value input, wrap it into a list
        try:
            len(knot)
        except TypeError:
            knot = [knot]
        if direction != 0 and direction != 1 and direction != 2:
            raise ValueError('direction must be 0, 1 or 2')

        (n1, n2, n3, dim) = self.controlpoints.shape
        if direction == 0:
            C = np.matrix(np.identity(n1))
            for k in knot:
                C = self.basis1.insert_knot(k) * C
            self.controlpoints = np.tensordot(C, self.controlpoints, axes=(1, 0))
        elif direction == 1:
            C = np.matrix(np.identity(n2))
            for k in knot:
                C = self.basis2.insert_knot(k) * C
            self.controlpoints = np.tensordot(C,
                                              self.controlpoints,
                                              axes=(1, 1)).transpose((1, 0, 2, 3))
        elif direction == 2:
            C = np.matrix(np.identity(n3))
            for k in knot:
                C = self.basis3.insert_knot(k) * C
            self.controlpoints = np.tensordot(C,
                                              self.controlpoints,
                                              axes=(1, 2)).transpose((1, 2, 0, 3))

    def split(self, direction, knots):
        """ Split a volume into two or more separate representations with C0
        continuity between them.
        @param direction: The parametric direction (u=0, v=1, w=2)
        @type  direction: Int
        @param knots    : splitting point(s)
        @type  knots    : Float or list of Floats
        @return         : The volume split into multiple pieces
        @rtype          : List of Volume
        """
        # for single-value input, wrap it into a list
        try:
            len(knots)
        except TypeError:
            knots = [knots]
        # error test input
        if direction != 0 and direction != 1 and direction != 2:
            raise ValueError('direction must be 0, 1 or 2')

        p = self.get_order()
        results = []
        splitting_vol = self.clone()
        basis = [self.basis1, self.basis2, self.basis3]
        # insert knots to produce C{-1} at all splitting points
        for k in knots:
            continuity = basis[direction].get_continuity(k)
            if continuity == np.inf:
                continuity = p[direction] - 1
            splitting_vol.insert_knot(direction, [k] * (continuity + 1))

        # everything is available now, just have to find the right index range
        # in the knot vector and controlpoints to store in each separate curve
        # piece
        last_cp_i = 0
        last_knot_i = 0
        (n1, n2, n3, dim) = splitting_vol.controlpoints.shape
        if direction == 0:
            for k in knots:
                mu = bisect_left(splitting_vol.basis1.knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[0], splitting_vol.basis1.knots[last_knot_i:mu + p[0]])
                cp = splitting_vol.controlpoints[last_cp_i:last_cp_i + n_cp, :, :, :]
                cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n2 * n3, dim))

                results.append(Volume(basis, self.basis2, self.basis3, cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp

            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[0], splitting_vol.basis1.knots[last_knot_i:])
            n_cp = len(basis)
            cp = splitting_vol.controlpoints[last_cp_i:, :, :, :]
            cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n2 * n3, dim))
            results.append(Volume(basis, self.basis2, self.basis3, cp, self.rational))
        elif direction == 1:
            for k in knots:
                mu = bisect_left(splitting_vol.basis2.knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[1], splitting_vol.basis2.knots[last_knot_i:mu + p[1]])
                cp = splitting_vol.controlpoints[:, last_cp_i:last_cp_i + n_cp, :, :]
                cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n3, dim))

                results.append(Volume(self.basis1, basis, self.basis3, cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp
            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[1], splitting_vol.basis2.knots[last_knot_i:])
            n_cp = len(basis)
            cp = splitting_vol.controlpoints[:, last_cp_i:, :, :]
            cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n3, dim))
            results.append(Volume(self.basis1, basis, self.basis3, cp, self.rational))
        else:
            for k in knots:
                mu = bisect_left(splitting_vol.basis3.knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[2], splitting_vol.basis3.knots[last_knot_i:mu + p[2]])
                cp = splitting_vol.controlpoints[:, :, last_cp_i:last_cp_i + n_cp, :]
                cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n2, dim))

                results.append(Volume(self.basis1, self.basis2, basis, cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp
            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[2], splitting_vol.basis3.knots[last_knot_i:])
            n_cp = len(basis)
            cp = splitting_vol.controlpoints[:, :, last_cp_i:, :]
            cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n2, dim))
            results.append(Volume(self.basis1, self.basis2, basis, cp, self.rational))

        return results

    def rebuild(self, p, n):
        """ Creates an approximation to this volume by resampling it using a
        uniform knot vector of order p and with n control points.
        @param p: Discretization order
        @type  p: Int or list of three int
        @param n: Number of control points
        @type  n: Int or list of three int
        @return : Approximation of this volume on a different basis
        @rtype  : Volume
        """
        try:
            len(p)
        except TypeError:
            p = [p, p, p]
        try:
            len(n)
        except TypeError:
            n = [n, n, n]

        old_basis = [self.basis1, self.basis2, self.basis3]
        basis = []
        u = []
        N = []
        # establish uniform open knot vectors
        for i in range(3):
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

        # find interpolation points as evaluation of existing volume
        x = self.evaluate(u[0], u[1], u[2])

        # solve interpolation problem
        cp = np.tensordot(np.linalg.inv(N[2]), x, axes=(1, 2))
        cp = np.tensordot(np.linalg.inv(N[1]), cp, axes=(1, 2))
        cp = np.tensordot(np.linalg.inv(N[0]), cp, axes=(1, 2))

        # re-order controlpoints so they match up with Volume constructor
        cp = cp.transpose((2, 1, 0, 3))
        cp = cp.reshape(n[0] * n[1] * n[2], cp.shape[3])

        # return new resampled curve
        return Volume(basis[0], basis[1], basis[2], cp)

    def write_g2(self, outfile):
        """write GoTools formatted SplineVolume to file"""
        outfile.write('700 1 0 0\n')  # volume header, gotools version 1.0.0
        outfile.write('%i %i\n' % (self.dimension, int(self.rational)))
        self.basis1.write_g2(outfile)
        self.basis2.write_g2(outfile)
        self.basis3.write_g2(outfile)

        (n1, n2, n3, n4) = self.controlpoints.shape
        for k in range(n3) + range(self.basis3.periodic + 1):
            for j in range(n2) + range(self.basis2.periodic + 1):
                for i in range(n1) + range(self.basis1.periodic + 1):
                    for d in range(n4):
                        outfile.write('%f ' % self.controlpoints[i, j, k, d])
                    outfile.write('\n')

    __call__ = evaluate

    def __len__(self):
        """return the number of control points (basis functions) for this volume"""
        return len(self.basis1) * len(self.basis2) * len(self.basis3)

    def __getitem__(self, i):
        (n1, n2, n3, dim) = self.controlpoints.shape
        i1 = int(i % n1)
        i2 = int(i / n1) % n2
        i3 = int(i / n1 / n2)
        return self.controlpoints[i1, i2, i3, :]

    def __setitem__(self, i, newCP):
        (n1, n2, n3, dim) = self.controlpoints.shape
        i1 = int(i % n1)
        i2 = int(i / n1) % n2
        i3 = int(i / n1 / n2)
        self.controlpoints[i1, i2, i3, :] = newCP
        return self

    def __repr__(self):
        result = str(self.basis1) + '\n'
        result += str(self.basis2) + '\n'
        result += str(self.basis3) + '\n'
        # print legacy controlpoint enumeration
        n1, n2, n3, dim = self.controlpoints.shape
        for k in range(n3):
            for j in range(n2):
                for i in range(n1):
                    result += str(self.controlpoints[i, j, k, :]) + '\n'
        return result
