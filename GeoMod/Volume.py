# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis, Surface
from GeoMod.SplineObject import SplineObject
from GeoMod.Utils import ensure_listlike
from bisect import bisect_left
import numpy as np

__all__ = ['Volume']


class Volume(SplineObject):
    """Volume()

    Represents a volume: an object with a three-dimensional parameter space."""

    def __init__(self, basis1=None, basis2=None, basis3=None, controlpoints=None, rational=False):
        """__init__([basis1=None], [basis2=None], [basis3=None], [controlpoints=None], [rational=False])

        Construct a volume with the given basis and control points.

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
        super(Volume, self).__init__([basis1, basis2, basis3], controlpoints, rational)

    def swap_parametrization(self, pardir1, pardir2):
        """Swaps two parameter directions.

        :param int pardir1: A direction
        :param int pardir2: Another direction
        :raises ValueError: If the parameter directions are not different and
            do not correspond to actual directions
        """
        if (pardir1 == 0 and pardir2 == 1) or (pardir1 == 1 and pardir2 == 0):
            self.controlpoints = self.controlpoints.transpose(
                (1, 0, 2, 3))  # re-order controlpoints
            self.bases[0], self.bases[1] = self.bases[1], self.bases[0]  # swap knot vectors
        elif (pardir1 == 0 and pardir2 == 2) or (pardir1 == 2 and pardir2 == 0):
            self.controlpoints = self.controlpoints.transpose(
                (2, 1, 0, 3))  # re-order controlpoints
            self.bases[0], self.bases[2] = self.bases[2], self.bases[0]  # swap knot vectors
        elif (pardir1 == 1 and pardir2 == 2) or (pardir1 == 2 and pardir2 == 1):
            self.controlpoints = self.controlpoints.transpose(
                (0, 2, 1, 3))  # re-order controlpoints
            self.bases[1], self.bases[2] = self.bases[2], self.bases[1]  # swap knot vectors
        else:
            raise ValueError(
                'pardir1 and pardir2 must be different from each other and either 0,1 or 2')

    def faces(self):
        """Return the six faces of this volume (with outward normal vectors) in
        order: umin, umax, vmin, vmax, wmin, wmax.

        :return: Boundary faces
        :rtype: (Surface)
        """
        (p1, p2, p3) = self.order()
        (n1, n2, n3, dim) = self.controlpoints.shape
        rat = self.rational
        umin = Surface(p3, p2, self.bases[2], self.bases[1], np.reshape(self.controlpoints[0, :, :, :],
                                                                    (n2 * n3, dim), rat))
        umax = Surface(p3, p2, self.bases[2], self.bases[1], np.reshape(self.controlpoints[-1, :, :, :],
                                                                    (n2 * n3, dim), rat))
        vmin = Surface(p3, p1, self.bases[2], self.bases[0], np.reshape(self.controlpoints[:, 0, :, :],
                                                                    (n1 * n3, dim), rat))
        vmax = Surface(p3, p1, self.bases[2], self.bases[0], np.reshape(self.controlpoints[:, -1, :, :],
                                                                    (n1 * n3, dim), rat))
        wmin = Surface(p2, p1, self.bases[1], self.bases[0], np.reshape(self.controlpoints[:, :, 0, :],
                                                                    (n1 * n2, dim), rat))
        wmax = Surface(p2, p1, self.bases[1], self.bases[0], np.reshape(self.controlpoints[:, :, -1, :],
                                                                    (n1 * n2, dim), rat))
        umax.swap_parametrization()
        vmax.swap_parametrization()
        wmax.swap_parametrization()
        return [umin, umax, vmin, vmax, wmin, wmax]

    def raise_order(self, raise_u, raise_v, raise_w):
        """Raise the order of the surface.

        :param int raise_u: Number of degrees to increase in the first direction
        :param int raise_v: Number of degrees to increase in the second direction
        :param int raise_w: Number of degrees to increase in the third direction
        """
        # create the new basis
        newBasis1 = self.bases[0].raise_order(raise_u)
        newBasis2 = self.bases[1].raise_order(raise_v)
        newBasis3 = self.bases[2].raise_order(raise_w)

        # set up an interpolation problem. This is in projective space, so no problems for rational cases
        interpolation_pts_u = newBasis1.greville()  # parametric interpolation points u
        interpolation_pts_v = newBasis2.greville()  # parametric interpolation points v
        interpolation_pts_w = newBasis3.greville()  # parametric interpolation points w
        N_u_old = self.bases[0].evaluate(interpolation_pts_u)
        N_u_new = newBasis1.evaluate(interpolation_pts_u)
        N_v_old = self.bases[1].evaluate(interpolation_pts_v)
        N_v_new = newBasis2.evaluate(interpolation_pts_v)
        N_w_old = self.bases[2].evaluate(interpolation_pts_w)
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
        self.bases = [newBasis1, newBasis2, newBasis3]

    def split(self, direction, knots):
        """Split a volume into two or more separate representations with C0
        continuity between them.

        :param int direction: The parametric direction to split in
        :param knots: The splitting points
        :type knots: float or [float]
        :return: The new volumes
        :rtype: [Volume]
        """
        # for single-value input, wrap it into a list
        knots = ensure_listlike(knots)
        # error test input
        if direction != 0 and direction != 1 and direction != 2:
            raise ValueError('direction must be 0, 1 or 2')

        p = self.order()
        results = []
        splitting_vol = self.clone()
        basis = [self.bases[0], self.bases[1], self.bases[2]]
        # insert knots to produce C{-1} at all splitting points
        for k in knots:
            continuity = basis[direction].continuity(k)
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
                mu = bisect_left(splitting_vol.bases[0].knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[0], splitting_vol.bases[0].knots[last_knot_i:mu + p[0]])
                cp = splitting_vol.controlpoints[last_cp_i:last_cp_i + n_cp, :, :, :]
                cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n2 * n3, dim))

                results.append(Volume(basis, self.bases[1], self.bases[2], cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp

            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[0], splitting_vol.bases[0].knots[last_knot_i:])
            n_cp = basis.num_functions()
            cp = splitting_vol.controlpoints[last_cp_i:, :, :, :]
            cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n2 * n3, dim))
            results.append(Volume(basis, self.bases[1], self.bases[2], cp, self.rational))
        elif direction == 1:
            for k in knots:
                mu = bisect_left(splitting_vol.bases[1].knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[1], splitting_vol.bases[1].knots[last_knot_i:mu + p[1]])
                cp = splitting_vol.controlpoints[:, last_cp_i:last_cp_i + n_cp, :, :]
                cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n3, dim))

                results.append(Volume(self.bases[0], basis, self.bases[2], cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp
            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[1], splitting_vol.bases[1].knots[last_knot_i:])
            n_cp = basis.num_functions()
            cp = splitting_vol.controlpoints[:, last_cp_i:, :, :]
            cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n3, dim))
            results.append(Volume(self.bases[0], basis, self.bases[2], cp, self.rational))
        else:
            for k in knots:
                mu = bisect_left(splitting_vol.bases[2].knots, k)
                n_cp = mu - last_knot_i
                basis = BSplineBasis(p[2], splitting_vol.bases[2].knots[last_knot_i:mu + p[2]])
                cp = splitting_vol.controlpoints[:, :, last_cp_i:last_cp_i + n_cp, :]
                cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n2, dim))

                results.append(Volume(self.bases[0], self.bases[1], basis, cp, self.rational))
                last_knot_i = mu
                last_cp_i += n_cp
            # with n splitting points, we're getting n+1 pieces. Add the final one:
            basis = BSplineBasis(p[2], splitting_vol.bases[2].knots[last_knot_i:])
            n_cp = basis.num_functions()
            cp = splitting_vol.controlpoints[:, :, last_cp_i:, :]
            cp = np.reshape(cp.transpose((2, 1, 0, 3)), (n_cp * n1 * n2, dim))
            results.append(Volume(self.bases[0], self.bases[1], basis, cp, self.rational))

        return results

    def rebuild(self, p, n):
        """Creates an approximation to this volume by resampling it using
        uniform knot vectors of order *p* with *n* control points.

        :param int p: Polynomial discretization order
        :param int n: Number of control points
        :return: A new approximate volume
        :rtype: Volume
        """
        p = ensure_listlike(p, dups=3)
        n = ensure_listlike(n, dups=3)

        old_basis = [self.bases[0], self.bases[1], self.bases[2]]
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
        """Write the volume in GoTools format.

        :param file-like outfile: The file to write to
        """
        outfile.write('700 1 0 0\n')  # volume header, gotools version 1.0.0
        outfile.write('%i %i\n' % (self.dimension, int(self.rational)))
        self.bases[0].write_g2(outfile)
        self.bases[1].write_g2(outfile)
        self.bases[2].write_g2(outfile)

        (n1, n2, n3, n4) = self.controlpoints.shape
        for k in range(n3) + range(self.bases[2].periodic + 1):
            for j in range(n2) + range(self.bases[1].periodic + 1):
                for i in range(n1) + range(self.bases[0].periodic + 1):
                    for d in range(n4):
                        outfile.write('%f ' % self.controlpoints[i, j, k, d])
                    outfile.write('\n')

    def __len__(self):
        """Return the number of control points (basis functions) for the surface."""
        return self.bases[0].num_functions() * self.bases[1].num_functions() * self.bases[2].num_functions()

    def __getitem__(self, i):
        """Get the control point at a given index.

        :rtype: numpy.array
        """
        (n1, n2, n3, dim) = self.controlpoints.shape
        i1 = int(i % n1)
        i2 = int(i / n1) % n2
        i3 = int(i / n1 / n2)
        return self.controlpoints[i1, i2, i3, :]

    def __setitem__(self, i, newCP):
        """Set the control point at a given index.

        :param int i: Index
        :param numpy.array newCP: New control point
        """
        (n1, n2, n3, dim) = self.controlpoints.shape
        i1 = int(i % n1)
        i2 = int(i / n1) % n2
        i3 = int(i / n1 / n2)
        self.controlpoints[i1, i2, i3, :] = newCP
        return self

    def __repr__(self):
        result = str(self.bases[0]) + '\n'
        result += str(self.bases[1]) + '\n'
        result += str(self.bases[2]) + '\n'
        # print legacy controlpoint enumeration
        n1, n2, n3, dim = self.controlpoints.shape
        for k in range(n3):
            for j in range(n2):
                for i in range(n1):
                    result += str(self.controlpoints[i, j, k, :]) + '\n'
        return result