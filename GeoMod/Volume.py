# -*- coding: utf-8 -*-

from GeoMod import BSplineBasis, Surface
from GeoMod.SplineObject import SplineObject
from GeoMod.Utils import ensure_listlike, check_direction
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

    def faces(self):
        """Return the six faces of this volume in order: umin, umax, vmin, vmax, wmin, wmax.

        :return: Boundary faces
        :rtype: (Surface)
        """
        (p1, p2, p3) = self.order()
        (n1, n2, n3, dim) = self.controlpoints.shape
        rat = self.rational
        umin = Surface(self.bases[2], self.bases[1], np.reshape(self.controlpoints[0, :, :, :],  (n2 * n3, dim), rat))
        umax = Surface(self.bases[2], self.bases[1], np.reshape(self.controlpoints[-1, :, :, :], (n2 * n3, dim), rat))
        vmin = Surface(self.bases[2], self.bases[0], np.reshape(self.controlpoints[:, 0, :, :],  (n1 * n3, dim), rat))
        vmax = Surface(self.bases[2], self.bases[0], np.reshape(self.controlpoints[:, -1, :, :], (n1 * n3, dim), rat))
        wmin = Surface(self.bases[1], self.bases[0], np.reshape(self.controlpoints[:, :, 0, :],  (n1 * n2, dim), rat))
        wmax = Surface(self.bases[1], self.bases[0], np.reshape(self.controlpoints[:, :, -1, :], (n1 * n2, dim), rat))
        result = [umin, umax, vmin, vmax, wmin, wmax]
        for s in result:
            s.swap()
        return result

    def corners(self):
        """Return the eight corner control points in (parametric) in order (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1),...

        :return: Corners
        :rtype: (np.ndarray)
        .. warning:: For rational splines, this will return the corners in projective coordinates, including weights.
        """
        (n1, n2, n3, dim) = self.controlpoints.shape
        return self.controlpoints[::n1-1, ::n2-1, ::n3-1, :].reshape((8,dim))

    def split(self, knots, direction):
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
        direction = check_direction(direction, self.pardim)

        p = self.order()
        results = []
        splitting_vol = self.clone()
        basis = [self.bases[0], self.bases[1], self.bases[2]]
        # insert knots to produce C{-1} at all splitting points
        for k in knots:
            continuity = basis[direction].continuity(k)
            if continuity == np.inf:
                continuity = p[direction] - 1
            splitting_vol.insert_knot([k] * (continuity + 1), direction)

        b = splitting_vol.bases[direction]
        if b.periodic > -1:
            mu = bisect_left(b.knots, knots[0])
            b.roll(mu)
            splitting_vol.controlpoints = np.roll(splitting_vol.controlpoints, -mu, direction)
            b.knots = b.knots[:-b.periodic-1]
            b.periodic = -1
            if len(knots) > 1:
                return splitting_vol.split(knots[1:], direction)
            else:
                return splitting_vol

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
        vol = self
        for i in range(self.pardim):
            if self.periodic(i):
                vol = vol.split(vol.start(i), i)
        outfile.write('700 1 0 0\n')  # volume header, gotools version 1.0.0
        outfile.write('%i %i\n' % (vol.dimension, int(vol.rational)))
        vol.bases[0].write_g2(outfile)
        vol.bases[1].write_g2(outfile)
        vol.bases[2].write_g2(outfile)

        (n1, n2, n3, n4) = vol.controlpoints.shape
        for k in chain(range(n3), range(vol.bases[2].periodic + 1)):
            for j in chain(range(n2), range(vol.bases[1].periodic + 1)):
                for i in chain(range(n1), range(vol.bases[0].periodic + 1)):
                    for d in range(n4):
                        outfile.write('%f ' % vol.controlpoints[i, j, k, d])
                    outfile.write('\n')

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
