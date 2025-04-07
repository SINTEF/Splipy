# -*- coding: utf-8 -*-

from math import pi

import numpy as np
from scipy.spatial import ConvexHull

from .basis import BSplineBasis
from .curve import Curve
from .surface import Surface
from .splineobject import SplineObject
from .utils import ensure_listlike, check_direction, sections
from . import state

__all__ = ['TrimmedSurface']


class TrimmedSurface(Surface):
    """TrimmedSurface()

    Represents a surface: an object with a two-dimensional parameter space and
    one or more interior closed trimming loops."""

    _intended_pardim = 2

    def __init__(self, basis1=None, basis2=None, controlpoints=None, rational=False, loops=None, **kwargs):
        """  Construct a surface with the given basis and control points.

        The default is to create a linear one-element mapping from and to the
        unit square.

        :param BSplineBasis basis1: The basis of the first parameter direction
        :param BSplineBasis basis2: The basis of the second parameter direction
        :param array-like controlpoints: An *n1* × *n2* × *d* matrix of control points
        :param bool rational: Whether the surface is rational (in which case the
            control points are interpreted as pre-multiplied with the weight,
        :param [[Curve]] Loops: One or more loops. Clockwise loops define holes,
            while counter-clockwise define surface interiors (materials always on the "left")
        :raises RuntimeError: If the loops are not contained to dimension 2 (parametric
            space), or if they are not closed, or if they are not looping properly
        """
        super(Surface, self).__init__([basis1, basis2], controlpoints, rational, **kwargs)
        # make sure to make deep copies of the loops so nothing bad happens
        self.boundaries = [[l.clone() for l in one_loop] for one_loop in loops]

        # error check input curves
        for one_loop in self.boundaries:
            for curve in one_loop:
                if not curve.dimension == 2:
                    raise RuntimeError('Boundary curves need to have dimension 2')
            for i in range(len(one_loop)):
                # print(state.parametric_absolute_tolerance)
                # print(one_loop[i-1][-1,:], ' ', one_loop[i][0,:])
                if not np.allclose(one_loop[i-1][-1,:], one_loop[i][0,:],
                                   rtol=state.parametric_relative_tolerance,
                                   atol=state.parametric_absolute_tolerance):
                    raise RuntimeError('Boundary curves not closed')

        self.__compute_convex_hulls()

    def edges(self):
        """Return the four edge curves in (parametric) order: umin, umax, vmin, vmax

        :return: Edge curves
        :rtype: (Curve)
        """
        return tuple(self.section(*args) for args in sections(2, 1))

    def __compute_convex_hulls(self):
        self.rotation   = []
        self.convexhull = []
        for loop in self.boundaries:
            # don't know if we really need the error test for boundary loops, but
            # I'll just leave it in, because I wrote i and kind of like it

            # get all controlpoints for all the curves, make sure not to double-
            # count the end and start of subsequent curve-pieces
            x = np.vstack([curve[1:,:2] for curve in loop])
            # compute vectors between the control points (velocity approximation?)
            dx = np.diff(np.append(x, [x[0,:]], axis=0), axis=0)
            # compute angle at all points
            theta = np.arctan2(dx[:,1], dx[:,0])
            # compute angle difference at all (control-)points
            dt = np.diff(np.append(theta, theta[0]))
            dt[np.where( dt<-pi) ] += 2*pi
            dt[np.where( dt>+pi) ] -= 2*pi

            total_rotation = np.sum(dt)
            if np.allclose(total_rotation, 2*pi):
                self.rotation.append('counterclockwise')
            elif np.allclose(total_rotation, -2*pi):
                self.rotation.append('clockwise')
            else:
                raise RuntimeError('Boundary loops does not loop exactly once')

            hull = ConvexHull(x)
            self.convexhull.append(x[hull.vertices,:])
            # print(self.convexhull[-1])

    def is_contained(self, u, v):
        """ Returns a boolean mask if the input points are inside (True) or 
        outside (False) of the trimming curves."""

        raise NotImplementedError('This has yet to be implemented')

        # do a quick test based on convex hull
        self.__is_contained_coarse(u,v)
        # for all points that are still undecided, do a fine-grained newton
        # iteration approach
        self.__is_contained_fine(u,v)

        return False

    def __is_contained_fine(self, u, v):
        """ Does a fine test based on parametric curve representation to see if
        points are inside or outside trimming domain. Trimming curves are high-
        polynomial representations, so figuring this out means newton iteration
        to locate nearest point on curve and decide if this is inside or outside
        domain."""
        return False

    def __is_contained_coarse(self, u, v):
        """ Does a course test based on control-grid to see if points are inside or
        outside domain. Inside control-grid means inside a trimming loop and outputs
        False. Outside the *convex hull* of a control-grid means"""
        return False

