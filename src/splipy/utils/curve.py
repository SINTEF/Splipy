from __future__ import annotations

__doc__ = "Implementation of various curve utilities"

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from splipy.curve import Curve
    from splipy.typing import FloatArray, Points


def curve_length_parametrization(pts: Points, normalize: bool = False, reps: int = 1) -> FloatArray:
    """Calculate knots corresponding to a curvelength parametrization of a set of
    points.

    :param numpy.array pts: A set of points
    :param bool normalize: Whether to normalize the parametrization
    :param reps int: How many repetitions of the first and last knot to return.
    :return: The parametrization
    :rtype: [float]
    """
    points = np.asarray(pts)
    npts = points.shape[0]
    knots = np.zeros((npts + (reps - 1) * 2,), dtype=float)

    distances = np.linalg.norm(points[1:, ...] - points[:-1, ...], axis=1)
    distances = np.cumsum(distances)
    knots[reps : reps - 1 + npts] = distances

    if reps > 1:
        knots[-reps + 1 :] = knots[-reps]

    if normalize:
        knots /= knots[-1]

    return knots


def get_curve_points(curve: Curve) -> FloatArray:
    """Evaluate the curve in all its knots.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :return: The curve points
    :type: numpy.array
    """
    return curve(curve.knots(0))
