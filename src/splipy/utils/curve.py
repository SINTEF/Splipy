__doc__ = 'Implementation of various curve utilities'

import numpy as np


def curve_length_parametrization(pts, normalize=False):
    """Calculate knots corresponding to a curvelength parametrization of a set of
    points.

    :param numpy.array pts: A set of points
    :param bool normalize: Whether to normalize the parametrization
    :return: The parametrization
    :rtype: [float]
    """
    knots = [0.0]
    for i in range(1, pts.shape[0]):
        knots.append(knots[-1] + np.linalg.norm(pts[i,:] - pts[i-1,:]))

    if normalize:
        length = knots[-1]
        knots = [k/length for k in knots]

    return knots


def get_curve_points(curve):
    """Evaluate the curve in all its knots.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :return: The curve points
    :type: numpy.array
    """
    return curve(curve.knots(0))
