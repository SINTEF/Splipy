from __future__ import annotations

from typing import Sequence, Optional, Union

import numpy as np

from ..types import Scalars, FArray

__doc__ = 'Implementation of various curve utilities'


def curve_length_parametrization(
    pts: Union[Sequence[Scalars], FArray],
    normalize: bool = False,
    buffer: Optional[FArray] = None,
) -> FArray:
    """Calculate knots corresponding to a curvelength parametrization of a set of
    points.

    :param numpy.array pts: A set of points
    :param bool normalize: Whether to normalize the parametrization
    :param numpy.array buffer: If given, the parametrization is stored in this array.
    :return: The parametrization
    :rtype: [float]
    """
    points = np.array(pts, dtype=float)

    knots = buffer if buffer is not None else np.empty((len(points) - 1,), dtype=float)
    knots[:] = np.cumsum(np.linalg.norm(points[1:] - points[:-1], axis=1))

    if normalize:
        knots /= knots[-1]

    return knots


def get_curve_points(curve):
    """Evaluate the curve in all its knots.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :return: The curve points
    :type: numpy.array
    """
    return curve(curve.knots(0))
