from __future__ import annotations

__doc__ = "Implementation of various curve utilities"

from dataclasses import dataclass
from typing import TYPE_CHECKING, cast

import numpy as np
from scipy.optimize import minimize_scalar, root
from scipy.spatial.distance import cdist

from . import rot

if TYPE_CHECKING:
    from splipy.curve import Curve
    from splipy.typing import FloatArray


def curve_length_parametrization(pts, normalize=False):
    """Calculate knots corresponding to a curvelength parametrization of a set of
    points.

    :param numpy.array pts: A set of points
    :param bool normalize: Whether to normalize the parametrization
    :return: The parametrization
    :rtype: [float]
    """
    knots = [0.0]
    for i in range(1, len(pts)):
        knots.append(knots[-1] + np.linalg.norm(pts[i] - pts[i - 1]))

    if normalize:
        length = knots[-1]
        knots = [k / length for k in knots]

    return knots


def get_curve_points(curve):
    """Evaluate the curve in all its knots.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :return: The curve points
    :type: numpy.array
    """
    return curve(curve.knots(0))


def normal(curve: Curve, t: float) -> FloatArray:
    """Evaluate the right unit normal of a 2D curve at a parameter value.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :param float t: Parameter value
    :return: The right normal vector, or the zero vector if the tangent vanishes
    :rtype: numpy.array
    """
    tangent = curve.derivative(t)
    n = rot(tangent, -np.pi / 2)
    norm = np.linalg.norm(n)

    if norm > 1e-10:
        return n / norm

    return np.zeros(2)


@dataclass
class Feature:
    """Describes a kink or high-curvature point of a curve, together with the
    parametric interval that should be edited to remedy it when thickening.

    :ivar float t: Parameter value where the feature occurs on the curve
    :ivar str kind: Feature type, "kink" or "curvature"
    :ivar T1: Left endpoint of the edit interval
    :vartype T1: float or None
    :ivar T2: Right endpoint of the edit interval
    :vartype T2: float or None
    :ivar side: Side of the curve that needs editing, "left" or "right"
    :vartype side: str or None
    :ivar point: Intersection point between the two offset curves
    :vartype point: FloatArray or None
    """

    t: float
    kind: str
    T1: float | None = None
    T2: float | None = None
    side: str | None = None
    point: FloatArray | None = None


def merge_edit_intervals(features: list[Feature]) -> list[Feature]:
    """Remove features whose edit intervals overlap, keeping the intervals disjoint.

    :param list[Feature] features: List of features
    :return: List of features with unnecessary features removed
    :rtype: list[Feature]
    """
    if not features:
        return []

    features = sorted(features, key=lambda f: cast("float", f.T1))
    result = [features[0]]

    for feature in features[1:]:
        T1 = cast("float", feature.T1)
        T2 = cast("float", feature.T2)
        prev = result[-1]
        P1 = cast("float", prev.T1)
        P2 = cast("float", prev.T2)

        # No overlap.
        if T1 >= P2:
            result.append(feature)
            continue

        # Feature contained in previous.
        if P1 <= T1 and T2 <= P2:
            continue

        # Previous contained in feature.
        if T1 <= P1 and P2 <= T2:
            result[-1] = feature
            continue

        raise RuntimeError(
            f"Partially overlapping edits: {prev.kind}[{P1}, {P2}] and {feature.kind}[{T1}, {T2}]"
        )

    return result


def offset_points(
    curve: Curve,
    t: FloatArray | float,
    amount: float,
    change_1: tuple[FloatArray, str] | None = None,
    change_2: tuple[FloatArray, str] | None = None,
) -> tuple[FloatArray, FloatArray]:
    """Compute left and right offset points for a 2D curve at given parameter values.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :param t: Parameter value(s) to offset
    :type t: float or FloatArray
    :param float amount: Offset amount
    :param change_1: Override for the first offset point, as (point, side)
    :type change_1: tuple[FloatArray, str] or None
    :param change_2: Override for the last offset point, as (point, side)
    :type change_2: tuple[FloatArray, str] or None
    :return: The left and right offset points
    :rtype: tuple[numpy.array, numpy.array]
    """
    t = np.atleast_1d(np.asarray(t, dtype=np.float64))
    n = len(t)
    if n == 0:
        return np.empty((0, 2)), np.empty((0, 2))

    x = curve.evaluate(t)
    v = curve.derivative(t)
    speed = np.linalg.norm(v, axis=1)

    if np.all(speed < 1e-13):
        raise ValueError("Curve derivative is zero everywhere")

    for i in range(n):
        if speed[i] < 1e-13:
            j = i - 1 if i > 0 else i + 1
            v[i] = v[j]
            speed[i] = np.linalg.norm(v[i])

    v /= speed[:, None]
    normals = np.column_stack((-v[:, 1], v[:, 0]))

    left_points = x + amount * normals
    right_points = x - amount * normals

    if change_1 is not None:
        point, side = change_1

        if side == "left":
            left_points[0] = point
        elif side == "right":
            right_points[0] = point
        else:
            raise ValueError(f"Unknown side '{side}'")

    if change_2 is not None:
        point, side = change_2

        if side == "left":
            left_points[-1] = point
        elif side == "right":
            right_points[-1] = point
        else:
            raise ValueError(f"Unknown side '{side}'")

    return left_points, right_points


def first_intersection(
    curve: Curve,
    feature_t: float,
    offset: float,
    side: int,
    samples: int = 150,
    root_tol: float = 1e-6,
) -> tuple[float, float, FloatArray] | None:
    """Find the first intersection between the two offset branches of a curve on
    either side of a feature.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :param float feature_t: Parameter value of the feature
    :param float offset: Offset distance from the curve
    :param int side: +1 for the left offset, -1 for the right offset
    :param int samples: Number of samples used to represent the curve
    :param float root_tol: Tolerance for what counts as an intersection
    :return: (T1, T2, point), or None if no intersection is found
    :rtype: tuple[float, float, numpy.array] or None
    """
    t_start = curve.start()[0]
    t_end = curve.end()[0]

    # Binary search for a reasonable search interval around the feature.
    target_length = 8.0 * abs(offset)

    if curve.length(t_start, feature_t) > target_length:
        upper, lower = t_start, feature_t
        for _ in range(10):
            mid = (lower + upper) / 2
            if curve.length(mid, feature_t) > target_length:
                lower = mid
            else:
                upper = mid
        t_start = upper

    if curve.length(feature_t, t_end) > target_length:
        lower, upper = feature_t, t_end
        for _ in range(10):
            mid = (lower + upper) / 2
            if curve.length(feature_t, mid) > target_length:
                upper = mid
            else:
                lower = mid
        t_end = lower

    if feature_t <= t_start or feature_t >= t_end:
        return None

    eps_split = 1e-8
    left_params = np.linspace(t_start, feature_t - eps_split, samples)
    right_params = np.linspace(feature_t + eps_split, t_end, samples)

    left_at_left, right_at_left = offset_points(curve, left_params, offset)
    left_points = left_at_left if side == 1 else right_at_left

    left_at_right, right_at_right = offset_points(curve, right_params, offset)
    right_points = left_at_right if side == 1 else right_at_right

    distances = cdist(left_points, right_points)
    n_guesses = min(15, distances.size)

    flat_idx = np.argpartition(distances.ravel(), n_guesses - 1)[:n_guesses]
    left_idx, right_idx = np.unravel_index(flat_idx, distances.shape)

    def residual(z: FloatArray) -> FloatArray:
        T1, T2 = z

        T1 = min(max(T1, t_start), feature_t)
        T2 = min(max(T2, feature_t), t_end)

        left_arr, right_arr = offset_points(curve, np.array([T1, T2]), offset)
        chosen = left_arr if side == 1 else right_arr

        return cast("FloatArray", chosen[0] - chosen[1])

    solutions: list[list[float]] = []

    for i, j in zip(left_idx, right_idx):
        guess = [left_params[i], right_params[j]]
        sol = root(residual, guess)

        if not sol.success:
            continue

        T1, T2 = sol.x
        if not (t_start <= T1 <= feature_t and feature_t <= T2 <= t_end):
            continue

        error = np.linalg.norm(residual(np.array([T1, T2])))
        if error < root_tol:
            solutions.append([T1, T2])

    if len(solutions) == 0:
        return None

    solutions_arr = np.unique(np.round(np.asarray(solutions), decimals=8), axis=0)

    # Choose largest T1 and smallest T2 to get the intersection closest to the feature.
    idx = np.lexsort((solutions_arr[:, 1], -solutions_arr[:, 0]))[0]
    T1, T2 = solutions_arr[idx]

    point_arr, other_arr = offset_points(curve, np.array([T1]), offset)
    point = (point_arr if side == 1 else other_arr)[0]

    return T1, T2, point


def find_corner_features(curve: Curve, angle_tol: float = np.deg2rad(5.0)) -> list[Feature]:
    """Find kinks of a curve that are sharp enough to require editing.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :param float angle_tol: Angle tolerance for what counts as a kink
    :return: List of valid corner features
    :rtype: list[Feature]
    """
    corners = []

    t_start = curve.start()[0]
    t_end = curve.end()[0]

    # Epsilon used to avoid evaluating the derivative at the kink itself, and
    # to skip kinks at the endpoints.
    eps = 1e-8

    for corner_t in curve.get_kinks():
        if abs(corner_t - t_start) < eps:
            continue
        if abs(corner_t - t_end) < eps:
            continue

        v_left = np.asarray(curve.derivative(corner_t - eps))
        v_right = np.asarray(curve.derivative(corner_t + eps))

        left_norm = np.linalg.norm(v_left)
        right_norm = np.linalg.norm(v_right)

        if left_norm < 1e-12 or right_norm < 1e-12:
            continue

        v_left /= left_norm
        v_right /= right_norm

        angle = np.arccos(np.clip(np.dot(v_left, v_right), -1.0, 1.0))

        if angle < angle_tol:
            continue

        corners.append(Feature(t=corner_t, kind="kink"))

    return corners


def find_high_curvature_features(curve: Curve, offset: float = 1.0, samples: int = 5000) -> list[Feature]:
    """Find curvature maxima inside regions where the curvature exceeds the
    offset limit.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :param float offset: Offset distance from the curve that is considered
    :param int samples: Number of samples used to evaluate the curve curvature
    :return: List of high curvature features
    :rtype: list[Feature]
    """
    t = np.linspace(curve.start()[0], curve.end()[0], samples)
    kappa = cast("FloatArray", curve.curvature(t))
    limit = 1.0 / offset

    mask = kappa >= limit

    # Iterate through samples and register high curvature regions, finding
    # the curvature maximum inside each region.
    features = []
    i = 0

    while i < samples:
        if not mask[i]:
            i += 1
            continue

        start = i
        while i + 1 < samples and mask[i + 1]:
            i += 1
        end = i

        a = t[start]
        b = t[end]

        if a == b:
            feature_t = a
        else:
            result = minimize_scalar(lambda T: -float(curve.curvature(T)), bounds=(a, b), method="bounded")
            feature_t = result.x

        features.append(Feature(t=feature_t, kind="curvature"))

        i += 1

    return features


def find_all_features(curve: Curve, offset: float) -> list[Feature]:
    """Find all features of a curve that require editing given an offset amount.

    :param curve: The curve
    :type curve: :class:`splipy.Curve`
    :param float offset: Amount of curve offset
    :return: List of features that require editing
    :rtype: list[Feature]
    """
    features = []

    features.extend(find_high_curvature_features(curve, offset))
    features.extend(find_corner_features(curve))

    valid_features = []
    for feature in features:
        left_result = first_intersection(curve, feature.t, offset, side=+1)
        right_result = first_intersection(curve, feature.t, offset, side=-1)

        if left_result is not None:
            feature.T1, feature.T2, feature.point = left_result
            feature.side = "left"
        elif right_result is not None:
            feature.T1, feature.T2, feature.point = right_result
            feature.side = "right"
        else:
            continue

        valid_features.append(feature)

    return merge_edit_intervals(valid_features)
