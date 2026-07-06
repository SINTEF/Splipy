from __future__ import annotations

from collections.abc import Sequence
from typing import Literal, TypedDict

import numpy as np
import numpy.typing as npt

# Anything that can be converted to an array. Use this type for parameters that
# are (multidimensional) lists of points, such as controlpoints. Generally, only
# use this for public-facing functions. Use np.asarray to convert to an actual
# array, and then use FloatArray internally in Splipy. Avoid using np.array,
# since that will make a copy of objects that are already arrays. If you need a
# copy to mutate, be explicit and use np.asarray().copy().
type ArrayLike = npt.ArrayLike

# Alias to signal that we expect control points - makes no difference to type
# checking but it's helpful for humans.
type ControlPoints = ArrayLike

# Alias to signal that we expect points (which are not precisely the same thing
# as control points) - makes no difference to type checking but it's helpful for
# humans.
type Points = ArrayLike

# Anything that acts and behaves as a float. Note that np.floating and
# np.integer encompasses floats of many different sizes.
type Scalar = float | np.floating | int | np.integer

type Int = int | np.integer

type FloatArray = npt.NDArray[np.floating]
type IntArray = npt.NDArray[np.integer]

# Again, these are identical, but we can signal intent to a human reader based
# on which we use.
type Knots = Sequence[Scalar] | FloatArray | IntArray
type Params = Sequence[Scalar] | FloatArray | IntArray
type Point = Sequence[Scalar] | FloatArray | IntArray

type Direction = Literal["u", "v", "w", "U", "V", "W"] | int
type SectionElement = Literal[-1, 0] | None
type Section = tuple[SectionElement, ...]


class SectionKwargs(TypedDict, total=False):
    u: SectionElement
    v: SectionElement
    w: SectionElement
