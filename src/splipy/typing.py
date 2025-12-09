from __future__ import annotations

from typing import Literal, SupportsFloat, TypedDict, TypeVar

import numpy as np
import numpy.typing as npt

B = TypeVar("B", bound=npt.NBitBase)

type ArrayLike = npt.ArrayLike
type Scalar = SupportsFloat
type FloatArray = npt.NDArray[np.floating]
type Direction = Literal["u", "v", "w", "U", "V", "W"] | int
type SectionElement = Literal[-1, 0] | None
type Section = tuple[SectionElement, ...]


class SectionKwargs(TypedDict, total=False):
    u: SectionElement
    v: SectionElement
    w: SectionElement
