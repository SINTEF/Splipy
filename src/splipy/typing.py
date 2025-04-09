from __future__ import annotations

from typing import Literal, SupportsFloat, TypeAlias, TypedDict, TypeVar

import numpy as np
import numpy.typing as npt

B = TypeVar("B", bound=npt.NBitBase)

ArrayLike: TypeAlias = npt.ArrayLike
ScalarLike: TypeAlias = SupportsFloat
FloatArray: TypeAlias = npt.NDArray[np.floating]
Direction: TypeAlias = Literal["u", "v", "w", "U", "V", "W"] | int
SectionElement: TypeAlias = Literal[-1, 0] | None
Section: TypeAlias = tuple[SectionElement, ...]


class SectionKwargs(TypedDict, total=False):
    u: SectionElement
    v: SectionElement
    w: SectionElement
