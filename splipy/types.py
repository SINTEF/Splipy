from collections.abc import Sequence
from typing import Literal, TypedDict, Union

from numpy import float64, int_
from numpy.typing import NDArray

Direction = Union[int, Literal["u", "U", "v", "V", "w", "W"]]

FArray = NDArray[float64]
IArray = NDArray[int_]

Scalar = Union[
    float,
    float64,
]

Scalars = Union[
    Sequence[Scalar],
    FArray,
]

ScalarOrScalars = Union[Scalar, Scalars]

SectionElt = Literal[-1, 0, None]
SectionLike = Sequence[SectionElt]
Section = tuple[SectionElt, ...]


class SectionKwargs(TypedDict, total=False):
    u: SectionElt
    v: SectionElt
    w: SectionElt
