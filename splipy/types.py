from collections.abc import Sequence
from typing import Union, Literal

from numpy.typing import NDArray
from numpy import float_, int_


Direction = Union[int, Literal['u', 'U', 'v', 'V', 'w', 'W']]

FArray = NDArray[float_]
IArray = NDArray[int_]

Scalar = Union[
    float,
    float_,
]

Scalars = Union[
    Sequence[Scalar],
    FArray,
]

ScalarOrScalars = Union[Scalar, Scalars]
