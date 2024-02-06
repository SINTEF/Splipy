from __future__ import annotations

from typing import SupportsFloat, Sequence, Protocol, SupportsIndex, Iterator, Any, TypeVar, List, overload, Union
import numpy as np
from numpy.typing import ArrayLike


# T_co = TypeVar("T_co", covariant=True)


# class Seq(Protocol[T_co]):
#     def __len__(self) -> int: ...
#     @overload
#     def __getitem__(self, index: SupportsIndex, /) -> T_co: ...
#     @overload
#     def __getitem__(self, index: slice, /) -> Seq[T_co]: ...
#     def __iter__(self) -> Iterator[T_co]: ...
#     def __contains__(self, x: Any) -> bool: ...
#     def __reversed__(self, /) -> Iterator[T_co]: ...
#     def count(self, value: Any, /) -> int: ...
#     def index(self, value: Any, /) -> int: ...


# b: Seq[SupportsFloat]

# a: ArrayLike = b

f: float
n: np.floating
u: Union[float, np.floating]

reveal_type(abs(f))
reveal_type(abs(n))
reveal_type(abs(u))
