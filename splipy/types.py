from typing import Literal, Union, Protocol, Any, TypeVar, SupportsIndex, Iterator
import numpy as np


Direction = Union[int, Literal['u', 'U', 'v', 'V', 'w', 'W']]

T_co = TypeVar("T_co", covariant=True)


class Float(Protocol):
    def __float__(self) -> float: ...
    def __lt__(self, other: Any) -> bool: ...
    def __gt__(self, other: Any) -> bool: ...


class Seq(Protocol[T_co]):
    def __len__(self) -> int: ...
    def __getitem__(self, index: SupportsIndex, /) -> T_co: ...
    def __iter__(self) -> Iterator[T_co]: ...


OneOrMore = Union[T_co, Seq[T_co]]
