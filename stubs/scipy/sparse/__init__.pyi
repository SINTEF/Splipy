from typing import overload, TypeVar

from numpy.typing import NDArray
from numpy import float_, int_, generic


T = TypeVar("T", covariant=True, bound=generic)


class csr_matrix:
    @overload
    def __init__(self, shape: tuple[int, int]) -> None:
        ...

    @overload
    def __init__(
        self,
        data: tuple[NDArray[float_], NDArray[int_], NDArray[int_]],
        shape: tuple[int, int],
    ) -> None:
        ...

    def toarray(self) -> NDArray[float_]:
        ...

    def __matmul__(self, other: NDArray[T]) -> NDArray[T]:
        ...
