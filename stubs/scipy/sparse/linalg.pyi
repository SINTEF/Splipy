from typing import TypeVar

from numpy.typing import NDArray
from numpy import floating

from . import csr_matrix


T = TypeVar("T", covariant=True, bound=floating)

def spsolve(mx: csr_matrix, rhs: NDArray[T]) -> NDArray[T]:
    ...
