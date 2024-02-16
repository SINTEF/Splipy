from typing import Sequence, Literal
from typing_extensions import Self

from numpy.typing import NDArray
from numpy import float_

from .function import Function
from .sample import Integral


class Matrix:
    def solve(
        self,
        rhs: NDArray[float_],
        constrain: NDArray[float_],
    ) -> NDArray[float_]:
        ...


class Domain:
    def basis(
        self,
        name: Literal["spline"],
        degree: Sequence[int],
        knotmultiplicities: Sequence[Sequence[int]]
    ) -> Function:
        ...

    def integrate(
        self,
        x: Function,
        ischeme: str,
    ) -> Matrix:
        ...

    def integral(
        self,
        x: Function,
        degree: int,
    ) -> Integral:
        ...


def rectilinear(args: Sequence[NDArray[float_]]) -> tuple[Domain, Function]:
    ...
