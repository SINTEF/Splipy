from numpy.typing import NDArray
from numpy import float_

from .function import Function
from .sample import Integral


class newton:
    def __init__(
        self,
        x: str,
        residual: Integral,
        constrain: NDArray[float_],
    ) -> None:
        ...

    def solve(self, tol: float, maxiter: int) -> NDArray[float_]:
        ...
