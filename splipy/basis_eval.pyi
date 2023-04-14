from typing import Tuple

from numpy.typing import NDArray
from numpy import floating, integer


def snap(knots: NDArray[floating], eval_pts: NDArray[floating], tolerance: float) -> None: ...

def evaluate(
    knots: NDArray[floating],
    order: int,
    eval_pts: NDArray[floating],
    periodic: int,
    tolerance: float,
    d: int,
    from_right: bool = True,
) -> Tuple[
    Tuple[
        NDArray[floating],
        NDArray[integer],
        NDArray[integer],
    ],
    Tuple[int, int],
]: ...
