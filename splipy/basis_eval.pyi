from numpy import float64, int_
from numpy.typing import NDArray

def snap(knots: NDArray[float64], eval_pts: NDArray[float64], tolerance: float) -> None: ...
def evaluate(
    knots: NDArray[float64],
    order: int,
    eval_pts: NDArray[float64],
    periodic: int,
    tolerance: float,
    d: int,
    from_right: bool = True,
) -> tuple[
    tuple[
        NDArray[float64],
        NDArray[int_],
        NDArray[int_],
    ],
    tuple[int, int],
]: ...
