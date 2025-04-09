from numpy import floating, int_
from numpy.typing import NDArray

def snap(knots: NDArray[floating], eval_pts: NDArray[floating], tolerance: float) -> None: ...
def evaluate(
    knots: NDArray[floating],
    order: int,
    eval_pts: NDArray[floating],
    periodic: int,
    tolerance: float,
    d: int,
    from_right: bool = True,
) -> tuple[
    tuple[
        NDArray[floating],
        NDArray[int_],
        NDArray[int_],
    ],
    tuple[int, int],
]: ...
