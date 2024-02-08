from numpy.typing import NDArray
from numpy import float_, int_


def snap(knots: NDArray[float_], eval_pts: NDArray[float_], tolerance: float) -> None: ...

def evaluate(
    knots: NDArray[float_],
    order: int,
    eval_pts: NDArray[float_],
    periodic: int,
    tolerance: float,
    d: int,
    from_right: bool = True,
) -> tuple[
    tuple[
        NDArray[float_],
        NDArray[int_],
        NDArray[int_],
    ],
    tuple[int, int],
]: ...
