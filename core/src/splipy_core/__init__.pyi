from numpy import floating, int_, integer
from numpy.typing import NDArray

type Scalar = float | floating | int | integer

def snap_point(knots: NDArray[floating], eval_pt: Scalar, tolerance: Scalar) -> Scalar: ...
def snap_points(knots: NDArray[floating], eval_pts: NDArray[floating], tolerance: Scalar) -> None: ...
def evaluate(
    knots: NDArray[floating],
    order: int,
    eval_pts: NDArray[floating],
    periodic: int,
    tolerance: Scalar,
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
