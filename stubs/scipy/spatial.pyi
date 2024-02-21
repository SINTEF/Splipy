from numpy.typing import NDArray
import numpy as np


class ConvexHull:
    vertices: NDArray[np.int_]

    def __init__(self, points: NDArray[np.float_]) -> None:
        ...
