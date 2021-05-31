from typing import List

import numpy as np


class Transform:
    """Superclass for control point transformations."""

    def __call__(self, cps: np.ndarray) -> np.ndarray:
        ...

    def __mul__(self, other: 'Transform') -> 'Transform':
        """Compose transforms in a sequence.

            (trf1 * trf2)(cps) == trf2(trf1(cps))
        """
        if not isinstance(other, Transform):
            return NotImplemented
        if isinstance(other, ChainedTransform):
            return ChainedTransform(self, *other.sequence)
        return ChainedTransform(self, other)


class ChainedTransform:

    sequence: List[Transform]

    def __init__(self, *sequence: Transform):
        self.sequence = list(sequence)

    def __call__(self, cps: np.ndarray) -> np.ndarray:
        for trf in self.sequence:
            cps = trf(cps)
        return cps


class SwapTransform:

    dir1: int
    dir2: int

    def __init__(self, dir1: int, dir2: int):
        self.dir1 = dir1
        self.dir2 = dir2

    def __call__(self, cps: np.ndarray) -> np.ndarray:
        new_directions = list(range(cps.ndim))
        new_directions[self.dir1] = self.dir2
        new_directions[self.dir2] = self.dir1
        return cps.transpose(new_directions)
