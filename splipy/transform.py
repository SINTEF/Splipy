from typing import List, Optional

import numpy as np


def evaluate(bases, cps, tensor=True):
    if tensor:
        idx = len(bases) - 1
        for N in bases[::-1]:
            cps = np.tensordot(N, cps, axes=(1, idx))
    else:
        cps = np.einsum('ij,j...->i...', bases[0], cps)
        for N in bases[1:]:
            cps = np.einsum('ij,ij...->i...', N, cps)
    return cps


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


class Evaluator:

    ns: List[np.ndarray]
    tensor: bool
    rational: bool
    squeeze: bool

    def __init__(self, ns: List[np.ndarray], tensor: bool, rational: bool, squeeze: bool):
        self.ns = ns
        self.tensor = tensor
        self.rational = rational
        self.squeeze = squeeze

    def __call__(self, cps: np.ndarray) -> np.ndarray:
        result = evaluate(self.ns, cps, self.tensor)

        if self.rational:
            for i in range(cps.shape[-1]):
                result[..., i] /= result[..., -1]
            result = np.delete(result, cps.shape[-1] - 1, -1)

        if self.squeeze:
            result = result.reshape(-1)

        return result


class DerivativeEvaluator:

    dns: List[np.ndarray]
    tensor: bool
    rational: Optional[List[np.ndarray]]
    squeeze: bool

    def __init__(self, dns: List[np.ndarray], tensor: bool, rational: Optional[List[np.ndarray]], squeeze: bool):
        self.dns = dns
        self.tensor = tensor
        self.rational = rational
        self.squeeze = squeeze

    def __call__(self, cps: np.ndarray) -> np.ndarray:
        result = evaluate(self.dns, cps, self.tensor)

        if self.rational is not None:
            non_derivative = evaluate(self.rational, cps, self.tensor)
            W = non_derivative[..., -1]  # W
            Wd = result[..., -1]         # W'
            for i in range(cps.shape[-1] - 1):
                result[..., i] = result[..., i] / W - non_derivative[..., i] * Wd / W / W
            result = np.delete(result, cps.shape[-1] - 1, -1)

        if self.squeeze:
            result = result.reshape(-1)

        return result
