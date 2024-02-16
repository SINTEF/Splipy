from __future__ import annotations

from typing import Any, Union

from numpy.typing import NDArray
from numpy import float_


class Function:
    def grad(self, other: Function) -> Function:
        ...

    def sum(self, axis: int) -> Function:
        ...

    def vector(self, n: int) -> Function:
        ...

    def __mul__(self, other: Any) -> Function:
        ...


class Namespace:
    def __setattr__(
        self,
        name: str,
        value: Union[
            Function,
            int,
            float,
            NDArray[float_],
            str,
        ]
    ) -> None:
        ...

    def eval_nm(self, code: str) -> Function:
        ...

    def eval_n(self, code: str) -> Function:
        ...


def outer(x: Function, y: Function) -> Function:
    ...

def J(x: Function) -> Function:
    ...

def matmat(x: Function, y: NDArray[float_]) -> Function:
    ...
