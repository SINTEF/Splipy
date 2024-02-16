from __future__ import annotations

from itertools import islice
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, Optional, Sequence, TextIO, Union

import numpy as np
from typing_extensions import Self

from splipy.basis import BSplineBasis
from splipy.splineobject import SplineObject

from .master import MasterIO

if TYPE_CHECKING:
    from types import TracebackType

    from splipy.splinemodel import SplineModel


class SPL(MasterIO):
    filename: Path
    fstream: TextIO

    def __init__(self, filename: Union[Path, str]) -> None:
        self.filename = Path(filename)

    def __enter__(self) -> Self:
        self.fstream = self.filename.open("r").__enter__()
        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        self.fstream.__exit__(exc_type, exc_val, exc_tb)

    def lines(self) -> Iterator[str]:
        for line in self.fstream:
            yield line.split("#", maxsplit=1)[0].strip()

    def write(self, obj: Union[SplineObject, SplineModel, Sequence[SplineObject]]) -> None:
        raise OSError("Writing to SPL not supported")

    def read(self) -> list[SplineObject]:
        lines = self.lines()

        version = next(lines).split()
        assert version[0] == "C"
        assert version[3] == "0"  # No support for rational SPL yet
        pardim = int(version[1])
        physdim = int(version[2])

        orders = [int(k) for k in islice(lines, pardim)]
        ncoeffs = [int(k) for k in islice(lines, pardim)]
        totcoeffs = int(np.prod(ncoeffs))
        nknots = [a + b for a, b in zip(orders, ncoeffs)]

        next(lines)  # Skip spline accuracy

        knots = [[float(k) for k in islice(lines, nkts)] for nkts in nknots]
        bases = [BSplineBasis(p, kts, -1) for p, kts in zip(orders, knots)]

        cpts = np.array([float(k) for k in islice(lines, totcoeffs * physdim)], dtype=float)
        cpts = cpts.reshape(physdim, *(ncoeffs[::-1])).transpose()

        if 1 <= pardim <= 3:
            patch = SplineObject.constructor(pardim)(bases, cpts, raw=True)
        else:
            patch = SplineObject(bases, controlpoints=cpts, raw=True)

        return [patch]
