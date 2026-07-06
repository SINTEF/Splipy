from __future__ import annotations

from itertools import islice
from pathlib import Path
from typing import TYPE_CHECKING, Self, TextIO

import numpy as np

from splipy.basis import BSplineBasis
from splipy.splineobject import SplineObject

from .master import MasterIO

if TYPE_CHECKING:
    from collections.abc import Iterator
    from types import TracebackType

    from splipy.curve import Curve


class SPL(MasterIO):
    fstream: TextIO
    filename: str
    trimming_curves: list[Curve]

    def __init__(self, filename: str) -> None:
        if not filename.endswith(".spl"):
            filename += ".spl"
        self.filename = filename
        self.trimming_curves = []

    def __enter__(self) -> Self:
        self.fstream = Path(self.filename).open()
        return self

    def lines(self) -> Iterator[str]:
        for line in self.fstream:
            yield line.split("#", maxsplit=1)[0].strip()

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

        cpts = np.array(list(islice(lines, totcoeffs * physdim)), dtype=float)
        cpts = cpts.reshape(physdim, *(ncoeffs[::-1])).transpose()

        obj = SplineObject.construct_subclass(bases, cpts, rational=False, raw=True)
        return [obj]

    def __exit__(
        self,
        exc_type: type[BaseException],
        exc_value: BaseException,
        traceback: TracebackType,
    ) -> None:
        self.fstream.close()
