from __future__ import annotations

from itertools import islice
from pathlib import Path

import numpy as np

from splipy.basis import BSplineBasis
from splipy.curve import Curve
from splipy.splineobject import SplineObject
from splipy.surface import Surface
from splipy.volume import Volume

from .master import MasterIO


class SPL(MasterIO):
    def __init__(self, filename):
        if not filename.endswith(".spl"):
            filename += ".spl"
        self.filename = filename
        self.trimming_curves = []

    def __enter__(self):
        self.fstream = Path(self.filename).open()
        return self

    def lines(self):
        for line in self.fstream:
            yield line.split("#", maxsplit=1)[0].strip()

    def read(self):
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

        cpts = np.array([float(k) for k in islice(lines, totcoeffs * physdim)])
        cpts = cpts.reshape(physdim, *(ncoeffs[::-1])).transpose()

        if pardim == 1:
            patch = Curve(*bases, controlpoints=cpts, raw=True)
        elif pardim == 2:
            patch = Surface(*bases, controlpoints=cpts, raw=True)
        elif pardim == 3:
            patch = Volume(*bases, controlpoints=cpts, raw=True)
        else:
            patch = SplineObject(bases, controlpoints=cpts, raw=True)

        return [patch]

    def __exit__(self, exc_type, exc_value, traceback):
        self.fstream.close()
