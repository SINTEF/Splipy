from __future__ import annotations

from .basis import BSplineBasis
from .curve import Curve
from .splinemodel import SplineModel
from .splineobject import SplineObject
from .surface import Surface
from .trimmedsurface import TrimmedSurface
from .volume import Volume

__all__ = ["BSplineBasis", "Curve", "Surface", "Volume", "SplineObject", "SplineModel", "TrimmedSurface"]

__version__ = "1.10.1"
