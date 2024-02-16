from .basis import BSplineBasis
from .curve import Curve
from .splinemodel import SplineModel
from .splineobject import SplineObject
from .surface import Surface
from .trimmedsurface import TrimmedSurface
from .volume import Volume

__version__ = "1.8.2"
__all__ = [
    "BSplineBasis",
    "SplineObject",
    "Curve",
    "Surface",
    "Volume",
    "TrimmedSurface",
    "SplineModel",
]
