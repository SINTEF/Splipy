from .g2 import G2
from .ofoam import OpenFOAM
from .spl import SPL
from .stl import STL
from .svg import SVG

__all__ = ["G2", "SVG", "SPL", "STL", "OpenFOAM"]

from importlib.util import find_spec

has_cv2 = find_spec("cv2")


# GRDECL depends on optional cv2
has_grdecl = has_cv2
if has_grdecl:
    from .grdecl import GRDECL  # noqa
    __all__.append("GRDECL")


# ThreeDM depends on optional rhino3dm
has_rhino = find_spec("rhino3dm")
if has_rhino:
    from .threedm import ThreeDM  # noqa
    __all__.append("ThreeDM")
