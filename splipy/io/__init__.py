from .g2  import G2
from .svg import SVG
from .spl import SPL
from .stl import STL
from .ofoam import OpenFOAM

__all__ = ['G2', 'SVG', 'SPL', 'STL', 'OpenFOAM']

from importlib.util import find_spec


has_cv2 = find_spec('cv2')


# GRDECL depends on optional cv2
has_grdecl = has_cv2
if has_grdecl:
    from .grdecl import GRDECL
    __all__.append('GRDECL')


# ThreeDM depends on optional rhino3dm
has_rhino = find_spec('rhino3dm')
if has_rhino:
    from .threedm import ThreeDM
    __all__.append('ThreeDM')
