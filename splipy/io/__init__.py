from .g2  import G2
from .svg import SVG
from .spl import SPL
from .stl import STL
from .ofoam import OpenFOAM


from importlib.util import find_spec
has_cv2 = find_spec('cv2')
has_rhino = find_spec('rhino3dm')


# GRDECL depends on optional cv2
has_grdecl = has_cv2
if has_cv2:
    pass


# ThreeDM depends on optional rhino3dm
if has_rhino:
    pass


__all__ = ['G2', 'SVG', 'SPL', 'STL', 'OpenFOAM']

if has_grdecl:
    __all__.append('GRDECL')

if has_rhino:
    __all__.append('ThreeDM')
