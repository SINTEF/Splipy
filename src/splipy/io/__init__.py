from __future__ import annotations

import importlib.util

from .g2 import G2
from .ofoam import OpenFOAM
from .spl import SPL
from .stl import STL
from .svg import SVG

has_grdecl = (
    bool(importlib.util.find_spec("h5py"))
    and bool(importlib.util.find_spec("cv2"))
    and bool(importlib.util.find_spec("tqdm"))
)
has_rhino = importlib.util.find_spec("rhino3dm")

# GRDECL depends on optional cv2, h5py, tqdm
if has_grdecl:
    from .grdecl import GRDECL  # noqa: F401

# ThreeDM depends on optional rhino3dm
if has_rhino:
    from .threedm import ThreeDM  # noqa: F401


__all__ = ["G2", "SVG", "SPL", "STL", "OpenFOAM"]

if has_grdecl:
    __all__.append("GRDECL")

if has_rhino:
    __all__.append("ThreeDM")
