# -*- coding: utf-8 -*-

from splipy.BSplineBasis import BSplineBasis
from splipy.SplineObject import SplineObject
from splipy.Curve import Curve
from splipy.Surface import Surface
from splipy.Volume import Volume
from splipy.TrimmedSurface import TrimmedSurface

# SplineModel imports io.G2 which imports TrimmedSurface,
# so we must import SplineModel after TrimmedSurface
from splipy.SplineModel import SplineModel

__version__ = '1.4.0'
__all__ = ['BSplineBasis',
           'SplineObject',
           'Curve',
           'Surface',
           'Volume',
           'SplineModel',
           'TrimmedSurface',
           'curve_factory',
           'surface_factory',
           'volume_factory',
           'utils',
           'io']
