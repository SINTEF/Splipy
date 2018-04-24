# -*- coding: utf-8 -*-

from splipy.BSplineBasis import BSplineBasis
from splipy.SplineObject import SplineObject
from splipy.Curve import Curve
from splipy.Surface import Surface
from splipy.Volume import Volume
from splipy.SplineModel import SplineModel
from splipy.TrimmedSurface import TrimmedSurface

__version__ = '1.3.1'
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
