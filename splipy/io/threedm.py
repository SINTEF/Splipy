from __future__ import annotations

from itertools import chain
from typing import TYPE_CHECKING, Iterable, Iterator, Optional, Protocol, Sequence, Union, cast

import numpy as np
import rhino3dm as rhino
from typing_extensions import Self

from splipy import curve_factory
from splipy.basis import BSplineBasis
from splipy.curve import Curve
from splipy.surface import Surface

from .master import MasterIO

if TYPE_CHECKING:
    from pathlib import Path
    from types import TracebackType

    from splipy.splinemodel import SplineModel
    from splipy.splineobject import SplineObject
    from splipy.types import FArray


# The rhino3dm library has incomplete type information. In particular, some
# classes which are iterable and support the sequence protocol don't advertise
# that fact. All the casting and the protocols here are just to fix that.
# TODO(Eivind): Remove all this cruft if Rhino ever fixes their types.
class Point:
    X: float
    Y: float
    Z: float
    W: float


class CurvePointList(Protocol):
    def __getitem__(self, i: int) -> Point:
        ...

    def __len__(self) -> int:
        ...


class SurfacePointList(Protocol):
    def __getitem__(self, i: tuple[int, int]) -> Point:
        ...


class ThreeDM(MasterIO):
    filename: str
    trimming_curves: list
    fstream: rhino.File3dm

    def __init__(self, filename: Union[str, Path]) -> None:
        self.filename = str(filename)
        self.trimming_curves = []

    def __enter__(self) -> Self:
        self.fstream = rhino.File3dm.Read(self.filename)
        return self

    def __exit__(
        self,
        exc_type: Optional[type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType],
    ) -> None:
        pass

    def write(self, obj: Union[SplineObject, SplineModel, Sequence[SplineObject]]) -> None:
        raise OSError("Writing to 3DM not supported")

    def read(self) -> list[SplineObject]:
        result: list[SplineObject] = []

        for obj in cast(Iterator[rhino.File3dmObject], self.fstream.Objects):
            geom: Union[rhino.GeometryBase, rhino.Polyline] = obj.Geometry

            if isinstance(geom, rhino.Extrusion):
                geom = geom.ToBrep(splitKinkyFaces=True)

            if isinstance(geom, rhino.Brep):
                faces = cast(Sequence[rhino.BrepFace], geom.Faces)
                for idx in range(len(faces)):
                    nsrf = faces[idx].UnderlyingSurface().ToNurbsSurface()
                    result.append(self.read_surface(nsrf))

            if isinstance(geom, rhino.Line):
                a = (geom.From.X, geom.From.Y, geom.From.Z)
                b = (geom.To.X, geom.To.Y, geom.To.Z)
                result.append(curve_factory.line(a, b))
                continue

            if isinstance(geom, rhino.PolylineCurve):
                geom = geom.ToPolyline()

            if isinstance(geom, (rhino.Polyline, rhino.Circle, rhino.Curve, rhino.BezierCurve, rhino.Arc)):
                geom = geom.ToNurbsCurve()

            if isinstance(geom, rhino.NurbsCurve):
                result.append(self.read_curve(geom))

            if isinstance(geom, (rhino.Cylinder, rhino.Sphere, rhino.Surface)):
                geom = geom.ToNurbsSurface()

            if isinstance(geom, rhino.NurbsSurface):
                result.append(self.read_surface(geom))

        return result

    def read_surface(self, nsrf: rhino.NurbsSurface) -> Surface:
        knotsu = np.fromiter(chain([0.0], cast(Iterable[float], nsrf.KnotsU), [0.0]), dtype=float)
        knotsu[0] = knotsu[1]
        knotsu[-1] = knotsu[-2]

        knotsv = np.fromiter(chain([0.0], cast(Iterable[float], nsrf.KnotsV), [0.0]), dtype=float)
        knotsv[0] = knotsv[1]
        knotsv[-1] = knotsv[-2]

        basisu = BSplineBasis(nsrf.OrderU, knotsu, -1)
        basisv = BSplineBasis(nsrf.OrderV, knotsv, -1)

        cpts: FArray = np.ndarray((nsrf.Points.CountU * nsrf.Points.CountV, 3 + nsrf.IsRational), dtype=float)
        points = cast(SurfacePointList, nsrf.Points)
        for v in range(nsrf.Points.CountV):
            count_u = nsrf.Points.CountU
            for u in range(count_u):
                cpts[u + v * count_u, 0] = points[u, v].X
                cpts[u + v * count_u, 1] = points[u, v].Y
                cpts[u + v * count_u, 2] = points[u, v].Z
                if nsrf.IsRational:
                    cpts[u + v * count_u, 3] = points[u, v].W

        return Surface(basisu, basisv, cpts, nsrf.IsRational)

    def read_curve(self, ncrv: rhino.NurbsCurve) -> Curve:
        knots = np.fromiter(chain([0.0], cast(Iterable[float], ncrv.Knots), [0.0]), dtype=float)
        knots[0] = knots[1]
        knots[-1] = knots[-2]

        basis = BSplineBasis(ncrv.Order, knots, -1)

        points = cast(CurvePointList, ncrv.Points)
        cpts: FArray = np.ndarray((len(points), ncrv.Dimension + ncrv.IsRational))
        for u in range(0, len(points)):
            cpts[u, 0] = points[u].X
            cpts[u, 1] = points[u].Y
            if ncrv.Dimension > 2:
                cpts[u, 2] = points[u].Z
            if ncrv.IsRational:
                cpts[u, 3] = points[u].W

        return Curve(basis, cpts, ncrv.IsRational)
