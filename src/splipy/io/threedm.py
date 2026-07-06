# Make type-checking work in cases where rhino3dm is not installed
# mypy: disable-error-code="no-any-unimported"

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol, Self, cast

import numpy as np
from rhino3dm import (  # type: ignore[import-not-found,unused-ignore]
    Arc,
    BezierCurve,
    Brep,
    BrepFace,
    Circle,
    Cylinder,
    Extrusion,
    File3dm,
    File3dmObject,
    GeometryBase,
    Line,
    NurbsCurve,
    NurbsSurface,
    Point4d,
    Polyline,
    PolylineCurve,
    Sphere,
)
from rhino3dm import Curve as threedmCurve  # name conflict with splipy
from rhino3dm import Surface as threedmSurface  # name conflict with splipy

from splipy import BSplineBasis, Curve, Surface, curve_factory

from .master import MasterIO

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from types import TracebackType

    from splipy.splinemodel import SplineModel
    from splipy.splineobject import SplineObject
    from splipy.typing import Point


# The rhino3dm type hints are incomplete, hence we have some shims.
class NurbsSurfacePointListShim(Protocol):
    def __getitem__(self, index: Any) -> Point4d: ...


class NurbsCurvePointListShim(Protocol):
    def __getitem__(self, index: Any) -> Point4d: ...
    def __len__(self) -> int: ...


class ThreeDM(MasterIO):
    filename: str
    trimming_curves: list[Curve]

    onlywrite: bool
    fstream: File3dm

    def __init__(self, filename: str) -> None:
        if filename[-4:] != ".3dm":
            filename += ".3dm"
        self.filename = filename
        self.trimming_curves = []

    def __enter__(self) -> Self:
        return self

    def write(self, _: SplineObject | Sequence[SplineObject] | SplineModel) -> None:
        raise OSError("Writing to 3DM not supported")

    def read(self) -> list[SplineObject]:
        if not hasattr(self, "fstream"):
            self.onlywrite = False
            self.fstream = File3dm.Read(self.filename)

        if self.onlywrite:
            raise OSError(f"Could not read from file {self.filename}")

        result: list[SplineObject] = []

        for obj in cast("Iterable[File3dmObject]", self.fstream.Objects):
            geom: GeometryBase | Polyline = obj.Geometry
            if type(geom) is Extrusion:
                geom = geom.ToBrep(splitKinkyFaces=True)
            if type(geom) is Brep:
                for face in cast("Iterable[BrepFace]", geom.Faces):
                    nsrf = face.UnderlyingSurface().ToNurbsSurface()
                    result.append(self.read_surface(nsrf))

            if type(geom) is Line:
                result.append(curve_factory.line(cast("Point", geom.From), cast("Point", geom.To)))
                continue
            if type(geom) is PolylineCurve:
                geom = geom.ToPolyline()
            if (
                type(geom) is Polyline
                or type(geom) is Circle
                or type(geom) is threedmCurve
                or type(geom) is BezierCurve
                or type(geom) is Arc
            ):
                geom = geom.ToNurbsCurve()

            if type(geom) is NurbsCurve:
                result.append(self.read_curve(geom))

            if type(geom) is Cylinder or type(geom) is Sphere or type(geom) is threedmSurface:
                geom = geom.ToNurbsSurface()
            if type(geom) is NurbsSurface:
                result.append(self.read_surface(geom))

        return result

    def read_surface(self, nsrf: NurbsSurface) -> Surface:
        knotsu: list[float] = [0]
        for i in cast("Iterable[float]", nsrf.KnotsU):
            knotsu.append(i)
        knotsu.append(knotsu[len(knotsu) - 1])
        knotsu[0] = knotsu[1]

        knotsv: list[float] = [0]
        for i in cast("Iterable[float]", nsrf.KnotsV):
            knotsv.append(i)
        knotsv.append(knotsv[len(knotsv) - 1])
        knotsv[0] = knotsv[1]

        basisu = BSplineBasis(nsrf.OrderU, knotsu, -1)
        basisv = BSplineBasis(nsrf.OrderV, knotsv, -1)

        cpts = np.ndarray((nsrf.Points.CountU * nsrf.Points.CountV, 3 + nsrf.IsRational))
        pts = cast("NurbsSurfacePointListShim", nsrf.Points)
        for v in range(0, nsrf.Points.CountV):
            for u in range(0, nsrf.Points.CountU):
                cpts[u + v * nsrf.Points.CountU, 0] = pts[u, v].X
                cpts[u + v * nsrf.Points.CountU, 1] = pts[u, v].Y
                cpts[u + v * nsrf.Points.CountU, 2] = pts[u, v].Z
                if nsrf.IsRational:
                    cpts[u + v * nsrf.Points.CountU, 3] = pts[u, v].W

        return Surface(basisu, basisv, cpts, nsrf.IsRational)

    def read_curve(self, ncrv: NurbsCurve) -> Curve:
        knots: list[float] = [0]
        for i in cast("Iterable[float]", ncrv.Knots):
            knots.append(i)
        knots[0] = knots[1]
        knots.append(knots[len(knots) - 1])
        basis = BSplineBasis(ncrv.Order, knots, -1)

        points = cast("NurbsCurvePointListShim", ncrv.Points)

        cpts = np.ndarray((len(points), ncrv.Dimension + ncrv.IsRational))
        for u in range(0, len(points)):
            print(type(points[u]))
            cpts[u, 0] = points[u].X
            cpts[u, 1] = points[u].Y
            if ncrv.Dimension > 2:
                cpts[u, 2] = points[u].Z
            if ncrv.IsRational:
                cpts[u, 3] = points[u].W

        return Curve(basis, cpts, ncrv.IsRational)

    def __exit__(
        self,
        exc_type: type[BaseException],
        exc_value: BaseException,
        traceback: TracebackType,
    ) -> None:
        # Apperently File3DM objects don't need to dedicated cleanup/close code
        pass
