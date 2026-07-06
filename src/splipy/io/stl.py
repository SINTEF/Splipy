from __future__ import annotations

import struct
from abc import ABC, abstractmethod
from itertools import chain
from pathlib import Path
from typing import TYPE_CHECKING, BinaryIO, Self, TextIO

import numpy as np

from splipy.splinemodel import SplineModel
from splipy.splineobject import SplineObject
from splipy.surface import Surface
from splipy.utils import ensure_listlike
from splipy.volume import Volume

from .master import MasterIO

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence
    from types import TracebackType

    from splipy.typing import FloatArray

ASCII_FACET = """facet normal 0 0 0
outer loop
vertex {face[0][0]:.4f} {face[0][1]:.4f} {face[0][2]:.4f}
vertex {face[1][0]:.4f} {face[1][1]:.4f} {face[1][2]:.4f}
vertex {face[2][0]:.4f} {face[2][1]:.4f} {face[2][2]:.4f}
endloop
endfacet
"""

BINARY_HEADER = "80sI"
BINARY_FACET = "12fH"


type Point = tuple[float, float, float]
type TriFace = tuple[Point, Point, Point]
type QuadFace = tuple[Point, Point, Point, Point]
type Face = TriFace | QuadFace


class STL_Writer[IO: (TextIO, BinaryIO)](ABC):
    fp: IO

    def __init__(self, stream: IO) -> None:
        self.fp = stream
        self._write_header()

    @abstractmethod
    def _write_header(self) -> None: ...

    @abstractmethod
    def close(self) -> None: ...

    @abstractmethod
    def _write(self, face: TriFace) -> None: ...

    def _split(self, face: QuadFace) -> tuple[TriFace, TriFace]:
        p1, p2, p3, p4 = face
        return (p1, p2, p3), (p3, p4, p1)

    def add_face(self, face: Face) -> None:
        """Add one face with 3 or 4 vertices."""
        if len(face) == 4:
            face1, face2 = self._split(face)
            self._write(face1)
            self._write(face2)
        elif len(face) == 3:
            self._write(face)
        else:
            raise ValueError("only 3 or 4 vertices for each face")

    def add_faces(self, faces: Iterable[Face]) -> None:
        """Add many faces."""
        for face in faces:
            self.add_face(face)


class ASCII_STL_Writer(STL_Writer[TextIO]):
    """Export 3D objects build of 3 or 4 vertices as ASCII STL file."""

    def _write_header(self) -> None:
        self.fp.write("solid python\n")

    def close(self) -> None:
        self.fp.write("endsolid python\n")

    def _write(self, face: TriFace) -> None:
        self.fp.write(ASCII_FACET.format(face=face))


class BINARY_STL_Writer(STL_Writer[BinaryIO]):
    """Export 3D objects build of 3 or 4 vertices as binary STL file."""

    counter: int

    def __init__(self, stream: BinaryIO) -> None:
        self.counter = 0
        super().__init__(stream)

    def _write_header(self) -> None:
        self.fp.seek(0)
        self.fp.write(struct.pack(BINARY_HEADER, b"Python Binary STL Writer", self.counter))

    def close(self) -> None:
        self._write_header()

    def _write(self, face: TriFace) -> None:
        self.counter += 1
        data = [
            0.0,
            0.0,
            0.0,
            face[0][0],
            face[0][1],
            face[0][2],
            face[1][0],
            face[1][1],
            face[1][2],
            face[2][0],
            face[2][1],
            face[2][2],
            0,
        ]
        self.fp.write(struct.pack(BINARY_FACET, *data))


class STL(MasterIO):
    filename: str
    binary: bool

    writer: STL_Writer

    def __init__(self, filename: str, binary: bool = True) -> None:
        if filename[-4:] != ".stl":
            filename += ".stl"
        self.filename = filename
        self.binary = binary

    def __enter__(self) -> Self:
        if self.binary:
            self.writer = BINARY_STL_Writer(Path(self.filename).open("wb"))
        else:
            self.writer = ASCII_STL_Writer(Path(self.filename).open("w"))
        return self

    def write(
        self, obj: SplineObject | Sequence[SplineObject] | SplineModel, n: int | Sequence[int] | None = None
    ) -> None:
        if isinstance(obj, SplineModel):
            if obj.pardim == 3:  # volume model
                for node in obj.boundary():
                    assert isinstance(node.obj, Surface)
                    self.write_surface(node.obj, n)
            elif obj.pardim == 2:  # surface model
                for comp in obj.objects():
                    assert isinstance(comp, Surface)
                    self.write_surface(comp, n)

        elif isinstance(obj, Volume):
            for surf in obj.faces():
                if surf is not None:  # happens with periodic volumes
                    self.write_surface(surf, n)

        elif isinstance(obj, Surface):
            self.write_surface(obj, n)

        elif isinstance(obj, SplineObject):
            raise ValueError("Unsopported object for STL format")

        else:
            for sub in obj:
                self.write(sub)

    def write_surface(self, surface: Surface, n: int | Sequence[int] | None = None) -> None:
        # choose evaluation points as one of three cases:
        #   1. specified with input
        #   2. linear splines, only picks knots
        #   3. general splines choose 2*order-1 per knot span
        if n is not None:
            n = ensure_listlike(n, 2)

        u: FloatArray
        if n is not None:
            u = np.linspace(surface.start(0), surface.end(0), n[0])
        elif surface.order(0) == 2:
            u = surface.knots(0)
        else:
            knots = surface.knots(0)
            p = surface.order(0)

            linear_knots = chain.from_iterable(
                np.linspace(k0, k1, 2 * p - 3, endpoint=False) for k0, k1 in zip(knots[:-1], knots[1:])
            )
            all_knots = chain(linear_knots, knots)
            u = np.sort(np.fromiter(all_knots, dtype=float))

        v: FloatArray
        if n is not None:
            v = np.linspace(surface.start(1), surface.end(1), n[1])
        elif surface.order(1) == 2:
            v = surface.knots(1)
        else:
            knots = surface.knots(1)
            p = surface.order(1)

            linear_knots = chain.from_iterable(
                np.linspace(k0, k1, 2 * p - 3, endpoint=False) for k0, k1 in zip(knots[:-1], knots[1:])
            )
            all_knots = chain(linear_knots, knots)
            v = np.sort(np.fromiter(all_knots, dtype=float))

        # perform evaluation and make sure that we have 3 components (in case of 2D geometries)
        x = surface(u, v)
        if x.shape[2] != 3:
            x.resize((x.shape[0], x.shape[1], 3))

        # compute tiny quad pieces
        faces = [
            (x[i, j], x[i, j + 1], x[i + 1, j + 1], x[i + 1, j])
            for i in range(x.shape[0] - 1)
            for j in range(x.shape[1] - 1)
        ]

        self.writer.add_faces(faces)

    def __exit__(
        self,
        exc_type: type[BaseException],
        exc_value: BaseException,
        traceback: TracebackType,
    ) -> None:
        self.writer.close()
        self.writer.fp.close()
