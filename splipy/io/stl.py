import struct
from abc import ABC, abstractmethod
from typing import Union, Optional, Type, Sequence, TextIO, BinaryIO, cast
from types import TracebackType
from pathlib import Path

import numpy as np
from typing_extensions import Self

from ..surface import Surface
from ..volume import Volume
from ..utils import ensure_listlike
from ..splinemodel import SplineModel
from ..splineobject import SplineObject
from ..types import Scalars

from .master import MasterIO


ASCII_FACET = """facet normal 0 0 0
outer loop
vertex {face[0][0]:.4f} {face[0][1]:.4f} {face[0][2]:.4f}
vertex {face[1][0]:.4f} {face[1][1]:.4f} {face[1][2]:.4f}
vertex {face[2][0]:.4f} {face[2][1]:.4f} {face[2][2]:.4f}
endloop
endfacet
"""

BINARY_HEADER ="80sI"
BINARY_FACET = "12fH"


Face = Sequence[Scalars]


class StlWriter(ABC):

    @abstractmethod
    def _write_header(self) -> None:
        ...

    @abstractmethod
    def _write(self, face: Face) -> None:
        ...

    @abstractmethod
    def close(self) -> None:
        ...

    def _split(self, face: Face) -> tuple[Face, Face]:
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
            raise ValueError('only 3 or 4 vertices for each face')

    def add_faces(self, faces: Sequence[Face]) -> None:
        """ Add many faces. """
        for face in faces:
            self.add_face(face)


class AsciiStlWriter(StlWriter):
    """Export 3D objects build of 3 or 4 vertices as ASCII STL file."""

    fp: TextIO

    def __init__(self, stream: TextIO) -> None:
        self.fp = stream
        self._write_header()

    def _write_header(self) -> None:
        self.fp.write("solid python\n")

    def close(self) -> None:
        self.fp.write("endsolid python\n")

    def _write(self, face: Face) -> None:
        self.fp.write(ASCII_FACET.format(face=face))


class BinaryStlWriter(StlWriter):
    """Export 3D objects build of 3 or 4 vertices as binary STL file."""

    counter: int
    fp: BinaryIO

    def __init__(self, stream: BinaryIO) -> None:
        self.counter = 0
        self.fp = stream
        self._write_header()

    def close(self) -> None:
        self._write_header()

    def _write_header(self) -> None:
        self.fp.seek(0)
        self.fp.write(struct.pack(BINARY_HEADER, b'Python Binary STL Writer', self.counter))

    def _write(self, face: Face) -> None:
        self.counter += 1
        data = [
            0., 0., 0.,
            face[0][0], face[0][1], face[0][2],
            face[1][0], face[1][1], face[1][2],
            face[2][0], face[2][1], face[2][2],
            0
        ]
        self.fp.write(struct.pack(BINARY_FACET, *data))


class STL(MasterIO):

    filename: str
    binary: bool

    writer: Union[BinaryStlWriter, AsciiStlWriter]

    def __init__(self, filename: Union[str, Path], binary: bool = True) -> None:
        self.filename = str(filename)
        self.binary = binary

    def __enter__(self) -> Self:
        if self.binary:
            self.writer = BinaryStlWriter(open(self.filename, 'wb'))
        else:
            self.writer = AsciiStlWriter(open(self.filename, 'w'))
        return self

    def __exit__(
        self,
        exc_type: Optional[Type[BaseException]],
        exc_val: Optional[BaseException],
        exc_tb: Optional[TracebackType]
    ) -> None:
        self.writer.close()
        self.writer.fp.__exit__(exc_type, exc_val, exc_tb)

    def write(
        self,
        obj: Union[SplineModel, Sequence[SplineObject], SplineObject],
        n: Optional[Union[int, Sequence[int]]] = None,
    ) -> None:
        if isinstance(obj, SplineModel):
            if obj.pardim == 3: # volume model
                for node in obj.boundary():
                    self.write_surface(cast(Surface, node.obj), n)
            elif obj.pardim == 2: # surface model
                for surface in obj:
                    self.write_surface(cast(Surface, surface), n)

        elif isinstance(obj, Volume):
            for face in obj.faces():
                if face is not None: # happens with periodic volumes
                    self.write_surface(face, n)

        elif isinstance(obj, Surface):
            self.write_surface(obj, n)

        elif isinstance(obj, Sequence):
            for o in obj:
                self.write(o)

        else:
            raise ValueError('Unsopported object for STL format')

    def write_surface(
        self,
        surface: Surface,
        n: Optional[Union[int, Sequence[int]]] = None,
    ) -> None:
        # choose evaluation points as one of three cases:
        #   1. specified with input
        #   2. linear splines, only picks knots
        #   3. general splines choose 2*order-1 per knot span
        if n is not None:
            n = ensure_listlike(n,2)

        if n is not None:
            u = np.linspace(surface.start(0), surface.end(0), n[0])
        elif surface.order(0) == 2:
            u = surface.knots(0)
        else:
            knots = surface.knots(0)
            p = surface.order(0)
            ut = [np.linspace(k0,k1, 2*p-3, endpoint=False) for (k0,k1) in zip(knots[:-1], knots[1:])]
            ut = [point for element in ut for point in element] + list(knots)
            u = np.sort(ut)

        if n is not None:
            v = np.linspace(surface.start(1), surface.end(1), n[1])
        elif surface.order(1) == 2:
            v = surface.knots(1)
        else:
            knots = surface.knots(1)
            p = surface.order(1)
            vt = [np.linspace(k0,k1, 2*p-3, endpoint=False) for (k0,k1) in zip(knots[:-1], knots[1:])]
            vt = [point for element in vt for point in element] + list(knots)
            v = np.sort(vt)

        # perform evaluation and make sure that we have 3 components (in case of 2D geometries)
        x = surface(u,v)
        if x.shape[2] != 3:
            x.resize((x.shape[0],x.shape[1],3))

        # compute tiny quad pieces
        faces = [[x[i,j], x[i,j+1], x[i+1,j+1], x[i+1,j]] for i in range(x.shape[0]-1) for j in range(x.shape[1]-1)]

        self.writer.add_faces(faces)
