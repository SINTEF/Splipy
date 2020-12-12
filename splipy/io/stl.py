#coding:utf-8

import struct

import numpy as np

from ..surface import Surface
from ..volume import Volume
from ..utils import ensure_listlike
from ..splinemodel import SplineModel

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


class ASCII_STL_Writer(object):
    """ Export 3D objects build of 3 or 4 vertices as ASCII STL file.
    """
    def __init__(self, stream):
        self.fp = stream
        self._write_header()

    def _write_header(self):
        self.fp.write("solid python\n")

    def close(self):
        self.fp.write("endsolid python\n")

    def _write(self, face):
        self.fp.write(ASCII_FACET.format(face=face))

    def _split(self, face):
        p1, p2, p3, p4 = face
        return (p1, p2, p3), (p3, p4, p1)

    def add_face(self, face):
        """ Add one face with 3 or 4 vertices. """
        if len(face) == 4:
            face1, face2 = self._split(face)
            self._write(face1)
            self._write(face2)
        elif len(face) == 3:
            self._write(face)
        else:
            raise ValueError('only 3 or 4 vertices for each face')

    def add_faces(self, faces):
        """ Add many faces. """
        for face in faces:
            self.add_face(face)

class BINARY_STL_Writer(ASCII_STL_Writer):
    """ Export 3D objects build of 3 or 4 vertices as binary STL file.
    """
    def __init__(self, stream):
        self.counter = 0
        #### new-style classes way of calling super constructor
        # super(Binary_STL_Writer, self).__init__(stream)

        #### old-style classes way of doing it
        ASCII_STL_Writer.__init__(self, stream)

    def close(self):
        self._write_header()

    def _write_header(self):
        self.fp.seek(0)
        self.fp.write(struct.pack(BINARY_HEADER, b'Python Binary STL Writer', self.counter))

    def _write(self, face):
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
    def __init__(self, filename, binary=True):
        if filename[-4:] != '.stl':
            filename += '.stl'
        self.filename = filename
        self.binary   = binary

    def __enter__(self):
        if self.binary:
            fp = open(self.filename, 'wb')
            self.writer = BINARY_STL_Writer(fp)
        else:
            fp = open(self.filename, 'w')
            self.writer = ASCII_STL_Writer(fp)
        return self

    def write(self, obj, n=None):
        if isinstance(obj, SplineModel):
            if obj.pardim == 3: # volume model
                for surface in obj.boundary():
                    self.write_surface(surface.obj,n)
            elif obj.pardim == 2: # surface model
                for surface in obj:
                    self.write_surface(surface, n)

        elif isinstance(obj, Volume):
            for surface in obj.faces():
                self.write_surface(surface, n)

        elif isinstance(obj, Surface):
            self.write_surface(obj, n)

        else:
            raise ValueError('Unsopported object for STL format')

    def write_surface(self, surface, n=None):
        # choose evaluation points as one of three cases:
        #   1. specified with input
        #   2. linear splines, only picks knots
        #   3. general splines choose 2*order-1 per knot span
        if n != None:
            n = ensure_listlike(n,2)

        if n != None:
            u = np.linspace(surface.start(0), surface.end(0), n[0])
        elif surface.order(0) == 2:
            u = surface.knots(0)
        else:
            knots = surface.knots(0)
            p = surface.order(0)
            u = [np.linspace(k0,k1, 2*p-3, endpoint=False) for (k0,k1) in zip(knots[:-1], knots[1:])]
            u = [point for element in u for point in element] + knots
            u = np.sort(u)

        if n != None:
            v = np.linspace(surface.start(1), surface.end(1), n[1])
        elif surface.order(1) == 2:
            v = surface.knots(1)
        else:
            knots = surface.knots(1)
            p = surface.order(1)
            v = [np.linspace(k0,k1, 2*p-3, endpoint=False) for (k0,k1) in zip(knots[:-1], knots[1:])]
            v = [point for element in v for point in element] + knots
            v = np.sort(v)

        # perform evaluation and make sure that we have 3 components (in case of 2D geometries)
        x = surface(u,v)
        if x.shape[2] != 3:
            x.resize((x.shape[0],x.shape[1],3))

        # compute tiny quad pieces
        faces = [[x[i,j], x[i,j+1], x[i+1,j+1], x[i+1,j]] for i in range(x.shape[0]-1) for j in range(x.shape[1]-1)]

        self.writer.add_faces(faces)

    def __exit__(self, exc_type, exc_value, traceback):
        self.writer.close()
        self.writer.fp.close()

