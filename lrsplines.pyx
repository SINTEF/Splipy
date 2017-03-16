# distutils: language = c++

from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref, preincrement as preinc

import numpy as np
from splipy.utils import check_direction


cdef extern from '<iostream>' namespace 'std':
    cdef cppclass istream:
        pass

cdef extern from '<fstream>' namespace 'std':
    cdef cppclass ifstream(istream):
        ifstream(const char *) except +
        bool is_open()
        void close()

cdef extern from 'HashSet.h':
    cdef cppclass HashSet_iterator[T]:
        T operator*()
        HashSet_iterator[T] operator++()
        bool equal(HashSet_iterator[T])

cdef extern from 'Basisfunction.h' namespace 'LR':
    cdef cppclass Basisfunction_ 'LR::Basisfunction':
        int getId()
        void getControlPoint(vector[double]&)

cdef extern from 'Element.h' namespace 'LR':
    cdef cppclass Element_ 'LR::Element':
        int getId()
        int getDim()
        double getParmin(int)
        double getParmax(int)
        HashSet_iterator[Basisfunction_*] supportBegin()
        HashSet_iterator[Basisfunction_*] supportEnd()

cdef extern from 'LRSpline.h' namespace 'LR':
    cdef cppclass LRSpline_ 'LR::LRSpline':
        int dimension()
        int nVariate()
        int nBasisFunctions()
        double startparam(int)
        double endparam(int)
        int order(int)
        vector[Element_*].iterator elementBegin()
        vector[Element_*].iterator elementEnd()
        HashSet_iterator[Basisfunction_*] basisBegin()
        HashSet_iterator[Basisfunction_*] basisEnd()

cdef extern from 'LRSplineSurface.h' namespace 'LR':
    cdef cppclass LRSplineSurface_ 'LR::LRSplineSurface' (LRSpline_):
        LRSplineSurface() except +
        void read(istream) except +


cdef class BasisFunction:

    cdef Basisfunction_* bf

    @property
    def id(self):
        return self.bf.getId()

    @property
    def controlpoint(self):
        cdef vector[double] data
        self.bf.getControlPoint(data)
        return list(data)


cdef class Element:

    cdef Element_* el

    @property
    def id(self):
        return self.el.getId()

    @property
    def pardim(self):
        return self.el.getDim()

    def start(self, direction=None):
        if direction is None:
            return tuple(self.el.getParmin(i) for i in range(self.pardim))
        direction = check_direction(direction, self.pardim)
        return self.el.getParmin(direction)

    def end(self, direction=None):
        if direction is None:
            return tuple(self.el.getParmax(i) for i in range(self.pardim))
        direction = check_direction(direction, self.pardim)
        return self.el.getParmax(direction)

    def basis_functions(self):
        cdef HashSet_iterator[Basisfunction_*] it = self.el.supportBegin()
        cdef HashSet_iterator[Basisfunction_*] end = self.el.supportEnd()
        while not it.equal(end):
            bf = BasisFunction()
            bf.bf = deref(it)
            yield bf
            preinc(it)


cdef class LRSplineObject:

    cdef LRSpline_* lr

    @property
    def dimension(self):
        return self.lr.dimension()

    @property
    def pardim(self):
        return self.lr.nVariate()

    @property
    def controlpoints(self):
        cps = np.empty((len(self), self.dimension))
        for i, bf in enumerate(self.basis_functions()):
            cps[i,:] = bf.controlpoint
        return cps

    def __len__(self):
        return self.lr.nBasisFunctions()

    def start(self, direction=None):
        if direction is None:
            return tuple(self.el.startparam(i) for i in range(self.pardim))
        direction = check_direction(direction, self.pardim)
        return self.el.startparam(direction)

    def end(self, direction=None):
        if direction is None:
            return tuple(self.lr.endparam(i) for i in range(self.pardim))
        direction = check_direction(direction, self.pardim)
        return self.lr.endparam(direction)

    def order(self, direction=None):
        if direction is None:
            return tuple(self.lr.order(i) for i in range(self.pardim))
        direction = check_direction(direction, self.pardim)
        return self.lr.order(direction)

    def elements(self):
        cdef vector[Element_*].iterator it = self.lr.elementBegin()
        cdef vector[Element_*].iterator end = self.lr.elementEnd()
        while it != end:
            el = Element()
            el.el = deref(it)
            yield el
            preinc(it)

    def basis_functions(self):
        cdef HashSet_iterator[Basisfunction_*] it = self.lr.basisBegin()
        cdef HashSet_iterator[Basisfunction_*] end = self.lr.basisEnd()
        while not it.equal(end):
            bf = BasisFunction()
            bf.bf = deref(it)
            yield bf
            preinc(it)


cdef class LRSurface(LRSplineObject):

    @staticmethod
    def from_file(str filename):
        cdef ifstream* stream
        cdef LRSplineSurface_* lr
        stream = new ifstream(filename.encode())
        lr = new LRSplineSurface_()
        surf = LRSurface()
        if stream.is_open():
            lr.read(deref(stream))
            surf.lr = lr
            stream.close()
            return surf
        raise FileNotFoundError()
