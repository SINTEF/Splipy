#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"

extern "C" {
  PyObject* Generate_Line(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_LineSegment(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_Circle(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_CircleSegment(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_Ellipse(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_EllipticSegment(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_Helix(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_InterpolateCurve(PyObject* self, PyObject* args, PyObject* kwds);
}
