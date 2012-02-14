#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"

extern "C" {
  PyObject* Generate_Box(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_Cylinder(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_LoftSurfaces(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_Parallelepiped(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_SweepSurfaceLinear(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_SweepSurfaceRotational(PyObject* self, PyObject* args, PyObject* kwds);
}
