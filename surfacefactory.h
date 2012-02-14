#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"

extern "C" {
  PyObject* Generate_Plane(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_CircularDisc(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_SphereSurface(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_CylinderSurface(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_ConeSurface(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_TorusSurface(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_Rectangle(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_SweepCurveLinear(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_SweepCurveRotational(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_TrimSurface(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_AddLoop(PyObject* self, PyObject* args, PyObject* kwds);
  PyObject* Generate_LoftCurves(PyObject* self, PyObject* args, PyObject* kwds);
}
