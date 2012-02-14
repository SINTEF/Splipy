#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"

extern "C" {
  PyObject* Generate_RegularizeSurface(PyObject* self, PyObject* args, PyObject* kwds);
}
