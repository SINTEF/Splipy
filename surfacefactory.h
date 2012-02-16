#pragma once

#include <Python.h>

extern "C" {
  void init_SurfaceFactory_Module();
  extern PyObject* SurfaceFactory_module;
}
