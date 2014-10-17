#pragma once

#include <Python.h>

namespace GeoModeller {

extern "C" {
  void init_SurfaceFactory_Module();
  extern PyObject* SurfaceFactory_module;
}

}
