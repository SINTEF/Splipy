#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"

namespace GeoModeller {

extern "C" {
  void init_SurfaceModelFactory_Module();
  extern PyObject* SurfaceModelFactory_module;
}

}
