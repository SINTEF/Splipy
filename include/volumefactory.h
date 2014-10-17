#pragma once

#include <Python.h>

namespace GeoModeller {

extern "C" {
  void init_VolumeFactory_Module();
  extern PyObject* VolumeFactory_module;
}

}
