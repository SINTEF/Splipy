#pragma once

#include <Python.h>

extern "C" {
  void init_VolumeFactory_Module();
  extern PyObject* VolumeFactory_module;
}
