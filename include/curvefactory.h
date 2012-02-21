#pragma once

#include <Python.h>

extern "C" {
  void init_CurveFactory_Module();
  extern PyObject* CurveFactory_module;
}
