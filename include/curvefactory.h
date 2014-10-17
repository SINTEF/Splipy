#pragma once

#include <Python.h>

namespace GeoModeller {

extern "C" {
  void init_CurveFactory_Module();
  extern PyObject* CurveFactory_module;
}

}
