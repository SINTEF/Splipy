#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"

namespace GeoModeller {

extern "C" {
  void init_Preprocess_Module();
  extern PyObject* Preprocess_module;
}

}
