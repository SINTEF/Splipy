#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/SplineVolume.h"

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::ParamVolume> data;
  } Volume;

  void init_Volume_Type();

  extern PyTypeObject Volume_Type;
}
shared_ptr<Go::SplineVolume> convertSplineVolume(shared_ptr<Go::ParamVolume> volume);
