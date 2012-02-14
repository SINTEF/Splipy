#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/compositemodel/SurfaceModel.h"

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::SurfaceModel> data;
  } SurfaceModel;

  void init_SurfaceModel_Type();

  extern PyTypeObject SurfaceModel_Type;
}
