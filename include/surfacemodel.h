#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/compositemodel/SurfaceModel.h"

namespace GeoModeller {

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::SurfaceModel> data;
  } SurfaceModel;

  void init_SurfaceModel_Type();

  extern PyTypeObject SurfaceModel_Type;
}

// helpers
void WriteSurfaceModelG2(std::ofstream& g2_file, SurfaceModel* model, bool convert);

}
