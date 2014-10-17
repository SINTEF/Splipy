#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/trivariatemodel/VolumeModel.h"

namespace GeoModeller {

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::VolumeModel> data;
  } VolumeModel;

  void init_VolumeModel_Type();

  extern PyTypeObject VolumeModel_Type;
}

// helpers
void WriteVolumeModelG2(std::ofstream& g2_file, VolumeModel* model, bool convert);

}
