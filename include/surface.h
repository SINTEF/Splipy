#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/geometry/ParamSurface.h"

namespace GeoModeller {

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::ParamSurface> data;
  } Surface;

  void init_Surface_Type();

  extern PyTypeObject Surface_Type;
}

// helper functions
shared_ptr<Go::SplineSurface> convertSplineSurface(shared_ptr<Go::ParamSurface> surface);
void printSurfaceToStream(std::ostream& str, shared_ptr<Go::ParamSurface> surf);
void WriteSurfaceG2(std::ofstream& g2_file, Surface* surface, bool convert);

}
