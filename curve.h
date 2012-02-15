#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/SplineCurve.h"

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::ParamCurve> data;
  } Curve;

  void init_Curve_Type();

  extern PyTypeObject Curve_Type;
}

// helpers
shared_ptr<Go::SplineCurve> convertSplineCurve(shared_ptr<Go::ParamCurve> curve);
void WriteCurveG2(std::ofstream& g2_file, Curve* curve, bool convert);
