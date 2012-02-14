#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/utils/Point.h"

extern "C" {
  typedef struct {
    PyObject_HEAD
    shared_ptr<Go::Point> data;
  } Point;

  void init_Point_Type();

  extern PyTypeObject Point_Type;
}
