#pragma once

#include "curve.h"
#include "surface.h"
#include "volume.h"

namespace GeoModeller {

void DoWriteHDF5Field(const std::string& fname, const std::string& fldname,
                      int patch, int level, PyObject* data, bool trunc);
void DoWriteHDF5Geometry(const std::string& fname, const std::string& fldname,
                         int patch, int level, PyObject* data, bool trunc);

PyObject* DoReadHDF5Field(const std::string& fname, const std::string& fldname,
                          int patch, int level);
PyObject* DoReadHDF5Geometry(const std::string& fname, const std::string& fldname,
                             int patch, int level);
}
