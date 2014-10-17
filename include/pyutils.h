#pragma once

#include <Python.h>
#include "GoTools/geometry/GoTools.h"

#include "point.h"
#include "curve.h"
#include "surface.h"
#include "volume.h"
#include "surfacemodel.h"
#include "volumemodel.h"

namespace GeoModeller {

void InitializeTypeObject(PyTypeObject* type_object);
shared_ptr<Go::Point> PyObject_AsGoPoint(PyObject* obj);
shared_ptr<Go::ParamCurve> PyObject_AsGoCurve(PyObject* obj);
shared_ptr<Go::ParamSurface> PyObject_AsGoSurface(PyObject* obj);
shared_ptr<Go::ParamVolume> PyObject_AsGoVolume(PyObject* obj);
shared_ptr<Go::SurfaceModel> PyObject_AsGoSurfaceModel(PyObject* obj);
shared_ptr<Go::VolumeModel> PyObject_AsGoVolumeModel(PyObject* obj);

void PyMethods_Append(std::vector<PyMethodDef>& defs, PyMethodDef* start);

}
