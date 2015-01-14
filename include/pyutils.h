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

//! \brief Convert a vector of doubles to a python point list
//! \param result The resulting PyList. Needs to be resized up front!
//! \param data The data
//! \param dim The dimensionality of the data
void VectorToPyPointList(PyObject* result,
                         const std::vector<double>& data, int dim);

//! \brief Tesselate a spline basis
//! \param begin Iterator to starting parameter
//! \param end Iterator past ending parameter
//! \param n Number of tesselation points per knot span
  template<typename RandomIterator>
std::vector<double> Tesselate(RandomIterator begin, RandomIterator end, int n)
{
  RandomIterator uit = begin;
  std::vector<double> gpar;
  double ucurr = 0.0, uprev = *(uit++);
  while (uit != end)
  {
    ucurr = *(uit++);
    if (ucurr > uprev) {
      if (n == 1)
        gpar.push_back(uprev);
      else for (int i = 0; i < n; i++)
      {
        double xg = (double)(2*i-n)/(double)n;
        gpar.push_back(0.5*(ucurr*(1.0+xg) + uprev*(1.0-xg)));
      }
    }
    uprev = ucurr;
  }

  if (ucurr > gpar.back())
    gpar.push_back(ucurr);

  return gpar;
}

}
