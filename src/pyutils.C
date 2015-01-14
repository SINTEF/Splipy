#include "pyutils.h"

namespace GeoModeller {

void InitializeTypeObject(PyTypeObject* type_object)
{ 
  static PyTypeObject py_type_object_header = { PyObject_HEAD_INIT(NULL) 0};
  int size = (long*)&(py_type_object_header.tp_name) - (long*)&py_type_object_header;

  memset(type_object, 0, sizeof(PyTypeObject));
  memcpy(type_object, &py_type_object_header, size);
}

shared_ptr<Go::Point> PyObject_AsGoPoint(PyObject* obj)
{
  if (!obj)
    return shared_ptr<Go::Point>();
  // Instance of a point
  if (PyObject_TypeCheck(obj,&Point_Type))
    return ((Point*)obj)->data;

  shared_ptr<Go::Point> result;

  // A python list
  if (PyObject_TypeCheck(obj,&PyList_Type)) {
    std::vector<double> vec;
    for (int i=0;i<PyList_Size(obj);++i) {
      PyObject* o = PyList_GetItem(obj,i);
      if (o)
        vec.push_back(PyFloat_AsDouble(o));
    }
    result.reset(new Go::Point(&vec[0],&vec[0]+vec.size(),true));
  }

  // A python tuple
  if (PyObject_TypeCheck(obj,&PyTuple_Type)) {
    std::vector<double> vec;
    int size = PyTuple_Size(obj);
    vec.resize(size);
    if (PyArg_ParseTuple(obj,"dd|d",&vec[0],&vec[1],&vec[2]))
      result.reset(new Go::Point(&vec[0],&vec[0]+vec.size(),true));
  }

  return result;
}

shared_ptr<Go::ParamCurve> PyObject_AsGoCurve(PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&Curve_Type))
    return ((Curve*)obj)->data;

  return shared_ptr<Go::ParamCurve>();
}

shared_ptr<Go::ParamSurface> PyObject_AsGoSurface(PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&Surface_Type))
    return ((Surface*)obj)->data;

  return shared_ptr<Go::ParamSurface>();
}

shared_ptr<Go::ParamVolume> PyObject_AsGoVolume(PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&Volume_Type))
    return ((Volume*)obj)->data;

  return shared_ptr<Go::ParamVolume>();
}

shared_ptr<Go::SurfaceModel> PyObject_AsGoSurfaceModel(PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&SurfaceModel_Type))
    return ((SurfaceModel*)obj)->data;

  return shared_ptr<Go::SurfaceModel>();
}

shared_ptr<Go::VolumeModel> PyObject_AsGoVolumeModel(PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&VolumeModel_Type))
    return ((VolumeModel*)obj)->data;

  return shared_ptr<Go::VolumeModel>();
}

void PyMethods_Append(std::vector<PyMethodDef>& defs, PyMethodDef* start)
{
  while (start && start->ml_name) {
    defs.push_back(*start);
    start++;
  }
}

void VectorToPyPointList(PyObject* result, const std::vector<double>& data, int dim)
{
  for (size_t i=0;i<data.size()/dim;++i) {
    if (dim == 1)
      PyList_SetItem(result, i, PyFloat_FromDouble(data[i]));
    else {
      Point* pt = (Point*)Point_Type.tp_alloc(&Point_Type, 0);
      pt->data.reset(new Go::Point(data.begin()+i*dim, data.begin()+(i+1)*dim));
      PyList_SetItem(result, i, (PyObject*)pt);
    }
  }
}

}
