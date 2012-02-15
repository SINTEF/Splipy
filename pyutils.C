#include "pyutils.h"

#ifdef HAS_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL GEOMOD_ARRAY_API
#define NO_IMPORT_ARRAY
#include <arrayobject.h>
#endif

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

  // A numpy array
#ifdef HAS_NUMPY
//  if (!result) {
//    PyArrayObject* vector = (PyArrayObject*)PyArray_FromAny(obj,NULL,0,0,0,NULL);
//    if (vector) {
//      int siz = vector->dimensions[0];
//      result.reset(new Go::Point((double*)vector->data,
//                                 ((double*)vector->data)+siz,true));
//    }
//  }
#endif
  
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
