#include "volume.h"

#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/trivariate/ElementaryVolume.h"

#include <fstream>
#include <sstream>

extern "C"
{
PyTypeObject Volume_Type;

PyObject* Volume_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Volume* self;
  self = (Volume*)type->tp_alloc(type,0);

  // only generators allocate the paramvolume!
  return (PyObject*)self;
}

void Volume_Dealloc(Volume* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* Volume_Str(Volume* self)
{
  std::stringstream str;
  if (self->data) {
    str << *self->data;
  } else
    str << "(empty)";
  return PyString_FromString(str.str().c_str());
}

PyObject* Volume_Add(PyObject* o1, PyObject* o2)
{
  Volume* vol = (Volume*)o1;
  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::Point> p = PyObject_AsGoPoint(o2);
  if (p) {
    result->data.reset(vol->data->clone());
    result->data->translate(*p);
  }

  return (PyObject*)result;
}

PyObject* Volume_Sub(PyObject* o1, PyObject* o2)
{
  Volume* vol = (Volume*)o1;
  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::Point> p = PyObject_AsGoPoint(o2);
  if (p) {
    result->data.reset(vol->data->clone());
    result->data->translate(-(*p));
  }

  return (PyObject*)result;
}

PyMethodDef Volume_methods[] = {
     {NULL,           NULL,                     0,            NULL}
   };

PyNumberMethods Volume_operators = {0};

PyDoc_STRVAR(volume__doc__, "A parametric description of a volume");
void init_Volume_Type()
{
  InitializeTypeObject(&Volume_Type);
  Volume_operators.nb_add      = Volume_Add;
  Volume_operators.nb_subtract = Volume_Sub;
  Volume_Type.tp_name = "GoTools.Volume";
  Volume_Type.tp_basicsize = sizeof(Volume);
  Volume_Type.tp_dealloc = (destructor)Volume_Dealloc;
  Volume_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES;
  Volume_Type.tp_doc = volume__doc__;
  Volume_Type.tp_methods = Volume_methods;
  Volume_Type.tp_as_number = &Volume_operators;
  Volume_Type.tp_base = 0;
  Volume_Type.tp_new = Volume_New;
  Volume_Type.tp_str = (reprfunc)Volume_Str;
  PyType_Ready(&Volume_Type);
}

}

shared_ptr<Go::SplineVolume> convertSplineVolume(shared_ptr<Go::ParamVolume> volume)
{
  if (volume->instanceType() == Go::Class_SplineVolume)
    return dynamic_pointer_cast<Go::SplineVolume, Go::ParamVolume>(volume);
  shared_ptr<Go::ElementaryVolume> e_volume = dynamic_pointer_cast<Go::ElementaryVolume, Go::ParamVolume>(volume);
  return shared_ptr<Go::SplineVolume>(e_volume->geometryVolume());
}

void WriteVolumeG2(std::ofstream& g2_file, Volume* volume, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume->data);
    if (!vol->isLeftHanded())
    {
      vol = shared_ptr<Go::SplineVolume>(vol->clone());
      vol->reverseParameterDirection(2);
    }
    vol->writeStandardHeader(g2_file);
    vol->write(g2_file);
  } else {
    volume->data->writeStandardHeader(g2_file);
    volume->data->write(g2_file);
  }
}
