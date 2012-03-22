#include "volume.h"

#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/trivariate/ElementaryVolume.h"
#include "GoTools/geometry/ParamSurface.h"

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

PyDoc_STRVAR(volume_clone__doc__,"Clone a volume\n"
                                 "@return: New copy of volume\n");
PyObject* Volume_Clone(PyObject* self, PyObject* args, PyObject* kwds)
{
  Volume* res = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::ParamVolume> vol = PyObject_AsGoVolume(self);
  res->data.reset(vol->clone());
 
  return (PyObject*)res;
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

PyDoc_STRVAR(volume_get_edges__doc__,"Return the six edge volumes\n"
                                      "@return: A list of the six edges");
PyObject* Volume_GetEdges(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;
  if (!volume->isSpline()) {
    Volume* vol = (Volume*)self;
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }

  shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume);
  std::vector<shared_ptr<Go::ParamSurface> > edges = vol->getAllBoundarySurfaces();

  Surface* surface0 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface1 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface2 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface3 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface4 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface5 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  surface0->data = shared_ptr<Go::ParamSurface>(edges.at(0));
  surface1->data = shared_ptr<Go::ParamSurface>(edges.at(1));
  surface2->data = shared_ptr<Go::ParamSurface>(edges.at(2));
  surface3->data = shared_ptr<Go::ParamSurface>(edges.at(3));
  surface4->data = shared_ptr<Go::ParamSurface>(edges.at(4));
  surface5->data = shared_ptr<Go::ParamSurface>(edges.at(5));

  PyObject* result = PyList_New(0);
  PyList_Append(result, (PyObject*) surface0);
  PyList_Append(result, (PyObject*) surface1);
  PyList_Append(result, (PyObject*) surface2);
  PyList_Append(result, (PyObject*) surface3);
  PyList_Append(result, (PyObject*) surface4);
  PyList_Append(result, (PyObject*) surface5);

  return result;
}

PyMethodDef Volume_methods[] = {
     {(char*)"Clone",     (PyCFunction)Volume_Clone,    METH_VARARGS|METH_KEYWORDS, volume_clone__doc__},
     {(char*)"GetEdges",  (PyCFunction)Volume_GetEdges, METH_VARARGS|METH_KEYWORDS, volume_get_edges__doc__},
     {NULL,             NULL,                       0,                          NULL}
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
