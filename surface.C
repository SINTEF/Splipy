#include "surface.h"
#include "pyutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/SplineSurface.h"

#include <sstream>

#ifdef HAS_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL GEOMOD_ARRAY_API
#define NO_IMPORT_ARRAY
#include <arrayobject.h>
#endif

extern "C"
{
PyTypeObject Surface_Type;

PyObject* Surface_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Surface* self;
  self = (Surface*)type->tp_alloc(type,0);

  // only generators allocate the paramsurface!
  return (PyObject*)self;
}

void Surface_Dealloc(Surface* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* Surface_Str(Surface* self)
{
  std::stringstream str;
  printSurfaceToStream(str,self->data);
  return PyString_FromString(str.str().c_str());
}

PyDoc_STRVAR(surface_raise_order__doc__,"Raise order of a spline surface");
PyObject* Surface_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"raise_u", "raise_v", NULL };
  int raise_u=0, raise_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords,&raise_u,&raise_v))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }
   static_pointer_cast<Go::SplineSurface>(surface)->raiseOrder(raise_u,raise_v);

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(surface_get_knots__doc__,"Return knots for a spline surface");
PyObject* Surface_GetKnots(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
  if (!surface->isSpline()) {
    Surface* surf = (Surface*)self;
    surf->data = convertSplineSurface(surface);
    surface = surf->data;
  }

  PyObject* result = PyTuple_New(2);
  for (int i=0;i<2;++i) {
    PyObject* list = PyList_New(0);
    std::vector<double> knots;
    static_pointer_cast<Go::SplineSurface>(surface)->basis(i).knotsSimple(knots);
    for (std::vector<double>::iterator it  = knots.begin();
                                       it != knots.end();++it) {
      PyList_Append(list,Py_BuildValue((char*)"d",*it));
    }
    PyTuple_SetItem(result,i,list);
  }

  return result;
}

PyDoc_STRVAR(surface_insert_knot__doc__,"Insert a knot in a parameter direction of spline surface");
PyObject* Surface_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", "knot", NULL };
  int direction=0;
  double knot;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"id",
                                   (char**)keyWords,&direction,&knot))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (!surface)
    return NULL;
   if (!surface->isSpline()) {
     Surface* surf = (Surface*)self;
     surf->data = convertSplineSurface(surface);
     surface = surf->data;
   }
   if (direction == 0)
     static_pointer_cast<Go::SplineSurface>(surface)->insertKnot_u(knot);
   else
     static_pointer_cast<Go::SplineSurface>(surface)->insertKnot_v(knot);

   Py_INCREF(Py_None);
   return Py_None;
}

PyMethodDef Surface_methods[] = {
     {(char*)"GetKnots",   (PyCFunction)Surface_GetKnots, METH_VARARGS|METH_KEYWORDS,surface_get_knots__doc__},
     {(char*)"InsertKnot", (PyCFunction)Surface_InsertKnot, METH_VARARGS|METH_KEYWORDS,surface_insert_knot__doc__},
     {(char*)"RaiseOrder", (PyCFunction)Surface_RaiseOrder, METH_VARARGS|METH_KEYWORDS,surface_raise_order__doc__},
     {NULL,           NULL,                     0,            NULL}
   };

PyDoc_STRVAR(surface__doc__, "A parametric description of a surface");
void init_Surface_Type()
{
  InitializeTypeObject(&Surface_Type);
  Surface_Type.tp_name = "GoTools.Surface";
  Surface_Type.tp_basicsize = sizeof(Surface);
  Surface_Type.tp_dealloc = (destructor)Surface_Dealloc;
  Surface_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  Surface_Type.tp_doc = surface__doc__;
  Surface_Type.tp_methods = Surface_methods;
  Surface_Type.tp_base = 0;
  Surface_Type.tp_new = Surface_New;
  Surface_Type.tp_str = (reprfunc)Surface_Str;
  PyType_Ready(&Surface_Type);
}

}

shared_ptr<Go::SplineSurface> convertSplineSurface(shared_ptr<Go::ParamSurface> surface)
{
  if (surface->instanceType() == Go::Class_SplineSurface)
    return dynamic_pointer_cast<Go::SplineSurface, Go::ParamSurface>(surface);
  shared_ptr<Go::ElementarySurface> e_surface = dynamic_pointer_cast<Go::ElementarySurface, Go::ParamSurface>(surface);
  return shared_ptr<Go::SplineSurface>(e_surface->geometrySurface());
}

void printSurfaceToStream(std::ostream& str, shared_ptr<Go::ParamSurface> surf)
{
  if (surf) {
    if (surf->instanceType() == Go::Class_Plane)
      str << "Plane:" << std::endl;
    if (surf->instanceType() == Go::Class_Sphere)
      str << "Sphere surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Cylinder)
      str << "Cylinder surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Cone)
      str << "Cone surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Torus)
      str << "Torus surface:" << std::endl;
    if (surf->instanceType() == Go::Class_Disc)
      str << "Circular disc:" << std::endl;
    str << *surf;
  } else
    str << "(empty)";
}

