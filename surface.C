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

PyMethodDef Surface_methods[] = {
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

