#include "curve.h"

#include "geomodeller.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"

#include <sstream>

#ifdef HAS_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL GEOMOD_ARRAY_API
#define NO_IMPORT_ARRAY
#include <arrayobject.h>
#endif

extern "C"
{
PyTypeObject Curve_Type;

PyObject* Curve_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Curve* self;
  self = (Curve*)type->tp_alloc(type,0);

  // only generators allocate the paramcurve!
  return (PyObject*)self;
}

void Curve_Dealloc(Curve* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* Curve_Str(Curve* self)
{
  std::stringstream str;
  if (self->data) {
    if (self->data->instanceType() == Go::Class_Line)
      str << "Line:" << std::endl;
    if (self->data->instanceType() == Go::Class_Circle)
      str << "Circle:" << std::endl;
    if (self->data->instanceType() == Go::Class_Ellipse)
      str << "Ellipse:" << std::endl;
    if (self->data->instanceType() == Go::Class_SplineCurve)
      str << "Spline curve:" << std::endl;
    str << *self->data;
  } else
    str << "(empty)";
  return PyString_FromString(str.str().c_str());
}

PyMethodDef Curve_methods[] = {
     {NULL,           NULL,                     0,            NULL}
   };

PyDoc_STRVAR(curve__doc__, "A parametric description of a curve");
void init_Curve_Type()
{
  InitializeTypeObject(&Curve_Type);
  Curve_Type.tp_name = "GoTools.Curve";
  Curve_Type.tp_basicsize = sizeof(Curve);
  Curve_Type.tp_dealloc = (destructor)Curve_Dealloc;
  Curve_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  Curve_Type.tp_doc = curve__doc__;
  Curve_Type.tp_methods = Curve_methods;
  Curve_Type.tp_base = 0;
  Curve_Type.tp_new = Curve_New;
  Curve_Type.tp_str = (reprfunc)Curve_Str;
  PyType_Ready(&Curve_Type);
}
}

shared_ptr<Go::SplineCurve> convertSplineCurve(shared_ptr<Go::ParamCurve> curve)
{
  if (curve->instanceType() == Go::Class_SplineCurve)
    return dynamic_pointer_cast<Go::SplineCurve, Go::ParamCurve>(curve);
  return shared_ptr<Go::SplineCurve>(curve->geometryCurve());
}
