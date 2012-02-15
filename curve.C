#include "curve.h"

#include "geomodeller.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"

#include <fstream>
#include <sstream>

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

PyDoc_STRVAR(curve_get_knots__doc__,"Get the unique knots of a spline curve");
PyObject* Curve_GetKnots(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  if (!curve)
    return NULL;
  Curve* crv = (Curve*)self;
  crv->data = convertSplineCurve(curve);
  PyObject* result = PyList_New(0);
  std::vector<double> knots;
  static_pointer_cast<Go::SplineCurve>(crv->data)->basis().knotsSimple(knots);
  for (std::vector<double>::iterator it  = knots.begin();
                                     it != knots.end();++it) {
    PyList_Append(result,Py_BuildValue((char*)"d",*it));
  }
                
  return result;
}

PyDoc_STRVAR(curve_insert_knot__doc__,"Insert a knot into a spline curve");
PyObject* Curve_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"knot", NULL };
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  double knot;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"d",
                                   (char**)keyWords,&knot) || !curve)
    return NULL;

  Curve* crv = (Curve*)self;
  crv->data = convertSplineCurve(curve);

  static_pointer_cast<Go::SplineCurve>(crv->data)->insertKnot(knot);

  Py_INCREF(Py_None);
  return Py_None;
}

PyMethodDef Curve_methods[] = {
     {(char*)"GetKnots",   (PyCFunction)Curve_GetKnots,   METH_VARARGS,               curve_get_knots__doc__},
     {(char*)"InsertKnot", (PyCFunction)Curve_InsertKnot, METH_VARARGS|METH_KEYWORDS, curve_insert_knot__doc__},
     {NULL,                NULL,                          0,                          NULL}
   };

PyDoc_STRVAR(curve__doc__, "A parametric description of a curve");
void init_Curve_Type()
{
  InitializeTypeObject(&Curve_Type);
  Curve_Type.tp_name = "GoTools.Curve";
  Curve_Type.tp_basicsize = sizeof(Curve);
  Curve_Type.tp_dealloc = (destructor)Curve_Dealloc;
  Curve_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES;
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

void WriteCurveG2(std::ofstream& g2_file, Curve* curve, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineCurve> crv = convertSplineCurve(curve->data);
    crv->writeStandardHeader(g2_file);
    crv->write(g2_file);
  } else {
    curve->data->writeStandardHeader(g2_file);
    curve->data->write(g2_file);
  }
}
