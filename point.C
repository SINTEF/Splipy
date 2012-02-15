#include "point.h"
#include "pyutils.h"
#include "geomodeller.h"

#include <sstream>

extern "C"
{
PyTypeObject Point_Type;

PyObject* Point_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Point* self;
  self = (Point*)type->tp_alloc(type,0);
  static const char* keyWords[] = {"x", "y", "z", "list", NULL };
  double x=0,y=0,z=0;
  PyObject* list=NULL;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|dddO",
                                   (char**)keyWords,&x,&y,&z,&list))
    return NULL;

  self->data = PyObject_AsGoPoint(list);
  if (!self->data) {
    if (modState.dim == 2)
      self->data.reset(new Go::Point(x,y));
    else
      self->data.reset(new Go::Point(x,y,z));
  }
  return (PyObject*)self;
}

void Point_Dealloc(Point* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* Point_Str(Point* self)
{
  std::stringstream str;
  str << self->data->size() << "D point: " << (*self->data)[0] << " " << (*self->data)[1];
  if (self->data->size() == 3)
    str << " " << (*self->data)[2];
  return PyString_FromString(str.str().c_str());
}

PyObject* Point_Add(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  Point* result = NULL;
  shared_ptr<Go::Point> p2 = PyObject_AsGoPoint(o2);
  if (p1 && p2) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point(*p1+*p2));
  }

  return (PyObject*)result;
}

PyObject* Point_Sub(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  Point* result = NULL;
  shared_ptr<Go::Point> p2 = PyObject_AsGoPoint(o2);
  if (p1 && p2) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point(*p1-*p2));
  }

  return (PyObject*)result;
}

PyObject* Point_Mul(PyObject* o1, PyObject* o2)
{
  PyObject* o1r;
  PyObject* o2r;
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  if (p1) {
    o1r = o1;
    o2r = o2;
  } else {
    p1 = PyObject_AsGoPoint(o2);
    o1r = o2;
    o2r = o1;
  }

  Point* result = NULL;
  if (p1 && PyObject_TypeCheck(o2r,&PyFloat_Type)) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point(*p1*PyFloat_AsDouble(o2r)));
  } else if (PyObject_TypeCheck(o2r,&Point_Type)) {
    shared_ptr<Go::Point> p2 = PyObject_AsGoPoint(o2r);
    double res = *p1*(*p2);
    result = (Point*)Py_BuildValue((char*)"d",res);
  }

  return (PyObject*)result;
}

PyObject* Point_Div(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  Point* result = NULL;
  if (PyObject_TypeCheck(o2,&PyFloat_Type)) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point(*p1/PyFloat_AsDouble(o2)));
  } 

  return (PyObject*)result;
}

PyObject* Point_Neg(PyObject* o1)
{
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(-(*p1)));

  return (PyObject*)result;
}

// the % operator
PyObject* Point_Mod(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  shared_ptr<Go::Point> p2 = PyObject_AsGoPoint(o2);
  Point* result = NULL;
  if (o1 && o2) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point((*p1) % (*p2)));
  }

  return (PyObject*)result;
}

PyObject* Point_Abs(PyObject* o1)
{
  shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o1);
  Point* result = NULL;
  if (p1)
    result = (Point*)Py_BuildValue((char*)"d",p1->length());

  return (PyObject*)result;
}

PyMethodDef Point_methods[] = {
     {NULL,           NULL,                     0,            NULL}
   };

PyNumberMethods Point_operators = {0};

PyDoc_STRVAR(point__doc__, "A point");
void init_Point_Type()
{
  InitializeTypeObject(&Point_Type);
  Point_operators.nb_add = Point_Add;
  Point_operators.nb_subtract = Point_Sub;
  Point_operators.nb_multiply = Point_Mul;
  Point_operators.nb_divide = Point_Div;
  Point_operators.nb_remainder = Point_Mod;
  Point_operators.nb_negative = Point_Neg;
  Point_operators.nb_absolute = Point_Abs;
  Point_Type.tp_name = "GoTools.Point";
  Point_Type.tp_basicsize = sizeof(Point);
  Point_Type.tp_dealloc = (destructor)Point_Dealloc;
  Point_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES;
  Point_Type.tp_doc = point__doc__;
  Point_Type.tp_methods = Point_methods;
  Point_Type.tp_as_number = &Point_operators;
  Point_Type.tp_base = 0;
  Point_Type.tp_new = Point_New;
  Point_Type.tp_str = (reprfunc)Point_Str;
  PyType_Ready(&Point_Type);
}

}
