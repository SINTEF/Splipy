#include "point.h"
#include "pyutils.h"
#include "geomodeller.h"
#include "GoTools/geometry/GeometryTools.h"

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

PyDoc_STRVAR(point_normalize__doc__,"Normalize a point (vector) to have eucledian length one\n"
                                    "@return: None");
PyObject* Point_Normalize(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(self);
  if (!point)
    return NULL;

  point->normalize();

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(point_clone__doc__,"Clone a point\n"
                                "@return: New copy of point\n"
                                "@rtype: Point\n");
PyObject* Point_Clone(PyObject* self, PyObject* args, PyObject* kwds)
{
  Point* res = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  shared_ptr<Go::Point> pt = PyObject_AsGoPoint(self);
  res->data.reset(new Go::Point(*pt));
 
  return (PyObject*)res;
}

PyDoc_STRVAR(point_rotate__doc__,"Rotate a point around an axis\n"
                                 "@param axis: The axis to rotate around\n"
                                 "@type axis: Point, list of floats or tuple of floats\n"
                                 "@param angle: Angle to rotate point with in radians\n"
                                 "@type angle: float\n"
                                 "@return: None");
PyObject* Point_Rotate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", "angle", NULL };
  PyObject* axiso;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&axiso,&angle))
    return NULL;

  shared_ptr<Go::Point> point = PyObject_AsGoPoint(self);
  shared_ptr<Go::Point> axis  = PyObject_AsGoPoint(axiso);
  if (!point || !axis)
    return NULL;

  // special case 4D points as rational control points with weights
  if(point->dimension() == 4) { 
    double w = (*point)[3];
    (*point)[0] /= w;
    (*point)[1] /= w;
    (*point)[2] /= w;
  } else if(point->dimension() != 3)
    return NULL;

  Go::GeometryTools::rotatePoint(*axis, angle, point->begin());

  if(point->dimension() == 4) { 
    // GeometryTools::rotatePoint only uses the three first components of double* argument
    // now projecting control points back into 4D space
    double w = (*point)[3];
    (*point)[0] *= w;
    (*point)[1] *= w;
    (*point)[2] *= w;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

PyObject* Point_Str(Point* self)
{
  std::stringstream str;
  int dim = self->data->size();
  str << dim << "D point:" ;
  for(int i=0; i<dim; ++i) 
    str << " " << (*self->data)[i];
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
  else if (PyObject_TypeCheck(o2r,&PyInt_Type)) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point(*p1*PyInt_AsLong(o2r)));
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
  } else if (PyObject_TypeCheck(o2r,&PyInt_Type)) {
    result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
    result->data.reset(new Go::Point(*p1/PyInt_AsLong(o2)));
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

Py_ssize_t Point_NmbComponent(PyObject* self)
{
  shared_ptr<Go::Point> sm = PyObject_AsGoPoint(self);
  if (!sm)
    return 0;

  return sm->size();
}

PyObject* Point_GetComponent(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::Point> sm = PyObject_AsGoPoint(self);
  if (!sm || i < 0 || i >= sm->size())
    return NULL;

  return PyFloat_FromDouble((*sm)[i]);
}

PyMethodDef Point_methods[] = {
     {(char*)"Clone",     (PyCFunction)Point_Clone,     METH_VARARGS,               point_clone__doc__},
     {(char*)"Normalize", (PyCFunction)Point_Normalize, METH_VARARGS,               point_rotate__doc__},
     {(char*)"Rotate",    (PyCFunction)Point_Rotate,    METH_VARARGS|METH_KEYWORDS, point_normalize__doc__},
     {NULL,               NULL,                         0,                          NULL}
   };

PyNumberMethods Point_operators = {0};
PySequenceMethods Point_seq_operators = {0};

PyDoc_STRVAR(point__doc__, "A point");
void init_Point_Type()
{
  Point_seq_operators.sq_item = Point_GetComponent;
  Point_seq_operators.sq_length = Point_NmbComponent;
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
  Point_Type.tp_as_sequence = &Point_seq_operators;
  Point_Type.tp_base = 0;
  Point_Type.tp_new = Point_New;
  Point_Type.tp_str = (reprfunc)Point_Str;
  PyType_Ready(&Point_Type);
}

}
