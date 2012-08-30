#include "curve.h"

#include "geomodeller.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/GeometryTools.h"

#include <fstream>
#include <sstream>

extern "C"
{
PyTypeObject Curve_Type;

PyObject* Curve_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Curve* self;
  self = (Curve*)type->tp_alloc(type,0);
  static const char* keyWords[] = {"order", "knots", "coefs", "rational", NULL };
  PyObject* knotso=0;
  PyObject* coefso=0;
  bool rational=false;
  int order=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|iOOb",
                                   (char**)keyWords, &order, &knotso, &coefso, &rational))
    return NULL;

  // curve is specified
  if (knotso && coefso) {
    std::vector<double> knots;
    std::vector<double> coefs;
    if (PyObject_TypeCheck(knotso,&PyList_Type)) {
      for (int i=0;i<PyList_Size(knotso);++i) {
        PyObject* o = PyList_GetItem(knotso,i);
        if (o && PyObject_TypeCheck(o,&PyFloat_Type))
          knots.push_back(PyFloat_AsDouble(o));
        else if (o && PyObject_TypeCheck(o,&PyInt_Type))
          knots.push_back(PyInt_AsLong(o));
      }
    }
    if (PyObject_TypeCheck(coefso,&PyList_Type)) {
      for (int i=0;i<PyList_Size(coefso);++i) {
        PyObject* o = PyList_GetItem(coefso,i);
        shared_ptr<Go::Point> p1 = PyObject_AsGoPoint(o); 
        if (p1) {
          for (double* p = p1->begin(); p != p1->end(); ++p)
            coefs.push_back(*p);
        } else if (o && PyObject_TypeCheck(o,&PyFloat_Type))
          coefs.push_back(PyFloat_AsDouble(o));
        else if (o && PyObject_TypeCheck(o,&PyInt_Type))
          coefs.push_back(PyInt_AsLong(o));
      }
    }
    ((Curve*)self)->data.reset(new Go::SplineCurve(knots.size()-order, order, knots.begin(),
                                                   coefs.begin(), modState.dim, rational));
  }

  return (PyObject*)self;
}

void Curve_Dealloc(Curve* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyDoc_STRVAR(curve_append_curve__doc__,"Merge another curve with this one, with possible reparametrization\n"
                                       "@param curve: The curve to append\n"
                                       "@type curve: Curve\n"
                                       "@param continuity: (optional) Required continuity\n"
                                       "@type continuity: int\n"
                                       "@param reparam: (optional) Specify whether or not there should be reparametrization\n"
                                       "@type reparam: bool\n"
                                       "@return: None");
PyObject* Curve_AppendCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curve", "continuity", "reparam", NULL };
  int continuity = 0;
  bool reparam   = true;
  PyObject *oCrv;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|ib",
                                   (char**)keyWords, &oCrv, &continuity, &reparam))
    return NULL;

  shared_ptr<Go::ParamCurve> crv2 = PyObject_AsGoCurve(self);
  shared_ptr<Go::ParamCurve> crv;
  if (crv2->instanceType() != Go::Class_SplineCurve) {
    std::cerr << "Converting param curve to a spline curve (append not implemented for param curves)" << std::endl;
    crv = ((Curve*)self)->data = convertSplineCurve(crv2);
  }
  else
    crv = dynamic_pointer_cast<Go::SplineCurve,Go::ParamCurve>(crv);

  shared_ptr<Go::SplineCurve> otherCrv = convertSplineCurve(PyObject_AsGoCurve(oCrv));
  if (!crv || !otherCrv)
    return NULL;

  double dist;
  crv->appendCurve(otherCrv.get(), continuity, dist, reparam);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(curve_clone__doc__,"Clone a curve\n"
                                "@return: New copy of curve\n"
                                "@rtype: Curve\n");
PyObject* Curve_Clone(PyObject* self, PyObject* args, PyObject* kwds)
{
  Curve* res = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  shared_ptr<Go::ParamCurve> crv = PyObject_AsGoCurve(self);
  res->data.reset(crv->clone());
 
  return (PyObject*)res;
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

PyDoc_STRVAR(curve_evaluate__doc__,"Evaluate curve at a parameter value\n"
                                   "@param value: The parameter value\n"
                                   "@type value: float\n"
                                   "@return: The value of the curve\n"
                                   "@rtype: Point");
PyObject* Curve_Evaluate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"value", NULL };
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  double value=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"d",
                                   (char**)keyWords,&value) || !curve)
    return NULL;

  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(curve->dimension()));
  curve->point(*result->data,value);

  return (PyObject*)result;
}

PyDoc_STRVAR(curve_flip_parametrization__doc__,"Flip curve parametrization\n"
                                               "@return: None");
PyObject* Curve_FlipParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  if (!curve)
    return NULL;

  curve->reverseParameterDirection();

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(curve_get_knots__doc__,"Get the unique knots of a spline curve\n"
                                    "@return: List with the knot values\n"
                                    "@rtype: List of float");
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

PyDoc_STRVAR(curve_get_order__doc__,"Get the curve order (polynomial degree + 1)\n"
                                    "@return: Order\n"
                                    "@rtype: int");
PyObject* Curve_GetOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  if (!curve)
    return NULL;

  shared_ptr<Go::SplineCurve> sc = convertSplineCurve(curve);
  if(!sc)
    return NULL;

  return Py_BuildValue((char*)"i",sc->order());
}

PyDoc_STRVAR(curve_insert_knot__doc__,"Insert a knot into a spline curve\n"
                                      "@param knot: The knot to insert\n"
                                      "@type knot: float\n"
                                      "@return: None");
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

PyDoc_STRVAR(curve_normalize__doc__,"Normalize a curve in the parameter domain\n"
                                    "@return: None");
PyObject* Curve_Normalize(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  if (!curve)
    return NULL;

  Curve* crv = (Curve*)self;
  crv->data = convertSplineCurve(curve);

  crv->data->setParameterInterval(0,1);

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(curve_project__doc__,"Project the curve onto an axis or plane along parallel to the cartesian coordinate system\n"
                                  "@param axis: The axis or plane to project onto (\"X\",\"Y\",\"Z\" or a comibation of these)\n"
                                  "@type axis: string\n"
                                  "@return: None");
PyObject* Curve_Project(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", NULL };
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  char *sAxis;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&sAxis) || !curve)
    return NULL;

  bool bAxis[3];
  bAxis[0] = true;
  bAxis[1] = true;
  bAxis[2] = true;
  while(*sAxis != 0) {
    if(*sAxis == 'x' || *sAxis == 'X')
      bAxis[0] = false;
    else if(*sAxis == 'y' || *sAxis == 'Y')
      bAxis[1] = false;
    else if(*sAxis == 'z' || *sAxis == 'Z')
      bAxis[2] = false;
    sAxis++;
  }

  Curve* crv = (Curve*)self;
  crv->data = convertSplineCurve(curve);

  shared_ptr<Go::SplineCurve>   scurve   = static_pointer_cast<Go::SplineCurve>(crv->data);
  bool                          rational = scurve->rational();
  int                           dim      = scurve->dimension();
  std::vector<double>::iterator coefs    = (rational) ? scurve->rcoefs_begin() : scurve->coefs_begin();
  std::vector<double>::iterator coefsEnd = (rational) ? scurve->rcoefs_end()   : scurve->coefs_end();
  while(coefs != coefsEnd) {
    if(bAxis[0] && dim>0)
      coefs[0] = 0.0;
    if(bAxis[1] && dim>1)
      coefs[1] = 0.0;
    if(bAxis[2] && dim>2)
      coefs[2] = 0.0;
    coefs += (dim+rational);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(curve_raise_order__doc__,"Raise the order of the curve's B-spline basis without changing the shape of the curve\n"
                                      "@param n: Specifies how many times the order will be raised\n"
                                      "@type knot: int\n"
                                      "@return: None");
PyObject* Curve_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"n", NULL };
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  int amount;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&amount) || !curve)
    return NULL;

  Curve* crv = (Curve*)self;
  crv->data = convertSplineCurve(curve);

  static_pointer_cast<Go::SplineCurve>(crv->data)->raiseOrder(amount);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(curve_rotate__doc__,"Rotate a curve around an axis\n"
                                 "@param axis: The axis to rotate around\n"
                                 "@type axis: Point, list of floats or tuple of floats\n"
                                 "@param angle: Angle to rotate curve with in radians\n"
                                 "@type angle: float\n"
                                 "@return: None");
PyObject* Curve_Rotate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", "angle", NULL };
  PyObject* axiso;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamCurve> crv = PyObject_AsGoCurve(self);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  if (!crv || !axis)
    return NULL;

   Curve* curve = (Curve*)self;
   curve->data = convertSplineCurve(crv);
   crv = curve->data;

   Go::GeometryTools::rotateSplineCurve(*axis, angle,
                        *static_pointer_cast<Go::SplineCurve>(crv));

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(curve_split__doc__, "Split the curve into segments\n"
                                 "@param params: The parameter values to split at\n"
                                 "@type params: Float or list of floats\n"
                                 "@return: The resulting curves\n"
                                 "@rtype: List of Curve");
PyObject* Curve_Split(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"params", NULL };
  PyObject* params;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&params) || !params)
    return NULL;

  std::vector<double> p;
  if (PyObject_TypeCheck(params,&PyList_Type)) {
    for (int i=0;i<PyList_Size(params);++i) {
      PyObject* o = PyList_GetItem(params,i);
      if (o)
        p.push_back(PyFloat_AsDouble(o));
    }
  } else if (PyObject_TypeCheck(params,&PyFloat_Type) || 
             PyObject_TypeCheck(params,&PyInt_Type))
    p.push_back(PyFloat_AsDouble(params));

  if (p.empty())
    return NULL;

  std::vector<shared_ptr<Go::ParamCurve> > curves;
  if (p.size() > 1) {
    shared_ptr<Go::SplineCurve> crv = convertSplineCurve(((Curve*)self)->data);
    std::vector<shared_ptr<Go::SplineCurve> > curves2 = crv->split(p);
    for (size_t i=0;i<curves2.size();++i)
      curves.push_back(static_pointer_cast<Go::ParamCurve,Go::SplineCurve>(curves2[i]));
  } else
    curves = ((Curve*)self)->data->split(p[0]);

  PyObject* result = PyList_New(0);
  for (size_t i=0;i<curves.size();++i) {
    Curve* crv = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
    crv->data = curves[i];
    PyList_Append(result,(PyObject*)crv);
  }

  return result;
}

PyDoc_STRVAR(curve_translate__doc__,"Translate a curve along a given vector\n"
                                    "@param vector: The vector to translate along\n"
                                    "@type axis: Point, list of floats or tuple of floats\n"
                                    "@return: None");
PyObject* Curve_Translate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"vector", NULL };
  PyObject* veco;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&veco))
    return NULL;

  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(self);
  shared_ptr<Go::Point>      vec   = PyObject_AsGoPoint(veco);
  if (!curve || !vec)
    return NULL;

   if (curve->geometryCurve() != NULL) {
     Curve* crv = (Curve*)self;
     crv->data = convertSplineCurve(curve);
     curve = crv->data;
   }

   Go::GeometryTools::translateSplineCurve(*vec, *static_pointer_cast<Go::SplineCurve>(curve));

   Py_INCREF(Py_None);
   return Py_None;
}

PyObject* Curve_Add(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::ParamCurve> curve = PyObject_AsGoCurve(o1);
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(o2);
  Curve* result = NULL;
  if (curve && point) {
    result = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
    if (curve->instanceType() == Go::Class_SplineCurve)
      result->data = shared_ptr<Go::ParamCurve>(
          new Go::SplineCurve(*static_pointer_cast<Go::SplineCurve>(curve)));
    else
      result->data = convertSplineCurve(curve);

    Go::GeometryTools::translateSplineCurve(*point, 
                         *static_pointer_cast<Go::SplineCurve>(result->data));
  }

  return (PyObject*) result;
}

Py_ssize_t Curve_NmbComponent(PyObject* self)
{
  shared_ptr<Go::ParamCurve> pc = PyObject_AsGoCurve(self);
  if (!pc)
    return 0;

  shared_ptr<Go::SplineCurve> sc = convertSplineCurve(pc);
  if(!sc)
    return 0;

  return sc->numCoefs();
}

PyObject* Curve_GetComponent(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::ParamCurve> pc = PyObject_AsGoCurve(self);
  if (!pc)
    return NULL;

  if(pc->dimension() != 3) 
    return NULL;

  shared_ptr<Go::SplineCurve> sc = convertSplineCurve(pc);
  if(!sc)
    return NULL;
  
  if(i < 0 || i >= sc->numCoefs())
    return NULL;

  double x,y,z,w;
  if(sc->rational()) {
    w = *(sc->rcoefs_begin() + i*(sc->dimension()+1)+3)    ;
    x = *(sc->rcoefs_begin() + i*(sc->dimension()+1)+0) / w;
    y = *(sc->rcoefs_begin() + i*(sc->dimension()+1)+1) / w;
    z = *(sc->rcoefs_begin() + i*(sc->dimension()+1)+2) / w;
  } else {
    x = *(sc->coefs_begin() + i*(sc->dimension() )+0) / w;
    y = *(sc->coefs_begin() + i*(sc->dimension() )+1) / w;
    z = *(sc->coefs_begin() + i*(sc->dimension() )+2) / w;
  }
  
  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(x,y,z));

  return (PyObject*) result;
}

PyMethodDef Curve_methods[] = {
     {(char*)"AppendCurve",         (PyCFunction)Curve_AppendCurve,         METH_VARARGS|METH_KEYWORDS, curve_append_curve__doc__},
     {(char*)"Clone",               (PyCFunction)Curve_Clone,               METH_VARARGS,               curve_clone__doc__},
     {(char*)"Evaluate",            (PyCFunction)Curve_Evaluate,            METH_VARARGS|METH_KEYWORDS, curve_evaluate__doc__},
     {(char*)"FlipParametrization", (PyCFunction)Curve_FlipParametrization, METH_VARARGS,               curve_flip_parametrization__doc__},
     {(char*)"GetKnots",            (PyCFunction)Curve_GetKnots,            METH_VARARGS,               curve_get_knots__doc__},
     {(char*)"GetOrder",            (PyCFunction)Curve_GetOrder,            METH_VARARGS,               curve_get_order__doc__},
     {(char*)"InsertKnot",          (PyCFunction)Curve_InsertKnot,          METH_VARARGS|METH_KEYWORDS, curve_insert_knot__doc__},
     {(char*)"Normalize",           (PyCFunction)Curve_Normalize,           METH_VARARGS,               curve_normalize__doc__},
     {(char*)"Project",             (PyCFunction)Curve_Project,             METH_VARARGS|METH_KEYWORDS, curve_project__doc__},
     {(char*)"RaiseOrder",          (PyCFunction)Curve_RaiseOrder,          METH_VARARGS|METH_KEYWORDS, curve_raise_order__doc__},
     {(char*)"Rotate",              (PyCFunction)Curve_Rotate,              METH_VARARGS|METH_KEYWORDS, curve_rotate__doc__},
     {(char*)"Split",               (PyCFunction)Curve_Split,               METH_VARARGS|METH_KEYWORDS, curve_split__doc__},
     {(char*)"Translate",           (PyCFunction)Curve_Translate,           METH_VARARGS|METH_KEYWORDS, curve_translate__doc__},
     {NULL,                         NULL,                                   0,                          NULL}
   };

PyNumberMethods Curve_operators = {0};
PySequenceMethods Curve_seq_operators = {0};

PyDoc_STRVAR(curve__doc__, "A parametric description of a curve");
void init_Curve_Type()
{
  Curve_seq_operators.sq_item = Curve_GetComponent;
  Curve_seq_operators.sq_length = Curve_NmbComponent;
  InitializeTypeObject(&Curve_Type);
  Curve_operators.nb_add = Curve_Add;
  Curve_Type.tp_name = "GoTools.Curve";
  Curve_Type.tp_basicsize = sizeof(Curve);
  Curve_Type.tp_dealloc = (destructor)Curve_Dealloc;
  Curve_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES;
  Curve_Type.tp_doc = curve__doc__;
  Curve_Type.tp_methods = Curve_methods;
  Curve_Type.tp_base = 0;
  Curve_Type.tp_new = Curve_New;
  Curve_Type.tp_str = (reprfunc)Curve_Str;
  Curve_Type.tp_as_number = &Curve_operators;
  Curve_Type.tp_as_sequence = &Curve_seq_operators;
  PyType_Ready(&Curve_Type);
}
}

shared_ptr<Go::SplineCurve> convertSplineCurve(shared_ptr<Go::ParamCurve> curve)
{
  if (!curve)
    return shared_ptr<Go::SplineCurve>();

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
