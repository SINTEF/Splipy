//===========================================================================
//                                                                           
// File: surface.C
//                                                                           
// Created: Mon Feb 20 13:22:00 2012                                         
//                                                                           
// Author: Arne Morten Kvarving <arne.morten.kvarving@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description: GoTools.Surface python class implementation
//                                                                           
//===========================================================================

#include "undef.h"
#include "surface.h"
#include "curve.h"
#include "pyutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"
#include "sislP.h"

#include <fstream>
#include <sstream>

namespace GeoModeller {

extern "C"
{
PyTypeObject Surface_Type;

PyObject* Surface_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  Surface* self;
  self = (Surface*)type->tp_alloc(type,0);
  static const char* keyWords[] = {"order1", "order2", "knots1", "knots2 ", "coefs", "rational", NULL };
  PyObject* knots1o=0;
  PyObject* knots2o=0;
  PyObject* coefso=0;
  bool rational=false;
  int order1=0;
  int order2=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|iiOOOb",
                                   (char**)keyWords, &order1, &order2, &knots1o, 
                                                  &knots2o, &coefso, &rational))
    return NULL;

  // surface is specified
  if (knots1o && knots2o && coefso) {
    std::vector<double> knots1;
    std::vector<double> knots2;
    std::vector<double> coefs;
    if (PyObject_TypeCheck(knots1o,&PyList_Type)) {
      for (int i=0;i<PyList_Size(knots1o);++i) {
        PyObject* o = PyList_GetItem(knots1o,i);
        if (o && PyObject_TypeCheck(o,&PyFloat_Type))
          knots1.push_back(PyFloat_AsDouble(o));
        else if (o && PyObject_TypeCheck(o,&PyInt_Type))
          knots1.push_back(PyInt_AsLong(o));
      }
    }
    if (PyObject_TypeCheck(knots2o,&PyList_Type)) {
      for (int i=0;i<PyList_Size(knots2o);++i) {
        PyObject* o = PyList_GetItem(knots2o,i);
        if (o && PyObject_TypeCheck(o,&PyFloat_Type))
          knots2.push_back(PyFloat_AsDouble(o));
        else if (o && PyObject_TypeCheck(o,&PyInt_Type))
          knots2.push_back(PyInt_AsLong(o));
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
        if (o && PyObject_TypeCheck(o,&PyInt_Type))
          coefs.push_back(PyInt_AsLong(o));
      }
    }
    int dim = coefs.size()/((knots1.size()-order1)*(knots2.size()-order2));
    if (rational)
      dim--;
    ((Surface*)self)->data.reset(new Go::SplineSurface(knots1.size()-order1,
                                                       knots2.size()-order2, 
                                                       order1, order2,
                                                       knots1.begin(), knots2.begin(),
                                                       coefs.begin(), dim, rational));
  }

  return (PyObject*)self;
}

void Surface_Dealloc(Surface* self)
{
  self->data.reset();
  self->ob_type->tp_free((PyObject*)self);
}

PyDoc_STRVAR(surface_clone__doc__,
             "Clone a surface\n"
             "@param coefs: Coefficients for new surface\n"
             "@type coefs: List of (list of float)\n"
             "@param ncomp: Number of components in coefficients\n"
             "@type ncomp: Integer\n"
             "@return: New (copy of) surface\n");
PyObject* Surface_Clone(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"coefs", "ncomp", NULL};
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);

  PyObject* coefso=NULL;
  int dim = 1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|Oi",
                                   (char**)keyWords,&coefso,&dim))
    return NULL;

  Surface* res = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  if (coefso) {
    shared_ptr<Go::SplineSurface> surf = convertSplineSurface(parSurf);
    int nCoefs = surf->numCoefs_u()*surf->numCoefs_v();
    if (PyList_Size(coefso)/nCoefs == 0) {
      PyErr_SetString(PyExc_ValueError, "Too few coefficients");
      return NULL;
    }

    std::vector<double> coefs;
    for (int i=0;i<PyList_Size(coefso);++i) {
      PyObject* o = PyList_GetItem(coefso,i);
      if (PyObject_TypeCheck(o,&PyList_Type)) {
        dim = PyList_Size(o);
        for (size_t l=0;l<PyList_Size(o);++l)
          coefs.push_back(PyFloat_AsDouble(PyList_GetItem(o,l)));
      } else
        coefs.push_back(PyFloat_AsDouble(o));
    }
    res->data.reset(new Go::SplineSurface(surf->basis_u(), surf->basis_v(),
                                          coefs.begin(), dim, false));
  } else
    res->data.reset(parSurf->clone());
 
  return (PyObject*)res;
}

PyDoc_STRVAR(surface_append__doc__,
             "Append another surface to this surface\n"
             "@return: The surface\n"
             "@rtype: Surface\n");
PyObject* Surface_Append(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"other", "dir", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(parSurf);
  PyObject* other;
  int dir=0; 
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|i",
                                   (char**)keyWords,&other) || !surf)
    return NULL;

  shared_ptr<Go::ParamSurface> otherPar = PyObject_AsGoSurface(other);
  if (!otherPar) {
    PyErr_SetString(PyExc_TypeError, "Expected surface");
    return NULL;
  }
  if (dir > 2 || dir < 0) {
    PyErr_SetString(PyExc_ValueError, "Invalid value for `dir`");
    return NULL;
  }

  double dist;
  surf->appendSurface(otherPar.get(), dir, 0, dist, false);

  Py_INCREF(self);
  return self;
}

PyObject* Surface_Str(Surface* self)
{
  std::stringstream str;
  printSurfaceToStream(str,self->data);
  return PyString_FromString(str.str().c_str());
}

PyDoc_STRVAR(surface_evaluate__doc__,
             "Evaluate surface at given parameter values\n"
             "@param value_u: The u parameter value\n"
             "@type value_u: float\n"
             "@param value_v: The v parameter value\n"
             "@type value_v: float\n"
             "@return: The value of the surface at the given parameters\n"
             "@rtype: Point");
PyObject* Surface_Evaluate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"value_u", "value_v", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  double value_u=0, value_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"dd",
                                   (char**)keyWords,&value_u,&value_v) || !parSurf)
    return NULL;

  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(parSurf->dimension()));
  parSurf->point(*result->data,value_u,value_v);

  return (PyObject*)result;
}

PyDoc_STRVAR(surface_evaluategrid__doc__,
             "Evaluate surface on a tensor-product grid\n"
             "@param param_u: The u parameters\n"
             "@type param_u: List of float\n"
             "@param param_v: The v parameters\n"
             "@type param_v: List of float\n"
             "@return: The values of the surface at the given parameters\n"
             "@rtype: List of Points (floats for scalars)");
PyObject* Surface_EvaluateGrid(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"param_u", "param_v", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(parSurf);
  PyObject* paramuo;
  PyObject* paramvo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords,&paramuo,&paramvo) || !parSurf)
    return NULL;

  if (!PyObject_TypeCheck(paramuo,&PyList_Type) ||
      !PyObject_TypeCheck(paramvo,&PyList_Type)) {
    PyErr_SetString(PyExc_TypeError, "Expected list");
    return NULL;
  }

  std::vector<double> paramu;
  for (int i=0;i<PyList_Size(paramuo);++i)
    paramu.push_back(PyFloat_AsDouble(PyList_GetItem(paramuo,i)));
  std::vector<double> paramv;
  for (int i=0;i<PyList_Size(paramvo);++i)
    paramv.push_back(PyFloat_AsDouble(PyList_GetItem(paramvo,i)));

  std::vector<double> res(paramu.size()*paramv.size()*surf->dimension());
  surf->gridEvaluator(res, paramu, paramv);

  PyObject* result = PyList_New(res.size()/surf->dimension());
  VectorToPyPointList(result, res, surf->dimension());

  return result;
}

PyDoc_STRVAR(surface_evaluate_normal__doc__,
             "Evaluate surface normal at given parameter values\n"
             "@param value_u: The u parameter value\n"
             "@type value_u: float\n"
             "@param value_v: The v parameter value\n"
             "@type value_v: float\n"
             "@return: The value of the surface normal at the given parameters\n"
             "@rtype: Point");
PyObject* Surface_EvaluateNormal(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"value_u", "value_v", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  double value_u=0, value_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"dd",
                                   (char**)keyWords,&value_u,&value_v) || !parSurf)
    return NULL;

  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(parSurf->dimension()));
  parSurf->normal(*result->data,value_u,value_v);

  return (PyObject*)result;
}

PyDoc_STRVAR(surface_evaluate_tangent__doc__,
             "Evaluate the two tangent vectors of the surface at given parameter values\n"
             "@param value_u: The u parameter value\n"
             "@type value_u: float\n"
             "@param value_v: The v parameter value\n"
             "@type value_v: float\n"
             "@param from_right_u: Evaluate u in the limit from positive direction (Optional)\n"
             "@type from_right_u: bool\n"
             "@param from_right_v: Evaluate v in the limit from positive direction (Optional)\n"
             "@type from_right_v: bool\n"
             "@return: The value of the two tangents at the given parameters\n"
             "@rtype: Tuple of two Point");
PyObject* Surface_EvaluateTangent(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"value_u", "value_v", "from_right_u", "from_right_v", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  bool from_right_u=true;
  bool from_right_v=true;
  double value_u=0, value_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"dd|bb",
                                   (char**)keyWords,&value_u,&value_v,
                                   &from_right_u, &from_right_v) || !parSurf)
    return NULL;

  // evaluate the surface and two tangent vectors
  std::vector<Go::Point> tangents(3);
  parSurf->point(tangents, value_u, value_v, from_right_u, from_right_v);

  // make the python Point objects
  Point* du = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  du->data.reset(new Go::Point(tangents[1]));

  Point* dv = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  dv->data.reset(new Go::Point(tangents[2]));

  // return a tuple of the two results
  PyObject* result = PyTuple_New(2);

  PyTuple_SetItem( result, 0, (PyObject*) du);
  PyTuple_SetItem( result, 1, (PyObject*) dv);

  return result;
}

PyDoc_STRVAR(surface_flip_parametrization__doc__,
             "Flip surface parametrization\n"
             "@param direction: The parametric direction to flip (0=u, 1=v)\n"
             "@type direction: int\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_FlipParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", NULL };
  PyObject* axiso;
  int direction = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&direction))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  Go::SplineSurface *spSurf = parSurf->asSplineSurface();
  spSurf->reverseParameterDirection(direction == 0);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_force_rational__doc__,
             "Enforce a rational representation of the spline surface\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_ForceRational(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  Go::SplineSurface *spSurf = parSurf->asSplineSurface();
  spSurf->representAsRational();

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_lower_order__doc__,
             "Lower the order of a spline surface (need full continuity)\n"
             "@param lower_u: Lower of order in u\n"
             "@type lower_u: int (>= 0)\n"
             "@param lower_v: Lower of order in v\n"
             "@type lower_v: int (>= 0)\n"
             "@returns: The surface\n"
             "@rtype: Surface");
PyObject* Surface_LowerOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"lower_u", "lower_v", NULL };
  int lower_u=0, lower_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords,&lower_u,&lower_v))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }
  shared_ptr<Go::SplineSurface> spSurf = 
                       static_pointer_cast<Go::SplineSurface>(parSurf);
  std::vector<double>::const_iterator first =  spSurf->basis(0).begin()+lower_u;
  std::vector<double>::const_iterator last  =  spSurf->basis(0).end()-lower_u;
  Go::BsplineBasis b1 = Go::BsplineBasis(spSurf->order_u()-lower_u,first,last);
  first =  spSurf->basis(1).begin()+lower_v;
  last  =  spSurf->basis(1).end()-lower_v;
  Go::BsplineBasis b2 = Go::BsplineBasis(spSurf->order_v()-lower_v,first,last);

  if (spSurf->rational()) {
    if (PyErr_WarnEx(PyExc_RuntimeWarning,
                     "The geometry basis is rational (using NURBS). "
                     "The basis for the unknown fields of one degree "
                     "higher will however be non-rational. "
                     "This may affect accuracy.", 1) == -1)
      return NULL;
  }

  std::vector<double> ug(b1.numCoefs()), vg(b2.numCoefs());
  for (size_t i = 0; i < ug.size(); i++)
    ug[i] = b1.grevilleParameter(i);
  for (size_t i = 0; i < vg.size(); i++)
    vg[i] = b2.grevilleParameter(i);

  // Evaluate the spline surface at all points
  std::vector<double> XYZ(spSurf->dimension()*ug.size()*vg.size());
  spSurf->gridEvaluator(XYZ,ug,vg);

  // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
  Surface* pySurf = (Surface*)self;
  pySurf->data = convertSplineSurface(parSurf);
  pySurf->data.reset(Go::SurfaceInterpolator::regularInterpolation(b1, b2, ug,
                                                                   vg, XYZ,
                                                                   spSurf->dimension(),
                                                                   false, XYZ));

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_raise_order__doc__,
             "Raise order of a spline surface\n"
             "@param raise_u: Raise of order in u\n"
             "@type raise_u: int (>= 0)\n"
             "@param raise_v: Raise of order in v\n"
             "@type raise_v: int (>= 0)\n"
             "@returns: The surface\n"
             "@rtype: Surface");
PyObject* Surface_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"raise_u", "raise_v", NULL };
  int raise_u=0, raise_v=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords,&raise_u,&raise_v))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }
  static_pointer_cast<Go::SplineSurface>(parSurf)->raiseOrder(raise_u,raise_v);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_get_const_par_curve__doc__,
             "Generate and return a SplineCurve that represents a constant parameter curve on the surface\n"
             "@param parameter: Value of the fixed parameter\n"
             "@type parameter: float\n"
             "@param pardir: 0 for constant u-parameter, 1 for v\n"
             "@type pardir: int\n"
             "@return A newly constructed SplineCurve representing the specific constant parameter curve\n"
             "@rtype: Curve");
PyObject* Surface_GetConstParCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  int parDir = 0;
  double parameter = 0;

  static const char* keyWords[] = {"parameter", "pardir", NULL};
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!PyArg_ParseTupleAndKeywords(args, kwds, (char*)"di",
                                   (char**)keyWords, &parameter, &parDir) || !parSurf)
    return NULL;

  if (!parSurf->isSpline()) {
    Surface *pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(parSurf);
  Go::SplineCurve* spCurve = spSurf->constParamCurve(parameter, parDir == 0);

  Curve* pyCurve = (Curve*)Curve_Type.tp_alloc(&Curve_Type, 0);
  pyCurve->data = shared_ptr<Go::ParamCurve>(spCurve);

  return (PyObject*)pyCurve;
}

PyDoc_STRVAR(surface_get_edges__doc__,
             "Return the four edge curves in (parametric) order: bottom, right, top, left\n"
             "@return: The four edges\n"
             "@rtype: List of Curve");
PyObject* Surface_GetEdges(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  Go::SplineSurface* spSurf = parSurf->asSplineSurface();
  Curve* curve0 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Curve* curve1 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Curve* curve2 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  Curve* curve3 = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  curve0->data = shared_ptr<Go::SplineCurve>(spSurf->edgeCurve(0));
  curve1->data = shared_ptr<Go::SplineCurve>(spSurf->edgeCurve(1));
  curve2->data = shared_ptr<Go::SplineCurve>(spSurf->edgeCurve(2));
  curve3->data = shared_ptr<Go::SplineCurve>(spSurf->edgeCurve(3));

  PyObject* result = PyList_New(0);
  PyList_Append(result, (PyObject*) curve0);
  PyList_Append(result, (PyObject*) curve1);
  PyList_Append(result, (PyObject*) curve2);
  PyList_Append(result, (PyObject*) curve3);

  return result;
}


PyDoc_STRVAR(surface_get_knots__doc__,
             "Return knots for a spline surface\n"
             "@param with_multiplicities: (optional) Set to true to obtain the knot vectors with multiplicities\n"
             "@type with_multiplicities: Boolean\n"
             "@return: List with the knot values\n"
             "@rtype: Tuple with List of float");
PyObject* Surface_GetKnots(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"with_multiplicities", NULL };
  bool withmult=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|b",
                                   (char**)keyWords, &withmult))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  PyObject* result = PyTuple_New(2);
  shared_ptr<Go::SplineSurface> spSurf = static_pointer_cast<Go::SplineSurface>(parSurf);
  for (int i=0;i<2;++i) {
    PyObject* list = PyList_New(0);
    std::vector<double> knots;
    if (withmult) {
      Go::BsplineBasis& basis = (i==0?spSurf->basis_u():spSurf->basis_v());
      knots = basis.getKnots();
    } else
      spSurf->basis(i).knotsSimple(knots);
    for (std::vector<double>::iterator it  = knots.begin();
                                       it != knots.end();++it) {
      PyList_Append(list,Py_BuildValue((char*)"d",*it));
    }
    PyTuple_SetItem(result,i,list);
  }

  return result;
}

PyDoc_STRVAR(surface_get_greville__doc__,
             "Return Greville points for a spline surface\n"
             "@return: List with the parameter values\n"
             "@rtype: Tuple with List of float");
PyObject* Surface_GetGreville(PyObject* self)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  PyObject* result = PyTuple_New(2);
  shared_ptr<Go::SplineSurface> spSurf = static_pointer_cast<Go::SplineSurface>(parSurf);
  for (int i=0;i<2;++i) {
    PyObject* list = PyList_New(0);
    for (size_t j=0;j<spSurf->basis(i).numCoefs();++j)
      PyList_Append(list,Py_BuildValue((char*)"d",
                     spSurf->basis(i).grevilleParameter(j)));
    PyTuple_SetItem(result,i,list);
  }

  return result;
}

PyDoc_STRVAR(surface_get_order__doc__,
             "Return spline surface order (polynomial degree + 1) in all parametric directions\n"
             "@return: B-Spline order\n"
             "@rtype: List of two integers");
PyObject* Surface_GetOrder(PyObject* self, PyObject* args)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  shared_ptr<Go::SplineSurface> spline = static_pointer_cast<Go::SplineSurface>(parSurf);
  if(!spline) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to coerce surface to spline form");
    return NULL;
  }

  PyObject* result = PyList_New(0);

  PyList_Append(result, Py_BuildValue((char*) "i", spline->order_u()));
  PyList_Append(result, Py_BuildValue((char*) "i", spline->order_v()));

  return result;
}

PyDoc_STRVAR(surface_get_sub_surf__doc__,
             "Get a Spline Surface which represent a part of 'this' Surface\n"
             "@param from_par: The parametric lower left corner of the sub surface\n"
             "@type  from_par: Point, list of floats or tuple of floats\n"
             "@param to_par:   The parametric upper right corner of the sub surface\n"
             "@type  to_par:   Point, list of floats or tuple of floats\n"
             "@return: A parametric rectangular subregion of this Surface represented as a Spline Surface\n"
             "@rtype: Surface");
PyObject* Surface_GetSubSurf(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  static const char* keyWords[] = {"from_par", "to_par", NULL };
  PyObject *lowerLefto, *upperRighto;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords, &lowerLefto, &upperRighto) || !parSurf)
    return NULL;

  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }
  shared_ptr<Go::SplineSurface> spSurf = static_pointer_cast<Go::SplineSurface>(parSurf);
  if(!spSurf) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineSurface");
    return NULL;
  }

  shared_ptr<Go::Point> lowerLeft  = PyObject_AsGoPoint(lowerLefto);
  shared_ptr<Go::Point> upperRight = PyObject_AsGoPoint(upperRighto);

  Go::SplineSurface *subSurf = spSurf->subSurface((*lowerLeft)[0], (*lowerLeft)[1], (*upperRight)[0], (*upperRight)[1]);

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data = shared_ptr<Go::ParamSurface>(subSurf);

  return (PyObject*) result;
}

PyDoc_STRVAR(surface_get_parameter_at_point__doc__,
             "Get surface parameter values at a geometric Point\n"
             "@param point: The geometric point to intersect \n"
             "@type point: Point\n"
             "@param tolerance: (optional) Tolerance in search\n"
             "@type tolerance: Float\n"
             "@return: Parameters, distance and actual point on geometry\n"
             "@rtype: Tuple of (list of Float, Float, Point)");
PyObject* Surface_GetParameterAtPoint(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  static const char* keyWords[] = {"point", "tolerance", NULL };
  PyObject* point=nullptr;
  double tolerance=1e-5;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|d",
                                   (char**)keyWords, &point, &tolerance) || !parSurf)
    return NULL;
  shared_ptr<Go::Point> pt = PyObject_AsGoPoint(point);
  std::vector<double> clo(2);
  double dist;
  Go::Point clopt(3);
  parSurf->closestPoint(*pt, clo[0], clo[1], clopt, dist, tolerance);
  Point* rpt = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  rpt->data.reset(new Go::Point(clopt));

  PyObject* list = PyList_New(2);
  VectorToPyPointList(list, clo, 1);

  PyObject* result = PyTuple_New(3);
  PyTuple_SetItem(result, 0, list);
  PyTuple_SetItem(result, 1, PyFloat_FromDouble(dist));
  PyTuple_SetItem(result, 2, (PyObject*)rpt);

  return result;
}

PyDoc_STRVAR(surface_insert_knot__doc__,
             "Insert a knot in a spline surface\n"
             "@param direction: Direction to insert knot in\n"
             "@type direction: int (0 or 1)\n"
             "@param knot: The knot to insert\n"
             "@type knot: float\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", "knot", NULL };
  int direction=0;
  double knot;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"id",
                                   (char**)keyWords,&direction,&knot))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }
  if (direction == 0)
    static_pointer_cast<Go::SplineSurface>(parSurf)->insertKnot_u(knot);
  else
    static_pointer_cast<Go::SplineSurface>(parSurf)->insertKnot_v(knot);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_interpolate__doc__,
             "Interpolate a field onto the surface' basis\n"
             "@param coefs: The values to interpolate (sample in the Greville points)\n"
             "@type coefs: List of Float"
             "@return: New basis coefficients\n"
             "@rtype: List of float");
PyObject* Surface_Interpolate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"coefs", NULL };
  PyObject* pycoefs=NULL;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&pycoefs))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;
  if (!PyObject_TypeCheck(pycoefs,&PyList_Type)) {
    PyErr_SetString(PyExc_TypeError, "Expected list");
    return NULL;
  }

  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }
  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(parSurf);

  std::vector<double> coefs;
  for (int i=0;i<PyList_Size(pycoefs);++i) {
    PyObject* o = PyList_GetItem(pycoefs,i);
    coefs.push_back(PyFloat_AsDouble(o));
  }

  std::vector<double> greville_u(surf->basis_u().numCoefs());
  std::vector<double> greville_v(surf->basis_v().numCoefs());
  for(int i=0; i<surf->basis_u().numCoefs(); i++)
    greville_u[i] = surf->basis_u().grevilleParameter(i);
  for(int i=0; i<surf->basis_v().numCoefs(); i++)
    greville_v[i] = surf->basis_v().grevilleParameter(i);
  std::vector<double> weights(0);

  int dim = coefs.size()/(greville_u.size()*greville_v.size());

  Go::SplineSurface* res =
        Go::SurfaceInterpolator::regularInterpolation(surf->basis_u(),
                                                      surf->basis_v(),
                                                      greville_u,
                                                      greville_v,
                                                      coefs,
                                                      dim,
                                                      false,
                                                      weights);

  PyObject* result = PyList_New(0);
  for (std::vector<double>::const_iterator it  = res->coefs_begin();
                                           it != res->coefs_end();++it)
    PyList_Append(result,Py_BuildValue((char*)"d",*it));

  delete res;

  return result;
}

PyDoc_STRVAR(surface_intersect__doc__,
             "Check if this surface intersects another surface or surface.\n"
             "Set tolerance 'gap' for intersection tolerance \n"
             "@param obj: The object to test against\n"
             "@type obj: Curve or Surface\n"
             "@return: True if the surface intersects obj\n"
             "@rtype: Bool");
PyObject* Surface_Intersect(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"obj", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  double knot;
  PyObject *pyObj;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&pyObj) || !parSurf)
    return NULL;

  if(PyObject_TypeCheck(pyObj, &Curve_Type) ) {
    shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(pyObj);
    if(!parCrv) {
      PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::ParamCurve");
      return NULL;
    }

    // input arguments
    double         geomRes       = modState.gapTolerance;
    double         startCrv      = parCrv->startparam();
    double         endCrv        = parCrv->endparam();
    Go::RectDomain domain        = parSurf->containingDomain();
    double         itStartCrv    = (startCrv + endCrv)/2.0;
    double         itStartSurf[] = { (domain.umin()+domain.umax())/2,
                                     (domain.vmin()+domain.vmax())/2 };
    //output arguments
    double         crvPos;
    double         surfPos[2];
    double         dist;
    Go::Point      crvPt;
    Go::Point      surfPt;
    // newton iteration search for closest point on surface
    Go::ClosestPoint::closestPtCurveSurf(parCrv.get(), parSurf.get(), geomRes,
                                         startCrv, endCrv,         // curve parameter range 
                                         &domain,                  // surface parameter range 
                                         itStartCrv, itStartSurf,  // iteration start
                                         crvPos, surfPos,          // [out] parametric closest points
                                         dist,                     // [out] geometric distance between points
                                         crvPt,  surfPt);          // [out] geometric closest points
    if(dist < modState.gapTolerance) {
      // Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
      // result->data.reset(new Go::Point(crvPt));
      // return (PyObject*) result;
      return Py_BuildValue((char*)"b",true);
    } else {
      // Py_INCREF(Py_None);
      // return Py_None;
      return Py_BuildValue((char*)"b",false);
    }
  } else if(PyObject_TypeCheck(pyObj, &Surface_Type) ) {
    shared_ptr<Go::ParamSurface> parSurf2 = PyObject_AsGoSurface(pyObj);
    if(!parSurf2) {
      PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::ParamSurface");
      return NULL;
    }
    // convert to spline surfaces
    shared_ptr<Go::SplineSurface> spSurf1 = convertSplineSurface(parSurf);
    shared_ptr<Go::SplineSurface> spSurf2 = convertSplineSurface(parSurf2);
    if(!spSurf1 || !spSurf2) {
      PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineSurface");
      return NULL;
    }
    // convert to SISL surfaces
    SISLSurf* sislSurf1 = GoSurf2SISL(*spSurf1, true);
    SISLSurf* sislSurf2 = GoSurf2SISL(*spSurf2, true);
    if(!sislSurf1 || !sislSurf2) {
      PyErr_SetString(PyExc_RuntimeError, "Unable to obtain SISLSurf");
      return NULL;
    }
    // setup parameters
    double epsco = 0.0;                   // Computational resolution (not used)
    double epsge = modState.gapTolerance; // Geometry resolution
    int    numPts;                        // number of single intersection points
    double *parPt1;                       // parameter values of the single intersection pts
    double *parPt2;                       // parameter values of the single intersection pts
    int    numCrv;                        // number of intersection curves
    SISLIntcurve **curves;                // intersection curves
    int    status;                        // warnings and error messages

    // actual SISL call
    s1859(sislSurf1, sislSurf2, epsco, epsge, // input arguments
          &numPts, &parPt1, &parPt2,          // output points
          &numCrv, &curves,                   // output curves
          &status);                           // output errors

    // error handling
    if (status != 0) {
      std::ostringstream ss;
      ss << "SISL returned " << (status > 0 ? "warning" : "error") << " code " << status;
      if (status > 0) {
        // Warnings may throw exceptions, depending on user settings,
        // in that case we are obliged to treat it as one
        if (PyErr_WarnEx(PyExc_RuntimeWarning, ss.str().c_str(), 1) == -1)
          return NULL;
      } else {
        PyErr_SetString(PyExc_RuntimeError, ss.str().c_str());
        return NULL;
      }
    }

    // return results
    if(numPts > 0 || numCrv > 0) {
      return Py_BuildValue((char*)"b",true);
    } else {
      return Py_BuildValue((char*)"b",false);
    }

  }
  return NULL;
}

PyDoc_STRVAR(surface_project__doc__,
             "Project the surface onto an axis or plane parallel to the cartesian coordinate system\n"
             "@param axis: The axis or plane to project onto (\"X\",\"Y\",\"Z\" or a combination of these)\n"
             "@type axis: string\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_Project(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", NULL };
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  char *sAxis;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&sAxis) || !parSurf) {
    Py_INCREF(self);
    return self;
  }

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

  Surface* pySurf = (Surface*)self;
  pySurf->data = convertSplineSurface(parSurf);

  shared_ptr<Go::SplineSurface> spSurf   = static_pointer_cast<Go::SplineSurface>(pySurf->data);
  bool                          rational = spSurf->rational();
  int                           dim      = spSurf->dimension();
  std::vector<double>::iterator coefs    = (rational) ? spSurf->rcoefs_begin() : spSurf->coefs_begin();
  std::vector<double>::iterator coefsEnd = (rational) ? spSurf->rcoefs_end()   : spSurf->coefs_end();
  while(coefs != coefsEnd) {
    if(bAxis[0] && dim>0)
      coefs[0] = 0.0;
    if(bAxis[1] && dim>1)
      coefs[1] = 0.0;
    if(bAxis[2] && dim>2)
      coefs[2] = 0.0;
    coefs += (dim+rational);
  }

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_reparametrize__doc__,
             "Re-parametrize a surface\n"
             "@param umin: The minimum u value\n"
             "@type umin: float\n"
             "@param umax: The maximum u value\n"
             "@type umax: float\n"
             "@param vmin: The minimum v value\n"
             "@type vmin: float\n"
             "@param vmax: The maximum v value\n"
             "@type vmax: float\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_ReParametrize(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"umin", "umax", "vmin", "vmax", NULL };
  double umin=0, umax=1, vmin=0, vmax=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|dddd",
                                   (char**)keyWords,&umin,&umax,&vmin,&vmax))
    return NULL;
  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(self);
  if (surface) {
    shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(surface);
    if (!surface->isSpline())
      ((Surface*)self)->data = spSurf;
    spSurf->setParameterDomain(umin, umax, vmin, vmax);
  }

  Py_INCREF(self);
  return self;
}


PyDoc_STRVAR(surface_rotate__doc__,
             "Rotate a surface around an axis\n"
             "@param axis: The axis to rotate around\n"
             "@type axis: Point, list of floats or tuple of floats\n"
             "@param angle: Angle to rotate surface with in radians\n"
             "@type angle: float\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_Rotate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", "angle", NULL };
  PyObject* axiso;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  if (!axis) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for axis: expected pointlike");
    return NULL;
  }

  if (parSurf->instanceType() == Go::Class_Cylinder) { 
    // since we can't move cylinders, we have to create a new one and trash the old one
    shared_ptr<Go::Cylinder> cyl = static_pointer_cast<Go::Cylinder>(parSurf);
    Go::Point x, y, z, c; // x,y,z-axis and center point
    double r = cyl->getRadius();
    cyl->getCoordinateAxes(x,y,z);
    c = cyl->getLocation();

    Go::GeometryTools::rotatePoint(*axis, angle, c.begin());
    Go::GeometryTools::rotatePoint(*axis, angle, x.begin());
    Go::GeometryTools::rotatePoint(*axis, angle, y.begin());
    Go::GeometryTools::rotatePoint(*axis, angle, z.begin());

    Go::Cylinder *result = new Go::Cylinder(r, c, z, x);
    if(cyl->isBounded()) {
      Go::RectDomain domain = cyl->containingDomain();
      result->setParameterBounds(domain.umin(), domain.vmin(), domain.umax(), domain.vmax());
    }
    ((Surface*)self)->data.reset(result);

  } else if (!parSurf->isSpline()) {
    // for all other elementary surfaces, we convert to SplineSurface and work on with that
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
    Go::GeometryTools::rotateSplineSurf(*axis, angle,
                       *static_pointer_cast<Go::SplineSurface>(parSurf));
  } else {
    Go::GeometryTools::rotateSplineSurf(*axis, angle,
                       *static_pointer_cast<Go::SplineSurface>(parSurf));
  }

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_swap_parametrization__doc__,
             "Swaps the two surface parameter directions\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_SwapParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  if (!parSurf->isSpline()) {
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
  }

  Go::SplineSurface *spSurf = parSurf->asSplineSurface();
  spSurf->swapParameterDirection();

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(surface_translate__doc__,
             "Translate a surface along a given vector\n"
             "@param vector: The vector to translate along\n"
             "@type vector: Point, list of floats or tuple of floats\n"
             "@return: The surface\n"
             "@rtype: Surface");
PyObject* Surface_Translate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"vector", NULL };
  PyObject* veco;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&veco))
    return NULL;

  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  shared_ptr<Go::Point> vec = PyObject_AsGoPoint(veco);
  if (!vec) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for vector: expected pointlike");
    return NULL;
  }

  if (parSurf->instanceType() == Go::Class_Cylinder) { 
    // since we can't move cylinders, we have to create a new one and trash the old one
    shared_ptr<Go::Cylinder> cyl = static_pointer_cast<Go::Cylinder>(parSurf);
    Go::Point x, y, z, c; // x,y,z-axis and center point
    double r = cyl->getRadius();
    cyl->getCoordinateAxes(x,y,z);
    c = cyl->getLocation();

    Go::Cylinder *result = new Go::Cylinder(r, c + (*vec), z, x);
    if(cyl->isBounded()) {
      Go::RectDomain domain = cyl->containingDomain();
      result->setParameterBounds(domain.umin(), domain.vmin(), domain.umax(), domain.vmax());
    }
    ((Surface*)self)->data.reset(result);

  } else if (!parSurf->isSpline()) {
    // for all other elementary surfaces, we convert to SplineSurface and work on with that
    Surface* pySurf = (Surface*)self;
    pySurf->data = convertSplineSurface(parSurf);
    parSurf = pySurf->data;
    Go::GeometryTools::translateSplineSurf(*vec, 
                              *static_pointer_cast<Go::SplineSurface>(parSurf));
  } else {
    Go::GeometryTools::translateSplineSurf(*vec, 
                              *static_pointer_cast<Go::SplineSurface>(parSurf));
  }


  Py_INCREF(self);
  return self;
}

PyObject* Surface_Reduce(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(parSurf);
  if(!spSurf) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineSurface");
    return NULL;
  }

  PyObject* knots_u = PyList_New(0);
  vector<double>::const_iterator kit;
  for (kit = spSurf->basis_u().begin(); kit != spSurf->basis_u().end(); ++kit)
    PyList_Append(knots_u, Py_BuildValue("d", *kit));

  PyObject* knots_v = PyList_New(0);
  for (kit = spSurf->basis_v().begin(); kit != spSurf->basis_v().end(); ++kit)
    PyList_Append(knots_v, Py_BuildValue("d", *kit));

  PyObject* coefs = PyList_New(0);
  vector<double>::const_iterator cit = spSurf->rational() ? spSurf->rcoefs_begin() : spSurf->coefs_begin();
  vector<double>::const_iterator cit_end = spSurf->rational() ? spSurf->rcoefs_end() : spSurf->coefs_end();
  for (; cit != cit_end; ++cit)
    PyList_Append(coefs, Py_BuildValue("d", *cit));
  
  return Py_BuildValue("(O(iiOOOO))", &Surface_Type, spSurf->order_u(), spSurf->order_v(),
                       knots_u, knots_v, coefs, spSurf->rational() ? Py_True : Py_False);
}

static std::array<std::vector<double>,2>
  getTesselationParams(shared_ptr<Go::SplineSurface>& surf, int n[2])
{
  // Grab parameter values in evaluation points
  std::array<std::vector<double>,2> gpar;
  gpar[0] = Tesselate(surf->basis(0).begin(), surf->basis(0).end(), n[0]);
  gpar[1] = Tesselate(surf->basis(1).begin(), surf->basis(1).end(), n[1]);

  return gpar;
}

PyDoc_STRVAR(surface_get_tesselationparams__doc__,
             "Obtain tesselation parameters for a surface\n"
             "@param n1: Number of tesselation points per knotspan in u\n"
             "@type n1: Int\n"
             "@param n2: Number of tesselation points per knotspan in v\n"
             "@type n2: Int\n"
             "@rtype: Tuple with (list of float, list of float)");
PyObject* Surface_GetTesselationParams(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(parSurf);
  static const char* keyWords[] = {"n1", "n2", NULL};
  int np[2] = {1,1};
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|ii",
                                   (char**)keyWords,np,np+1))
    return NULL;

  // Grab parameter values in evaluation points
  std::array<std::vector<double>,2> gpar = getTesselationParams(surf, np);

  size_t nx = gpar[0].size();
  size_t ny = gpar[1].size();

  PyObject* g1 = PyList_New(nx);
  for (size_t i=0;i<nx;++i)
    PyList_SetItem(g1, i, Py_BuildValue((char*)"d",gpar[0][i]));
  PyObject* g2 = PyList_New(ny);
  for (size_t i=0;i<ny;++i)
    PyList_SetItem(g2, i, Py_BuildValue((char*)"d",gpar[1][i]));

  PyObject* result = PyTuple_New(2);
  PyTuple_SetItem(result, 0, g1);
  PyTuple_SetItem(result, 1, g2);

  return result;
}

PyDoc_STRVAR(surface_tesselate__doc__,
             "Tesselate a surface\n"
             "@param n1: Number of tesselation points per knotspan in u\n"
             "@type n1: Int\n"
             "@param n2: Number of tesselation points per knotspan in v\n"
             "@type n2: Int\n"
             "@rtype: Tuple with (list of nodes, list of elements)");
PyObject* Surface_Tesselate(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  shared_ptr<Go::SplineSurface> surf = convertSplineSurface(parSurf);
  static const char* keyWords[] = {"n1", "n2", NULL};
  int np[2] = {1,1};
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|ii",
                                   (char**)keyWords,np,np+1))
    return NULL;

  // Grab parameter values in evaluation points
  std::array<std::vector<double>,2> gpar = getTesselationParams(surf, np);

  // Evaluate the surface at all points
  size_t nx = gpar[0].size();
  size_t ny = gpar[1].size();
  std::vector<double> XYZ(surf->dimension()*nx*ny);
  surf->gridEvaluator(XYZ,gpar[0],gpar[1]);

  // Establish the block grid coordinates
  PyObject* gc = PyList_New(XYZ.size());
  for (int i=0;i<XYZ.size();++i)
    PyList_SetItem(gc, i, Py_BuildValue((char*)"d", XYZ[i]));

  PyObject* ge = PyList_New(0);
  // Establish the block grid topology
  int ie, nse1 = np[0];
  int je, nse2 = np[1];
  int nel1 = (nx-1)/nse1;
  int n[4], ip = 0;
  size_t i, j, l;
  for (j = je = 1, n[1] = 0; j < ny; j++)
  {
    n[0] = n[1];
    n[1] = n[0] + 1;
    n[2] = n[1] + nx;
    n[3] = n[1] + nx-1;
    for (i = ie = 1; i < nx; i++)
    {
      for (l = 0; l < 4; l++)
        PyList_Append(ge, Py_BuildValue((char*)"i",n[l]++));
      if (i%nse1 == 0) ie++;
    }
    if (j%nse2 == 0) je++;
  }

  PyObject* result = PyTuple_New(2);
  PyTuple_SetItem(result, 0, gc);
  PyTuple_SetItem(result, 1, ge);

  return result;
}

PyObject* Surface_Add(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(o1);
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(o2);
  Surface* result = NULL;
  if (parSurf && point) {
    result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    if (parSurf->isSpline())
      result->data = shared_ptr<Go::ParamSurface>(
          new Go::SplineSurface(*static_pointer_cast<Go::SplineSurface>(parSurf)));
    else
      result->data = convertSplineSurface(parSurf);

    Go::GeometryTools::translateSplineSurf(*point, 
                         *static_pointer_cast<Go::SplineSurface>(result->data));
  }

  if (!result) {
    PyErr_SetString(PyExc_TypeError, "Expected Surface + Point");
  }

  return (PyObject*) result;
}

PyObject* Surface_Scale(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(o1);
  double scale = 1.0;
  if (PyObject_TypeCheck(o2,&PyFloat_Type))
    scale = PyFloat_AsDouble(o2);
  if (PyObject_TypeCheck(o2,&PyInt_Type))
    scale = PyInt_AsLong(o2);
  if (parSurf && scale > 0) {
    if (!parSurf->isSpline())
      ((Surface*)o1)->data = convertSplineSurface(parSurf);
    shared_ptr<Go::SplineSurface> spSurf = 
      dynamic_pointer_cast<Go::SplineSurface,
                           Go::ParamSurface>(((Surface*)o1)->data);
    for (std::vector<double>::iterator it  = spSurf->ctrl_begin(); 
                                       it != spSurf->ctrl_end(); ) {
      for (int i=0;i<parSurf->dimension();++i)
        *it++ *= scale;
      if (spSurf->rational())
        it++;
    }
  }

  return o1;
}

PyMethodDef Surface_methods[] = {
     {(char*)"Append",              (PyCFunction)Surface_Append,                METH_VARARGS|METH_KEYWORDS, surface_clone__doc__},
     {(char*)"Clone",               (PyCFunction)Surface_Clone,                 METH_VARARGS|METH_KEYWORDS, surface_clone__doc__},
     {(char*)"Evaluate",            (PyCFunction)Surface_Evaluate,              METH_VARARGS|METH_KEYWORDS, surface_evaluate__doc__},
     {(char*)"EvaluateGrid",        (PyCFunction)Surface_EvaluateGrid,          METH_VARARGS|METH_KEYWORDS, surface_evaluategrid__doc__},
     {(char*)"EvaluateNormal",      (PyCFunction)Surface_EvaluateNormal,        METH_VARARGS|METH_KEYWORDS, surface_evaluate_normal__doc__},
     {(char*)"EvaluateTangent",     (PyCFunction)Surface_EvaluateTangent,       METH_VARARGS|METH_KEYWORDS, surface_evaluate_tangent__doc__},
     {(char*)"FlipParametrization", (PyCFunction)Surface_FlipParametrization,   METH_VARARGS|METH_KEYWORDS, surface_flip_parametrization__doc__},
     {(char*)"ForceRational",       (PyCFunction)Surface_ForceRational,         METH_VARARGS,               surface_force_rational__doc__},
     {(char*)"GetConstParCurve",    (PyCFunction)Surface_GetConstParCurve,      METH_VARARGS|METH_KEYWORDS, surface_get_const_par_curve__doc__},
     {(char*)"GetEdges",            (PyCFunction)Surface_GetEdges,              METH_VARARGS|METH_KEYWORDS, surface_get_edges__doc__},
     {(char*)"GetGreville",         (PyCFunction)Surface_GetGreville,           0,                          surface_get_greville__doc__},
     {(char*)"GetKnots",            (PyCFunction)Surface_GetKnots,              METH_VARARGS|METH_KEYWORDS, surface_get_knots__doc__},
     {(char*)"GetOrder",            (PyCFunction)Surface_GetOrder,              METH_VARARGS              , surface_get_order__doc__},
     {(char*)"GetSubSurf",          (PyCFunction)Surface_GetSubSurf,            METH_VARARGS|METH_KEYWORDS, surface_get_sub_surf__doc__},
     {(char*)"GetTesselationParams",(PyCFunction)Surface_GetTesselationParams,  METH_VARARGS|METH_KEYWORDS, surface_get_tesselationparams__doc__},
     {(char*)"GetParameterAtPoint", (PyCFunction)Surface_GetParameterAtPoint,   METH_VARARGS|METH_KEYWORDS, surface_get_parameter_at_point__doc__},
     {(char*)"InsertKnot",          (PyCFunction)Surface_InsertKnot,            METH_VARARGS|METH_KEYWORDS, surface_insert_knot__doc__},
     {(char*)"Interpolate",         (PyCFunction)Surface_Interpolate,           METH_VARARGS|METH_KEYWORDS, surface_interpolate__doc__},
     {(char*)"Intersect",           (PyCFunction)Surface_Intersect,             METH_VARARGS|METH_KEYWORDS, surface_intersect__doc__},
     {(char*)"LowerOrder",          (PyCFunction)Surface_LowerOrder,            METH_VARARGS|METH_KEYWORDS, surface_lower_order__doc__},
     {(char*)"Project",             (PyCFunction)Surface_Project,               METH_VARARGS|METH_KEYWORDS, surface_project__doc__},
     {(char*)"RaiseOrder",          (PyCFunction)Surface_RaiseOrder,            METH_VARARGS|METH_KEYWORDS, surface_raise_order__doc__},
     {(char*)"ReParametrize",       (PyCFunction)Surface_ReParametrize,         METH_VARARGS|METH_KEYWORDS, surface_reparametrize__doc__},
     {(char*)"Rotate",              (PyCFunction)Surface_Rotate,                METH_VARARGS|METH_KEYWORDS, surface_rotate__doc__},
     {(char*)"SwapParametrization", (PyCFunction)Surface_SwapParametrization,   METH_VARARGS|METH_KEYWORDS, surface_swap_parametrization__doc__},
     {(char*)"Tesselate",           (PyCFunction)Surface_Tesselate,             METH_VARARGS|METH_KEYWORDS, surface_tesselate__doc__},
     {(char*)"Translate",           (PyCFunction)Surface_Translate,             METH_VARARGS|METH_KEYWORDS, surface_translate__doc__},
     {(char*)"__reduce__",          (PyCFunction)Surface_Reduce,                0,                          NULL},
     {NULL,           NULL,                     0,            NULL}
   };

Py_ssize_t Surface_NmbComponent(PyObject* self)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return 0;

  shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(parSurf);
  if(!spSurf)
    return 0;

  return spSurf->numCoefs_u() * spSurf->numCoefs_v();
}

PyObject* Surface_GetComponent(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(self);
  if (!parSurf)
    return NULL;

  if(parSurf->dimension() != 3) {
    PyErr_SetString(PyExc_ValueError, "Not a three-dimensional surface");
    return NULL;
  }

  shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(parSurf);
  if(!spSurf) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineSurface");
    return NULL;
  }
  
  int nCoefs = spSurf->numCoefs_u() * spSurf->numCoefs_v();
  if(i < 0 || i >= nCoefs)
    return NULL;

  double x,y,z,w;
  int dim = spSurf->dimension();
  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  if(spSurf->rational()) {
    vector<double>::const_iterator cp = spSurf->rcoefs_begin() + i*(dim+1);
    result->data.reset(new Go::Point(cp, cp+(dim+1)));
  } else {
    x = *(spSurf->coefs_begin() + i*(spSurf->dimension() )+0) ;
    y = *(spSurf->coefs_begin() + i*(spSurf->dimension() )+1) ;
    z = *(spSurf->coefs_begin() + i*(spSurf->dimension() )+2) ;
    result->data.reset(new Go::Point(x,y,z));
  }

  return (PyObject*) result;
}

PyNumberMethods Surface_operators = {0};
PySequenceMethods Surface_seq_operators = {0};

PyDoc_STRVAR(surface__doc__,
             "A parametric description of a surface");
void init_Surface_Type()
{
  Surface_seq_operators.sq_item = Surface_GetComponent;
  Surface_seq_operators.sq_length = Surface_NmbComponent;
  InitializeTypeObject(&Surface_Type);
  Surface_operators.nb_add = Surface_Add;
  Surface_operators.nb_inplace_multiply = Surface_Scale;
  Surface_Type.tp_name = "GoTools.Surface";
  Surface_Type.tp_basicsize = sizeof(Surface);
  Surface_Type.tp_dealloc = (destructor)Surface_Dealloc;
  Surface_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  Surface_Type.tp_doc = surface__doc__;
  Surface_Type.tp_methods = Surface_methods;
  Surface_Type.tp_base = 0;
  Surface_Type.tp_new = Surface_New;
  Surface_Type.tp_str = (reprfunc)Surface_Str;
  Surface_Type.tp_as_number = &Surface_operators;
  Surface_Type.tp_as_sequence = &Surface_seq_operators;
  PyType_Ready(&Surface_Type);
}

}

//! \brief Check if a surface is a spline surface
//! \return The surface itself if it s a spline surface,
//          the surface converted to a spline surface otherwise.
shared_ptr<Go::SplineSurface> convertSplineSurface(shared_ptr<Go::ParamSurface> parSurf)
{
  if (!parSurf)
    return shared_ptr<Go::SplineSurface>();

  if (parSurf->instanceType() == Go::Class_SplineSurface)
    return dynamic_pointer_cast<Go::SplineSurface, Go::ParamSurface>(parSurf);
  if (parSurf->instanceType() == Go::Class_SurfaceOnVolume)
    return shared_ptr<Go::SplineSurface>(parSurf->asSplineSurface());
  shared_ptr<Go::ElementarySurface> e_surface = dynamic_pointer_cast<Go::ElementarySurface, Go::ParamSurface>(parSurf);
  if (e_surface)
    return shared_ptr<Go::SplineSurface>(e_surface->geometrySurface());
  return shared_ptr<Go::SplineSurface>();
}

void printSurfaceToStream(std::ostream& str, shared_ptr<Go::ParamSurface> parSurf)
{
  if (parSurf) {
    // fetch GoTools class type 
    std::string       sClassName  = Go::GoTools::className(parSurf->instanceType());
    std::stringstream ssClassName ;
    ssClassName << sClassName;
    if(sClassName.compare("Unknown") == 0)
      ssClassName << "(" << parSurf->instanceType() << ")";
    
    // print class type
    str << "Type: " << ssClassName.str() << std::endl;
    // print full raw data
    str << *parSurf;
  } else
    str << "(empty)";
}

void WriteSurfaceG2(std::ofstream& g2_file, Surface* pySurf, bool convert)
{
  if (!pySurf->data)
    return;
  if (convert) {
    shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(pySurf->data);
    if (spSurf) {
      spSurf->writeStandardHeader(g2_file);
      spSurf->write(g2_file);
    }
  } else {
    pySurf->data->writeStandardHeader(g2_file);
    pySurf->data->write(g2_file);
  }
}

}
