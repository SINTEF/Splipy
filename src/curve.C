#include "undef.h"
#include "curve.h"

#include "geomodeller.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/utils/LUDecomp.h"
#include "GoTools/utils/LUDecomp_implementation.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

#include <fstream>
#include <sstream>

namespace GeoModeller {

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

    int dim = coefs.size()/(knots.size()-order);
    if (rational)
      dim--;
    ((Curve*)self)->data.reset(new Go::SplineCurve(knots.size()-order, order,
                                                   knots.begin(), coefs.begin(),
                                                   dim, rational));
  }

  Py_INCREF(self);
  return (PyObject*)self;
}

void Curve_Dealloc(Curve* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyDoc_STRVAR(curve_append_curve__doc__,
             "Merge another curve with this one, with possible reparametrization\n"
             "@param curve: The curve to append\n"
             "@type curve: Curve\n"
             "@param continuity: (optional) Required continuity\n"
             "@type continuity: int\n"
             "@param reparam: (optional) Specify whether or not there should be reparametrization\n"
             "@type reparam: bool\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_AppendCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"curve", "continuity", "reparam", NULL };
  int continuity = 0;
  bool reparam   = true;
  PyObject *oCrv;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|ib",
                                   (char**)keyWords, &oCrv, &continuity, &reparam))
    return NULL;

  shared_ptr<Go::ParamCurve> parCrv2 = PyObject_AsGoCurve(self);
  shared_ptr<Go::ParamCurve> parCrv;
  if (parCrv2->instanceType() != Go::Class_SplineCurve) {
    if (!PyErr_WarnEx(PyExc_RuntimeWarning, "Converting parametric curve to spline curve", 1) == -1)
      return NULL;
    parCrv = ((Curve*)self)->data = convertSplineCurve(parCrv2);
  }
  else
    parCrv = dynamic_pointer_cast<Go::SplineCurve,Go::ParamCurve>(parCrv2);

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(PyObject_AsGoCurve(oCrv));
  if (!parCrv || !spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }

  double dist;
  parCrv->appendCurve(spCrv.get(), continuity, dist, reparam);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_clone__doc__,
             "Clone a curve\n"
             "@param coefs: Coefficients for new curve\n"
             "@type coefs: List of (list of float)\n"
             "@param ncomp: Number of components in coefficients\n"
             "@type ncomp: Integer\n"
             "@return: New (copy of) curve\n"
             "@rtype: Curve\n");
PyObject* Curve_Clone(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"coefs", "ncomp", NULL};
  Curve* res = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);

  PyObject* coefso=NULL;
  int dim = 1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|Oi",
                                   (char**)keyWords,&coefso,&dim))
    return NULL;

  if (coefso) {
    shared_ptr<Go::SplineCurve> curv = convertSplineCurve(parCrv);
    int nCoefs = curv->numCoefs();
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
    res->data.reset(new Go::SplineCurve(curv->basis(), coefs.begin(), dim, false));
  } else
    res->data.reset(parCrv->clone());
 
  return (PyObject*)res;
}

PyObject* Curve_Str(Curve* self)
{
  std::stringstream str;
  if (self->data) {
    // fetch GoTools class type 
    std::string       sClassName  = Go::GoTools::className(self->data->instanceType());
    std::stringstream ssClassName ;
    ssClassName << sClassName;
    if(sClassName.compare("Unknown") == 0)
      ssClassName << "(" << self->data->instanceType() << ")";
    
    // print class type
    str << "Type: " << ssClassName.str() << std::endl;
    // print full raw data
    str << *self->data;
  } else
    str << "(empty)";
  return PyString_FromString(str.str().c_str());
}

PyDoc_STRVAR(curve_evaluate__doc__,
             "Evaluate curve at a parameter value\n"
             "@param value: The parameter value\n"
             "@type value: float\n"
             "@param derivatives: The number of derivatives to obtain\n"
             "@type derivatives: Integer\n"
             "@return: The value of the curve\n"
             "@rtype: Point or tuple of Points if derivs > 0");
PyObject* Curve_Evaluate(PyObject* self, PyObject* args, PyObject* kwds)
{
  try {
    static const char* keyWords[] = {"value", "derivatives", NULL };
    shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
    double value=0;
    int derivs=0;
    if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"d|i",
                                     (char**)keyWords,&value,&derivs) || !parCrv)
      return NULL;

    if (derivs == 0) {
      Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
      result->data.reset(new Go::Point(parCrv->dimension()));
      parCrv->point(*result->data,value);
      return (PyObject*)result;
    } else {
      std::vector<Go::Point> pts((derivs+1));
      parCrv->point(pts, value, derivs);
      PyObject* result = PyTuple_New(pts.size());
      for (size_t i=0;i<derivs+1;++i) {
        Point* pt = (Point*)Point_Type.tp_alloc(&Point_Type,0);
        if(parCrv->dimension() == 2)
          pt->data.reset(new Go::Point(pts[i][0], pts[i][1]));
        else
          pt->data.reset(new Go::Point(pts[i][0], pts[i][1], pts[i][2]));
        PyTuple_SetItem(result, i, (PyObject*)pt);
      }
      return result;
    }

    Py_INCREF(Py_None);
    return Py_None;
  } catch(std::exception e) {
    PyErr_SetString(PyExc_Exception, e.what());
    return NULL;
  }
}

PyDoc_STRVAR(curve_evaluategrid__doc__,
             "Evaluate curve on a grid\n"
             "@param params: The parameters\n"
             "@type params: List of float\n"
             "@return: The values of the curve at the given parameters\n"
             "@rtype: List of floats");
PyObject* Curve_EvaluateGrid(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"params", NULL};
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  shared_ptr<Go::SplineCurve> crv = convertSplineCurve(parCrv);
  PyObject* paramso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&paramso) || !parCrv)
    return NULL;

  if (!PyObject_TypeCheck(paramso,&PyList_Type)) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for params: expected list");
    return NULL;
  }

  std::vector<double> params;
  for (int i=0;i<PyList_Size(paramso);++i)
    params.push_back(PyFloat_AsDouble(PyList_GetItem(paramso,i)));

  std::vector<double> res(params.size()*crv->dimension());
  crv->gridEvaluator(res, params);

  PyObject* result = PyList_New(res.size()/crv->dimension());
  VectorToPyPointList(result, res, crv->dimension());

  return result;
}

PyDoc_STRVAR(curve_evaluate_tangent__doc__,
             "Evaluate the curve tangent at a parameter value\n"
             "@param value: The parameter value\n"
             "@type value: float\n"
             "@return: The tangent of the curve\n"
             "@rtype: Point");
PyObject* Curve_EvaluateTangent(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"value", NULL };
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  double value=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"d",
                                   (char**)keyWords,&value) || !parCrv)
    return NULL;

  std::vector<Go::Point> pts(2);
  parCrv->point(pts, value, 1);

  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(pts[1]));

  return (PyObject*)result;
}

PyDoc_STRVAR(curve_flip_parametrization__doc__,
             "Flip curve parametrization\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_FlipParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  parCrv->reverseParameterDirection();

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_force_rational__doc__,
             "Enforce a rational representation of the spline curve\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_ForceRational(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if(!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }
  ((Curve*)self)->data = spCrv;

  spCrv->representAsRational();

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_lower_order__doc__,
             "Lower the order of a spline curve (need full continuity)\n"
             "@param lower: Lower of order\n"
             "@type lower: int (>= 0)\n"
             "@returns The curve\n"
             "@rtype: Curve");
PyObject* Curve_LowerOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"lower", NULL};
  int lower = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwds, (char*)"i", (char**)keyWords, &lower))
    return NULL;

  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;
  if (!parCrv->geometryCurve()) {
    Curve* pyCrv = (Curve*)self;
    pyCrv->data = convertSplineCurve(parCrv);
    parCrv = pyCrv->data;
  }

  shared_ptr<Go::SplineCurve> spCrv = static_pointer_cast<Go::SplineCurve>(parCrv);
  std::vector<double>::const_iterator first = spCrv->basis().begin() + lower;
  std::vector<double>::const_iterator last  = spCrv->basis().end() - lower;
  Go::BsplineBasis basis = Go::BsplineBasis(spCrv->order() - lower, first, last);

  if (spCrv->rational()) {
    if (PyErr_WarnEx(PyExc_RuntimeWarning,
                     "The geometry basis is rational (using NURBS). "
                     "The basis for the unknown fields of one degree "
                     "higher will however be non-rational. "
                     "This may affect accuracy.", 1) == -1)
      return NULL;
  }

  std::vector<double> greville(basis.numCoefs());
  for (size_t i = 0; i < greville.size(); i++)
    greville[i] = basis.grevilleParameter(i);

  // Evaluate the spline curve at all points
  std::vector<double> XYZ(spCrv->dimension() * greville.size());
  spCrv->gridEvaluator(XYZ, greville);
  
  // Project the coordinates onto the new basis (the second XYZ is dummy here)
  Curve* pyCrv = (Curve*)self;
  pyCrv->data = convertSplineCurve(parCrv);
  pyCrv->data.reset(Go::CurveInterpolator::regularInterpolation(basis, greville, XYZ,
                                                                spCrv->dimension(), false, XYZ));

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_get_greville__doc__,
             "Return Greville points for a spline curve\n"
             "@return: List with the parameter values\n"
             "@rtype: List of float");
PyObject* Curve_GetGreville(PyObject* self)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;
  if (!parCrv->geometryCurve()) {
    Curve* pyCrv = (Curve*)self;
    pyCrv->data = convertSplineCurve(parCrv);
    parCrv = pyCrv->data;
  }

  shared_ptr<Go::SplineCurve> spCrv = static_pointer_cast<Go::SplineCurve>(parCrv);
  PyObject* result = PyList_New(0);
  for (size_t j=0;j<spCrv->basis().numCoefs();++j)
    PyList_Append(result,Py_BuildValue((char*)"d",
                   spCrv->basis().grevilleParameter(j)));

  return result;
}

PyDoc_STRVAR(curve_get_kinks__doc__,
             "Get the parametric value of the curve kinks (where the knot vector has multiplicity p and the curve is C0). Set 'knot' tolerance for more control.\n"
             "@return: List of the knot values\n"
             "@rtype: List of float");
PyObject* Curve_GetKinks(PyObject* self)
{

  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;
  Curve* pyCrv = (Curve*)self;
  pyCrv->data = convertSplineCurve(parCrv);
  PyObject* result = PyList_New(0);
  std::vector<double> knots;
  shared_ptr<Go::SplineCurve> spCrv = static_pointer_cast<Go::SplineCurve>(pyCrv->data);
  knots = spCrv->basis().getKnots();
  int p = spCrv->basis().order();
  int n = knots.size();
  int multiplicity = 0;
  double prevKnot = knots[0];
  for(int i=1; i<knots.size()-1; i++) {
    if(fabs(knots[i]-prevKnot) < modState.knotTolerance) {
      multiplicity++;
    } else {
      multiplicity = 1;
      prevKnot = knots[i];
    }
    if(multiplicity == p-1)
      PyList_Append(result,Py_BuildValue((char*)"d",prevKnot));
  }
  return result;
}

PyDoc_STRVAR(curve_get_knots__doc__,
             "Get the knots of a spline curve\n"
             "@param with_multiplicities: (optional) Set to true to obtain the knot vector with multiplicities\n"
             "@type with_multiplicities: Boolean\n"
             "@return: List with the knot values\n"
             "@rtype: List of float");
PyObject* Curve_GetKnots(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"with_multiplicities", NULL };

  bool withmult=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|b",
                                   (char**)keyWords, &withmult))
    return NULL;
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;
  Curve* pyCrv = (Curve*)self;
  pyCrv->data = convertSplineCurve(parCrv);
  PyObject* result = PyList_New(0);
  std::vector<double> knots;
  shared_ptr<Go::SplineCurve> spCrv = static_pointer_cast<Go::SplineCurve>(pyCrv->data);
  if (withmult) 
    knots = spCrv->basis().getKnots();
  else
    spCrv->basis().knotsSimple(knots);
  for (std::vector<double>::iterator it  = knots.begin();
                                     it != knots.end();++it) {
    PyList_Append(result,Py_BuildValue((char*)"d",*it));
  }
                
  return result;
}

PyDoc_STRVAR(curve_get_order__doc__,
             "Get the curve order (polynomial degree + 1)\n"
             "@return: Order\n"
             "@rtype: Integer");
PyObject* Curve_GetOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if(!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }

  return Py_BuildValue((char*)"i",spCrv->order());
}

PyDoc_STRVAR(curve_get_parameter_at_point__doc__,
             "Get all Curve parameter values at a geometric Point\n"
             "@param point: The geometric point to intersect \n"
             "@type point: Point, list of floats or tuple of floats\n"
             "@return: Parameter values of all intersection points\n"
             "@rtype: List of float");
PyObject* Curve_GetParameterAtPoint(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"point", NULL };
  PyObject* pointo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,
                                   &pointo))
    return NULL;

  // get Point from Python object
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(pointo);
  if(!point) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for point: expected pointlike");
    return NULL;
  }

  // get ParamCurve from Python object
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  // get SplineCurve from ParamCurve
  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if (!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }

  // get SISL curve from SplineCurve
  SISLCurve *sislCrv = Curve2SISL(*spCrv, true);
  if (!sislCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain SISLCurve");
    return NULL;
  }

  // setup SISL parameters
  double *pt      = &((*point)[0]);        // pointer to the data of the geometric point
  int    pointDim = point->dimension();    // dimension of the geometric point
  double epsge    = modState.gapTolerance; // Geometry resolution
  int    numPts;          // number of intersection points
  int    numCrv;          // number of intersection curves (Wait, what? intersection CURVES?)
  double *parPts;         // array of parameter-values at intersection points
  SISLIntcurve **intCrvs; // array of intersection curves  (I have no idea how you get these guys.
                          // We're silently ignoring and hoping for the best)
  int    status;          // errors/warnings from evaluation


  s1871(sislCrv, pt, pointDim, epsge, // input arguments
        &numPts, &parPts,         // output points
        &numCrv, &intCrvs,        // output curves
        &status);                 // output errors

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
  PyObject* result = PyList_New(0);
  std::vector<double> knots;
  for (int i=0; i<numPts; i++) {
    PyList_Append(result,Py_BuildValue((char*)"d",parPts[i]));
  }
                
  return result;
}

PyDoc_STRVAR(curve_get_tesselationparams__doc__,
             "Obtain tesselation parameters for a curve\n"
             "@param n: Number of tesselation points per knotspan\n"
             "@type n: Int\n"
             "@rtype: List of float");
PyObject* Curve_GetTesselationParams(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  shared_ptr<Go::SplineCurve> crv = convertSplineCurve(parCrv);
  static const char* keyWords[] = {"n", NULL};
  int np = 1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|i",
                                   (char**)keyWords,&np))
    return NULL;

  // Grab parameter values in evaluation points
  std::vector<double> gpar = Tesselate(crv->basis().begin(), crv->basis().end(), np);

  size_t nx = gpar.size();
  PyObject* result = PyList_New(nx);
  for (size_t i=0;i<nx;++i)
    PyList_SetItem(result, i, Py_BuildValue((char*)"d",gpar[i]));

  return result;
}

PyDoc_STRVAR(curve_insert_knot__doc__,
             "Insert a knot into a spline curve\n"
             "@param knot: The knot to insert\n"
             "@type knot: float\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"knot", NULL };
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  double knot;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"d",
                                   (char**)keyWords,&knot) || !parCrv)
    return NULL;

  Curve* pyCrv = (Curve*)self;
  pyCrv->data = convertSplineCurve(parCrv);

  static_pointer_cast<Go::SplineCurve>(pyCrv->data)->insertKnot(knot);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_interpolate__doc__,
             "Interpolate a field onto the curve's basis\n"
             "@param coefs: The values to interpolate (sample in the Greville points)\n"
             "@type coefs: List of Float"
             "@return: New basis coefficients\n"
             "@rtype: List of float");
PyObject* Curve_Interpolate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"coefs", NULL };
  PyObject* pycoefs=NULL;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&pycoefs))
    return NULL;

  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  if (!PyObject_TypeCheck(pycoefs,&PyList_Type)) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for coefs: expected list");
    return NULL;
  }

  if (!parCrv->geometryCurve()) {
    Curve* pyCrv = (Curve*)self;
    pyCrv->data = convertSplineCurve(parCrv);
    parCrv = pyCrv->data;
  }
  shared_ptr<Go::SplineCurve> crv = convertSplineCurve(parCrv);

  std::pair<std::vector<double>, int> coefs = PyPointListToVector(pycoefs);

  std::vector<double> greville_u(crv->basis().numCoefs());
  for(int i=0; i<crv->basis().numCoefs(); i++)
    greville_u[i] = crv->basis().grevilleParameter(i);
  std::vector<double> weights(0);

  Go::SplineCurve* res =
        Go::CurveInterpolator::regularInterpolation(crv->basis(),
                                                    greville_u,
                                                    coefs.first,
                                                    coefs.second,
                                                    false,
                                                    weights);

  PyObject* result = PyList_New(0);
  for (std::vector<double>::const_iterator it  = res->coefs_begin();
                                           it != res->coefs_end();++it)
    PyList_Append(result,Py_BuildValue((char*)"d",*it));

  delete res;

  return result;
}

PyDoc_STRVAR(curve_intersect__doc__,
             "Check if this curve intersects another curve or surface.\n"
             "Set tolerance 'gap' for intersection tolerance \n"
             "@param obj: The object to test against\n"
             "@type obj: Curve or Surface\n"
             "@return: One intersection point (if any)\n"
             "@rtype: Point or None");
PyObject* Curve_Intersect(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"obj", NULL };
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  double knot;
  PyObject *pyObj;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&pyObj) || !parCrv)
    return NULL;

  if(PyObject_TypeCheck(pyObj, &Curve_Type) ) {
    shared_ptr<Go::ParamCurve> parCrv2 = PyObject_AsGoCurve(pyObj);
    if(!parCrv2) {
      PyErr_SetString(PyExc_TypeError, "Unable to treat obj as curve");
      return NULL;
    }
    double par1, par2, dist;
    Go::Point ptc1, ptc2;
    Go::ClosestPoint::closestPtCurves(parCrv.get(), parCrv2.get(), par1, par2, dist, ptc1, ptc2);
    if(dist < modState.gapTolerance) {
      Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
      result->data.reset(new Go::Point(ptc1));
      return (PyObject*) result;
    } else {
      Py_INCREF(Py_None);
      return Py_None;
    }
  } else if(PyObject_TypeCheck(pyObj, &Surface_Type) ) {
    shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(pyObj);
    if(!parSurf) {
      PyErr_SetString(PyExc_TypeError, "Unable to treat obj as surface");
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
      Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
      result->data.reset(new Go::Point(crvPt));
      return (PyObject*) result;
    } else {
      Py_INCREF(Py_None);
      return Py_None;
    }
  }

  PyErr_SetString(PyExc_TypeError, "Invalid type for obj: expected Curve or Surface");
  return NULL;
}

PyDoc_STRVAR(curve_normalize__doc__,
             "Normalize a curve in the parameter domain\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_Normalize(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  Curve* pyCrv = (Curve*)self;
  pyCrv->data = convertSplineCurve(parCrv);

  pyCrv->data->setParameterInterval(0,1);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_project__doc__,
             "Project the curve onto an axis or plane along parallel to the cartesian coordinate system\n"
             "@param axis: The axis or plane to project onto (\"X\",\"Y\",\"Z\" or a comibation of these)\n"
             "@type axis: string\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_Project(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", NULL };
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  char *sAxis;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&sAxis) || !parCrv)
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

  Curve* pyCrv = (Curve*)self;
  pyCrv->data = convertSplineCurve(parCrv);

  shared_ptr<Go::SplineCurve>   spCrv    = static_pointer_cast<Go::SplineCurve>(pyCrv->data);
  bool                          rational = spCrv->rational();
  int                           dim      = spCrv->dimension();
  std::vector<double>::iterator coefs    = (rational) ? spCrv->rcoefs_begin() : spCrv->coefs_begin();
  std::vector<double>::iterator coefsEnd = (rational) ? spCrv->rcoefs_end()   : spCrv->coefs_end();
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

PyDoc_STRVAR(curve_raise_order__doc__,
             "Raise the order of the curve's B-spline basis without changing the shape of the curve\n"
             "@param n: Specifies how many times the order will be raised\n"
             "@type n: int\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  try {
    static const char* keyWords[] = {"n", NULL };
    shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
    int amount;

    if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                     (char**)keyWords,&amount) || !parCrv)
      return NULL;

    Curve* pyCrv = (Curve*)self;
    pyCrv->data = convertSplineCurve(parCrv);

    static_pointer_cast<Go::SplineCurve>(pyCrv->data)->raiseOrder(amount);

    Py_INCREF(self);
    return self;
  } catch(std::exception e) {
    PyErr_SetString(PyExc_Exception, e.what());
    return NULL;
  }
}

PyDoc_STRVAR(curve_reparametrize__doc__,
             "Re-parametrize a curve\n"
             "@param umin: The minimum u value\n"
             "@type umin: float\n"
             "@param umax: The maximum u value\n"
             "@type umax: float\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_ReParametrize(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"umin", "umax", NULL };
  double umin=0, umax=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|dd",
                                   (char**)keyWords,&umin,&umax))
    return NULL;
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if(!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }
  if(umin >= umax) {
    PyErr_SetString(PyExc_ValueError, "Invalid parameter range: requires umin < umax");
    return NULL;
  }

  spCrv->setParameterInterval(umin, umax);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_rebuild__doc__,
             "Rebuild the curve by resampling it to a given order and number of control points\n"
             "The rebuilt curve will match the old one at the end-points, the end tangents and \n"
             "n-4 internal uniformly distributed points \n"
             "@param n: Specifies how many control points the resulting curve should have\n"
             "@type n: int\n"
             "@param p: Specifies the resulting order (polynomial degree +1) of the curve\n"
             "@type p: int\n"
             "@return: Rebuilt curve\n"
             "@rtype: Curve");
PyObject* Curve_Rebuild(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"n", "p", NULL };
  int n,p;
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  int amount;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords, &n, &p))
    return NULL;

  if (n < p || p < 1 || n-4 < 0) {
    PyErr_SetString(PyExc_ValueError, "Invalid argument values (requires n >= max(p,4), p > 0)");
    return NULL;
  }

  // get a few needed curve values
  int    dim = parCrv->dimension();
  double u0  = parCrv->startparam();
  double u1  = parCrv->endparam();

  // set boundary conditions (match tangent at start/end)
  std::vector<Go::Point> ptsStart(2);
  std::vector<Go::Point> ptsStop(2);
  parCrv->point(ptsStart, u0, 1, true );
  parCrv->point(ptsStop,  u1, 1, false);
  shared_ptr<Go::Point> p0(new Go::Point(ptsStart[1]));
  shared_ptr<Go::Point> p1(new Go::Point(ptsStop[1]));
  *p0 *= (u1-u0)/(n-p+1); // scale tangents to new parameterization
  *p1 *= (u1-u0)/(n-p+1);


  // build the uniform knot vector
  std::vector<double> knots(n+p);
  int k = 0;
  for(int i=0; i<p; i++) 
    knots[k++] = 0;
  for(int i=0; i<n-p; i++) 
    knots[k++] = i+1;
  for(int i=0; i<p; i++) 
    knots[k++] = n-p+1;
  Go::BsplineBasis basis(n, p, knots.begin());

  // generate the parametric interpolation points
  int nPts = n;
  if(p0 != NULL) nPts--; // keeping backward compatible if one doesn't want 
  if(p1 != NULL) nPts--; // derivatives specified at start/end
  std::vector<double> paramOld(nPts);  // parameter values as seen from *curve
  std::vector<double> paramNew(nPts);  // interpolation points of rebuilt curve 
  for(int i=0; i<nPts; i++) {
    paramOld[i] = u0 + i*(u1-u0)/(nPts-1);
    paramNew[i] = i*((double) n-p+1)/(nPts-1);
  }

  // set up to solve the interpolation problem
  std::vector<std::vector<double> > A(n, std::vector<double>(n, 0));
  std::vector<std::vector<double> > b(n, std::vector<double>(dim, 0));
  
  // assembling up interpolation matrix A
  std::vector<double> tmp(2*p);
  k=0; // row iterator over matrix A
  if(p0 != NULL) {
      basis.computeBasisValues(paramNew[0], &tmp[0], 1);
      for (int j = 0; j < p; ++j)
        A[k][j] = tmp[2*j+1];
      k++;
  }
  for (int i = 0; i < nPts; ++i, ++k) {
      double par = paramNew[i];
      int    ki  = basis.knotIntervalFuzzy(par); // knot-interval of param.
      basis.computeBasisValues(paramNew[i], &tmp[0], 0);

      for (int j = 0; j < p; ++j)
          if ((ki-p+1+j>=0) && (ki-p+1+j<=n))
              A[k][ki-p+1+j] = tmp[j];
  }
  if(p1 != NULL) {
      basis.computeBasisValues(paramNew.back(), &tmp[0], 1);
      for (int j = 0; j < p; ++j)
        A[k][n-p+j] = tmp[2*j+1];
  }

  // generating right-hand side b by sampling interpolation points
  Go::Point evalPt;
  k = 0;
  if(p0 != NULL) {
    for(int j=0; j<dim; j++)
      b[k][j] = (*p0)[j];
    k++;
  }
  for(int i=0; i<nPts; i++, k++) {
    evalPt = parCrv->point(paramOld[i]);
    for(int j=0; j<dim; j++)
      b[k][j] = evalPt[j];
  }
  if(p1 != NULL)
    for(int j=0; j<dim; j++)
      b[k][j] = (*p1)[j];


  // Now we are ready to solve Ac = b.  b will be overwritten by solution
  Go::LUsolveSystem(A, n, &b[0]);

  // copy results
  std::vector<double> coefs(n*dim);
  for (int i = 0; i < n; ++i) {
      std::copy(b[i].begin(), b[i].end(), &coefs[i * dim]);
  }

  // wrap results in a python object and return
  Curve* pyCrv = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
  pyCrv->data.reset(new Go::SplineCurve(basis, coefs.begin(), dim));

  return (PyObject*)pyCrv;
}

PyDoc_STRVAR(curve_rotate__doc__,
             "Rotate a curve around an axis\n"
             "@param axis: The axis to rotate around\n"
             "@type axis: Point, list of floats or tuple of floats\n"
             "@param angle: Angle to rotate curve with in radians\n"
             "@type angle: float\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_Rotate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"axis", "angle", NULL };
  PyObject* axiso;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  if (!axis) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for axis: expected pointlike");
    return NULL;
  }

   Curve* pyCurve = (Curve*)self;
   pyCurve->data = convertSplineCurve(parCrv);
   parCrv = pyCurve->data;

   Go::GeometryTools::rotateSplineCurve(*axis, angle,
                        *static_pointer_cast<Go::SplineCurve>(parCrv));

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(curve_split__doc__,
             "Split the curve into segments\n"
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

  if (p.empty()) {
    PyErr_SetString(PyExc_ValueError, "No valid parameter values found");
    return NULL;
  }

  std::vector<shared_ptr<Go::ParamCurve> > curves;
  if (p.size() > 1) {
    shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(((Curve*)self)->data);
    std::vector<shared_ptr<Go::SplineCurve> > curves2 = spCrv->split(p);
    for (size_t i=0;i<curves2.size();++i)
      curves.push_back(static_pointer_cast<Go::ParamCurve,Go::SplineCurve>(curves2[i]));
  } else
    curves = ((Curve*)self)->data->split(p[0]);

  PyObject* result = PyList_New(0);
  for (size_t i=0;i<curves.size();++i) {
    Curve* pyCrv = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
    pyCrv->data = curves[i];
    PyList_Append(result,(PyObject*)pyCrv);
  }

  return result;
}

PyDoc_STRVAR(curve_get_sub_curve__doc__,
             "Get a Curve which represent a part of 'this' Curve\n"
             "@param from_par: The parametric left end of the sub curve\n"
             "@type  from_par: float\n"
             "@param to_par:   The upper right end of the sub curve\n"
             "@type  to_par:   float\n"
             "@return: A parametric subcurve of this Curve represented as a Spline Curve\n"
             "@rtype: Curve");
PyObject* Curve_GetSubCurve(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  static const char* keyWords[] = {"from_par", "to_par", NULL };
  double Left=0.0, Right=1.0;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"dd",
                                   (char**)keyWords, &Left, &Right) || !parCrv)
    return NULL;

  if (!parCrv->geometryCurve()) {
    Curve* pyCrv = (Curve*)self;
    pyCrv->data = convertSplineCurve(parCrv);
    parCrv = pyCrv->data;
  }
  shared_ptr<Go::SplineCurve> spCrv = static_pointer_cast<Go::SplineCurve>(parCrv);
  if(!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }

  Go::SplineCurve *subCrv = spCrv->subCurve(Left, Right);

  Curve* result = (Curve*)Surface_Type.tp_alloc(&Curve_Type,0);
  result->data = shared_ptr<Go::ParamCurve>(subCrv);

  return (PyObject*) result;
}

PyDoc_STRVAR(curve_tesselate__doc__,
             "Tesselate a curve\n"
             "@param n: Number of tesselation points per knotspan\n"
             "@type n: Int\n"
             "@rtype: Tuple with (list of nodes, list of elements)");
PyObject* Curve_Tesselate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"n", NULL};
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  shared_ptr<Go::SplineCurve> crv = convertSplineCurve(parCrv);
  int np = 1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|i",
                                   (char**)keyWords,&np))
    return NULL;

  // Grab parameter values in evaluation points
  std::vector<double> gpar = Tesselate(crv->basis().begin(), crv->basis().end(), np);

  // Evaluate the surface at all points
  size_t nx = gpar.size();
  std::vector<double> XYZ(crv->dimension()*nx);
  crv->gridEvaluator(XYZ,gpar);

  // Establish the block grid coordinates
  PyObject* gc = PyList_New(XYZ.size());
  for (int i=0;i<XYZ.size();++i)
    PyList_SetItem(gc, i, Py_BuildValue((char*)"d", XYZ[i]));


  // Establish the block grid topology
  PyObject* ge = PyList_New(0);
  int nse1 = np;
  int n[2], ie = 1, ip = 0;
  n[0] = 0;
  n[1] = n[0] + 1;
  size_t i, l;
  for (i = 1; i < nx; i++)
  {
    for (l = 0; l < 2; l++)
      PyList_Append(ge, Py_BuildValue((char*)"i",n[l]++));
    if (i%nse1 == 0) ie++;
  }

  PyObject* result = PyTuple_New(2);
  PyTuple_SetItem(result, 0, gc);
  PyTuple_SetItem(result, 1, ge);

  return result;
}


PyDoc_STRVAR(curve_translate__doc__,
             "Translate a curve along a given vector\n"
             "@param vector: The vector to translate along\n"
             "@type axis: Point, list of floats or tuple of floats\n"
             "@return: The curve\n"
             "@rtype: Curve");
PyObject* Curve_Translate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"vector", NULL };
  PyObject* veco;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&veco))
    return NULL;

  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  shared_ptr<Go::Point> vec = PyObject_AsGoPoint(veco);
  if (!vec) {
    PyErr_SetString(PyExc_TypeError, "Invalid type for vector: expected pointlike");
    return NULL;
  }

  if (parCrv->geometryCurve() != NULL) {
    Curve* pyCrv = (Curve*)self;
    pyCrv->data = convertSplineCurve(parCrv);
    parCrv = pyCrv->data;
  }

  Go::GeometryTools::translateSplineCurve(*vec, *static_pointer_cast<Go::SplineCurve>(parCrv));

  Py_INCREF(self);
  return self;
}

PyObject* Curve_Reduce(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if (!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::SplineCurve");
    return NULL;
  }

  PyObject* knots = PyList_New(0);
  vector<double>::const_iterator kit;
  for (kit = spCrv->knotsBegin(); kit != spCrv->knotsEnd(); ++kit)
    PyList_Append(knots, Py_BuildValue("d", *kit));

  PyObject* coefs = PyList_New(0);
  vector<double>::const_iterator cit = spCrv->rational() ? spCrv->rcoefs_begin() : spCrv->coefs_begin();
  vector<double>::const_iterator cit_end = spCrv->rational() ? spCrv->rcoefs_end() : spCrv->coefs_end();
  for (; cit != cit_end; ++cit)
    PyList_Append(coefs, Py_BuildValue("d", *cit));
  
  return Py_BuildValue("(O(iOOO))",
                       &Curve_Type, spCrv->order(), knots, coefs,
                       spCrv->rational() ? Py_True : Py_False);
}

PyObject* Curve_Add(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(o1);
  shared_ptr<Go::Point> point = PyObject_AsGoPoint(o2);
  Curve* pyCrv = NULL;

  if (parCrv && point) {
    pyCrv = (Curve*)Curve_Type.tp_alloc(&Curve_Type,0);
    if (parCrv->instanceType() == Go::Class_SplineCurve)
      pyCrv->data = shared_ptr<Go::ParamCurve>(
          new Go::SplineCurve(*static_pointer_cast<Go::SplineCurve>(parCrv)));
    else
      pyCrv->data = convertSplineCurve(parCrv);

    Go::GeometryTools::translateSplineCurve(*point, 
                         *static_pointer_cast<Go::SplineCurve>(pyCrv->data));
  }

  if (!pyCrv) {
    PyErr_SetString(PyExc_TypeError, "Expected Curve + Point");
  }

  return (PyObject*) pyCrv;
}

Py_ssize_t Curve_NmbComponent(PyObject* self)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return 0;

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if(!spCrv)
    return 0;

  return spCrv->numCoefs();
}

PyObject* Curve_GetComponent(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::ParamCurve> parCrv = PyObject_AsGoCurve(self);
  if (!parCrv)
    return NULL;

  if(parCrv->dimension() != 3) {
    PyErr_SetString(PyExc_ValueError, "Not a three-dimensional curve");
    return NULL;
  }

  shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(parCrv);
  if(!spCrv) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtian Go::SplineCurve");
    return NULL;
  }

  if(i < 0 || i >= spCrv->numCoefs()) {
    PyErr_SetString(PyExc_IndexError, "Index out of bounds");
    return NULL;
  }

  int dim = spCrv->dimension();
  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  if (spCrv->rational()) {
    vector<double>::const_iterator cp = spCrv->rcoefs_begin() + i*(dim+1);
    result->data.reset(new Go::Point(cp, cp+(dim+1)));
  } else {
    double x = *(spCrv->coefs_begin() + i*dim + 0);
    double y = *(spCrv->coefs_begin() + i*dim + 1);
    double z = *(spCrv->coefs_begin() + i*dim + 2);
    result->data.reset(new Go::Point(x,y,z));
  }

  return (PyObject*) result;
}

PyMethodDef Curve_methods[] = {
     {(char*)"AppendCurve",         (PyCFunction)Curve_AppendCurve,         METH_VARARGS|METH_KEYWORDS, curve_append_curve__doc__},
     {(char*)"Clone",               (PyCFunction)Curve_Clone,               METH_VARARGS|METH_KEYWORDS, curve_clone__doc__},
     {(char*)"Evaluate",            (PyCFunction)Curve_Evaluate,            METH_VARARGS|METH_KEYWORDS, curve_evaluate__doc__},
     {(char*)"EvaluateGrid",        (PyCFunction)Curve_EvaluateGrid,        METH_VARARGS|METH_KEYWORDS, curve_evaluategrid__doc__},
     {(char*)"EvaluateTangent",     (PyCFunction)Curve_EvaluateTangent,     METH_VARARGS|METH_KEYWORDS, curve_evaluate_tangent__doc__},
     {(char*)"FlipParametrization", (PyCFunction)Curve_FlipParametrization, METH_VARARGS,               curve_flip_parametrization__doc__},
     {(char*)"ForceRational",       (PyCFunction)Curve_ForceRational,       METH_VARARGS,               curve_force_rational__doc__},
     {(char*)"GetGreville",         (PyCFunction)Curve_GetGreville,         0,                          curve_get_greville__doc__},
     {(char*)"GetKinks",            (PyCFunction)Curve_GetKinks,            0,                          curve_get_kinks__doc__},
     {(char*)"GetKnots",            (PyCFunction)Curve_GetKnots,            METH_VARARGS|METH_KEYWORDS, curve_get_knots__doc__},
     {(char*)"GetOrder",            (PyCFunction)Curve_GetOrder,            METH_VARARGS,               curve_get_order__doc__},
     {(char*)"GetParameterAtPoint", (PyCFunction)Curve_GetParameterAtPoint, METH_VARARGS|METH_KEYWORDS, curve_get_parameter_at_point__doc__},
     {(char*)"GetTesselationParams",(PyCFunction)Curve_GetTesselationParams,METH_VARARGS|METH_KEYWORDS, curve_get_tesselationparams__doc__},
     {(char*)"GetSubCurve",         (PyCFunction)Curve_GetSubCurve,         METH_VARARGS|METH_KEYWORDS, curve_get_sub_curve__doc__},
     {(char*)"InsertKnot",          (PyCFunction)Curve_InsertKnot,          METH_VARARGS|METH_KEYWORDS, curve_insert_knot__doc__},
     {(char*)"Interpolate",         (PyCFunction)Curve_Interpolate,         METH_VARARGS|METH_KEYWORDS, curve_interpolate__doc__},
     {(char*)"Intersect",           (PyCFunction)Curve_Intersect,           METH_VARARGS|METH_KEYWORDS, curve_intersect__doc__},
     {(char*)"Normalize",           (PyCFunction)Curve_Normalize,           METH_VARARGS,               curve_normalize__doc__},
     {(char*)"Project",             (PyCFunction)Curve_Project,             METH_VARARGS|METH_KEYWORDS, curve_project__doc__},
     {(char*)"RaiseOrder",          (PyCFunction)Curve_RaiseOrder,          METH_VARARGS|METH_KEYWORDS, curve_raise_order__doc__},
     {(char*)"LowerOrder",          (PyCFunction)Curve_LowerOrder,          METH_VARARGS|METH_KEYWORDS, curve_lower_order__doc__},
     {(char*)"ReParametrize",       (PyCFunction)Curve_ReParametrize,       METH_VARARGS|METH_KEYWORDS, curve_reparametrize__doc__},
     {(char*)"Rebuild",             (PyCFunction)Curve_Rebuild,             METH_VARARGS|METH_KEYWORDS, curve_rebuild__doc__},
     {(char*)"Rotate",              (PyCFunction)Curve_Rotate,              METH_VARARGS|METH_KEYWORDS, curve_rotate__doc__},
     {(char*)"Split",               (PyCFunction)Curve_Split,               METH_VARARGS|METH_KEYWORDS, curve_split__doc__},
     {(char*)"Tesselate",           (PyCFunction)Curve_Tesselate,           METH_VARARGS|METH_KEYWORDS, curve_tesselate__doc__},
     {(char*)"Translate",           (PyCFunction)Curve_Translate,           METH_VARARGS|METH_KEYWORDS, curve_translate__doc__},
     {(char*)"__reduce__",          (PyCFunction)Curve_Reduce,              0,                          NULL},
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

shared_ptr<Go::SplineCurve> convertSplineCurve(shared_ptr<Go::ParamCurve> parCrv)
{
  if (!parCrv)
    return shared_ptr<Go::SplineCurve>();

  if (parCrv->instanceType() == Go::Class_SplineCurve)
    return dynamic_pointer_cast<Go::SplineCurve, Go::ParamCurve>(parCrv);
  return shared_ptr<Go::SplineCurve>(parCrv->geometryCurve());
}

void WriteCurveG2(std::ofstream& g2_file, Curve* pyCurve, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineCurve> spCrv = convertSplineCurve(pyCurve->data);
    spCrv->writeStandardHeader(g2_file);
    spCrv->write(g2_file);
  } else {
    pyCurve->data->writeStandardHeader(g2_file);
    pyCurve->data->write(g2_file);
  }
}

}
