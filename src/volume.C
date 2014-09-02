#include "volume.h"

#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/trivariate/ElementaryVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
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

  static const char* keyWords[] = {"order1", "order2", "order3", "knots1", "knots2 ", "knots3", "coefs", "rational", NULL };
  PyObject* knots1o=0;
  PyObject* knots2o=0;
  PyObject* knots3o=0;
  PyObject* coefso=0;
  bool rational=false;
  int order1=0;
  int order2=0;
  int order3=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|iiiOOOOb",
                                   (char**)keyWords, &order1, &order2, &order3,
                                           &knots1o, &knots2o, &knots3o,
                                           &coefso, &rational))
    return NULL;

  // volume is specified
  if (knots1o && knots2o && knots3o && coefso) {
    std::vector<double> knots1;
    std::vector<double> knots2;
    std::vector<double> knots3;
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
    if (PyObject_TypeCheck(knots3o,&PyList_Type)) {
      for (int i=0;i<PyList_Size(knots3o);++i) {
        PyObject* o = PyList_GetItem(knots3o,i);
        if (o && PyObject_TypeCheck(o,&PyFloat_Type))
          knots3.push_back(PyFloat_AsDouble(o));
        else if (o && PyObject_TypeCheck(o,&PyInt_Type))
          knots3.push_back(PyInt_AsLong(o));
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
    int dim = coefs.size()/((knots1.size()-order1)*(knots2.size()-order2)*(knots3.size()-order3));
    if (rational)
      dim--;
    ((Volume*)self)->data.reset(new Go::SplineVolume(knots1.size()-order1,
                                                     knots2.size()-order2, 
                                                     knots3.size()-order3, 
                                                     order1, order2, order3,
                                                     knots1.begin(),
                                                     knots2.begin(),
                                                     knots3.begin(),
                                                     coefs.begin(), dim, rational));
  }

  return (PyObject*)self;
}

void Volume_Dealloc(Volume* self)
{
  self->data.reset();
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
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  res->data.reset(parVol->clone());
 
  return (PyObject*)res;
}

PyObject* Volume_Add(PyObject* o1, PyObject* o2)
{
  Volume* pyVol = (Volume*)o1;
  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::Point> p = PyObject_AsGoPoint(o2);
  if (p) {
    result->data.reset(pyVol->data->clone());
    result->data->translate(*p);
  }

  return (PyObject*)result;
}

PyObject* Volume_Sub(PyObject* o1, PyObject* o2)
{
  Volume* pyVol = (Volume*)o1;
  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::Point> p = PyObject_AsGoPoint(o2);
  if (p) {
    result->data.reset(pyVol->data->clone());
    result->data->translate(-(*p));
  }

  return (PyObject*)result;
}

PyDoc_STRVAR(volume_evaluate__doc__,"Evaluate volume at given parameter values\n"
                                     "@param value_u: The u parameter value\n"
                                     "@type value_u: float\n"
                                     "@param value_v: The v parameter value\n"
                                     "@type value_v: float\n"
                                     "@param value_w: The w parameter value\n"
                                     "@type value_w: float\n"
                                     "@return: The value of the volume\n"
                                     "@rtype: Point");
PyObject* Volume_Evaluate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"value_u", "value_v", "value_w", NULL };
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  double value_u=0, value_v=0, value_w=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ddd",
                                   (char**)keyWords,&value_u,&value_v,
                                                    &value_w) || !parVol)
    return NULL;

  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(parVol->dimension()));
  parVol->point(*result->data, value_u, value_v, value_w);

  return (PyObject*)result;
}

PyDoc_STRVAR(volume_flip_parametrization__doc__,"Flip (or reverse) volume parametrization\n"
                                                "@param direction: The parametric direction to flip (0=u, 1=v, 2=w)\n"
                                                "@type direction: int\n"
                                                "@return: The volume\n"
                                                "@rtype: Volume");
PyObject* Volume_FlipParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", NULL };
  PyObject* axiso;
  int direction = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&direction))
    return NULL;

  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;

  Volume* pyVol = (Volume*)self;
  if (!parVol->isSpline()) {
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }

  pyVol->data->reverseParameterDirection(direction);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(volume_get_bounding_box__doc__,"Generate and return the Spline Volume bounding box\n"
                                            "@return: 6 numbers representing the bounding box in order xmin,xmax,ymin,ymax,..."
                                            "@rtype: List of floats");
PyObject* Volume_GetBoundingBox(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }

  shared_ptr<Go::SplineVolume> spVol = static_pointer_cast<Go::SplineVolume>(parVol);
  Go::BoundingBox box = spVol->boundingBox();
  Go::Point low  = box.low();
  Go::Point high = box.high();

  PyObject* result = PyList_New(0);

  PyList_Append(result, Py_BuildValue((char*) "d", low [0]) );
  PyList_Append(result, Py_BuildValue((char*) "d", high[0]) );
  PyList_Append(result, Py_BuildValue((char*) "d", low [1]) );
  PyList_Append(result, Py_BuildValue((char*) "d", high[1]) );
  PyList_Append(result, Py_BuildValue((char*) "d", low [2]) );
  PyList_Append(result, Py_BuildValue((char*) "d", high[2]) );

  return result;
}

PyDoc_STRVAR(volume_get_const_par_surf__doc__,"Generate and return a SplineSurface that represents a constant parameter surface on the volume\n"
                                              "@param parameter: Value of the fixed parameter\n"
                                              "@type parameter: float\n"
                                              "@param pardir: 0 for constant u-parameter, 1 for v and 2 for w-parameter\n"
                                              "@type pardir: int\n"
                                              "@return: A newly constructed SplineSurface representing the specified constant parameter surface\n"
                                              "@rtype: Surface");
PyObject* Volume_GetConstParSurf(PyObject* self, PyObject* args, PyObject* kwds)
{
  int    parDir    = 0;
  double parameter = 0;

  static const char* keyWords[] = {"parameter", "pardir", NULL };
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"di",
                                   (char**)keyWords,&parameter,&parDir) || !parVol)
    return NULL;
                                                    
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }

  shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
  Go::SplineSurface *spSurf = spVol->constParamSurface(parameter, parDir);

  Surface* pySurf = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  pySurf->data = shared_ptr<Go::ParamSurface>(spSurf);

  return (PyObject*) pySurf;
}

PyDoc_STRVAR(volume_get_faces__doc__,"Return the six faces\n"
                                     "@return: A list of the six faces\n"
                                     "@rtype: List of Surface");
PyObject* Volume_GetFaces(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }

  shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
  std::vector<shared_ptr<Go::ParamSurface> > edges = spVol->getAllBoundarySurfaces();

  std::vector<Surface*> vec(6);
  PyObject* result = PyList_New(0);
  for (int i=0;i<6;++i) {
    Surface* srf = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    srf->data = shared_ptr<Go::ParamSurface>(edges.at(i));
    PyList_Append(result, (PyObject*) srf);
  }

  return result;
}

PyDoc_STRVAR(volume_get_sub_vol__doc__,"Get a Spline Volume which represent a part of 'this' Volume\n"
	                               "@param from_par: The parametric lower left corner of the sub volume\n"
                                       "@type  from_par: Point, list of floats or triple of floats\n"
                                       "@param to_par:   The parametric upper right corner of the sub volume\n"
                                       "@type  to_par:   Point, list of floats or triple of floats\n"
                                       "@return: A parametric rectangular subregion of this Volume represented as a Spline Volume\n"
                                         "@rtype: Volume");
PyObject* Volume_GetSubVol(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  static const char* keyWords[] = {"from_par", "to_par", NULL };
  PyObject *lowerLefto, *upperRighto;

  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO",
                                   (char**)keyWords, &lowerLefto, &upperRighto) || !parVol)
    return NULL;

  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }
  shared_ptr<Go::SplineVolume> spVol = static_pointer_cast<Go::SplineVolume>(parVol);
  if(!spVol)
    return NULL;

  shared_ptr<Go::Point> lowerLeft  = PyObject_AsGoPoint(lowerLefto);
  shared_ptr<Go::Point> upperRight = PyObject_AsGoPoint(upperRighto);

  Go::SplineVolume *subVol = spVol->subVolume((*lowerLeft)[0], (*lowerLeft)[1], (*lowerLeft)[2],(*upperRight)[0], (*upperRight)[1],(*upperRight)[2]);

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data = shared_ptr<Go::ParamVolume>(subVol);

  return (PyObject*) result;
}


PyDoc_STRVAR(volume_get_knots__doc__,"Return knots for a spline volume\n"
                                     "@param with_multiplicities: (optional) Set to true to obtain the knot vectors with multiplicities\n"
                                     "@type with_multiplicities: Boolean\n"
                                     "@return: List with the knot values\n"
                                     "@rtype: Tuple with List of float");
PyObject* Volume_GetKnots(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"with_multiplicities", NULL };
  bool withmult=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|b",
                                   (char**)keyWords, &withmult))
    return NULL;
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }
  PyObject* result = PyTuple_New(3);

  shared_ptr<Go::SplineVolume> spVol = static_pointer_cast<Go::SplineVolume>(parVol);
  for (int i=0;i<3;++i) {
    PyObject* list = PyList_New(0);
    std::vector<double> knots;
    if (withmult)
      knots = const_cast<Go::BsplineBasis&>(spVol->basis(i)).getKnots();
    else
      spVol->basis(i).knotsSimple(knots);

    for (std::vector<double>::iterator it  = knots.begin();
                                       it != knots.end();++it) {
      PyList_Append(list,Py_BuildValue((char*)"d",*it));
    }
    PyTuple_SetItem(result,i,list);
  }

  return result;
}

PyDoc_STRVAR(volume_get_order__doc__,"Return spline volume order (polynomial degree + 1) in all parametric directions\n"
                                     "@return: B-Spline order \n"
                                     "@rtype: List of three integers");
PyObject* Volume_GetOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }

  shared_ptr<Go::SplineVolume> spVol = static_pointer_cast<Go::SplineVolume>(parVol);

  PyObject* result = PyList_New(0);

  PyList_Append(result, Py_BuildValue((char*) "i", spVol->order(0)));
  PyList_Append(result, Py_BuildValue((char*) "i", spVol->order(1)));
  PyList_Append(result, Py_BuildValue((char*) "i", spVol->order(2)));

  return result;
}

PyDoc_STRVAR(volume_insert_knot__doc__,"Insert a knot in a spline volume\n"
                                       "@param direction: Direction to insert knot in\n"
                                       "@type direction: int (0, 1 or 2)\n"
                                       "@param knot: The knot to insert\n"
                                       "@type knot: float\n"
                                       "@return: The volume\n"
                                       "@rtype: Volume");
PyObject* Volume_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", "knot", NULL };
  int direction=0;
  double knot;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"id",
                                   (char**)keyWords,&direction,&knot))
    return NULL;

  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }
  static_pointer_cast<Go::SplineVolume>(parVol)->insertKnot(direction, knot);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(volume_make_rhs__doc__,"Make sure volume has a right-hand coordinate system\n"
                                    "@return: The volume\n"
                                    "@rtype: Volume");
PyObject* Volume_MakeRHS(PyObject* self, PyObject* args)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;

  // evaluate jacobian in mid point
  Go::Array<double,6> params = parVol->parameterSpan();
  double u = (params[0] + params[1]) / 2;
  double v = (params[2] + params[3]) / 2;
  double w = (params[4] + params[5]) / 2;
  vector<Go::Point> results(4); // one position and three derivatives
  parVol->point(results, u, v, w, 1);
  double jacobian = (results[1] % results[2])*results[3];
  if (jacobian < 0)
    parVol->reverseParameterDirection(2);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(volume_split__doc__, "Split the volume into segments\n"
                                  "@param params: The parameter values to split at\n"
                                  "@type params: Float or list of floats\n"
                                  "@param pardir: The direction to split along (0,1 or 2)\n"
                                  "@type pardir: int\n"
                                  "@return: The resulting volumes\n"
                                  "@rtype: List of Volume");
PyObject* Volume_Split(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"params", "pardir", NULL };
  PyObject* params;
  int pardir=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Oi",
                                   (char**)keyWords,&params,&pardir) || !params)
    return NULL;

  if (pardir < 0 || pardir > 2)
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

  shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(((Volume*)self)->data);
  std::vector<shared_ptr<Go::SplineVolume> > volumes = spVol->split(p,pardir);

  PyObject* result = PyList_New(0);
  for (size_t i=0;i<volumes.size();++i) {
    Volume* pyVol = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
    pyVol->data = volumes[i];
    PyList_Append(result,(PyObject*)pyVol);
  }

  return result;
}

PyDoc_STRVAR(volume_raise_order__doc__,"Raise order of a spline volume\n"
                                       "@param raise_u: Raise of order in u\n"
                                       "@type raise_u: int (>= 0)\n"
                                       "@param raise_v: Raise of order in v\n"
                                       "@type raise_v: int (>= 0)\n"
                                       "@param raise_w: Raise of order in w\n"
                                       "@type raise_w: int (>= 0)\n"
                                       "@returns: The volume\n"
                                       "@rtype: Volume");
PyObject* Volume_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"raise_u", "raise_v", "raise_w", NULL };
  int raise_u=0, raise_v=0, raise_w;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"iii",
                                   (char**)keyWords,&raise_u,&raise_v,&raise_w))
    return NULL;

  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }
  static_pointer_cast<Go::SplineVolume>(parVol)->raiseOrder(raise_u,raise_v,raise_w);

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(volume_lower_order__doc__,"Lower order of a spline volume\n"
                                       "@param lower_u: Lower of order in u\n"
                                       "@type lower_u: int (>= 0)\n"
                                       "@param Lower_v: Lower of order in v\n"
                                       "@type lower_v: int (>= 0)\n"
                                       "@param Lower_w: Lower of order in w\n"
                                       "@type lower_w: int (>= 0)\n"
                                       "@returns: The volume\n"
                                       "@rtype: Volume");
PyObject* Volume_LowerOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"lower_u", "lower_v", "lower_w", NULL };
  int lower_u=0, lower_v=0, lower_w=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"iii",
                                   (char**)keyWords,&lower_u,&lower_v,&lower_w))
    return NULL;

  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;
  if (!parVol->isSpline()) {
    Volume* pyVol = (Volume*)self;
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }
  shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
  std::vector<double>::const_iterator first =  spVol->basis(0).begin()+lower_u;
  std::vector<double>::const_iterator last  =  spVol->basis(0).end()-lower_u;
  Go::BsplineBasis b1 = Go::BsplineBasis(spVol->order(0)-lower_u,first,last);
  first =  spVol->basis(1).begin()+lower_v;
  last  =  spVol->basis(1).end()-lower_v;
  Go::BsplineBasis b2 = Go::BsplineBasis(spVol->order(1)-lower_v,first,last);
  first =  spVol->basis(2).begin()+lower_w;
  last  =  spVol->basis(2).end()-lower_w;
  Go::BsplineBasis b3 = Go::BsplineBasis(spVol->order(2)-lower_w,first,last);

  if (spVol->rational())
    std::cout << "WARNING: The geometry basis is rational (using NURBS)\n."
              << "         The basis for the unknown fields of one degree"
              << "         lower will however be non-rational.\n"
              << "         This may affect accuracy.\n"<< std::endl;

  std::vector<double> ug(b1.numCoefs()), vg(b2.numCoefs()), wg(b3.numCoefs());
  for (size_t i = 0; i < ug.size(); i++)
    ug[i] = b1.grevilleParameter(i);
  for (size_t i = 0; i < vg.size(); i++)
    vg[i] = b2.grevilleParameter(i);
  for (size_t i = 0; i < wg.size(); i++)
    wg[i] = b3.grevilleParameter(i);

  // Evaluate the spline surface at all points
  std::vector<double> XYZ(spVol->dimension()*ug.size()*vg.size()*wg.size());
  spVol->gridEvaluator(XYZ,ug,vg,wg);

  // Project the coordinates onto the new basis (the 2nd XYZ is dummy here)
  Volume* pySurf = (Volume*)self;
  pySurf->data.reset(Go::VolumeInterpolator::regularInterpolation(b1, b2, b3, ug,
                                                                  vg, wg, XYZ,
                                                                  spVol->dimension(),
                                                                  false, XYZ));

  return self;
}

PyDoc_STRVAR(volume_reparametrize__doc__,"Re-parametrize a volume\n"
                                         "@param umin: The minimum u value\n"
                                         "@type umin: float\n"
                                         "@param umax: The maximum u value\n"
                                         "@type umax: float\n"
                                         "@param vmin: The minimum v value\n"
                                         "@type vmin: float\n"
                                         "@param vmax: The maximum v value\n"
                                         "@type vmax: float\n"
                                         "@param wmin: The minimum w value\n"
                                         "@type wmin: float\n"
                                         "@param wmax: The maximum w value\n"
                                         "@type wmax: float\n"
                                         "@return: The volume\n"
                                         "@rtype: Volume");
PyObject* Volume_ReParametrize(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"umin", "umax", "vmin", "vmax", "wmin", "wmax", NULL };
  double umin=0, umax=1, vmin=0, vmax=1, wmin=0, wmax=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|dddddd",
                                   (char**)keyWords,&umin,&umax,&vmin,&vmax,
                                                    &wmin,&wmax))
    return NULL;
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (parVol) {
    shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
    if (!parVol->isSpline())
      ((Volume*)self)->data = spVol;
    spVol->setParameterDomain(umin, umax, vmin, vmax, wmin, wmax);
  }

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(volume_swap_parametrization__doc__,"Swaps two of the parametric directions\n"
                                                "@param pardir1: The first parametric direction \n"
                                                "@type pardir1: int\n"
                                                "@param pardir2: The second parametric direction \n"
                                                "@type pardir2: int\n"
                                                "@return: The volume\n"
                                                "@rtype: Volume");
PyObject* Volume_SwapParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"pardir1", "pardir2", NULL };
  PyObject* axiso;
  int pardir1 = 0;
  int pardir2 = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords,&pardir1, &pardir2))
    return NULL;

  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;

  Volume* pyVol = (Volume*)self;
  if (!parVol->isSpline()) {
    pyVol->data = convertSplineVolume(parVol);
    parVol = pyVol->data;
  }

  pyVol->data->swapParameterDirection(pardir1, pardir2);

  Py_INCREF(self);
  return self;
}

PyObject* Volume_Scale(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::ParamVolume> vol = PyObject_AsGoVolume(o1);
  double scale = PyFloat_AsDouble(o2);
  if (vol && scale > 0) {
    if (!vol->isSpline())
      ((Volume*)o1)->data = convertSplineVolume(vol);
    shared_ptr<Go::SplineVolume> svol = 
      dynamic_pointer_cast<Go::SplineVolume,
                           Go::ParamVolume>(((Volume*)o1)->data);
    for (std::vector<double>::iterator it  = svol->ctrl_begin(); 
                                       it != svol->ctrl_end(); ) {
      for (int i=0;i<vol->dimension();++i)
        *it++ *= scale;
      if (svol->rational())
        it++;
    }
  }

  return o1;
}

PyDoc_STRVAR(volume_translate__doc__,"Translate a volume along a given vector\n"
                                     "@param vector: The vector to translate along\n"
                                     "@type vector: Point, list of floats or tuple of floats\n"
                                     "@return: The volume\n"
                                     "@rtype: Volume");
PyObject* Volume_Translate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"vector", NULL };
  PyObject* veco;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&veco))
    return NULL;

  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  shared_ptr<Go::Point>       vec = PyObject_AsGoPoint(veco);
  if (!parVol || !vec)
    return NULL;

  if (!parVol->isSpline()) {
    Volume* volum = (Volume*)self;
    volum->data = convertSplineVolume(parVol);
    parVol = volum->data;
  }

  parVol->translate(*vec);

  Py_INCREF(self);
  return self;
}

PyMethodDef Volume_methods[] = {
     {(char*)"Clone",               (PyCFunction)Volume_Clone,                 METH_VARARGS|METH_KEYWORDS, volume_clone__doc__},
     {(char*)"Evaluate",            (PyCFunction)Volume_Evaluate,              METH_VARARGS|METH_KEYWORDS, volume_evaluate__doc__},
     {(char*)"FlipParametrization", (PyCFunction)Volume_FlipParametrization,   METH_VARARGS|METH_KEYWORDS, volume_flip_parametrization__doc__},
     {(char*)"GetBoundingBox",      (PyCFunction)Volume_GetBoundingBox,        METH_VARARGS              , volume_get_bounding_box__doc__},
     {(char*)"GetConstParSurf",     (PyCFunction)Volume_GetConstParSurf,       METH_VARARGS|METH_KEYWORDS, volume_get_const_par_surf__doc__},
     {(char*)"GetSubVol",           (PyCFunction)Volume_GetSubVol,             METH_VARARGS|METH_KEYWORDS, volume_get_sub_vol__doc__},
     {(char*)"GetFaces",            (PyCFunction)Volume_GetFaces,              METH_VARARGS|METH_KEYWORDS, volume_get_faces__doc__},
     {(char*)"GetKnots",            (PyCFunction)Volume_GetKnots,              METH_VARARGS|METH_KEYWORDS, volume_get_knots__doc__},
     {(char*)"GetOrder",            (PyCFunction)Volume_GetOrder,              METH_VARARGS,               volume_get_order__doc__},
     {(char*)"InsertKnot",          (PyCFunction)Volume_InsertKnot,            METH_VARARGS|METH_KEYWORDS, volume_insert_knot__doc__},
     {(char*)"LowerOrder",          (PyCFunction)Volume_LowerOrder,            METH_VARARGS|METH_KEYWORDS, volume_lower_order__doc__},
     {(char*)"MakeRHS",             (PyCFunction)Volume_MakeRHS,               METH_VARARGS,               volume_make_rhs__doc__},
     {(char*)"RaiseOrder",          (PyCFunction)Volume_RaiseOrder,            METH_VARARGS|METH_KEYWORDS, volume_raise_order__doc__},
     {(char*)"ReParametrize",       (PyCFunction)Volume_ReParametrize,         METH_VARARGS|METH_KEYWORDS, volume_reparametrize__doc__},
     {(char*)"Split",               (PyCFunction)Volume_Split,                 METH_VARARGS|METH_KEYWORDS, volume_split__doc__},
     {(char*)"SwapParametrization", (PyCFunction)Volume_SwapParametrization,   METH_VARARGS|METH_KEYWORDS, volume_swap_parametrization__doc__},
     {(char*)"Translate",           (PyCFunction)Volume_Translate,             METH_VARARGS|METH_KEYWORDS, volume_translate__doc__},
     {NULL,                         NULL,                                      0,                          NULL}
   };

 
Py_ssize_t Volume_NmbComponent(PyObject* self)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return 0;

  shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
  if(!spVol)
    return 0;

  return spVol->numCoefs(0)*spVol->numCoefs(1)*spVol->numCoefs(2);
}

PyObject* Volume_GetComponent(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(self);
  if (!parVol)
    return NULL;

  if(parVol->dimension() != 3) 
    return NULL;

  shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
  if(!spVol)
    return NULL;
  
  int nComp =  spVol->numCoefs(0)*spVol->numCoefs(1)*spVol->numCoefs(2);
  if(i < 0 || i >= nComp) 
    return NULL;

  double x,y,z,w;
  int dim = spVol->dimension();
  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  if(spVol->rational()) {
    vector<double>::const_iterator cp = spVol->rcoefs_begin() + i*(dim+1);
    result->data.reset(new Go::Point(cp, cp+(dim+1)));
  } else {
    x = *(spVol->coefs_begin() + i*(spVol->dimension() )+0) ;
    y = *(spVol->coefs_begin() + i*(spVol->dimension() )+1) ;
    z = *(spVol->coefs_begin() + i*(spVol->dimension() )+2) ;
    result->data.reset(new Go::Point(x,y,z));
  }

  return (PyObject*) result;
}

PyNumberMethods Volume_operators = {0};
PySequenceMethods Volume_seq_operators = {0};

PyDoc_STRVAR(volume__doc__, "A parametric description of a volume");
void init_Volume_Type()
{
  Volume_seq_operators.sq_item = Volume_GetComponent;
  Volume_seq_operators.sq_length = Volume_NmbComponent;
  InitializeTypeObject(&Volume_Type);
  Volume_operators.nb_add      = Volume_Add;
  Volume_operators.nb_subtract = Volume_Sub;
  Volume_operators.nb_inplace_multiply = Volume_Scale;
  Volume_Type.tp_name = "GoTools.Volume";
  Volume_Type.tp_basicsize = sizeof(Volume);
  Volume_Type.tp_dealloc = (destructor)Volume_Dealloc;
  Volume_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES;
  Volume_Type.tp_doc = volume__doc__;
  Volume_Type.tp_methods = Volume_methods;
  Volume_Type.tp_as_number = &Volume_operators;
  Volume_Type.tp_as_sequence = &Volume_seq_operators;
  Volume_Type.tp_base = 0;
  Volume_Type.tp_new = Volume_New;
  Volume_Type.tp_str = (reprfunc)Volume_Str;
  PyType_Ready(&Volume_Type);
}

}

shared_ptr<Go::SplineVolume> convertSplineVolume(shared_ptr<Go::ParamVolume> parVol)
{
  if (parVol->instanceType() == Go::Class_SplineVolume)
    return dynamic_pointer_cast<Go::SplineVolume, Go::ParamVolume>(parVol);
  shared_ptr<Go::ElementaryVolume> e_volume = dynamic_pointer_cast<Go::ElementaryVolume, Go::ParamVolume>(parVol);
  return shared_ptr<Go::SplineVolume>(e_volume->geometryVolume());
}

void printVolumeToStream(std::ostream& str, shared_ptr<Go::ParamVolume> parVol)
{
  if (parVol) {
    // fetch GoTools class type 
    std::string       sClassName  = Go::GoTools::className(parVol->instanceType());
    std::stringstream ssClassName ;
    ssClassName << sClassName;
    if(sClassName.compare("Unknown") == 0)
      ssClassName << "(" << parVol->instanceType() << ")";
    
    // print class type
    str << "Type: " << ssClassName.str() << std::endl;
    // print full raw data
    str << *parVol;
  } else
    str << "(empty)";
}

void WriteVolumeG2(std::ofstream& g2_file, Volume* pyVol, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(pyVol->data);
    if (spVol->isLeftHanded()) {
      spVol = shared_ptr<Go::SplineVolume>(spVol->clone());
      spVol->reverseParameterDirection(2);
    }
    spVol->writeStandardHeader(g2_file);
    spVol->write(g2_file);
  } else {
    pyVol->data->writeStandardHeader(g2_file);
    pyVol->data->write(g2_file);
  }
}

