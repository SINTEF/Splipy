#include "volume.h"

#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/trivariate/ElementaryVolume.h"
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
    ((Volume*)self)->data.reset(new Go::SplineVolume(knots1.size()-order1,
                                                     knots2.size()-order2, 
                                                     knots3.size()-order3, 
                                                     order1, order2, order3,
                                                     knots1.begin(),
                                                     knots2.begin(),
                                                     knots3.begin(),
                                                     coefs.begin(), modState.dim, rational));
  }

  return (PyObject*)self;
}

void Volume_Dealloc(Volume* self)
{
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
  shared_ptr<Go::ParamVolume> vol = PyObject_AsGoVolume(self);
  res->data.reset(vol->clone());
 
  return (PyObject*)res;
}

PyObject* Volume_Add(PyObject* o1, PyObject* o2)
{
  Volume* vol = (Volume*)o1;
  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::Point> p = PyObject_AsGoPoint(o2);
  if (p) {
    result->data.reset(vol->data->clone());
    result->data->translate(*p);
  }

  return (PyObject*)result;
}

PyObject* Volume_Sub(PyObject* o1, PyObject* o2)
{
  Volume* vol = (Volume*)o1;
  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  shared_ptr<Go::Point> p = PyObject_AsGoPoint(o2);
  if (p) {
    result->data.reset(vol->data->clone());
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
  shared_ptr<Go::ParamVolume> vol = PyObject_AsGoVolume(self);
  double value_u=0, value_v=0, value_w=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ddd",
                                   (char**)keyWords,&value_u,&value_v,
                                                    &value_w) || !vol)
    return NULL;

  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  result->data.reset(new Go::Point(vol->dimension()));
  vol->point(*result->data, value_u, value_v, value_w);

  return (PyObject*)result;
}

PyDoc_STRVAR(volume_flip_parametrization__doc__,"Flip (or reverse) volume parametrization\n"
                                                "@param direction: The parametric direction to flip (0=u, 1=v, 2=w)\n"
                                                "@type direction: int\n"
                                                "@return: None");
PyObject* Volume_FlipParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", NULL };
  PyObject* axiso;
  int direction = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&direction))
    return NULL;

  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;

   Volume* vol = (Volume*)self;
   if (!volume->isSpline()) {
     vol->data = convertSplineVolume(volume);
     volume = vol->data;
   }

   vol->data->reverseParameterDirection(direction);

   Py_INCREF(Py_None);
   return Py_None;
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
  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"di",
                                   (char**)keyWords,&parameter,&parDir) || !volume)
    return NULL;
                                                    
  if (!volume->isSpline()) {
    Volume* vol = (Volume*)self;
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }

  shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume);
  Go::SplineSurface *surf = vol->constParamSurface(parameter, parDir);

  Surface* surface = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  surface->data = shared_ptr<Go::ParamSurface>(surf);

  return (PyObject*) surface;
}

PyDoc_STRVAR(volume_get_edges__doc__,"Return the six edge volumes\n"
                                     "@return: A list of the six edges\n"
                                     "@rtype: List of Surface");
PyObject* Volume_GetEdges(PyObject* self, PyObject* args, PyObject* kwds)
{
  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;
  if (!volume->isSpline()) {
    Volume* vol = (Volume*)self;
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }

  shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume);
  std::vector<shared_ptr<Go::ParamSurface> > edges = vol->getAllBoundarySurfaces();

  Surface* surface0 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface1 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface2 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface3 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface4 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  Surface* surface5 = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  surface0->data = shared_ptr<Go::ParamSurface>(edges.at(0));
  surface1->data = shared_ptr<Go::ParamSurface>(edges.at(1));
  surface2->data = shared_ptr<Go::ParamSurface>(edges.at(2));
  surface3->data = shared_ptr<Go::ParamSurface>(edges.at(3));
  surface4->data = shared_ptr<Go::ParamSurface>(edges.at(4));
  surface5->data = shared_ptr<Go::ParamSurface>(edges.at(5));

  PyObject* result = PyList_New(0);
  PyList_Append(result, (PyObject*) surface0);
  PyList_Append(result, (PyObject*) surface1);
  PyList_Append(result, (PyObject*) surface2);
  PyList_Append(result, (PyObject*) surface3);
  PyList_Append(result, (PyObject*) surface4);
  PyList_Append(result, (PyObject*) surface5);

  return result;
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
  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;
  if (!volume->isSpline()) {
    Volume* vol = (Volume*)self;
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }
  PyObject* result = PyTuple_New(3);

  shared_ptr<Go::SplineVolume> vol = static_pointer_cast<Go::SplineVolume>(volume);
  for (int i=0;i<3;++i) {
    PyObject* list = PyList_New(0);
    std::vector<double> knots;
    if (withmult)
      knots = const_cast<Go::BsplineBasis&>(vol->basis(i)).getKnots();
    else
      vol->basis(i).knotsSimple(knots);

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
  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;
  if (!volume->isSpline()) {
    Volume* vol = (Volume*)self;
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }

  shared_ptr<Go::SplineVolume> spline = static_pointer_cast<Go::SplineVolume>(volume);

  PyObject* result = PyList_New(0);

  PyList_Append(result, Py_BuildValue((char*) "i", spline->order(0)));
  PyList_Append(result, Py_BuildValue((char*) "i", spline->order(1)));
  PyList_Append(result, Py_BuildValue((char*) "i", spline->order(2)));

  return result;
}

PyDoc_STRVAR(volume_insert_knot__doc__,"Insert a knot in a spline volume\n"
                                       "@param direction: Direction to insert knot in\n"
                                       "@type direction: int (0, 1 or 2)\n"
                                       "@param knot: The knot to insert\n"
                                       "@type knot: float\n"
                                       "@return: None");
PyObject* Volume_InsertKnot(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"direction", "knot", NULL };
  int direction=0;
  double knot;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"id",
                                   (char**)keyWords,&direction,&knot))
    return NULL;

  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;
  if (!volume->isSpline()) {
    Volume* vol = (Volume*)self;
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }
  static_pointer_cast<Go::SplineVolume>(volume)->insertKnot(direction, knot);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(volume_make_rhs__doc__,"Make sure volume has a right-hand coordinate system\n"
                                    "@return: None");
PyObject* Volume_MakeRHS(PyObject* self, PyObject* args)
{
  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;

  // evaluate jacobian in mid point
  Go::Array<double,6> params = volume->parameterSpan();
  double u = (params[0] + params[1]) / 2;
  double v = (params[2] + params[3]) / 2;
  double w = (params[4] + params[5]) / 2;
  vector<Go::Point> results(4); // one position and three derivatives
  volume->point(results, u, v, w, 1);
  double jacobian = (results[1] % results[2])*results[3];
  if (jacobian < 0)
    volume->reverseParameterDirection(2);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(volume_split__doc__, "Split the volume into segments\n"
                                  "@param params: The parameter values to split at\n"
                                  "@type params: Float or list of floats\n"
                                  "@param pardir: The direction to split along\n"
                                  "@type pardir: integer (0,1,2)\n"
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

  shared_ptr<Go::SplineVolume> vol = convertSplineVolume(((Volume*)self)->data);
  std::vector<shared_ptr<Go::SplineVolume> > volumes = vol->split(p,pardir);

  PyObject* result = PyList_New(0);
  for (size_t i=0;i<volumes.size();++i) {
    Volume* part = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
    part->data = volumes[i];
    PyList_Append(result,(PyObject*)part);
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
                                       "@returns: None");
PyObject* Volume_RaiseOrder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"raise_u", "raise_v", "raise_w", NULL };
  int raise_u=0, raise_v=0, raise_w;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"iii",
                                   (char**)keyWords,&raise_u,&raise_v,&raise_w))
    return NULL;

  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;
   if (!volume->isSpline()) {
     Volume* vol = (Volume*)self;
     vol->data = convertSplineVolume(volume);
     volume = vol->data;
   }
   static_pointer_cast<Go::SplineVolume>(volume)->raiseOrder(raise_u,raise_v,raise_w);

   Py_INCREF(Py_None);
   return Py_None;
}

PyDoc_STRVAR(volume_swap_parametrization__doc__,"Swaps two of the parametric directions\n"
                                                "@param pardir1: The first parametric direction \n"
                                                "@type pardir1: int\n"
                                                "@param pardir2: The second parametric direction \n"
                                                "@type pardir2: int\n"
                                                "@return: None");
PyObject* Volume_SwapParametrization(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"pardir1", "pardir2", NULL };
  PyObject* axiso;
  int pardir1 = 0;
  int pardir2 = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ii",
                                   (char**)keyWords,&pardir1, &pardir2))
    return NULL;

  shared_ptr<Go::ParamVolume> volume = PyObject_AsGoVolume(self);
  if (!volume)
    return NULL;

  Volume* vol = (Volume*)self;
  if (!volume->isSpline()) {
    vol->data = convertSplineVolume(volume);
    volume = vol->data;
  }

  vol->data->swapParameterDirection(pardir1, pardir2);

  Py_INCREF(Py_None);
  return Py_None;
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
                                     "@return: None");
PyObject* Volume_Translate(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"vector", NULL };
  PyObject* veco;
  double angle=0.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&veco))
    return NULL;

  shared_ptr<Go::ParamVolume> vol = PyObject_AsGoVolume(self);
  shared_ptr<Go::Point>       vec = PyObject_AsGoPoint(veco);
  if (!vol || !vec)
    return NULL;

   if (!vol->isSpline()) {
     Volume* volum = (Volume*)self;
     volum->data = convertSplineVolume(vol);
     vol = volum->data;
   }

   vol->translate(*vec);

   Py_INCREF(Py_None);
   return Py_None;
}

PyMethodDef Volume_methods[] = {
     {(char*)"Clone",               (PyCFunction)Volume_Clone,                 METH_VARARGS|METH_KEYWORDS, volume_clone__doc__},
     {(char*)"Evaluate",            (PyCFunction)Volume_Evaluate,              METH_VARARGS|METH_KEYWORDS, volume_evaluate__doc__},
     {(char*)"FlipParametrization", (PyCFunction)Volume_FlipParametrization,   METH_VARARGS|METH_KEYWORDS, volume_flip_parametrization__doc__},
     {(char*)"GetConstParSurf",     (PyCFunction)Volume_GetConstParSurf,       METH_VARARGS|METH_KEYWORDS, volume_get_const_par_surf__doc__},
     {(char*)"GetEdges",            (PyCFunction)Volume_GetEdges,              METH_VARARGS|METH_KEYWORDS, volume_get_edges__doc__},
     {(char*)"GetKnots",            (PyCFunction)Volume_GetKnots,              METH_VARARGS|METH_KEYWORDS, volume_get_knots__doc__},
     {(char*)"GetOrder",            (PyCFunction)Volume_GetOrder,              METH_VARARGS,               volume_get_order__doc__},
     {(char*)"InsertKnot",          (PyCFunction)Volume_InsertKnot,            METH_VARARGS|METH_KEYWORDS, volume_insert_knot__doc__},
     {(char*)"MakeRHS",             (PyCFunction)Volume_MakeRHS,               METH_VARARGS,               volume_make_rhs__doc__},
     {(char*)"RaiseOrder",          (PyCFunction)Volume_RaiseOrder,            METH_VARARGS|METH_KEYWORDS, volume_raise_order__doc__},
     {(char*)"Split",               (PyCFunction)Volume_Split,                 METH_VARARGS|METH_KEYWORDS, volume_split__doc__},
     {(char*)"SwapParametrization", (PyCFunction)Volume_SwapParametrization,   METH_VARARGS|METH_KEYWORDS, volume_swap_parametrization__doc__},
     {(char*)"Translate",           (PyCFunction)Volume_Translate,             METH_VARARGS|METH_KEYWORDS, volume_translate__doc__},
     {NULL,                         NULL,                                      0,                          NULL}
   };

 
Py_ssize_t Volume_NmbComponent(PyObject* self)
{
  shared_ptr<Go::ParamVolume> pv = PyObject_AsGoVolume(self);
  if (!pv)
    return 0;

  shared_ptr<Go::SplineVolume> sv = convertSplineVolume(pv);
  if(!sv)
    return 0;

  return sv->numCoefs(0)*sv->numCoefs(1)*sv->numCoefs(2);
}

PyObject* Volume_GetComponent(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::ParamVolume> pv = PyObject_AsGoVolume(self);
  if (!pv)
    return NULL;

  if(pv->dimension() != 3) 
    return NULL;

  shared_ptr<Go::SplineVolume> sv = convertSplineVolume(pv);
  if(!sv)
    return NULL;
  
  int nComp =  sv->numCoefs(0)*sv->numCoefs(1)*sv->numCoefs(2);
  if(i < 0 || i >= nComp) 
    return NULL;

  double x,y,z,w;
  int dim = sv->dimension();
  Point* result = (Point*)Point_Type.tp_alloc(&Point_Type,0);
  if(sv->rational()) {
    vector<double>::const_iterator cp = sv->rcoefs_begin() + i*(dim+1);
    result->data.reset(new Go::Point(cp, cp+(dim+1)));
  } else {
    x = *(sv->coefs_begin() + i*(sv->dimension() )+0) ;
    y = *(sv->coefs_begin() + i*(sv->dimension() )+1) ;
    z = *(sv->coefs_begin() + i*(sv->dimension() )+2) ;
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

shared_ptr<Go::SplineVolume> convertSplineVolume(shared_ptr<Go::ParamVolume> volume)
{
  if (volume->instanceType() == Go::Class_SplineVolume)
    return dynamic_pointer_cast<Go::SplineVolume, Go::ParamVolume>(volume);
  shared_ptr<Go::ElementaryVolume> e_volume = dynamic_pointer_cast<Go::ElementaryVolume, Go::ParamVolume>(volume);
  return shared_ptr<Go::SplineVolume>(e_volume->geometryVolume());
}

void WriteVolumeG2(std::ofstream& g2_file, Volume* volume, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume->data);
    if (vol->isLeftHanded()) {
      vol = shared_ptr<Go::SplineVolume>(vol->clone());
      vol->reverseParameterDirection(2);
    }
    vol->writeStandardHeader(g2_file);
    vol->write(g2_file);
  } else {
    volume->data->writeStandardHeader(g2_file);
    volume->data->write(g2_file);
  }
}

int WriteVolumeSTL(std::ofstream& stl_file, Volume* volume, bool ascii, int res[3])
{
  int nTriangles = 0;
  shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume->data);
  if (vol->isLeftHanded()) {
    vol = shared_ptr<Go::SplineVolume>(vol->clone());
    vol->reverseParameterDirection(2);
  }

  std::vector<shared_ptr<Go::ParamSurface> > edges = vol->getAllBoundarySurfaces();

  edges[1]->reverseParameterDirection(true);
  edges[2]->reverseParameterDirection(true);
  edges[5]->reverseParameterDirection(true);
  
  for(int i=0; i<6; i++) {
    Surface* surface = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    surface->data = shared_ptr<Go::ParamSurface>(edges.at(i));

    int surfRes[] = { res[i/2==0], res[2-(i/2==2)] };
    nTriangles += WriteSurfaceSTL(stl_file, surface, ascii, surfRes);
  }
  return nTriangles;
}
