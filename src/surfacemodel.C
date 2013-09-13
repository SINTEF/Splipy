#include "surfacemodel.h"
#include "pyutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/ElementarySurface.h"
#include "GoTools/geometry/SplineSurface.h"

#include <fstream>
#include <sstream>

extern "C"
{
PyTypeObject SurfaceModel_Type;

PyObject* SurfaceModel_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  SurfaceModel* self;
  self = (SurfaceModel*)type->tp_alloc(type,0);

  static const char* keyWords[] = {"surfaces", NULL };
  PyObject* surfaceso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|O",
                                   (char**)keyWords,&surfaceso))
    return NULL;

  if (PyObject_TypeCheck(surfaceso,&PyList_Type)) {
    std::vector<shared_ptr<Go::ftSurface> > surfaces;
    for (int i=0; i < PyList_Size(surfaceso); ++i) {
      PyObject* pySurf = PyList_GetItem(surfaceso,i);
      shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(pySurf);
      shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(parSurf);
      if (!spSurf)
        continue;
      surfaces.push_back(shared_ptr<Go::ftSurface>(new Go::ftSurface(spSurf,-1)));
    }
    if (surfaces.size() < 1)
      return NULL;

    self->data.reset(new Go::SurfaceModel(modState.approxTolerance,
                                          modState.gapTolerance,
                                          modState.neighbourTolerance,
                                          modState.kinkTolerance,
                                          modState.bendTolerance,
                                          surfaces));
  }

  return (PyObject*)self;
}

void SurfaceModel_Dealloc(SurfaceModel* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* SurfaceModel_Str(SurfaceModel* self)
{
  std::stringstream str;
  if (self->data) {
    std::cout << "Surface model with " << self->data->nmbEntities() << " surfaces" << std::endl;
    for (int i=0;i<self->data->nmbEntities();++i) {
      shared_ptr<Go::ParamSurface> surf = self->data->getSurface(i);
      printSurfaceToStream(str,surf);
    }
  } else
    str << "(empty)";
  return PyString_FromString(str.str().c_str());
}

PyObject* SurfaceModel_Append(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::SurfaceModel> model = PyObject_AsGoSurfaceModel(o1);
  shared_ptr<Go::ParamSurface> surf  = PyObject_AsGoSurface(o2);

  if (!model || !surf)
    return NULL;

  model->append(shared_ptr<Go::ftSurface>(new Go::ftSurface(surf,-1)));

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(surfacemodel_ctoc__doc__,"Returns whether or not SurfaceModel has a Corner to Corner configuration\n"
                                      "@return: True/False\n"
                                      "@rtype: Bool");
PyObject* SurfaceModel_CtoC(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::SurfaceModel> sm = PyObject_AsGoSurfaceModel(self);
  if (!sm)
    return NULL;

  return Py_BuildValue((char*)"b", sm->isCornerToCorner());
}

PyDoc_STRVAR(surfacemodel_get__doc__,"Returns the i'th Surface of this SurfaceModel\n"
                                     "@param i: index of the surface to return\n"
                                     "@type i: int\n"
                                     "@return: The i'th Surface\n"
                                     "@rtype: Surface");
PyObject* SurfaceModel_Get(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::SurfaceModel> sm = PyObject_AsGoSurfaceModel(self);
  if (!sm || i < 0 || i >= sm->nmbEntities())
    return NULL;

  Surface* result = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
  result->data = sm->getFace(i)->getUntrimmed(modState.gapTolerance,
                                              modState.neighbourTolerance,
                                              modState.kinkTolerance);

  return (PyObject*) result;
}

PyDoc_STRVAR(surfacemodel_get_bounding_box__doc__,"Generate and return the Surfacemodel bounding box\n"
                                                  "@return: 6 numbers representing the bounding box in order xmin,xmax,ymin,ymax,...\n"
                                                  "@rtype: List of floats");
PyObject* SurfaceModel_GetBoundingBox(PyObject* self, PyObject* args)
{
  shared_ptr<Go::SurfaceModel> surfacemodel = PyObject_AsGoSurfaceModel(self);
  if (!surfacemodel)
    return NULL;

  Go::BoundingBox box = surfacemodel->boundingBox();
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

PyDoc_STRVAR(surfacemodel_make_common_spline__doc__,"Make sure all surfaces in model share the same spline space\n"
                                                    "@return: None");
PyObject* SurfaceModel_MakeCommonSpline(PyObject* self, PyObject* args)
{
  shared_ptr<Go::SurfaceModel> sm = PyObject_AsGoSurfaceModel(self);
  if (!sm)
    return NULL;

  sm->makeCommonSplineSpaces();

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(surfacemodel_make_ctoc__doc__,"Force model into a corner to corner configuration\n"
                                           "@return: None");
PyObject* SurfaceModel_MakeCtoC(PyObject* self, PyObject* args)
{
  shared_ptr<Go::SurfaceModel> sm = PyObject_AsGoSurfaceModel(self);
  if (!sm)
    return NULL;

  sm->makeCornerToCorner();

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(surfacemodel_nmb_faces__doc__,"Returns the number of simple entities (Surfaces) in this model\n"
                                           "@return: The number of faces in this model\n"
                                           "@rtype: integer");
Py_ssize_t SurfaceModel_NmbFaces(PyObject* self)
{
  shared_ptr<Go::SurfaceModel> sm = PyObject_AsGoSurfaceModel(self);
  if (!sm)
    return 0;

  return sm->nmbEntities();
}

PySequenceMethods SurfaceModel_seq_operators = {0};

PyMethodDef SurfaceModel_methods[] = {
     {"GetBoundingBox",        (PyCFunction)SurfaceModel_GetBoundingBox,     METH_VARARGS,               surfacemodel_get_bounding_box__doc__},
     {"IsCornerToCorner",      (PyCFunction)SurfaceModel_CtoC,               METH_VARARGS,               surfacemodel_ctoc__doc__},
     {"MakeCommonSplineSpace", (PyCFunction)SurfaceModel_MakeCommonSpline,   METH_VARARGS,               surfacemodel_make_common_spline__doc__},
     {"MakeCornerToCorner",    (PyCFunction)SurfaceModel_MakeCtoC,           METH_VARARGS,               surfacemodel_make_ctoc__doc__},
     {NULL,                    NULL,                                         0,                          NULL}
   };

PyDoc_STRVAR(surface_model__doc__, "A collection of parametric surfaces");
void init_SurfaceModel_Type()
{
  SurfaceModel_seq_operators.sq_inplace_concat = SurfaceModel_Append;
  SurfaceModel_seq_operators.sq_item           = SurfaceModel_Get;
  SurfaceModel_seq_operators.sq_length         = SurfaceModel_NmbFaces;
  InitializeTypeObject(&SurfaceModel_Type);
  SurfaceModel_Type.tp_name = "GoTools.SurfaceModel";
  SurfaceModel_Type.tp_basicsize = sizeof(SurfaceModel);
  SurfaceModel_Type.tp_dealloc = (destructor)SurfaceModel_Dealloc;
  SurfaceModel_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  SurfaceModel_Type.tp_doc = surface_model__doc__;
  SurfaceModel_Type.tp_methods = SurfaceModel_methods;
  SurfaceModel_Type.tp_base = 0;
  SurfaceModel_Type.tp_new = SurfaceModel_New;
  SurfaceModel_Type.tp_str = (reprfunc)SurfaceModel_Str;
  SurfaceModel_Type.tp_as_sequence = &SurfaceModel_seq_operators;
  PyType_Ready(&SurfaceModel_Type);
}

}

void WriteSurfaceModelG2(std::ofstream& g2_file, SurfaceModel* model, bool convert)
{
  if (!model->data)
    return;
  for (int i=0;i<model->data->nmbEntities();++i) {
    if (convert) {
      shared_ptr<Go::SplineSurface> surf = model->data->getSplineSurface(i);
      if (surf) {
        surf->writeStandardHeader(g2_file);
        surf->write(g2_file);
      }
    } else {
      shared_ptr<Go::ParamSurface> surf = model->data->getSurface(i);
      if (surf) {
        surf->writeStandardHeader(g2_file);
        surf->write(g2_file);
      }
    }
  }
}
