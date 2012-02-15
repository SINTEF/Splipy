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
      PyObject* surfaceo = PyList_GetItem(surfaceso,i);
      shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(surfaceo);
      if (!surface)
        continue;
      surfaces.push_back(shared_ptr<Go::ftSurface>(new Go::ftSurface(surface,-1)));
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

PyMethodDef SurfaceModel_methods[] = {
     {NULL,           NULL,                     0,            NULL}
   };

PyDoc_STRVAR(surface_model__doc__, "A collection of parametric surfaces");
void init_SurfaceModel_Type()
{
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
      surf->writeStandardHeader(g2_file);
      surf->write(g2_file);
    } else {
      shared_ptr<Go::ParamSurface> surf = model->data->getSurface(i);
      surf->writeStandardHeader(g2_file);
      surf->write(g2_file);
    }
  }
}
