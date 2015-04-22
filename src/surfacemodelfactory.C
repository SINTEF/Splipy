#include "surfacemodelfactory.h"

#include "geomodeller.h"
#include "pyutils.h"
#include "surface.h"
#include "surfacemodel.h"

#include "GoTools/compositemodel/RegularizeFace.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ParamSurface.h"

namespace GeoModeller {

extern "C" {
PyObject* SurfaceModelFactory_module;

PyDoc_STRVAR(generate_regularize_surface__doc__,
             "Regularize a surface\n"
             "@param surface: Surface to regularize\n"
             "@type surface: Surface\n"
             "@return: The resulting multi-surface model\n"
             "@rtype: SurfaceModel");
PyObject* Generate_RegularizeSurface(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", NULL };
  PyObject* surfaceo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&surfaceo))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(surfaceo);
  if (!surface) {
    PyErr_SetString(PyExc_RuntimeError, "Unable to obtain Go::ParamSurface");
    return NULL;
  }

  SurfaceModel* result = (SurfaceModel*)Surface_Type.tp_alloc(&SurfaceModel_Type,0);
  shared_ptr<Go::BoundedSurface> par_surf_bd = 
        dynamic_pointer_cast<Go::BoundedSurface, Go::ParamSurface>(surface);
  if (par_surf_bd && !par_surf_bd->allIsSpline())
    surface = shared_ptr<Go::ParamSurface>(par_surf_bd->allSplineCopy());

  shared_ptr<Go::ftSurface> ft_surf(new Go::ftSurface(surface,-1));
  ft_surf->createInitialEdges(modState.gapTolerance,modState.kinkTolerance);
  shared_ptr<Go::RegularizeFace> reg_face(new Go::RegularizeFace(ft_surf,
                                                                 modState.gapTolerance,
                                                                 modState.kinkTolerance,
                                                                 modState.neighbourTolerance));

  std::vector<shared_ptr<Go::ftSurface> > sub_faces = reg_face->getRegularFaces();
  result->data.reset(new Go::SurfaceModel(modState.approxTolerance,
                                          modState.gapTolerance,
                                          modState.neighbourTolerance,
                                          modState.kinkTolerance,
                                          modState.bendTolerance,
                                          sub_faces));

  return (PyObject*)result;
}

PyMethodDef SurfaceModelFactory_methods[] = {
     {(char*)"RegularizeSurface", (PyCFunction)Generate_RegularizeSurface, METH_VARARGS|METH_KEYWORDS, generate_regularize_surface__doc__},
     {NULL,                       0,                                       0,                          NULL}
  };

PyDoc_STRVAR(surface_model_factory__doc__,
             "A module with methods for generating surface models");

PyMODINIT_FUNC
init_SurfaceModelFactory_Module()
{
  SurfaceModelFactory_module = Py_InitModule3((char*)"GoTools.SurfaceModelFactory",SurfaceModelFactory_methods,surface_model_factory__doc__);
  modState.addInfo(SurfaceModelFactory_module);
}

}

}
