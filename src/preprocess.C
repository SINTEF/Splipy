#include "preprocess.h"
#include "curve.h"
#include "geomodeller.h"
#include "pyutils.h"
#include "surface.h"
#include "volume.h"

#if ENABLE_GPM
#include <GPM/SplineModel.h>
#endif
namespace GeoModeller {

  template<class T>
static std::vector<std::vector<int> >
        GenerateNumbers(T& data, std::vector<bool>& per)
{
  std::vector< std::vector<int> > numbers;
#if ENABLE_GPM
  SplineModel model(data);

  model.buildTopology(&per);

  model.generateGlobalNumbersPETSc();
  model.getGlobalNaturalNumbering(numbers);
#endif

  return numbers;
}

extern "C" {
PyObject* Preprocess_module;

PyDoc_STRVAR(preprocess_natural_node_numbers__doc__,"Generate natural node numbers for a model\n"
                                                    "@param patches: The patches to number\n"
                                                    "@type patches: List of (Surface or Volume)\n"
                                                    "@return: Node numbers\n"
                                                    "@rtype: List of (list of integer)");
PyObject* Preprocess_NaturalNodeNumbers(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_GPM
  std::cerr << "Compiled without GPM support - no numbering generated" << std::endl;
  Py_INCREF(Py_None);
  return Py_None;
#else

  static const char* keyWords[] = {"patches", "perX", "perY", "perZ", NULL };
  bool periodic[3] = {false, false, false};
  PyObject* patcheso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|bbb",
                                   (char**)keyWords,&patcheso,
                                                    &periodic[0],
                                                    &periodic[1],
                                                    &periodic[2]))
    return NULL;

  std::vector<bool> per;
  per.assign(periodic, periodic+3);

  std::vector<std::vector<int> > numbers;

  // list of surfaces or volumes
  if (PyObject_TypeCheck(patcheso,&PyList_Type)) {
    PyObject* obj1 = PyList_GetItem(patcheso, 0);
    if (PyObject_TypeCheck(obj1,&Surface_Type)) {
      std::vector< shared_ptr<Go::SplineSurface> > data;
      for (int i=0; i < PyList_Size(patcheso); ++i) {
        PyObject* obj = PyList_GetItem(patcheso, i);
        shared_ptr<Go::ParamSurface> parSurf = PyObject_AsGoSurface(obj);
        shared_ptr<Go::SplineSurface> spSurf = convertSplineSurface(parSurf);
        if (!spSurf)
          continue;
        data.push_back(spSurf);
      }
      numbers = GenerateNumbers(data, per);
    }
    if (PyObject_TypeCheck(obj1,&Volume_Type)) {
      std::vector< shared_ptr<Go::SplineVolume> > data;
      for (int i=0; i < PyList_Size(patcheso); ++i) {
        PyObject* obj = PyList_GetItem(patcheso, i);
        shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(obj);
        shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
        if (!spVol)
          continue;
        data.push_back(spVol);
      }
      numbers = GenerateNumbers(data, per);
    }
  }
  if (PyObject_TypeCheck(patcheso,&SurfaceModel_Type)) {
    SurfaceModel* sm = (SurfaceModel*)patcheso;
    std::vector< shared_ptr<Go::SplineSurface> > data;
    for (int i=0;i<sm->data->nmbEntities();++i)
      data.push_back(sm->data->getSplineSurface(i));
    numbers = GenerateNumbers(data, per);
  }
  if (PyObject_TypeCheck(patcheso,&VolumeModel_Type)) {
    VolumeModel* sm = (VolumeModel*)patcheso;
    std::vector< shared_ptr<Go::SplineVolume> > data;
    for (int i=0;i<sm->data->nmbEntities();++i)
      data.push_back(sm->data->getSplineVolume(i));
    numbers = GenerateNumbers(data, per);
  }

  PyObject* result = PyList_New(numbers.size());
  for (size_t i = 0; i < numbers.size(); ++i) {
    PyObject* lst = PyList_New(numbers[i].size());
    for (size_t j=0; j < numbers[i].size(); ++j)
      PyList_SetItem(lst, j, Py_BuildValue((char*) "i", numbers[i][j]));
    PyList_SetItem(result, i, lst);
  }

  return (PyObject*)result;
#endif
}

PyMethodDef Preprocess_methods[] = {
     {(char*)"NaturalNodeNumbers",    (PyCFunction)Preprocess_NaturalNodeNumbers, METH_VARARGS|METH_KEYWORDS, preprocess_natural_node_numbers__doc__},
     {NULL,                           NULL,                                       0,                          NULL}
  };

PyDoc_STRVAR(preprocess__doc__,"A module with methods for preprocessing models");

PyMODINIT_FUNC
init_Preprocess_Module()
{
  Preprocess_module = Py_InitModule3((char*)"GoTools.Preprocess",Preprocess_methods,preprocess__doc__);
  modState.addInfo(Preprocess_module);
}

}

}
