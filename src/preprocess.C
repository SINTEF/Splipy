#include "preprocess.h"
#include "curve.h"
#include "geomodeller.h"
#include "pyutils.h"
#include "surface.h"
#include "volume.h"

#if ENABLE_GPM
#include <GPM/SplineModel.h>
#include <GPM/TopologySet.h>
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

static size_t GetIdx(int face_idx, size_t i, size_t j,
                     bool reverse_u, bool reverse_v, bool uv_flip,
                     shared_ptr<Go::SplineVolume>& vol)
{
  size_t start, incu, incv, i2, j2, idxu, idxv;
  if (face_idx == 0) {
    start = 0;
    incu = vol->numCoefs(0);
    incv = vol->numCoefs(0)*vol->numCoefs(1);
    idxu = uv_flip?2:1;
    idxv = uv_flip?1:2;
  } else if (face_idx == 1) {
    start = vol->numCoefs(0)-1;
    incu = vol->numCoefs(0);
    incv = vol->numCoefs(0)*vol->numCoefs(1);
    idxu = uv_flip?2:1;
    idxv = uv_flip?1:2;
  } else if (face_idx == 2) {
    start = 0;
    incu = 1;
    incv = vol->numCoefs(0)*vol->numCoefs(1);
    idxu = uv_flip?2:0;
    idxv = uv_flip?0:2;
  } else if (face_idx == 3) {
    start = vol->numCoefs(0)*(vol->numCoefs(1)-1);
    incu = 1;
    incv = vol->numCoefs(0)*vol->numCoefs(1);
    idxu = uv_flip?2:0;
    idxv = uv_flip?0:2;
  } else if (face_idx == 4) {
    start = 0;
    incu = 1;
    incv = vol->numCoefs(0);
    idxu = uv_flip?1:0;
    idxv = uv_flip?0:1;
  } else { // face_idx == 5
    start = vol->numCoefs(0)*vol->numCoefs(1)*(vol->numCoefs(2)-1);
    incu = 1;
    incv = vol->numCoefs(0);
    idxu = uv_flip?1:0;
    idxv = uv_flip?0:1;
  }

  i2 = reverse_u?vol->numCoefs(idxu)-i-1:i;
  j2 = reverse_v?vol->numCoefs(idxv)-j-1:j;
  if (uv_flip)
    std::swap(i2, j2);

  return (start+i2*incu+j2*incv)*3;
}

static std::pair<size_t,size_t> GetSizes(int face_idx, shared_ptr<Go::SplineVolume>& vol)
{
  static const int sizes[][2] = {{1,2}, {0,2}, {0,1}};

  return std::make_pair(sizes[face_idx/2][0], sizes[face_idx/2][1]);
}

PyDoc_STRVAR(preprocess_average_faces__doc__,"Try to make a model water tight by averaging over matching faces\n"
                                             "@param patches: The patches to number\n"
                                             "@type patches: List of (Surface or Volume)\n"
                                             "@return: None\n");
PyObject* Preprocess_AverageFaces(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_GPM
  std::cerr << "Compiled without GPM support - no averaging performed" << std::endl;
#else
  static const char* keyWords[] = {"patches", NULL };
  PyObject* patcheso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&patcheso))
    return NULL;

  // list of surfaces or volumes
  if (PyObject_TypeCheck(patcheso,&PyList_Type)) {
    PyObject* obj1 = PyList_GetItem(patcheso, 0);
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
      SplineModel model(data);
      double difference=0.0;
      for (auto it  = model.getTopology()->face_begin();
                it != model.getTopology()->face_end();++it) {
        if ((*it)->volume.size() > 1) {
          shared_ptr<Go::SplineVolume> vol1 = data[(*it)->volume[0]->id];
          shared_ptr<Go::SplineVolume> vol2 = data[(*it)->volume[1]->id];
          pair<size_t, size_t> sizes = GetSizes((*it)->face[0], vol1);
          for (size_t j=0;j<sizes.second;++j) {
            for (size_t i=0;i<sizes.first;++i) {
              int idx1 = GetIdx((*it)->face[0], i, j, (*it)->u_reverse[0], (*it)->v_reverse[0], (*it)->uv_flip[0], vol1);
              int idx2 = GetIdx((*it)->face[1], i, j, (*it)->u_reverse[1], (*it)->v_reverse[1], (*it)->uv_flip[1], vol2);
              Go::Point p1(*(vol1->coefs_begin()+idx1), *(vol1->coefs_begin()+idx1+1), *(vol1->coefs_begin()+idx1+2));
              Go::Point p2(*(vol2->coefs_begin()+idx2), *(vol2->coefs_begin()+idx2+1), *(vol2->coefs_begin()+idx2+2));
              difference = std::max(difference, (p1-p2).length());
              Go::Point avg = (p1 + p2)*0.5;
              vol1->replaceCoefficient(idx1/3, avg);
              vol2->replaceCoefficient(idx2/3, avg);
            }
          }
        }
      }
      std::cout << "Maximum pointwise difference found: " << difference << std::endl;
    }
  }
#endif
  Py_INCREF(Py_None);
  return Py_None;
}

PyMethodDef Preprocess_methods[] = {
     {(char*)"NaturalNodeNumbers",    (PyCFunction)Preprocess_NaturalNodeNumbers, METH_VARARGS|METH_KEYWORDS, preprocess_natural_node_numbers__doc__},
     {(char*)"AverageFaces",          (PyCFunction)Preprocess_AverageFaces, METH_VARARGS|METH_KEYWORDS, preprocess_average_faces__doc__},
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
