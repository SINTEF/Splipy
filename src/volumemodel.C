#include "volumemodel.h"
#include "pyutils.h"
#include "geomodeller.h"

#include "GoTools/geometry/ClassType.h"
#include "GoTools/trivariate/ElementaryVolume.h"
#include "GoTools/trivariate/SplineVolume.h"

#include <fstream>
#include <sstream>

extern "C"
{
PyTypeObject VolumeModel_Type;

PyObject* VolumeModel_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  VolumeModel* self;
  self = (VolumeModel*)type->tp_alloc(type,0);

  static const char* keyWords[] = {"volumes", NULL };
  PyObject* volumeso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"|O",
                                   (char**)keyWords,&volumeso))
    return NULL;

  if (PyObject_TypeCheck(volumeso,&PyList_Type)) {
    std::vector<shared_ptr<Go::ftVolume> > volumes;
    for (int i=0; i < PyList_Size(volumeso); ++i) {
      PyObject* volumeo = PyList_GetItem(volumeso,i);
      shared_ptr<Go::ParamVolume> parVol = PyObject_AsGoVolume(volumeo);
      if (!parVol)
        continue;
      shared_ptr<Go::SplineVolume> spVol = convertSplineVolume(parVol);
      if (!spVol)
        continue;
      volumes.push_back(shared_ptr<Go::ftVolume>(new Go::ftVolume(static_pointer_cast<Go::ParamVolume>(spVol),-1)));
    }
    if (volumes.size() < 1)
      return NULL;

    self->data.reset(new Go::VolumeModel(volumes,
                                         modState.gapTolerance,
                                         modState.neighbourTolerance,
                                         modState.kinkTolerance,
                                         modState.bendTolerance));
  }

  return (PyObject*)self;
}

void VolumeModel_Dealloc(VolumeModel* self)
{
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* VolumeModel_Str(VolumeModel* self)
{
  std::stringstream str;
  if (self->data) {
    std::cout << "Volume model with " << self->data->nmbEntities() << " volumes" << std::endl;
    for (int i=0;i<self->data->nmbEntities();++i) {
      shared_ptr<Go::ParamVolume> surf = self->data->getVolume(i);
      printVolumeToStream(str,surf);
    }
  } else
    str << "(empty)";
  return PyString_FromString(str.str().c_str());
}

PyObject* VolumeModel_Append(PyObject* o1, PyObject* o2)
{
  shared_ptr<Go::VolumeModel> model = PyObject_AsGoVolumeModel(o1);
  shared_ptr<Go::ParamVolume> surf  = PyObject_AsGoVolume(o2);

  if (!model || !surf)
    return NULL;

  model->append(shared_ptr<Go::ftVolume>(new Go::ftVolume(surf,-1)));

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(volumemodel_ctoc__doc__,"Returns whether or not VolumeModel has a Corner to Corner configuration\n"
                                     "@return: True/False\n"
                                     "@rtype: Bool");
PyObject* VolumeModel_CtoC(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::VolumeModel> vm = PyObject_AsGoVolumeModel(self);
  if (!vm)
    return NULL;

  return Py_BuildValue((char*)"b", vm->isCornerToCorner());
}

PyDoc_STRVAR(volumemodel_get__doc__,"Returns the i'th Volume of this VolumeModel\n"
                                    "@param i: index of the volume to return\n"
                                    "@type i: int\n"
                                    "@return: The i'th Volume\n"
                                    "@rtype: Volume");
PyObject* VolumeModel_Get(PyObject* self, Py_ssize_t i)
{
  shared_ptr<Go::VolumeModel> vm = PyObject_AsGoVolumeModel(self);
  if (!vm || i < 0 || i >= vm->nmbEntities())
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data   = vm->getSplineVolume(i);

  return (PyObject*) result;
}

PyDoc_STRVAR(volumemodel_get_bounding_box__doc__,"Generate and return the Volumemodel bounding box\n"
                                                 "@return: 6 numbers representing the bounding box in order xmin,xmax,ymin,ymax,...\n"
                                                 "@rtype: List of floats");
PyObject* VolumeModel_GetBoundingBox(PyObject* self, PyObject* args)
{
  shared_ptr<Go::VolumeModel> volumemodel = PyObject_AsGoVolumeModel(self);
  if (!volumemodel)
    return NULL;

  Go::BoundingBox box = volumemodel->boundingBox();
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

PyDoc_STRVAR(volumemodel_get_shell__doc__,"Generate and return all boundary faces of the volume model\n"
                                          "@return: The enclosing shell representation\n"
                                          "@rtype: List of Surfaces");
PyObject* VolumeModel_GetShell(PyObject* self, PyObject* args)
{
  shared_ptr<Go::VolumeModel> volumemodel = PyObject_AsGoVolumeModel(self);
  if (!volumemodel)
    return NULL;

  std::vector<shared_ptr<Go::ftSurface> > shell = volumemodel->getBoundaryFaces();

  PyObject* result = PyList_New(0);

  for(size_t i=0; i<shell.size(); i++) {
    Surface* pySurf = (Surface*)Surface_Type.tp_alloc(&Surface_Type,0);
    pySurf->data = shell[i]->surface();
    PyList_Append(result, (PyObject*) pySurf);
  }

  return result;
}

PyDoc_STRVAR(volumemodel_get_boundary_ids__doc__,"Return all boundary face ids of a volume model with a shoebox topology\n"
                                                 "@return: Lists of face ids\n"
                                                 "@rtype: List of (List of (tuple of integers))");
PyObject* VolumeModel_GetBoundaryIds(PyObject* self, PyObject* args)
{
  shared_ptr<Go::VolumeModel> volumemodel = PyObject_AsGoVolumeModel(self);
  if (!volumemodel)
    return NULL;

  PyObject* result = PyList_New(6);
  double face_loc, bnd_loc;
  bool face_in_plane, align;

  Go::BoundingBox bb = volumemodel->boundingBox();
  // loop over the six physical faces
  for (size_t i=0;i<6;++i) {
    PyObject* r = PyList_New(0);
    // loop over all volumes
    for (int j=0;j<volumemodel->nmbEntities();++j) {
      // loop over all faces of the box and check coordinates
      std::vector<shared_ptr<Go::SplineSurface> >
        srfs = volumemodel->getSplineVolume(j)->getBoundarySurfaces();
      for (int k=0;k<6;++k) {
        Go::BoundingBox fb = srfs[k]->boundingBox();
        face_in_plane = fabs(fb.low()[i/2]-fb.high()[i/2]) < modState.approxTolerance; // is face in correct plane?
        face_loc = (i%2 == 0?fb.low()[i/2]:fb.high()[i/2]); // face location
        bnd_loc = (i%2 == 0?bb.low()[i/2]:bb.high()[i/2]); // boundary location
        align = fabs(face_loc-bnd_loc) < modState.approxTolerance; // is face in same position as boundary?
        if (face_in_plane && align) {
          PyObject* pt = PyTuple_New(2);
          PyTuple_SetItem(pt, 0, Py_BuildValue((char*)"i", j+1));
          PyTuple_SetItem(pt, 1, Py_BuildValue((char*)"i", k+1));
          PyList_Append(r, pt);
        }
      }
    }
    PyList_SetItem(result, i, r);
  }

  return result;
}

PyDoc_STRVAR(volumemodel_make_common_spline__doc__,"Make sure all volumes in model share the same spline space\n"
                                                   "@return: None");
PyObject* VolumeModel_MakeCommonSpline(PyObject* self, PyObject* args)
{
  shared_ptr<Go::VolumeModel> vm = PyObject_AsGoVolumeModel(self);
  if (!vm)
    return NULL;

  vm->makeCommonSplineSpaces();

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(volumemodel_make_ctoc__doc__,"Force model into a corner to corner configuration\n"
                                          "@return: None");
PyObject* VolumeModel_MakeCtoC(PyObject* self, PyObject* args)
{
  shared_ptr<Go::VolumeModel> vm = PyObject_AsGoVolumeModel(self);
  if (!vm)
    return NULL;

  vm->makeCornerToCorner();

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(volumemodel_nmb_faces__doc__,"Returns the number of simple entities (Volumes) in this model\n"
                                          "@return: The number of faces in this model\n"
                                          "@rtype: integer");
Py_ssize_t VolumeModel_NmbFaces(PyObject* self)
{
  shared_ptr<Go::VolumeModel> vm = PyObject_AsGoVolumeModel(self);
  if (!vm)
    return 0;

  return vm->nmbEntities();
}

PySequenceMethods VolumeModel_seq_operators = {0};

PyMethodDef VolumeModel_methods[] = {
     {"GetBoundingBox",        (PyCFunction)VolumeModel_GetBoundingBox,     METH_VARARGS,               volumemodel_get_bounding_box__doc__},
     {"GetShell",              (PyCFunction)VolumeModel_GetShell,           METH_VARARGS,               volumemodel_get_shell__doc__},
     {"GetBoundaryIds",        (PyCFunction)VolumeModel_GetBoundaryIds,     METH_VARARGS,               volumemodel_get_boundary_ids__doc__},
     {"IsCornerToCorner",      (PyCFunction)VolumeModel_CtoC,               METH_VARARGS,               volumemodel_ctoc__doc__},
     {"MakeCommonSplineSpace", (PyCFunction)VolumeModel_MakeCommonSpline,   METH_VARARGS,               volumemodel_make_common_spline__doc__},
     {"MakeCornerToCorner",    (PyCFunction)VolumeModel_MakeCtoC,           METH_VARARGS,               volumemodel_make_ctoc__doc__},
     {NULL,                    NULL,                                        0,                          NULL}
   };

PyDoc_STRVAR(volume_model__doc__, "A collection of parametric volumes");
void init_VolumeModel_Type()
{
  VolumeModel_seq_operators.sq_inplace_concat = VolumeModel_Append;
  VolumeModel_seq_operators.sq_item           = VolumeModel_Get;
  VolumeModel_seq_operators.sq_length         = VolumeModel_NmbFaces;
  InitializeTypeObject(&VolumeModel_Type);
  VolumeModel_Type.tp_name = "GoTools.VolumeModel";
  VolumeModel_Type.tp_basicsize = sizeof(VolumeModel);
  VolumeModel_Type.tp_dealloc = (destructor)VolumeModel_Dealloc;
  VolumeModel_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  VolumeModel_Type.tp_doc = volume_model__doc__;
  VolumeModel_Type.tp_methods = VolumeModel_methods;
  VolumeModel_Type.tp_base = 0;
  VolumeModel_Type.tp_new = VolumeModel_New;
  VolumeModel_Type.tp_str = (reprfunc)VolumeModel_Str;
  VolumeModel_Type.tp_as_sequence = &VolumeModel_seq_operators;
  PyType_Ready(&VolumeModel_Type);
}

}

void WriteVolumeModelG2(std::ofstream& g2_file, VolumeModel* model, bool convert)
{
  if (!model->data)
    return;
  for (int i=0;i<model->data->nmbEntities();++i) {
    if (convert) {
      shared_ptr<Go::SplineVolume> vol = model->data->getSplineVolume(i);
      if (vol) {
        vol->writeStandardHeader(g2_file);
        vol->write(g2_file);
      }
    } else {
      shared_ptr<Go::ParamVolume> vol = model->data->getVolume(i);
      if (vol) {
        vol->writeStandardHeader(g2_file);
        vol->write(g2_file);
      }
    }
  }
}
