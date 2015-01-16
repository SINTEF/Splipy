#include "vtf.h"

#include "geomodeller.h"
#include "pyutils.h"

#include <fstream>
#include <sstream>

namespace GeoModeller {

extern "C"
{
PyTypeObject VTF_Type;

PyObject* VTF_New(PyTypeObject* type, PyObject* args, PyObject* kwds)
{
  PyVTF* self;
  self = (PyVTF*)type->tp_alloc(type,0);
  self->nBlock = 0;
#ifdef ENABLE_VTF
  self->step = 1;
  self->state = NULL;
  self->dBlocks = new std::map<std::string,VTFADisplacementBlock*>();
  self->vBlocks = new std::map<std::string,VTFAVectorBlock*>();
  self->sBlocks = new std::map<std::string,VTFAScalarBlock*>();
  self->sets = new std::map<std::string,VTFASetBlock*>();

  static const char* keyWords[] = {"filename", "binary", NULL};
  const char* fname=NULL;
  bool binary=true;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s|b",
                                   (char**)keyWords,&fname,&binary) || !fname)
    return NULL;

  self->data.reset(new VTFAFile());
  self->data->CreateVTFFile(fname,binary);

#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  return (PyObject*)self;
}

void VTF_Dealloc(PyVTF* self)
{
#ifdef ENABLE_VTF
  for (auto it = self->dBlocks->begin(); it != self->dBlocks->end();++it) {
    self->data->WriteBlock(it->second);
    delete it->second;
  }
  delete self->dBlocks;
  for (auto it = self->vBlocks->begin(); it != self->vBlocks->end();++it) {
    self->data->WriteBlock(it->second);
    delete it->second;
  }
  delete self->vBlocks;
  for (auto it = self->sBlocks->begin(); it != self->sBlocks->end();++it) {
    self->data->WriteBlock(it->second);
    delete it->second;
  }
  delete self->sBlocks;
  for (auto it = self->sets->begin(); it != self->sets->end();++it) {
    self->data->WriteBlock(it->second);
    delete it->second;
  }
  delete self->sets;

  self->data->WriteBlock(self->state);
  delete self->state;
  self->data->CloseFile();
  self->data.reset();
#endif
  self->ob_type->tp_free((PyObject*)self);
}

PyObject* VTF_Str(PyVTF* self)
{
  std::stringstream str;
#ifdef ENABLE_VTF
  if (self->data) {
    str << "opened";
  } else
    str << "not opened";
#endif
  return PyString_FromString(str.str().c_str());
}

PyDoc_STRVAR(vtf_addgeometryblock__doc__,"Add a geometry block to the VTF file\n"
                                         "@param nodes: Nodal coordinates (interleaved format)\n"
                                         "@type nodes: List of floats\n"
                                         "@param elements: Element nodal connectivities\n"
                                         "@type elements: List of int\n"
                                         "@return: Self");
PyObject* VTF_AddGeometryBlock(PyObject* self, PyObject* args, PyObject* kwds)
{
  PyVTF* vtf = (PyVTF*)self;
#ifdef ENABLE_VTF
  static const char* keyWords[] = {"nodes", "elements", "blockid", "dim", NULL };
  PyObject* nodeso=NULL;
  PyObject* elementso=NULL;
  int blockId=1;
  int dim=3;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OO|ii",
                                   (char**)keyWords,&nodeso,&elementso,&blockId,&dim))
    return NULL;

  if (!PyObject_TypeCheck(nodeso,&PyList_Type) ||
      !PyObject_TypeCheck(elementso,&PyList_Type)) {
    std::cerr << "Invalid parameters in VTF.AddGeometryBlock()" << std::endl;
    return NULL;
  }

  VTFANodeBlock block(blockId,0);
  block.SetNumNodes(PyList_Size(nodeso)/modState.dim);
  int j=0;
  for (int i=0;i<PyList_Size(nodeso)/modState.dim;++i) {
    double coord[3] = {0};
    for (size_t k=0;k<modState.dim;++k)
      coord[k] = PyFloat_AsDouble(PyList_GetItem(nodeso,j++));
    block.AddNode(coord[0],coord[1],coord[2]);
  }
  vtf->data->WriteBlock(&block);

  std::vector<int> mnpc;
  for (int i=0;i<PyList_Size(elementso);++i)
    mnpc.push_back(PyInt_AsLong(PyList_GetItem(elementso,i)));

  VTFAElementBlock eBlock(blockId,0,0);

  if (dim == 1) {
    int nel = mnpc.size()/2;
    eBlock.AddElements(VTFA_BEAMS,&mnpc.front(),nel);
  } else if (dim == 2) {
    int nel = mnpc.size()/4;
    eBlock.AddElements(VTFA_QUADS,&mnpc.front(),nel);
  } else { // dim == 3
    int nel = mnpc.size()/8;
    eBlock.AddElements(VTFA_HEXAHEDRONS,&mnpc.front(),nel);
  }
  eBlock.SetPartID(blockId);
  std::stringstream str;
  str << "Patch " << blockId;
  eBlock.SetPartName(str.str().c_str());
  eBlock.SetNodeBlockID(blockId);
  vtf->data->WriteBlock(&eBlock);
#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(vtf_addgeometrydescriptor__doc__,"Add a block connecting geometry blocks\n"
                                              "@param size: Number of blocks in part\n"
                                              "@type size: Int\n"
                                              "@param blockid: ID of block being added\n"
                                              "@type blockid: Integer\n"
                                              "@param dim: Dimensionality of geometry block\n"
                                              "@type dim: Integer\n"
                                              "@return: Self\n"
                                              "@rtype: VTFFile");
PyObject* VTF_AddGeometryDescriptor(PyObject* self, PyObject* args, PyObject* kwds)
{
  PyVTF* vtf = (PyVTF*)self;
#ifdef ENABLE_VTF
  static const char* keyWords[] = {"size", NULL };
  int size=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&size))
    return NULL;

  std::vector<int> blocks;
  for (size_t i=0;i<size;++i)
    blocks.push_back(i+1);

  VTFAGeometryBlock gBlock;
  gBlock.SetGeometryElementBlocks(&blocks.front(), blocks.size(), vtf->step);
  vtf->data->WriteBlock(&gBlock);
#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  Py_INCREF(self);
  return self;
}

PyDoc_STRVAR(vtf_addgeometryset__doc__,"Add a set describing a geometry subpart\n"
                                       "@param elems: Geometry element numbers\n"
                                       "@type elems: List of Int\n"
                                       "@param block: Geometry block elements belong to"
                                       "@type block: Int"
                                       "@param name: Name of set\n"
                                       "@type name: String\n"
                                       "@return: None\n");
PyObject* VTF_AddGeometrySet(PyObject* self, PyObject* args, PyObject* kwds)
{
  PyVTF* vtf = (PyVTF*)self;
#ifdef ENABLE_VTF
  static const char* keyWords[] = {"elements", "block", "name", NULL};
  PyObject* elemso=NULL;
  int block=0;
  const char* sname=NULL;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Ois",
                                   (char**)keyWords,&elemso,&block,&sname))
    return NULL;

  if (!PyObject_TypeCheck(elemso,&PyList_Type) || !sname)
    return NULL;

  std::vector<int> elemIds;
  for (int i=0;i<PyList_Size(elemso);++i)
    elemIds.push_back(PyInt_AsLong(PyList_GetItem(elemso,i)));

  std::string name(sname);

  if (vtf->sets->find(sname) == vtf->sets->end())
    vtf->sets->insert(std::make_pair(name, new VTFASetBlock(vtf->sets->size(), 0)));

  (*vtf->sets)[name]->AddItems(&elemIds.front(), elemIds.size(), block);
  (*vtf->sets)[name]->SetName(sname);
#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(vtf_addfield__doc__,"Add a field\n"
                                 "@param coefs: Field data\n"
                                 "@type coefs: List of float\n"
                                 "@param block: Block to map field to\n"
                                 "@type block: Int\n"
                                 "@return: Block ID\n"
                                 "@rtype: Int");
PyObject* VTF_AddField(PyObject* self, PyObject* args, PyObject* kwds)
{
  PyVTF* vtf = (PyVTF*)self;
#ifdef ENABLE_VTF
  static const char* keyWords[] = {"coefs", "block", NULL};
  PyObject* coefso=NULL;
  int block=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Oi|i",
                                   (char**)keyWords,&coefso,&block))
    return NULL;

  if (!PyObject_TypeCheck(coefso,&PyList_Type))
    return NULL;

  std::pair<std::vector<float>, int> resVec = PyPointListToVector<float>(coefso);

  if (resVec.second > 1) {
    VTFAResultBlock dBlock(++vtf->nBlock,VTFA_DIM_VECTOR,VTFA_RESMAP_NODE,0);
    dBlock.SetResults3D(&resVec.first.front(),resVec.first.size()/3);
    dBlock.SetMapToBlockID(block);
    vtf->data->WriteBlock(&dBlock);
  } else {
    VTFAResultBlock dBlock(++vtf->nBlock,VTFA_DIM_SCALAR,VTFA_RESMAP_NODE,0);
    dBlock.SetResults1D(&resVec.first.front(),resVec.first.size());
    dBlock.SetMapToBlockID(block);
    vtf->data->WriteBlock(&dBlock);
  }

#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  return Py_BuildValue((char*)"i",vtf->nBlock);
}

PyDoc_STRVAR(vtf_addfieldblocks__doc__,"Add a description block for a field\n"
                                       "@param blocks: Field block numbers\n"
                                       "@type coefs: List of Int\n"
                                       "@param name: Name of field\n"
                                       "@type name: String\n"
                                       "@param comp: Number of components in field\n"
                                       "@type comp: Int=1\n"
                                       "@param displacement: Write as a displacement field\n"
                                       "@type displacement: Boolean\n"
                                       "@return: None\n");
PyObject* VTF_AddFieldBlocks(PyObject* self, PyObject* args, PyObject* kwds)
{
  PyVTF* vtf = (PyVTF*)self;
#ifdef ENABLE_VTF
  static const char* keyWords[] = {"blocks", "name", "comp", "displacement", NULL};
  PyObject* blockso=NULL;
  int comp=1;
  const char* sname=NULL;
  bool displacement=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Os|ib",
                                   (char**)keyWords,&blockso,&sname,&comp,&displacement))
    return NULL;

  if (!PyObject_TypeCheck(blockso,&PyList_Type) || !sname)
    return NULL;

  std::string name(sname);

  std::vector<int> blockIds;
  for (int i=0;i<PyList_Size(blockso);++i)
    blockIds.push_back(PyInt_AsLong(PyList_GetItem(blockso,i)));

  if (comp > 1) {
    if (displacement) {
      if (vtf->dBlocks->find(name) == vtf->dBlocks->end())
        vtf->dBlocks->insert(std::make_pair(name, new VTFADisplacementBlock(++vtf->nBlock)));
      (*vtf->dBlocks)[name]->SetName(sname);
      (*vtf->dBlocks)[name]->SetRelativeDisplacementResults(1);
      (*vtf->dBlocks)[name]->SetResultBlocks(&blockIds.front(),blockIds.size(),vtf->step);
    } else {
      if (vtf->vBlocks->find(name) == vtf->vBlocks->end())
        vtf->vBlocks->insert(std::make_pair(name, new VTFAVectorBlock(++vtf->nBlock)));
      (*vtf->vBlocks)[name]->SetName(sname);
      (*vtf->vBlocks)[name]->SetResultBlocks(&blockIds.front(),blockIds.size(),vtf->step);
    }
  } else {
    if (vtf->sBlocks->find(name) == vtf->sBlocks->end())
      vtf->sBlocks->insert(std::make_pair(name, new VTFAScalarBlock(++vtf->nBlock)));
    (*vtf->sBlocks)[name]->SetName(sname);
    (*vtf->sBlocks)[name]->SetResultBlocks(&blockIds.front(),blockIds.size(),vtf->step);
  }
#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(vtf_addstate__doc__,"Add a state (timestep) descriptor\n"
                                 "@param timelevel: Time level\n"
                                 "@type timelevel: Float\n"
                                 "@return: None\n");
PyObject* VTF_AddState(PyObject* self, PyObject* args, PyObject* kwds)
{
  PyVTF* vtf = (PyVTF*)self;
#ifdef ENABLE_VTF
  static const char* keyWords[] = {"timelevel", NULL};
  double time=0.0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"d",
                                   (char**)keyWords,&time))
    return NULL;

  if (!vtf->state)
    vtf->state = new VTFAStateInfoBlock();
  std::stringstream str;
  str << "Time " << time;
  vtf->state->SetStepData(vtf->step++,str.str().c_str(),time,0);
#else
  std::cerr << "VTF support disabled. No data written" << std::endl;
#endif

  Py_INCREF(Py_None);
  return Py_None;
}

PyMethodDef VTF_methods[] = {
     {(char*)"AddGeometryBlock",      (PyCFunction)VTF_AddGeometryBlock,      METH_VARARGS|METH_KEYWORDS, vtf_addgeometryblock__doc__},
     {(char*)"AddGeometryDescriptor", (PyCFunction)VTF_AddGeometryDescriptor, METH_VARARGS|METH_KEYWORDS, vtf_addgeometrydescriptor__doc__},
     {(char*)"AddGeometrySet",        (PyCFunction)VTF_AddGeometrySet,        METH_VARARGS|METH_KEYWORDS, vtf_addgeometryset__doc__},
     {(char*)"AddField",              (PyCFunction)VTF_AddField,              METH_VARARGS|METH_KEYWORDS, vtf_addfield__doc__},
     {(char*)"AddFieldBlocks",        (PyCFunction)VTF_AddFieldBlocks,        METH_VARARGS|METH_KEYWORDS, vtf_addfieldblocks__doc__},
     {(char*)"AddState",              (PyCFunction)VTF_AddState,              METH_VARARGS|METH_KEYWORDS, vtf_addstate__doc__},
     {NULL,                           NULL,                                   0,                          NULL}
   };

PyNumberMethods VTF_operators = {0};
PySequenceMethods VTF_seq_operators = {0};

PyDoc_STRVAR(VTF__doc__, "A VTF file");
void init_PyVTF_Type()
{
  InitializeTypeObject(&VTF_Type);

  VTF_Type.tp_name = "GoTools.VTF";
  VTF_Type.tp_basicsize = sizeof(PyVTF);
  VTF_Type.tp_dealloc = (destructor)VTF_Dealloc;
  VTF_Type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_CHECKTYPES;
  VTF_Type.tp_doc = VTF__doc__;
  VTF_Type.tp_methods = VTF_methods;
  VTF_Type.tp_as_number = &VTF_operators;
  VTF_Type.tp_as_sequence = &VTF_seq_operators;
  VTF_Type.tp_base = 0;
  VTF_Type.tp_new = VTF_New;
  VTF_Type.tp_str = (reprfunc)VTF_Str;
  PyType_Ready(&VTF_Type);
}

}
}
