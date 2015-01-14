#pragma once

#include <Python.h>
#include "GoTools/utils/config.h" // for shared_ptr
#ifdef ENABLE_VTF
#include "VTFAPI.h"
#include "VTOAPIPropertyIDs.h"
#endif

#include <map>

namespace GeoModeller {

extern "C" {
  typedef struct {
    PyObject_HEAD
    int nBlock; //!< Running block counter
#ifdef ENABLE_VTF
    shared_ptr<VTFAFile> data;
    int step;   //!< Current time step
    VTFAStateInfoBlock* state;
    std::map<std::string, VTFADisplacementBlock*>* dBlocks; //!< displacement field blocks
    std::map<std::string, VTFAVectorBlock*>* vBlocks;       //!< vector field blocks
    std::map<std::string, VTFAScalarBlock*>* sBlocks;       //!< scalar field blocks
    std::map<std::string, VTFASetBlock*>* sets;
#endif
  } PyVTF;

  void init_PyVTF_Type();

  extern PyTypeObject VTF_Type;
}

}
