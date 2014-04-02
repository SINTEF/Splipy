#include "curve.h"
#include "surface.h"
#include "pyutils.h"

#include <iostream>
#include <fstream>

#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/igeslib/IGESconverter.h>


PyObject* DoReadIGES(const std::string& fname, const std::string& type)
{
  std::fstream in(fname);
  if(!in.good()) {
    std::cerr << "Cannot open file " << fname << std::endl;
    return NULL;
  }

  Go::IGESconverter conv;
  conv.readIGES(in);
  in.close();

  bool surfaces(true);
  bool curves(true);
  if (!type.empty() && !type.compare("curves"))
    surfaces = false;
  if (!type.empty() && !type.compare("surfaces"))
    curves = false;

  PyObject* result=PyList_New(0);
  std::vector<shared_ptr<Go::GeomObject> > entities = conv.getGoGeom();
  for(size_t i=0; i<entities.size(); i++) {
    if(entities[i]->instanceType() == Go::Class_SplineCurve && curves) {
      Curve* curr = (Curve*)Curve_Type.tp_alloc(&Curve_Type, 0);
      curr->data.reset((Go::ParamCurve*) entities[i]->clone());
      PyList_Append(result, (PyObject*) curr);
    } else if(entities[i]->instanceType() == Go::Class_SplineSurface && surfaces) {
      Surface* curr = (Surface*) Surface_Type.tp_alloc(&Surface_Type, 0);
      curr->data.reset((Go::ParamSurface*) entities[i]->clone());
      PyList_Append(result, (PyObject*) curr);
    } else {
      // std::cout << "Element not read: " << entities[i]->instanceType() << std::endl;
    }
  }

  return result;
}

