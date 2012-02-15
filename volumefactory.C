#include "volumefactory.h"

#include "curve.h"
#include "geomodeller.h"
#include "point.h"
#include "pyutils.h"
#include "surface.h"
#include "volume.h"

#include "GoTools/trivariate/ConeVolume.h"
#include "GoTools/trivariate/CylinderVolume.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/trivariate/Parallelepiped.h"
#include "GoTools/trivariate/SphereVolume.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/TorusVolume.h"
#include "GoTools/trivariate/SweepVolumeCreator.h"

#ifdef HAS_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL GEOMOD_ARRAY_API
#define NO_IMPORT_ARRAY
#include <arrayobject.h>
#endif

PyObject* Generate_Box(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"corner", "axis_x", "axis_y", "length_x", "length_y", "length_z", NULL };
  PyObject* cornero;
  PyObject* axisxo;
  PyObject* axisyo;
  double length_x, length_y, length_z;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOddd",
                                   (char**)keyWords,&cornero,&axisxo,
                                   &axisyo,&length_x,&length_y,&length_z))
    return NULL;
       
  shared_ptr<Go::Point> corner = PyObject_AsGoPoint(cornero);
  shared_ptr<Go::Point> axisx  = PyObject_AsGoPoint(axisxo);
  shared_ptr<Go::Point> axisy  = PyObject_AsGoPoint(axisyo);
  if (!corner || !axisx || !axisy)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  Go::Point dir_v_ortog = *axisy-(*axisx*(*axisy))/axisx->length2()*(*axisx);
  result->data.reset(new Go::Parallelepiped(*corner,*axisx, dir_v_ortog,
                                            *axisx % dir_v_ortog, 
                                            length_x, length_y, length_z));

  return (PyObject*)result;
}

PyObject* Generate_Cone(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"apex", "axis", "angle", "height", "radius", NULL };
  PyObject* apexo;
  PyObject* axiso;
  double angle, height = std::numeric_limits<double>::infinity(), radius=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOd|dd",
                                   (char**)keyWords,&apexo,&axiso,
                                                    &angle,&height,&radius))
    return NULL;
       
  shared_ptr<Go::Point> apex = PyObject_AsGoPoint(apexo);
  shared_ptr<Go::Point> axis  = PyObject_AsGoPoint(axiso);
  if (!apex || !axis)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  result->data.reset(new Go::ConeVolume(radius, *apex, *axis,
                                        someNormal(*axis), angle));
  static_pointer_cast<Go::ConeVolume>(result->data)->setParameters(0.0,height,2);

  return (PyObject*)result;
}

PyObject* Generate_Cylinder(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "boundary_point", "normal", "height", NULL };
  PyObject* centero;
  PyObject* bpointo;
  PyObject* normalo;
  double height = -1.f;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO|d",
                                   (char**)keyWords,&centero,&bpointo,
                                                    &normalo,&height))
    return NULL;
       
  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> bpoint = PyObject_AsGoPoint(bpointo);
  shared_ptr<Go::Point> normal = PyObject_AsGoPoint(normalo);
  if (!center || !bpoint || !normal)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  Go::Point v = *bpoint-*center;
  Go::Point x_axis = v- ((v*(*normal))/normal->length2())*(*normal);
  double rad = x_axis.length();
  result->data.reset(new Go::CylinderVolume(*center,rad,*normal,x_axis));
  if (height > 0)
    static_pointer_cast<Go::CylinderVolume>(result->data)->setParameters(0,height,2);

  return (PyObject*)result;
}

PyObject* Generate_Parallelepiped(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"corner", "axis_x", "axis_y", "axis_z", "length_x", "length_y", "length_z", NULL };
  PyObject* cornero;
  PyObject* axisxo;
  PyObject* axisyo;
  PyObject* axiszo;
  double length_x, length_y, length_z;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOOOddd",
                                   (char**)keyWords,&cornero,&axisxo,
                                                    &axisyo,&axiszo,&length_x,
                                                    &length_y,&length_z))
    return NULL;
       
  shared_ptr<Go::Point> corner = PyObject_AsGoPoint(cornero);
  shared_ptr<Go::Point> axisx  = PyObject_AsGoPoint(axisxo);
  shared_ptr<Go::Point> axisy  = PyObject_AsGoPoint(axisyo);
  shared_ptr<Go::Point> axisz  = PyObject_AsGoPoint(axiszo);
  if (!corner || !axisx || !axisy || !axisz)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  
  result->data.reset(new Go::Parallelepiped(*corner,*axisx, *axisy, *axisz, 
                                            length_x, length_y, length_z));

  return (PyObject*)result;
}

PyObject* Generate_Sphere(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "radius", NULL };
  PyObject* centero;
  double radius;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"Od",
                                   (char**)keyWords,&centero,&radius))
    return NULL;
       
  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  if (!center)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);

  Go::Point x_axis(1.f,0.f,0.f);
  Go::Point z_axis(0.f,0.f,1.f);
  
  result->data.reset(new Go::SphereVolume(radius,*center,x_axis,z_axis));

  return (PyObject*)result;
}

PyObject* Generate_Torus(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"center", "axis", "major_radius", "minor_radius", NULL };
  PyObject* centero;
  PyObject* axiso;
  double major_radius, minor_radius;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOdd",
                                   (char**)keyWords,&centero,&axiso,
                                                    &major_radius,&minor_radius))
    return NULL;
       
  shared_ptr<Go::Point> center = PyObject_AsGoPoint(centero);
  shared_ptr<Go::Point> axis   = PyObject_AsGoPoint(axiso);
  if (!center || !axis)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);

  result->data.reset(new Go::TorusVolume(major_radius,minor_radius, *center,
                                         *axis, someNormal(*axis)));

  return (PyObject*)result;
}

PyObject* Generate_LoftSurfaces(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surfaces", NULL };
  PyObject* surfaceso;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O",
                                   (char**)keyWords,&surfaceso))
    return NULL;

  if (!PyObject_TypeCheck(surfaceso,&PyList_Type))
    return NULL;

  std::vector<shared_ptr<Go::SplineSurface> > surfaces;
  for (int i=0; i < PyList_Size(surfaceso); ++i) {
    PyObject* surfo = PyList_GetItem(surfaceso,i);
    shared_ptr<Go::ParamSurface> surf = PyObject_AsGoSurface(surfo);
    if (!surf)
      continue;
    surfaces.push_back(convertSplineSurface(surf));
  }
  if (surfaces.size() < 2)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data.reset(Go::LoftVolumeCreator::loftVolume(surfaces.begin(),surfaces.size()));

  return (PyObject*)result;
}

PyObject* Generate_SweepSurfaceLinear(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "curve", "point", NULL };
  PyObject* o1;
  PyObject* o2;
  PyObject* pointo;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO",
                                   (char**)keyWords,&o1,&o2,&pointo))
    return NULL;

  shared_ptr<Go::Point> point = PyObject_AsGoPoint(pointo);
  if (!point)
    return NULL;

  Volume* result = NULL;
  shared_ptr<Go::ParamSurface> s = PyObject_AsGoSurface(o1);
  if (!s)
    s = PyObject_AsGoSurface(o2);
  shared_ptr<Go::ParamCurve> c = PyObject_AsGoCurve(o2);
  if (!c)
    c = PyObject_AsGoCurve(o1);

  if (s && c) {
    result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
    result->data.reset(Go::SweepVolumeCreator::linearSweptVolume(*convertSplineSurface(s),
                                                                   *convertSplineCurve(c),*point));
  }

  return (PyObject*)result;
}

PyObject* Generate_SweepSurfaceRotational(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"surface", "position", "axis", "angle", NULL };
  PyObject* surfaceo;
  PyObject* poso;
  PyObject* axiso;
  double angle=2*M_PI;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"OOO|d",
                                   (char**)keyWords,&surfaceo,&poso,
                                                    &axiso,&angle))
    return NULL;

  shared_ptr<Go::ParamSurface> surface = PyObject_AsGoSurface(surfaceo);
  shared_ptr<Go::Point> axis = PyObject_AsGoPoint(axiso);
  shared_ptr<Go::Point> pos = PyObject_AsGoPoint(poso);
  if (!surface || !axis || !pos)
    return NULL;

  Volume* result = (Volume*)Volume_Type.tp_alloc(&Volume_Type,0);
  result->data.reset(Go::SweepVolumeCreator::rotationalSweptVolume(*convertSplineSurface(surface),
                                                                   angle,*pos,*axis));

  return (PyObject*)result;
}
