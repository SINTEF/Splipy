#include "geomodeller.h"

#include "curve.h"
#include "point.h"
#include "surface.h"
#include "surfacemodel.h"
#include "volume.h"

#include "GoTools/geometry/GoTools.h"

#include <fstream>

#ifdef HAS_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL GEOMOD_ARRAY_API
#include <arrayobject.h>
#endif

GeoModellerState modState;

// Documentation for modules stored here so we don't have to extern them

#include "curvefactory.h"
PyDoc_STRVAR(generate_line__doc__, "Generate a line");
PyDoc_STRVAR(generate_line_segment__doc__, "Generate a line segment");
PyDoc_STRVAR(generate_circle__doc__, "Generate a circle");
PyDoc_STRVAR(generate_circle_segment__doc__, "Generate a circle segment");
PyDoc_STRVAR(generate_ellipse__doc__, "Generate an ellipse");
PyDoc_STRVAR(generate_elliptic_segment__doc__, "Generate an elliptic segment");
PyDoc_STRVAR(generate_helix__doc__, "Generate a helix");
PyDoc_STRVAR(generate_interpolate_curve__doc__, "Generate a curve from a point cloud");

#include "surfacefactory.h"
PyDoc_STRVAR(generate_plane__doc__,"Generate an infinite plane");
PyDoc_STRVAR(generate_sphere_surface__doc__,"Generate a sphere surface");
PyDoc_STRVAR(generate_cylinder_surface__doc__,"Generate a cylinder surface");
PyDoc_STRVAR(generate_cone_surface__doc__,"Generate a cone surface");
PyDoc_STRVAR(generate_torus_surface__doc__,"Generate a torus surface");
PyDoc_STRVAR(generate_circular_disc__doc__,"Generate a circular disc");
PyDoc_STRVAR(generate_sweep_curve_rotational__doc__,"Generate a surface by rotationally sweeping a curve");
PyDoc_STRVAR(generate_sweep_curve_linear__doc__,"Generate a surface by linearly sweeping a curve along a curve");
PyDoc_STRVAR(generate_rectangle__doc__,"Generate a rectangle");
PyDoc_STRVAR(generate_trim_surface__doc__,"Generate a surface by trimming one surface with another");
PyDoc_STRVAR(generate_addloop__doc__,"Generate a surface by adding a loop");
PyDoc_STRVAR(generate_loft_curves__doc__,"Generate a surface by lofting curves");

#include "volumefactory.h"
PyDoc_STRVAR(generate_box__doc__, "Generate a box");
PyDoc_STRVAR(generate_cone__doc__, "Generate a cone");
PyDoc_STRVAR(generate_cylinder__doc__, "Generate a cylinder");
PyDoc_STRVAR(generate_parallelepiped__doc__, "Generate a parallelepiped");
PyDoc_STRVAR(generate_sphere__doc__, "Generate a sphere");
PyDoc_STRVAR(generate_torus__doc__, "Generate a torus");
PyDoc_STRVAR(generate_torus_segment__doc__, "Generate a torus segment");
PyDoc_STRVAR(generate_loft_surfaces__doc__, "Generate a volume by lofting surfaces");
PyDoc_STRVAR(generate_sweep_surface_linear__doc__, "Generate a volume by sweeping a surface along a curve, or curve along a surface");
PyDoc_STRVAR(generate_sweep_surface_rotational__doc__, "Generate a volume by sweeping a rotated surface");

#include "surfacemodelfactory.h"
PyDoc_STRVAR(generate_regularize_surface__doc__, "Regularize a surface");

extern "C"
{

PyDoc_STRVAR(get_version__doc__, "Get GeoModeller module version");
PyObject* GeoMod_GetVersion(PyObject* self)
{
  char tmp[24];
  sprintf(tmp,"%i.%i.%i",GEOMODELLER_VERSION_MAJOR,
                         GEOMODELLER_VERSION_MINOR,
                         GEOMODELLER_VERSION_PATCH);
  return Py_BuildValue((char*)"s", tmp);
}

PyDoc_STRVAR(get_go_version__doc__, "Get GoTools version");
PyObject* GeoMod_GetGoToolsVersion(PyObject* self)
{
  char tmp[24];
  sprintf(tmp,"%i.%i.%i",GO_VERSION_MAJOR,
                         GO_VERSION_MINOR,
                         GO_VERSION_PATCH);
  return Py_BuildValue((char*)"s", tmp);
}

PyDoc_STRVAR(set_dim__doc__, "Set geometry dimensionality");
PyObject* GeoMod_SetDimension(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"dim", NULL };
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&modState.dim))
    return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(set_tolerance__doc__, "Set a tolerance");
PyObject* GeoMod_SetTolerance(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"type", "tolerance", NULL };
  const char* typec;
  double tolerance;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"sd",
                                   (char**)keyWords,&typec, &tolerance))
    return NULL;

  std::string type(typec);

  if (type == "gap" && modState.gapTolerance >= 0)
    modState.gapTolerance = tolerance;

  if (type == "approx" && modState.approxTolerance >= 0)
    modState.approxTolerance = tolerance;

  if (type == "neighbour" && modState.neighbourTolerance >= 0)
    modState.neighbourTolerance = tolerance;

  if (type == "kink" && modState.kinkTolerance >= 0)
    modState.kinkTolerance = tolerance;

  if (type == "bend" && modState.kinkTolerance >= 0)
    modState.bendTolerance = tolerance;

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(set_final_output__doc__, "Set final output");
PyObject* GeoMod_SetFinalOutput(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"file", NULL };
  char* str;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&str))
    return NULL;

  if (str && modState.finalOutput.empty())
    modState.finalOutput = str;

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(set_debug_level__doc__, "Set debug level");
PyObject* GeoMod_SetDebugLevel(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"level", NULL };
  int level;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&level))
    return NULL;

  // a negative debugLevel means it has been overriden from command line
  if (level >= 1 && modState.debugLevel > -1)
    modState.debugLevel = level;

  Py_INCREF(Py_None);
  return Py_None;
}

static void WriteCurve(std::ofstream& g2_file, Curve* curve, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineCurve> crv = convertSplineCurve(curve->data);
    crv->writeStandardHeader(g2_file);
    crv->write(g2_file);
  } else {
    curve->data->writeStandardHeader(g2_file);
    curve->data->write(g2_file);
  }
}

static void WriteSurface(std::ofstream& g2_file, Surface* surface, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineSurface> srf = convertSplineSurface(surface->data);
    srf->writeStandardHeader(g2_file);
    srf->write(g2_file);
  } else {
    surface->data->writeStandardHeader(g2_file);
    surface->data->write(g2_file);
  }
}

static void WriteVolume(std::ofstream& g2_file, Volume* volume, bool convert)
{
  if (convert) {
    shared_ptr<Go::SplineVolume> vol = convertSplineVolume(volume->data);
    if (!vol->isLeftHanded())
    {
      vol = shared_ptr<Go::SplineVolume>(vol->clone());
      vol->reverseParameterDirection(2);
    }
    vol->writeStandardHeader(g2_file);
    vol->write(g2_file);
  } else {
    volume->data->writeStandardHeader(g2_file);
    volume->data->write(g2_file);
  }
}

static void WriteSurfaceModel(std::ofstream& g2_file, SurfaceModel* model, bool convert)
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

static void WriteEntity(std::ofstream& g2_file, PyObject* obj, bool convert)
{
  if (PyObject_TypeCheck(obj,&Curve_Type))
    WriteCurve(g2_file,(Curve*)obj,convert);
  if (PyObject_TypeCheck(obj,&Surface_Type))
    WriteSurface(g2_file,(Surface*)obj,convert);
  if (PyObject_TypeCheck(obj,&Volume_Type))
    WriteVolume(g2_file,(Volume*)obj,convert);
  if (PyObject_TypeCheck(obj,&SurfaceModel_Type))
    WriteSurfaceModel(g2_file,(SurfaceModel*)obj,convert);
}

static void DoWrite(const std::string& fname, PyObject* objectso, bool convert)
{
  std::ofstream g2_file(fname);
  if (PyObject_TypeCheck(objectso,&PyList_Type)) {
    for (int i=0; i < PyList_Size(objectso); ++i) {
      PyObject* obj = PyList_GetItem(objectso,i);
      WriteEntity(g2_file,obj,convert);
    }
  } else if (PyObject_TypeCheck(objectso,&PyTuple_Type)) {
    for (int i=0;i < PyTuple_Size(objectso); ++i) {
      PyObject* obj = PyTuple_GetItem(objectso,i);
      WriteEntity(g2_file,obj,convert);
    }
    PyTuple_GetItem(objectso,0);
  } else
    WriteEntity(g2_file,objectso,convert);
}

PyDoc_STRVAR(writeg2__doc__,"Write results to G2 file");
PyObject* GeoMod_WriteG2(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"filename", "objects", "convert", "level", NULL };
  PyObject* objectso;
  char* fname = 0;  
  bool convert=false;
  int level=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"sO|bi",
                             (char**)keyWords,&fname,&objectso,&convert,&level))
    return NULL;

  if (level <= abs(modState.debugLevel))
    DoWrite(fname,objectso,convert);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(final_output__doc__,"Write final results to G2 file");
PyObject* GeoMod_FinalOutput(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"objects", "convert", NULL };
  PyObject* objectso;
  bool convert=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|b",
                                   (char**)keyWords,&objectso,&convert))
    return NULL;
  if (modState.finalOutput.empty())
    modState.finalOutput = "result.g2";

  DoWrite(modState.finalOutput,objectso,convert);

  Py_INCREF(Py_None);
  return Py_None;
}


PyMethodDef GeoMod_methods[] = {
     // setup and general methods
     {(char*)"GetVersion",            (PyCFunction)GeoMod_GetVersion,        METH_VARARGS, get_version__doc__},
     {(char*)"SetTolerance",          (PyCFunction)GeoMod_SetTolerance,      METH_VARARGS, set_tolerance__doc__},
     {(char*)"GetGoToolsVersion",     (PyCFunction)GeoMod_GetGoToolsVersion, METH_VARARGS, get_go_version__doc__},
     {(char*)"SetDimension",          (PyCFunction)GeoMod_SetDimension, METH_VARARGS|METH_KEYWORDS, set_dim__doc__},
     {(char*)"SetFinalOutput",        (PyCFunction)GeoMod_SetFinalOutput, METH_VARARGS|METH_KEYWORDS, set_final_output__doc__},

     // curve generators
     {(char*)"Circle",                (PyCFunction)Generate_Circle, METH_VARARGS|METH_KEYWORDS, generate_circle__doc__},
     {(char*)"CircleSegment",         (PyCFunction)Generate_CircleSegment, METH_VARARGS|METH_KEYWORDS, generate_circle_segment__doc__},
     {(char*)"Ellipse",               (PyCFunction)Generate_Ellipse, METH_VARARGS|METH_KEYWORDS, generate_ellipse__doc__},
     {(char*)"EllipticSegment",       (PyCFunction)Generate_EllipticSegment, METH_VARARGS|METH_KEYWORDS, generate_elliptic_segment__doc__},
     {(char*)"InterpolateCurve",      (PyCFunction)Generate_InterpolateCurve, METH_VARARGS|METH_KEYWORDS, generate_interpolate_curve__doc__},
     {(char*)"Helix",                 (PyCFunction)Generate_Helix, METH_VARARGS|METH_KEYWORDS, generate_helix__doc__},
     {(char*)"Line",                  (PyCFunction)Generate_Line, METH_VARARGS|METH_KEYWORDS, generate_line__doc__},
     {(char*)"LineSegment",           (PyCFunction)Generate_LineSegment, METH_VARARGS|METH_KEYWORDS, generate_line_segment__doc__},

     // surface generators
     {(char*)"Plane",                 (PyCFunction)Generate_Plane, METH_VARARGS|METH_KEYWORDS, generate_plane__doc__},
     {(char*)"CircularDisc",          (PyCFunction)Generate_CircularDisc, METH_VARARGS|METH_KEYWORDS, generate_plane__doc__},
     {(char*)"SphereSurface",         (PyCFunction)Generate_SphereSurface, METH_VARARGS|METH_KEYWORDS, generate_sphere_surface__doc__},
     {(char*)"CylinderSurface",       (PyCFunction)Generate_CylinderSurface, METH_VARARGS|METH_KEYWORDS, generate_sphere_surface__doc__},
     {(char*)"ConeSurface",           (PyCFunction)Generate_ConeSurface, METH_VARARGS|METH_KEYWORDS, generate_cone_surface__doc__},
     {(char*)"TorusSurface",          (PyCFunction)Generate_TorusSurface, METH_VARARGS|METH_KEYWORDS, generate_torus_surface__doc__},
     {(char*)"Rectangle",             (PyCFunction)Generate_Rectangle, METH_VARARGS|METH_KEYWORDS, generate_rectangle__doc__},
     {(char*)"RotationalCurveSweep",  (PyCFunction)Generate_SweepCurveRotational, METH_VARARGS|METH_KEYWORDS, generate_sweep_curve_rotational__doc__},
     {(char*)"LinearCurveSweep",      (PyCFunction)Generate_SweepCurveLinear, METH_VARARGS|METH_KEYWORDS, generate_sweep_curve_linear__doc__},
     {(char*)"TrimSurface",           (PyCFunction)Generate_TrimSurface, METH_VARARGS|METH_KEYWORDS, generate_trim_surface__doc__},
     {(char*)"AddLoop",               (PyCFunction)Generate_AddLoop, METH_VARARGS|METH_KEYWORDS, generate_addloop__doc__},
     {(char*)"LoftCurves",            (PyCFunction)Generate_LoftCurves, METH_VARARGS|METH_KEYWORDS, generate_loft_curves__doc__},

     // volume generators
     {(char*)"Box",                   (PyCFunction)Generate_Box, METH_VARARGS|METH_KEYWORDS, generate_box__doc__},
     {(char*)"Cone",                  (PyCFunction)Generate_Cone, METH_VARARGS|METH_KEYWORDS, generate_cone__doc__},
     {(char*)"Cylinder",              (PyCFunction)Generate_Cylinder, METH_VARARGS|METH_KEYWORDS, generate_cylinder__doc__},
     {(char*)"LoftSurfaces",          (PyCFunction)Generate_LoftSurfaces, METH_VARARGS|METH_KEYWORDS, generate_loft_surfaces__doc__},
     {(char*)"LinearSurfaceSweep",    (PyCFunction)Generate_SweepSurfaceLinear, METH_VARARGS|METH_KEYWORDS, generate_sweep_surface_linear__doc__},
     {(char*)"RotationalSurfaceSweep",(PyCFunction)Generate_SweepSurfaceRotational, METH_VARARGS|METH_KEYWORDS, generate_sweep_surface_rotational__doc__},
     {(char*)"Parallelepiped",        (PyCFunction)Generate_Parallelepiped, METH_VARARGS|METH_KEYWORDS, generate_parallelepiped__doc__},
     {(char*)"Sphere",                (PyCFunction)Generate_Sphere, METH_VARARGS|METH_KEYWORDS, generate_sphere__doc__},
     {(char*)"Torus",                 (PyCFunction)Generate_Torus, METH_VARARGS|METH_KEYWORDS, generate_torus__doc__},
     {(char*)"TorusSegment",          (PyCFunction)Generate_Torus, METH_VARARGS|METH_KEYWORDS, generate_torus_segment__doc__},

     // surface model generators
     {(char*)"RegularizeSurface",     (PyCFunction)Generate_RegularizeSurface, METH_VARARGS|METH_KEYWORDS, generate_regularize_surface__doc__},

     // I/O
     {(char*)"WriteG2",             (PyCFunction)GeoMod_WriteG2, METH_VARARGS|METH_KEYWORDS, writeg2__doc__},
     {(char*)"FinalOutput",         (PyCFunction)GeoMod_FinalOutput, METH_VARARGS|METH_KEYWORDS, final_output__doc__},
     {(char*)"SetDebugLevel",       (PyCFunction)GeoMod_SetDebugLevel, METH_VARARGS|METH_KEYWORDS, set_debug_level__doc__},
     
     // done - need a null entry for termination
     {NULL,                       NULL,                                  0,            NULL}
   };

PyMODINIT_FUNC
InitGeoModTypes()
{
  init_Curve_Type();
  init_Point_Type();
  init_Surface_Type();
  init_SurfaceModel_Type();
  init_Volume_Type();
}

PyMODINIT_FUNC
registerPythonTypes()
{
  InitGeoModTypes();
  Py_INCREF(&Point_Type);
  PyObject* geoModule;
  geoModule = Py_InitModule((char*)"GoTools", GeoMod_methods);
  if (!geoModule)
    return;
#ifdef USE_NUMPY
  import_array();
#endif
  PyModule_AddObject(geoModule,(char*)"Curve",(PyObject*)&Curve_Type);
  PyModule_AddObject(geoModule,(char*)"Point",(PyObject*)&Point_Type);
  PyModule_AddObject(geoModule,(char*)"Surface",(PyObject*)&Surface_Type);
  PyModule_AddObject(geoModule,(char*)"SurfaceModel",(PyObject*)&SurfaceModel_Type);
  PyModule_AddObject(geoModule,(char*)"Volume",(PyObject*)&Volume_Type);
}

}

// helper functions
Go::Point someNormal(const Go::Point& vec)
{
  if (vec.dimension() == 2)
    return Go::Point(-vec[1], vec[0]);

  Go::Point unit;
  if (fabs(vec[0]) < fabs(vec[1]) && fabs(vec[0]) < fabs(vec[2]))
    unit = Go::Point(1.0, 0.0, 0.0);
  else if (fabs(vec[1]) < fabs(vec[2]))
    unit = Go::Point(0.0, 1.0, 0.0);
  else
    unit = Go::Point(0.0, 0.0, 1.0);

  return vec % unit;
}
