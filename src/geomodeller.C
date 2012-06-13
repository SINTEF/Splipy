#include "geomodeller.h"

#include "curve.h"
#include "point.h"
#include "pyutils.h"
#include "surface.h"
#include "surfacemodel.h"
#include "volume.h"

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <fstream>

GeoModellerState modState;

#include "curvefactory.h"
#include "surfacefactory.h"
#include "volumefactory.h"
#include "surfacemodelfactory.h"

#ifdef ENABLE_OPENNURBS
#include "3dmutils.h"
#endif

extern "C"
{

PyDoc_STRVAR(get_go_version__doc__, "Get GoTools version\n"
                                    "@return: String with version in format major.minor.patch");
PyObject* GeoMod_GetVersion(PyObject* self)
{
  char tmp[24];
  sprintf(tmp,"%i.%i.%i",GO_VERSION_MAJOR,
                         GO_VERSION_MINOR,
                         GO_VERSION_PATCH);
  return Py_BuildValue((char*)"s", tmp);
}

PyDoc_STRVAR(set_dim__doc__, "Set geometry dimensionality\n"
                             "@param dimension: Dimensionality to set\n"
                             "@type dimension: int (2 or 3)\n"
                             "@return: None");
PyObject* GeoMod_SetDimension(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"dimension", NULL };
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&modState.dim))
    return NULL;

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(set_tolerance__doc__, "Set a tolerance\n"
                                   "@param type: Tolerance that should be changed\n"
                                   "@type type: 'gap', 'approx', 'neighbour', 'kink', or 'bend'\n"
                                   "@param tolerance: Tolerance to set\n"
                                   "@type tolerance: float\n"
                                   "@return: None");
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

PyDoc_STRVAR(set_final_output__doc__, "Set final output file\n"
                                      "@param file: The output file\n"
                                      "@type file:  string\n"
                                      "@return:     None");
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

PyDoc_STRVAR(set_debug_level__doc__, "Set debug level\n"
                                     "@param level: The requested debug level\n"
                                     "@type  level: int >= 1\n"
                                     "@return:      None\n"
                                     "\n");
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

// these statics are needed to handle the dynamic type of the input parameters
// we may be given a list of objects, or a single object
static void WriteEntity(std::ofstream& g2_file, PyObject* obj, bool convert)
{
  if (PyObject_TypeCheck(obj,&Curve_Type))
    WriteCurveG2(g2_file,(Curve*)obj,convert);
  if (PyObject_TypeCheck(obj,&Surface_Type))
    WriteSurfaceG2(g2_file,(Surface*)obj,convert);
  if (PyObject_TypeCheck(obj,&Volume_Type))
    WriteVolumeG2(g2_file,(Volume*)obj,convert);
  if (PyObject_TypeCheck(obj,&SurfaceModel_Type))
    WriteSurfaceModelG2(g2_file,(SurfaceModel*)obj,convert);
  g2_file << std::endl;
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

static void WriteEntity3DM(ONX_Model& model, PyObject* obj)
{
  if (PyObject_TypeCheck(obj,&Curve_Type)) {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = GoCurveToONCurve(obj);
    mo.m_bDeleteObject = true;
    mo.m_attributes.m_layer_index = 0;
  }
  if (PyObject_TypeCheck(obj,&Surface_Type)) {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = GoSurfaceToONSurface(obj);
    mo.m_bDeleteObject = true;
    mo.m_attributes.m_layer_index = 0;
  }
}

static void DoWrite3DM(const std::string& fname, PyObject* objectso)
{
  ONX_Model model; 
  model.m_properties.m_RevisionHistory.NewRevision();
  model.m_properties.m_Application.m_application_name = "GeoModeller";
  model.m_properties.m_Application.m_application_URL = "http://sintef.no";
  model.m_properties.m_Application.m_application_details = "The SINTEF ICT GoTools/OpenNURBS based python bindings";
  
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0;
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01;

  if (PyObject_TypeCheck(objectso,&PyList_Type)) {
    for (int i=0; i < PyList_Size(objectso); ++i) {
      PyObject* obj = PyList_GetItem(objectso,i);
      WriteEntity3DM(model,obj);
    }
  } else if (PyObject_TypeCheck(objectso,&PyTuple_Type)) {
    for (int i=0;i < PyTuple_Size(objectso); ++i) {
      PyObject* obj = PyTuple_GetItem(objectso,i);
      WriteEntity3DM(model,obj);
    }
    PyTuple_GetItem(objectso,0);
  } else
    WriteEntity3DM(model,objectso);

  FILE* fp = ON::OpenFile(fname.c_str(), "wb");
  ON_BinaryFile archive(ON::write3dm, fp);
  ON_TextLog error_log;
  model.Polish();
  model.Write(archive,4,"",&error_log);
  ON::CloseFile(fp);
}

PyDoc_STRVAR(write3dm__doc__,"Write entities to 3DM file\n"
                             "@param filename: The file to write\n"
                             "@type  filename: string\n"
                             "@param entities: The entities to write to file\n"
                             "@type  entities: Curve, Surface, Volume, SurfaceModel or a list of these\n"
                             "@param    level: Debug level to store file at\n"
                             "@type     level: int >= 1\n"
                             "@return: None");
PyObject* GeoMod_Write3DM(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"filename", "entities", "level", NULL };
  PyObject* objectso;
  char* fname = 0;  
  int level=1;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"sO|i",
                             (char**)keyWords,&fname,&objectso,&level))
    return NULL;

  if (level <= abs(modState.debugLevel))
    DoWrite3DM(fname,objectso);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(writeg2__doc__,"Write entities to G2 file\n"
                            "@param filename: The file to write\n"
                            "@type  filename: string\n"
                            "@param entities: The entities to write to file\n"
                            "@type  entities: Curve, Surface, Volume, SurfaceModel or a list of these\n"
                            "@param  convert: Convert to a spline entity before saving\n"
                            "@type   convert: bool\n"
                            "@param    level: Debug level to store file at\n"
                            "@type     level: int >= 1\n"
                            "@return: None");
PyObject* GeoMod_WriteG2(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"filename", "entities", "convert", "level", NULL };
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

PyDoc_STRVAR(readg2__doc__,"Read entities from G2 file\n"
                            "@param filename: The file to read\n"
                            "@type  filename: string\n"
                            "@return: Curve, Surface, Volume, SurfaceModel or a list of these");
PyObject* GeoMod_ReadG2(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"filename", NULL };
  PyObject* objectso;
  char* fname = 0;  
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&fname))
    return NULL;

  std::ifstream f(fname);
  PyObject* result=NULL;
  PyObject* curr=NULL;
  while (!f.eof() && f.good()) {
    if (curr) {
      if (!result)
        result = PyList_New(0);
      PyList_Append(result,curr);
    }
    Go::ObjectHeader header;
    f >> header;
    if (header.classType() == Go::SplineCurve::classType())
      curr = ReadG2<Curve, Go::SplineCurve>(f,Curve_Type);
    if (header.classType() == Go::SplineSurface::classType())
      curr = ReadG2<Surface, Go::SplineSurface>(f,Surface_Type);
    if (header.classType() == Go::SplineVolume::classType())
      curr = ReadG2<Volume, Go::SplineVolume>(f,Volume_Type);
    char c;
    while (f.get(c)) {
      if (!isspace(c)) {
        f.putback(c);
        break;
      }
    }
    f.peek(); // and peek next byte to update flags
  }
  if (!result)
    result = curr;
  if (result && PyObject_TypeCheck(result,&PyList_Type))
    PyList_Append(result,curr);

  if (!result)
    std::cerr << "Failed to load " << fname << std::endl;

  return result;
}

PyDoc_STRVAR(read3dm__doc__,"Read entities from a 3DM file\n"
                            "@param filename: The file to read\n"
                            "@type  filename: string\n"
                            "@param type: Only read objects of this type\n"
                            "@type type: string (Curves or Surfaces)\n"
                            "@return: Curve, Surface or a list of these");
PyObject* GeoMod_Read3DM(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_OPENNURBS
  std::cerr << "Compiled without OpenNURBS support, no data read" << std::endl;
  Py_INCREF(Py_None);
  return Py_None;
#else
  static const char* keyWords[] = {"filename", "type", NULL };
  PyObject* objectso;
  char* fname = 0;
  char* type = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s|s",
                                   (char**)keyWords,&fname,&type))
    return NULL;

  FILE* archive_fp = ON::OpenFile(fname,"rb");
  if (!archive_fp)
    return NULL;

  ON_BinaryFile archive(ON::read3dm, archive_fp);
  ONX_Model model;

  ON_TextLog dump_to_stdout;
  if (!model.Read(archive,&dump_to_stdout)) {
    std::cerr << "Failed to read 3DM file " << fname << std::endl;
    return NULL;
  }
  ON::CloseFile(archive_fp);

  bool surfaces(true);
  bool curves(true);
  if (type && !strcasecmp(type,"curves"))
    surfaces = false;
  if (type && !strcasecmp(type,"surfaces"))
    curves = false;

  PyObject* result=NULL;
  PyObject* curr=NULL;
  for (int i=0;i<model.m_object_table.Count();++i) {
    if (curr) {
      if (!result)
        result = PyList_New(0);
      PyList_Append(result,curr);
      curr = NULL;
    }
    if (curves && model.m_object_table[i].m_object->ObjectType() == ON::curve_object)
      curr = (PyObject*)ONCurveToGoCurve((const ON_Curve*)model.m_object_table[i].m_object);
    if (surfaces && model.m_object_table[i].m_object->ObjectType() == ON::brep_object) {
      ON_Brep *brep_obj = (ON_Brep*) model.m_object_table[i].m_object;
      int sc=0;
      while (brep_obj->FaceIsSurface(sc)) {
        if (curr) {
          if (!result)
            result = PyList_New(0);
          PyList_Append(result,curr);
        }
        ON_Surface *surf = brep_obj->m_S[sc++];
        curr = (PyObject*)ONSurfaceToGoSurface(surf);
      }
    }
  }
  model.Destroy();
  if (!result)
    result = curr;
  if (PyObject_TypeCheck(result,&PyList_Type))
    PyList_Append(result,curr);

  return result;
#endif
}

PyDoc_STRVAR(final_output__doc__,"Write final entities to G2 file\n"
                                 "@param entities: The entities to write to file\n"
                                 "@type  entities: Curve, Surface, Volume, SurfaceModel or a list of these\n"
                                 "@param  convert: Convert to a spline entity before saving\n"
                                 "@type   convert: bool\n"
                                 "@return: None");
PyObject* GeoMod_FinalOutput(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"entities", "convert", NULL };
  PyObject* objectso;
  bool convert=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"O|b",
                                   (char**)keyWords,&objectso,&convert))
    return NULL;
  if (modState.finalOutput.empty())
    modState.finalOutput = "result.g2";

  if (modState.convertSpline > 0)
    convert = modState.convertSpline==1?false:true;

  DoWrite(modState.finalOutput,objectso,convert);

  Py_INCREF(Py_None);
  return Py_None;
}


PyMethodDef GeoMod_methods[] = {
     // setup and general methods
     {(char*)"GetVersion",            (PyCFunction)GeoMod_GetVersion,        METH_VARARGS,               get_go_version__doc__},
     {(char*)"SetTolerance",          (PyCFunction)GeoMod_SetTolerance,      METH_VARARGS|METH_KEYWORDS, set_tolerance__doc__},
     {(char*)"SetDimension",          (PyCFunction)GeoMod_SetDimension,      METH_VARARGS|METH_KEYWORDS, set_dim__doc__},
     {(char*)"SetFinalOutput",        (PyCFunction)GeoMod_SetFinalOutput,    METH_VARARGS|METH_KEYWORDS, set_final_output__doc__},

     // I/O
     {(char*)"FinalOutput",           (PyCFunction)GeoMod_FinalOutput,       METH_VARARGS|METH_KEYWORDS, final_output__doc__},
     {(char*)"Read3DM",               (PyCFunction)GeoMod_Read3DM,           METH_VARARGS|METH_KEYWORDS, read3dm__doc__},
     {(char*)"ReadG2",                (PyCFunction)GeoMod_ReadG2,            METH_VARARGS|METH_KEYWORDS, readg2__doc__},
     {(char*)"SetDebugLevel",         (PyCFunction)GeoMod_SetDebugLevel,     METH_VARARGS|METH_KEYWORDS, set_debug_level__doc__},
     {(char*)"WriteG2",               (PyCFunction)GeoMod_WriteG2,           METH_VARARGS|METH_KEYWORDS, writeg2__doc__},
     {(char*)"Write3DM",              (PyCFunction)GeoMod_Write3DM,          METH_VARARGS|METH_KEYWORDS, write3dm__doc__},
     
     // done - need a null entry for termination
     {NULL,                       NULL,                                  0,            NULL}
   };

void InitGoToolsTypes()
{
  init_Curve_Type();
  init_CurveFactory_Module();
  init_Point_Type();
  init_Surface_Type();
  init_SurfaceFactory_Module();
  init_SurfaceModel_Type();
  init_SurfaceModelFactory_Module();
  init_Volume_Type();
  init_VolumeFactory_Module();
}

PyDoc_STRVAR(gotools_module__doc__,"GoTools bindings");

PyMODINIT_FUNC
initGoTools()
{
  InitGoToolsTypes();
  PyObject* geoModule;
  geoModule = Py_InitModule3((char*)"GoTools", GeoMod_methods,gotools_module__doc__);
  if (!geoModule)
    return;
  PyModule_AddObject(geoModule,(char*)"Curve",(PyObject*)&Curve_Type);
  PyModule_AddObject(geoModule,(char*)"CurveFactory",CurveFactory_module);
  PyModule_AddObject(geoModule,(char*)"Point",(PyObject*)&Point_Type);
  PyModule_AddObject(geoModule,(char*)"Surface",(PyObject*)&Surface_Type);
  PyModule_AddObject(geoModule,(char*)"SurfaceFactory",SurfaceFactory_module);
  PyModule_AddObject(geoModule,(char*)"SurfaceModel",(PyObject*)&SurfaceModel_Type);
  PyModule_AddObject(geoModule,(char*)"SurfaceModelFactory",SurfaceModelFactory_module);
  PyModule_AddObject(geoModule,(char*)"Volume",(PyObject*)&Volume_Type);
  PyModule_AddObject(geoModule,(char*)"VolumeFactory",VolumeFactory_module);
  modState.addInfo(geoModule);
}

}

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
