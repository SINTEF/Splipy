#include "geomodeller.h"

#include "curve.h"
#include "point.h"
#include "preprocess.h"
#include "pyutils.h"
#include "surface.h"
#include "surfacemodel.h"
#include "volumemodel.h"
#include "stlutils.h"
#include "volume.h"

#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <fstream>
#include <iostream>

#include "curvefactory.h"
#include "surfacefactory.h"
#include "volumefactory.h"
#include "surfacemodelfactory.h"
#include "igesutils.h"

#ifdef ENABLE_HDF5
#include "hdf5utils.h"
#endif

#ifdef ENABLE_OPENNURBS
#include "3dmutils.h"
#endif

namespace GeoModeller {

GeoModellerState modState;

extern "C"
{

PyDoc_STRVAR(get_go_version__doc__, "Get GoTools version\n"
                                    "@return: Version in format major.minor.patch\n"
                                    "@rtype: string");
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
                                   "@type type: 'gap', 'approx', 'neighbour', 'kink', 'bend' or 'refine'\n"
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

  if (type == "refine" && modState.refineTolerance >= 0)
    modState.refineTolerance = tolerance;

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(get_tolerance__doc__, "Get a tolerance\n"
                                   "@param type: Requested tolerance\n"
                                   "@type type: 'gap', 'approx', 'neighbour', 'kink', 'bend' or 'refine'\n"
                                   "@return: Tolerance\n"
                                   "@rtype: Float");
PyObject* GeoMod_GetTolerance(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"type", NULL };
  const char* typec;
  double tolerance;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s",
                                   (char**)keyWords,&typec, &tolerance))
    return NULL;

  std::string type(typec);

  double result=0.0;

  if (type == "gap")
    result = modState.gapTolerance;

  if (type == "approx")
    result = modState.approxTolerance;

  if (type == "neighbour")
    result = modState.neighbourTolerance;

  if (type == "kink")
    result = modState.kinkTolerance;

  if (type == "bend")
    result = modState.bendTolerance;

  if (type == "refine")
    result = modState.refineTolerance;


  return Py_BuildValue((char*)"d", result);
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

PyDoc_STRVAR(get_processor_count__doc__, "Get processor count\n"
                                         "@return: Number of processes\n"
                                         "@rtype: Integer");
PyObject* GeoMod_GetProcessorCount(PyObject* self)
{
  return Py_BuildValue((char*)"i", modState.procCount);
}

PyDoc_STRVAR(set_processor_count__doc__, "Set processor count\n"
                                         "@param procs: The number of processors\n"
                                         "@type  procs: int >= 1\n"
                                         "@return:      None\n"
                                         "\n");
PyObject* GeoMod_SetProcessorCount(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"procs", NULL };
  int procs;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"i",
                                   (char**)keyWords,&procs))
    return NULL;

  if (procs >= 1)
    modState.procCount = procs;

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
  if (PyObject_TypeCheck(obj,&VolumeModel_Type))
    WriteVolumeModelG2(g2_file,(VolumeModel*)obj,convert);
  g2_file << std::endl;
}

static void DoWrite(std::string fname, PyObject* objectso, bool convert)
{
  std::ios_base::openmode mode = std::ios_base::out;
  if (fname == "splinegui") {
    fname = "/tmp/sgui.fifo";
    mode = std::ios_base::out|std::ios_base::app;
  }
  std::ofstream g2_file(fname.c_str(),mode);
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

static int WriteSTLEntity(std::ofstream& stl_file, bool ascii, PyObject* obj, int resolution[3])
{
  if (PyObject_TypeCheck(obj,&Surface_Type))
    return WriteSurfaceSTL(stl_file, (Surface*)obj, ascii, resolution);
  if (PyObject_TypeCheck(obj,&Volume_Type))
    return WriteVolumeSTL(stl_file, (Volume*)obj, ascii, resolution);

  return 0;
}

static void DoWriteSTL(std::string fname, bool ascii, PyObject* objectso, int resolution[3])
{
  // open file for writing
  std::ios_base::openmode mode = std::ios_base::out;
  if (!ascii)
    mode = std::ios_base::out|std::ios_base::binary;
  std::ofstream stl_file(fname.c_str(),mode);
  unsigned int nTriangles = 0;

  // write header
  if(ascii) {
    stl_file << "solid MadeByGeomodeller" << std::endl;
    stl_file << std::setprecision(15)     << std::endl;
  } else {
    char header[80];
    for(int i=0; i<80; i++) header[i] = 0;
    sprintf(header, "BINARY STL BY GEOMODELER");
    stl_file.write(header, 80);
    // make room for nTriangles here, insert it later
    stl_file.write((char*) &nTriangles, sizeof(unsigned int)) ; 
  }

  if (PyObject_TypeCheck(objectso,&PyList_Type)) {
    for (int i=0; i < PyList_Size(objectso); ++i) {
      PyObject* obj = PyList_GetItem(objectso,i);
      nTriangles += WriteSTLEntity(stl_file, ascii, obj, resolution);
    }
  } else if (PyObject_TypeCheck(objectso,&PyTuple_Type)) {
    for (int i=0;i < PyTuple_Size(objectso); ++i) {
      PyObject* obj = PyTuple_GetItem(objectso,i);
      nTriangles += WriteSTLEntity(stl_file, ascii, obj, resolution);
    }
    PyTuple_GetItem(objectso,0);
  } else
    nTriangles += WriteSTLEntity(stl_file, ascii, objectso, resolution);

  // fill the last part of the header
  if(!ascii) {
    stl_file.seekp(80, std::ios_base::beg);
    stl_file.write((char*) &nTriangles, sizeof(unsigned int)) ; 
  }
}

PyDoc_STRVAR(write3dm__doc__,"Write entities to native Rhino 3DM file\n"
                             "@param filename: The file to write\n"
                             "@type  filename: String\n"
                             "@param entities: The entities to write to file\n"
                             "@type  entities: Curve, Surface or a list of these\n"
                             "@param    level: Debug level to store file at\n"
                             "@type     level: Int >= 1\n"
                             "@return: None");
PyObject* GeoMod_Write3DM(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_OPENNURBS
  std::cerr << "Compiled without OpenNURBS support, no data written" << std::endl;
  Py_INCREF(Py_None);
  return Py_None;
#else
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
#endif
}

PyDoc_STRVAR(writeg2__doc__,"Write entities to native GoTools G2 file\n"
                            "@param filename: The file to write\n"
                            "@type  filename: String\n"
                            "@param entities: The entities to write to file\n"
                            "@type  entities: Curve, Surface, Volume, SurfaceModel, VolumeModel or a list of these\n"
                            "@param  convert: Convert to a spline entity before saving\n"
                            "@type   convert: Boolean\n"
                            "@param    level: Debug level to store file at\n"
                            "@type     level: Int >= 1\n"
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

PyDoc_STRVAR(writestl__doc__,"Write entities to STL file\n"
                             "@param filename   : The file to write\n"
                             "@type  filename   : String\n"
                             "@param entities   : The entities to write to file\n"
                             "@type  entities   : Surface, Volume or a list of these\n"
                             "@param resolution : Number of tesselation points used for triangles\n"
                             "@type  resolution : Int or list of Int\n"
                             "@param ascii      : True to store ascii file instead of binary\n"
                             "@type  ascii      : Boolean\n"
                             "@return: None");
PyObject* GeoMod_WriteSTL(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"filename", "entities", "resolution", "ascii", NULL };
  PyObject* objectso, *objectRes;
  char* fname = 0;  
  bool ascii=false;
  int resolution[] = {10,10,10};
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"sO|Ob",
                             (char**)keyWords,&fname,&objectso,&objectRes,&ascii))
    return NULL;
  if (PyObject_TypeCheck(objectRes, &PyList_Type)) {
    for (int i=0; i < PyList_Size(objectRes); ++i) {
      PyObject* obj = PyList_GetItem(objectRes,i);
      if(PyObject_TypeCheck(obj, &PyInt_Type) ) {
        resolution[i] = PyInt_AsLong(obj);
      } else if(PyObject_TypeCheck(obj, &PyFloat_Type) ) {
        resolution[i] = PyFloat_AsDouble(obj);
      } else {
        std::cerr << "Invalid resolution. No file written" << std::endl;
        return Py_None;
      }
    }
  } else if(PyObject_TypeCheck(objectRes, &PyInt_Type) ) {
    for(int i=0; i<3; i++)
      resolution[i] = PyInt_AsLong(objectRes);
  } else if(PyObject_TypeCheck(objectRes, &PyFloat_Type) ) {
    for(int i=0; i<3; i++)
      resolution[i] = PyFloat_AsDouble(objectRes);
  } else {
    std::cerr << "Invalid resolution. No file written" << std::endl;
    return Py_None;
  }

  DoWriteSTL(fname,ascii,objectso,resolution);

  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(readg2__doc__,"Read entities from G2 file\n"
                            "@param filename: The file to read\n"
                            "@type  filename: string\n"
                            "@return: The requested entities\n"
                            "@rtype: Curve, Surface, Volume, SurfaceModel, VolumeModel or a list of these");
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

PyDoc_STRVAR(readiges__doc__,"Read entities from an IGES file\n"
                             "@param filename: The file to read\n"
                             "@type  filename: string\n"
                             "@param type: Only read objects of this type\n"
                             "@type type: string (Curves or Surfaces)\n"
                             "@return: The requested entities\n"
                             "@rtype: Curve, Surface or a list of these");
PyObject* GeoMod_ReadIGES(PyObject* self, PyObject* args, PyObject* kwds)
{
  static const char* keyWords[] = {"filename", "type", NULL };
  PyObject* objectso;
  char* fname = 0;
  char* type = 0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"s|s",
                                   (char**)keyWords,&fname,&type))
    return NULL;

  return DoReadIGES(fname,type?type:"");
}

PyDoc_STRVAR(read3dm__doc__,"Read entities from a 3DM file\n"
                            "@param filename: The file to read\n"
                            "@type  filename: string\n"
                            "@param type: Only read objects of this type\n"
                            "@type type: string (Curves or Surfaces)\n"
                            "@return: The requested entities\n"
                            "@rtype: Curve, Surface or a list of these");
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

  return DoRead3DM(fname,type?type:"");
#endif
}

PyDoc_STRVAR(read_hdf5field__doc__,"Read a field from a HDF5 file\n"
                                   "@param filename: The file to read\n"
                                   "@type  filename: string\n"
                                   "@param fieldname: Field name\n"
                                   "@type fieldname: string \n"
                                   "@param groupname: Group to read dataset from\n"
                                   "@type: groupname: string\n"
                                   "@return: The requested field\n"
                                   "@rtype: List of float");
PyObject* GeoMod_ReadHDF5Field(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_HDF5
  std::cerr << "Compiled without HDF5 support, no data read" << std::endl;
  Py_INCREF(Py_None);
  return Py_None;
#else
  static const char* keyWords[] = {"filename", "fieldname", "group", "ints", NULL };
  char* fname = 0;
  char* fldname = 0;
  char* group = 0;
  bool ints=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"sss|b",
                                   (char**)keyWords,&fname,&fldname,&group,&ints))
    return NULL;

  if (!fname || !fldname || !group)
    return NULL;

  return DoReadHDF5Field(fname, fldname, group, ints);
#endif
}

PyDoc_STRVAR(read_hdf5geometry__doc__,"Read a geometry from a HDF5 file\n"
                                      "@param filename: The file to read\n"
                                      "@type  filename: string\n"
                                      "@param fieldname: Geometry name\n"
                                      "@type fieldname: string \n"
                                      "@param patch: patch number\n"
                                      "@type patch: integer (>= 1)\n"
                                      "@param level: time level\n"
                                      "@type level: integer\n"
                                      "@return: The requested geometry\n"
                                      "@rtype: Curve, Surface or Volume");
PyObject* GeoMod_ReadHDF5Geometry(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_HDF5
  std::cerr << "Compiled without HDF5 support, no data read" << std::endl;
  Py_INCREF(Py_None);
  return Py_None;
#else
  static const char* keyWords[] = {"filename", "geometry", "patch", "level", NULL };
  char* fname = 0;
  char* geoname = 0;
  int patch=1;
  int level=0;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ssii",
                                   (char**)keyWords,&fname,&geoname,&patch,&level))
    return NULL;

  if (!fname || !geoname || patch < 1)
    return NULL;

  return DoReadHDF5Geometry(fname, geoname, patch, level);
#endif
}

PyDoc_STRVAR(write_hdf5field__doc__,"Write a field to a HDF5 file\n"
                                    "@param filename: The file to read\n"
                                    "@type  filename: String\n"
                                    "@param fieldname: Field name\n"
                                    "@type fieldname: String \n"
                                    "@param patch: patch number\n"
                                    "@type patch: Integer (>= 1)\n"
                                    "@param level: time level\n"
                                    "@type level: Integer\n"
                                    "@param coefs: Field coefficients\n"
                                    "@type coefs: List of float\n"
                                    "@param truncate: Truncate hdf5 file?\n"
                                    "@type truncate: Boolean\n"
                                    "@return: None\n");
PyObject* GeoMod_WriteHDF5Field(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_HDF5
  std::cerr << "Compiled without HDF5 support, no data written" << std::endl;
#else
  static const char* keyWords[] = {"filename", "fieldname", "patch", "level", "coefs", "truncate", NULL };
  char* fname = 0;
  char* fldname = 0;
  int patch=1;
  int level=0;
  PyObject* coefso;
  bool truncate=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ssiiO|b",
                                   (char**)keyWords,&fname,&fldname,&patch,
                                                    &level,&coefso,&truncate))
    return NULL;

  if (!coefso || !PyObject_TypeCheck(coefso, &PyList_Type) ||
      !fname || !fldname || patch < 1)
    return NULL;

  DoWriteHDF5Field(fname, fldname, patch, level, coefso, truncate);
#endif
  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(write_hdf5geometry__doc__,"Write a geometry to a HDF5 file\n"
                                       "@param filename: The file to read\n"
                                       "@type  filename: String\n"
                                       "@param geometry: Geometry name\n"
                                       "@type geometry: String \n"
                                       "@param patch: patch number\n"
                                       "@type patch: Integer (>= 1)\n"
                                       "@param level: time level\n"
                                       "@type level: Integer\n"
                                       "@param Object: Object\n"
                                       "@type Object: Curve, Surface or Volume\n"
                                       "@param truncate: Truncate hdf5 file?\n"
                                       "@type truncate: Boolean\n"
                                       "@return: None\n");
PyObject* GeoMod_WriteHDF5Geometry(PyObject* self, PyObject* args, PyObject* kwds)
{
#ifndef ENABLE_HDF5
  std::cerr << "Compiled without HDF5 support, no data written" << std::endl;
#else
  static const char* keyWords[] = {"filename", "geometry", "patch", "level", "object", "truncate", NULL };
  char* fname = 0;
  char* geoname = 0;
  int patch=1;
  int level=0;
  PyObject* objecto;
  bool truncate=false;
  if (!PyArg_ParseTupleAndKeywords(args,kwds,(char*)"ssiiO|b",
                                   (char**)keyWords,&fname,&geoname,&patch,
                                                    &level,&objecto,&truncate))
    return NULL;

  if (!objecto || !fname || !geoname || patch < 1)
    return NULL;

  DoWriteHDF5Geometry(fname, geoname, patch, level, objecto, truncate);
#endif
  Py_INCREF(Py_None);
  return Py_None;
}

PyDoc_STRVAR(final_output__doc__,"Write final entities to G2 file\n"
                                 "@param entities: The entities to write to file\n"
                                 "@type  entities: Curve, Surface, Volume, SurfaceModel, VolumeModel or a list of these\n"
                                 "@param  convert: Convert to a spline entity before saving\n"
                                 "@type   convert: Boolean\n"
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
     {(char*)"GetTolerance",          (PyCFunction)GeoMod_GetTolerance,      METH_VARARGS|METH_KEYWORDS, get_tolerance__doc__},
     {(char*)"SetTolerance",          (PyCFunction)GeoMod_SetTolerance,      METH_VARARGS|METH_KEYWORDS, set_tolerance__doc__},
     {(char*)"SetDimension",          (PyCFunction)GeoMod_SetDimension,      METH_VARARGS|METH_KEYWORDS, set_dim__doc__},
     {(char*)"SetFinalOutput",        (PyCFunction)GeoMod_SetFinalOutput,    METH_VARARGS|METH_KEYWORDS, set_final_output__doc__},
     {(char*)"GetProcessorCount",     (PyCFunction)GeoMod_GetProcessorCount, 0,                          get_processor_count__doc__},
     {(char*)"SetProcessorCount",     (PyCFunction)GeoMod_SetProcessorCount, METH_VARARGS|METH_KEYWORDS, set_processor_count__doc__},

     // I/O
     {(char*)"FinalOutput",           (PyCFunction)GeoMod_FinalOutput,       METH_VARARGS|METH_KEYWORDS, final_output__doc__},
     {(char*)"ReadIGES",              (PyCFunction)GeoMod_ReadIGES,          METH_VARARGS|METH_KEYWORDS, readiges__doc__},
     {(char*)"Read3DM",               (PyCFunction)GeoMod_Read3DM,           METH_VARARGS|METH_KEYWORDS, read3dm__doc__},
     {(char*)"ReadG2",                (PyCFunction)GeoMod_ReadG2,            METH_VARARGS|METH_KEYWORDS, readg2__doc__},
     {(char*)"ReadHDF5Field",         (PyCFunction)GeoMod_ReadHDF5Field,     METH_VARARGS|METH_KEYWORDS, read_hdf5field__doc__},
     {(char*)"ReadHDF5Geometry",      (PyCFunction)GeoMod_ReadHDF5Geometry,  METH_VARARGS|METH_KEYWORDS, read_hdf5geometry__doc__},
     {(char*)"SetDebugLevel",         (PyCFunction)GeoMod_SetDebugLevel,     METH_VARARGS|METH_KEYWORDS, set_debug_level__doc__},
     {(char*)"WriteG2",               (PyCFunction)GeoMod_WriteG2,           METH_VARARGS|METH_KEYWORDS, writeg2__doc__},
     {(char*)"WriteSTL",              (PyCFunction)GeoMod_WriteSTL,          METH_VARARGS|METH_KEYWORDS, writestl__doc__},
     {(char*)"Write3DM",              (PyCFunction)GeoMod_Write3DM,          METH_VARARGS|METH_KEYWORDS, write3dm__doc__},
     {(char*)"WriteHDF5Field",        (PyCFunction)GeoMod_WriteHDF5Field,    METH_VARARGS|METH_KEYWORDS, write_hdf5field__doc__},
     {(char*)"WriteHDF5Geometry",     (PyCFunction)GeoMod_WriteHDF5Geometry, METH_VARARGS|METH_KEYWORDS, write_hdf5geometry__doc__},
     
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
  init_VolumeModel_Type();
  init_SurfaceModelFactory_Module();
  init_Volume_Type();
  init_VolumeFactory_Module();
  init_Preprocess_Module();
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
  PyModule_AddObject(geoModule,(char*)"VolumeModel",(PyObject*)&VolumeModel_Type);
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

}
