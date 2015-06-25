#include "undef.h"
#include <Python.h>

#include "geomodeller.h"

#ifdef ENABLE_OPENNURBS
#include <opennurbs.h>
#endif

using namespace GeoModeller;

void syntax()
{
  std::cout << "Syntax: geoModeler <script> [finaloutput=] [debuglevel=] [tospline=]" << std::endl;
}

std::string getexepath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

void addPath()
{
  std::string path = getexepath();
  path.erase(path.rfind('/'));
  path += ':';
  // we want to use sys.path so it includes site-packages
  // if this fails, default to using Py_GetPath
  PyObject *sysMod(PyImport_ImportModule((char*)"sys")); // must call Py_DECREF 
  PyObject *sysModDict(PyModule_GetDict(sysMod)); // borrowed ref, no need to de
  PyObject *pathObj(PyDict_GetItemString(sysModDict, "path")); // borrowed ref, 

  if( pathObj && PyList_Check(pathObj) )
  {
    for( int i = 0; i < PyList_Size(pathObj); i++ )
    {
      PyObject *e = PyList_GetItem(pathObj, i); // borrowed ref, no need to dele
      if( e && PyString_Check(e) )
      {
          path += PyString_AsString(e); // returns internal data, don't delete o
          path += ':';
      }
    }
  }
  else
  {
    path += Py_GetPath();
  }
  Py_DECREF(sysMod); // release ref to sysMod

  PySys_SetPath((char *)path.c_str());
}

int main(int argc, char** argv)
{
  const char* file=NULL;
  std::vector<char*> argv_;
  argv_.reserve(argc);
  argv_.push_back(NULL);
  for (int i=1;i<argc;++i) {
    if (!strncasecmp(argv[i],"finaloutput=",12))
      modState.finalOutput = argv[i]+12;
    else if (!strncasecmp(argv[i],"debuglevel=",11))
      modState.debugLevel = -atoi(argv[i]+11);
    else if (!strncasecmp(argv[i],"tospline=",9))
      modState.convertSpline = 1+atoi(argv[i]+9);
    else if (!file)
      file = argv[i];
    else
      argv_.push_back(argv[i]);
  }
  argv_[0] = (char*)file;
  if (!file) {
    syntax();
    exit(1);
  }

  FILE* f = fopen(file,"rb");
  if (!f) {
    std::cerr << "Error opening script file" << std::endl;
    exit(2);
  }
  Py_Initialize();
#ifndef __APPLE__
  addPath();
#endif
#ifdef ENABLE_OPENNURBS
  ON::Begin();
#endif
  initGoTools();
  PySys_SetArgv(argv_.size(), &argv_[0]);
  int result=PyRun_SimpleFile(f,file);
  fclose(f);
#ifndef __APPLE__
  Py_Finalize();
#endif
#ifdef ENABLE_OPENNURBS
  ON::End();
#endif
  return result;
}
