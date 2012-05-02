#include <Python.h>

#include "geomodeller.h"

#ifdef ENABLE_OPENNURBS
#include <opennurbs.h>
#endif

void syntax()
{
  std::cout << "Syntax: geoModeler <script> [finaloutput=] [debuglevel=] [tospline=]" << std::endl;
}

int main(int argc, char** argv)
{
  const char* file=NULL;
  for (int i=1;i<argc;++i) {
    if (!strncasecmp(argv[i],"finaloutput=",12))
      modState.finalOutput = argv[i]+12;
    else if (!strncasecmp(argv[i],"debuglevel=",11))
      modState.debugLevel = -atoi(argv[i]+11);
    else if (!strncasecmp(argv[i],"tospline=",9))
      modState.convertSpline = 1+atoi(argv[i]+9);
    else
      file = argv[i];
  }
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
#ifdef ENABLE_OPENNURBS
  ON::Begin();
#endif
  initGoTools();
  int result=PyRun_SimpleFile(f,file);
  fclose(f);
  Py_Finalize();
#ifdef ENABLE_OPENNURBS
  ON::End();
#endif
  return result;
}
