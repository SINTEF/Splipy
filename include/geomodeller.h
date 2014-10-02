#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/utils/Point.h"

#include <thread>
#include <sstream>

#include "geomodversion.h"

class GeoModellerState {
  public:
    GeoModellerState()
    { 
      std::stringstream str;
      str << GEOMODELLER_VERSION_MAJOR << "." << GEOMODELLER_VERSION_MINOR << "." << GEOMODELLER_VERSION_PATCH;
      version = str.str();
      author = "Arne Morten Kvarving, Kjetil A. Johannessen, Timo van Opstahl, Knut Nordanger, Eivind Fonn";
      date = __DATE__;
      credits = "The GoTools authors";

      dim = 3;
      convertSpline = 0;
      // Negative values means it's overridden from command line
      gapTolerance = 1.e-4;
      approxTolerance = 1.e-3;
      neighbourTolerance = 1.e-2;
      kinkTolerance = 1.e-2;
      bendTolerance = 1.e-1;
      refineTolerance = 1.e-3;
      debugLevel = 1;
      procCount = std::thread::hardware_concurrency();
    }

    void addInfo(PyObject* module)
    {
      PyModule_AddStringConstant(module,(char*)"__author__",  author.c_str());
      PyModule_AddStringConstant(module,(char*)"__credits__", credits.c_str());
      PyModule_AddStringConstant(module,(char*)"__date__",    date.c_str());
      PyModule_AddStringConstant(module,(char*)"__package__", (char*)"GoTools");
      PyModule_AddStringConstant(module,(char*)"__version__", version.c_str());
    }

    int dim; //!< The dimension of the geometries we will create
    int convertSpline; //!< Make output convert to splines by default
    int debugLevel; //!< The current debug level
    int procCount; //!< Processor count
    double gapTolerance; //!< The tolerance of water-tight models
    double approxTolerance; //!< The tolerance of some kind of approximation :)
    double neighbourTolerance; //!< How much crap we can stomach from our neighbour
    double kinkTolerance; //!< Kink tolerance
    double bendTolerance; //!< Bend tolerance
    double refineTolerance; //!< Tolerance for skipping refinement-inserted knots as they already exist
    std::string finalOutput; //!< The final output file

    //!< Python binding metadata
    std::string author; //!!< The binding authors
    std::string credits; //!!< The binding credits
    std::string date; //!!< The compilation date
    std::string version; //!< The binding version
};

extern GeoModellerState modState;

extern "C"
{
  void initGoTools();
}

// helper functions
Go::Point someNormal(const Go::Point& vec);

  template<class PythonClass, class GoType>
PyObject* ReadG2(std::istream& str, PyTypeObject& type)
{
  PythonClass* result = (PythonClass*)type.tp_alloc(&type,0);
  result->data.reset(new GoType());
  result->data->read(str);

  return (PyObject*)result;
}
