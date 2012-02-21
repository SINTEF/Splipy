#pragma once

#include <Python.h>
#include "GoTools/utils/config.h"
#include "GoTools/utils/Point.h"

#define GEOMODELLER_VERSION_MAJOR 0
#define GEOMODELLER_VERSION_MINOR 1
#define GEOMODELLER_VERSION_PATCH 0

class GeoModellerState {
  public:
    GeoModellerState()
    { 
      dim = 3;
      // Negative values means it's overridden from command line
      gapTolerance = 1.e-4;
      approxTolerance = 1.e-3;
      neighbourTolerance = 1.e-2;
      kinkTolerance = 1.e-2;
      bendTolerance = 1.e-1;
      debugLevel = 1;
    }

    int dim; //!< The dimension of the geometries we will create
    int debugLevel; //!< The current debug level
    double gapTolerance; //!< The tolerance of water-tight models
    double approxTolerance; //!< The tolerance of some kind of approximation :)
    double neighbourTolerance; //!< How much crap we can stomach from our neighbour
    double kinkTolerance; //!< Kink tolerance
    double bendTolerance; //!< Bend tolerance
    std::string finalOutput; //!< The final output file
};

extern GeoModellerState modState;

extern "C"
{
  void initGoTools();
}

// helper functions
Go::Point someNormal(const Go::Point& vec);

  template<class PythonClass, class GoType>
PyObject* ReadG2(std::ifstream& str, PyTypeObject& type)
{
  PythonClass* result = (PythonClass*)type.tp_alloc(&type,0);
  result->data.reset(new GoType());
  result->data->read(str);

  return (PyObject*)result;
}
