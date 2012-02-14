# - Tries to find the GoTools Parametrization library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(Parametrization_INCLUDE_DIRS ${parametrization_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Parametrization header files")
  # Library
  SET(Parametrization_LIBRARIES parametrization
    CACHE FILE "GoTools Parametrization library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(Parametrization_INCLUDE_DIRS 
  "GoTools/parametrization/PrOrganizedPoints.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/include"
)

# Find library
FIND_LIBRARY(Parametrization_LIBRARIES
  NAMES parametrization
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/lib"
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(Parametrization_FOUND FALSE)
IF(Parametrization_INCLUDE_DIRS AND Parametrization_LIBRARIES)
  SET(Parametrization_FOUND TRUE)
ENDIF(Parametrization_INCLUDE_DIRS AND Parametrization_LIBRARIES)
