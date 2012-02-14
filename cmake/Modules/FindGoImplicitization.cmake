# - Tries to find the GoTools Implicitization library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoImplicitization_INCLUDE_DIRS ${GoImplicitization_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Implicitization header files")
  # Library
  SET(GoImplicitization_LIBRARIES GoImplicitization
    CACHE FILE "GoTools Implicitization library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoImplicitization_INCLUDE_DIRS 
  "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/include"
)

# Find library
FIND_LIBRARY(GoImplicitization_LIBRARIES
  NAMES GoImplicitization
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/lib"
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(GoImplicitization_FOUND FALSE)
IF(GoImplicitization_INCLUDE_DIRS AND GoImplicitization_LIBRARIES)
  SET(GoImplicitization_FOUND TRUE)
ENDIF(GoImplicitization_INCLUDE_DIRS AND GoImplicitization_LIBRARIES)
