# - Tries to find the SISL library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(SISL_INCLUDE_DIRS ${sisl_SOURCE_DIR}/include
    CACHE PATH "Path to SISL header files")
  # Library
  SET(SISL_LIBRARIES sisl
    CACHE FILE "SISL library")
ENDIF(GoTools_BUILD_ALL)

# Find header files
FIND_PATH(SISL_INCLUDE_DIRS "sisl.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/sisl/include"
  "$ENV{PROGRAMFILES}/SINTEF/sisl/include"
)

# Find library
FIND_LIBRARY(SISL_LIBRARIES
  NAMES sisl 
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/sisl/lib"
  "$ENV{PROGRAMFILES}/SINTEF/sisl/lib"
)

# Check that we have found everything
SET(SISL_FOUND FALSE)
IF(SISL_INCLUDE_DIRS AND SISL_LIBRARIES)
  SET(SISL_FOUND TRUE)
ENDIF(SISL_INCLUDE_DIRS AND SISL_LIBRARIES)
