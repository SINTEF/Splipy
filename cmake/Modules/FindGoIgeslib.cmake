# - Tries to find the GoTools Igeslib library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoIgeslib_INCLUDE_DIRS ${GoIgeslib_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Igeslib header files")
  # Library
  SET(GoIgeslib_LIBRARIES GoIgeslib
    CACHE FILE "GoTools Igeslib library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoIgeslib_INCLUDE_DIRS "GoTools/igeslib/IGESconverter.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/include"
)

# Find library
FIND_LIBRARY(GoIgeslib_LIBRARIES
  NAMES GoIgeslib
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/lib"
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(GoIgeslib_FOUND FALSE)
IF(GoIgeslib_INCLUDE_DIRS AND GoIgeslib_LIBRARIES)
  SET(GoIgeslib_FOUND TRUE)
ENDIF(GoIgeslib_INCLUDE_DIRS AND GoIgeslib_LIBRARIES)
