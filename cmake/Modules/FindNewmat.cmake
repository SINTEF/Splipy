# - Tries to find the Newmat library
#
# Written by: jan.b.thomassen@sintef.no
#


# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(Newmat_INCLUDE_DIRS ${newmat_SOURCE_DIR}/include
    CACHE PATH "Path to Newmat header files")
  # Library
  SET(Newmat_LIBRARIES newmat
    CACHE FILE "Newamt library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(Newmat_INCLUDE_DIRS "newmat.h"
  "/usr/include/newmat"
  "/usr/local/include/newmat"
  "$ENV{HOME}/include/newmat"
  "$ENV{HOME}/install/include/newmat"
  "$ENV{PROGRAMFILES}/newmat/include/newmat"
  "$ENV{PROGRAMFILES}/SINTEF/newmat/include/newmat"
)

# Find library
FIND_LIBRARY(Newmat_LIBRARIES
  NAMES newmat
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "$ENV{PROGRAMFILES}/newmat/lib"
  "$ENV{PROGRAMFILES}/SINTEF/newmat/lib"
)

# Check that we have found everything
SET(Newmat_FOUND FALSE)
IF(Newmat_INCLUDE_DIRS AND Newmat_LIBRARIES)
  SET(Newmat_FOUND TRUE)
ENDIF(Newmat_INCLUDE_DIRS AND Newmat_LIBRARIES)
