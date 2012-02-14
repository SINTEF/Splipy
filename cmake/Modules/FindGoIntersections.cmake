# - Tries to find the GoTools Intersections library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoIntersections_INCLUDE_DIRS ${GoIntersections_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Intersections header files")
  # Library
  SET(GoIntersections_LIBRARIES GoIntersections
    CACHE FILE "GoTools Intersections library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoIntersections_INCLUDE_DIRS 
  "GoTools/intersections/SfSfIntersector.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/include"
)

# Find library
FIND_LIBRARY(GoIntersections_LIBRARIES
  NAMES GoIntersections
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/lib"
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(GoIntersections_FOUND FALSE)
IF(GoIntersections_INCLUDE_DIRS AND GoIntersections_LIBRARIES)
  SET(GoIntersections_FOUND TRUE)
ENDIF(GoIntersections_INCLUDE_DIRS AND GoIntersections_LIBRARIES)
