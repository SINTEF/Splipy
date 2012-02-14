# - Tries to find the GoTools Topology library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoTopology_INCLUDE_DIRS ${GoTopology_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools Topology header files")
  # Library
  SET(GoTopology_LIBRARIES GoTopology
    CACHE FILE "GoTools Topology library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoTopology_INCLUDE_DIRS 
  "GoTools/topology/tpTopologyTable.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/include"
)

# Find library
FIND_LIBRARY(GoTopology_LIBRARIES
  NAMES GoTopology
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/lib"
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(GoTopology_FOUND FALSE)
IF(GoTopology_INCLUDE_DIRS AND GoTopology_LIBRARIES)
  SET(GoTopology_FOUND TRUE)
ENDIF(GoTopology_INCLUDE_DIRS AND GoTopology_LIBRARIES)
