# - Tries to find the GoTools CompositeModel library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(GoCompositeModel_INCLUDE_DIRS ${GoCompositeModel_SOURCE_DIR}/include
    CACHE PATH "Path to GoTools CompositeModel header files")
  # Library
  SET(GoCompositeModel_LIBRARIES GoCompositeModel
    CACHE FILE "GoTools CompositeModel library")
ENDIF(GoTools_BUILD_ALL)


# Find header files
FIND_PATH(GoCompositeModel_INCLUDE_DIRS 
  "GoTools/compositemodel/CompositeModel.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/GoTools/include"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/include"
)

# Find library
FIND_LIBRARY(GoCompositeModel_LIBRARIES
  NAMES GoCompositeModel
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/GoTools/lib"
  "$ENV{PROGRAMFILES}/SINTEF/GoTools/lib"
  PATH_SUFFIXES GoTools
)

# Check that we have found everything
SET(GoCompositeModel_FOUND FALSE)
IF(GoCompositeModel_INCLUDE_DIRS AND GoCompositeModel_LIBRARIES)
  SET(GoCompositeModel_FOUND TRUE)
ENDIF(GoCompositeModel_INCLUDE_DIRS AND GoCompositeModel_LIBRARIES)
