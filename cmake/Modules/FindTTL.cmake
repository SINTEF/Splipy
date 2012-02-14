# - Tries to find the TTL library
#
# Written by: jan.b.thomassen@sintef.no
#

# 'GoTools_BUILD_ALL' will be defined in the top-level CMakeLists.txt
# file if we are building all of GoTools in one project.
IF(GoTools_BUILD_ALL)
  # Header files
  SET(TTL_INCLUDE_DIRS ${ttl_SOURCE_DIR}/include
    CACHE PATH "Path to TTL header files")
  # Library
  SET(TTL_LIBRARIES ttl
    CACHE FILE "TTL library")
ENDIF(GoTools_BUILD_ALL)

# Find header files
FIND_PATH(TTL_INCLUDE_DIRS 
  "ttl/ttl.h"
  "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  "C:/Program Files (x86)/ttl/include"
  "$ENV{PROGRAMFILES}/SINTEF/ttl/include/newmat"
)

# Find library
FIND_LIBRARY(TTL_LIBRARIES
  NAMES ttl
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  "C:/Program Files (x86)/ttl/lib"
  "$ENV{PROGRAMFILES}/SINTEF/ttl/lib"
)

# Check that we have found everything
SET(TTL_FOUND FALSE)
IF(TTL_INCLUDE_DIRS AND TTL_LIBRARIES)
  SET(TTL_FOUND TRUE)
ENDIF(TTL_INCLUDE_DIRS AND TTL_LIBRARIES)
