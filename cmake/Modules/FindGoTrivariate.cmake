IF(GoTrivariate_INCLUDE_DIRS AND GoTrivariate_LIBRARIES)
  SET(GoTrivariate_FIND_QUIETLY TRUE)
ENDIF(GoTrivariate_INCLUDE_DIRS AND GoTrivariate_LIBRARIES)

FIND_PATH(GoTrivariate_INCLUDE_DIRS
  NAMES GoTools/trivariate/SplineVolume.h
  PATHS "$ENV{HOME}/include"
  "$ENV{HOME}/install/include"
  /sima/libs/GoTools/include
)

FIND_LIBRARY(GoTrivariate_LIBRARIES
  NAMES GoTrivariate
  PATHS "$ENV{HOME}/lib"
  "$ENV{HOME}/install/lib"
  /sima/libs/GoTools/lib
  PATH_SUFFIXES GoTools
)

INCLUDE(FindPackageHandleStandardArgs)
IF(GoTrivariate_LIBRARIES)
  find_package_handle_standard_args(GoTrivariate DEFAULT_MSG
                                    GoTrivariate_INCLUDE_DIRS GoTrivariate_LIBRARIES)
ENDIF(GoTrivariate_LIBRARIES)
