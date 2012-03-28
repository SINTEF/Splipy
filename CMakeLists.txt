PROJECT(GoTools-python-bindings)

SET(GEOMODELLER_VERSION_MAJOR 0)
SET(GEOMODELLER_VERSION_MINOR 1)
SET(GEOMODELLER_VERSION_PATCH 0)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)

# Generate header with version info
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/include/geomodversion.h.in
               ${PROJECT_SOURCE_DIR}/include/geomodversion.h)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_SOURCE_DIR}/cmake/Modules)

FIND_PACKAGE(PythonLibs REQUIRED)
FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(GoTools REQUIRED)
FIND_PACKAGE(GoTrivariate REQUIRED)
FIND_PACKAGE(GoCompositeModel REQUIRED)
FIND_PACKAGE(SISL REQUIRED)
FIND_PACKAGE(Parametrization REQUIRED)
FIND_PACKAGE(TTL REQUIRED)
FIND_PACKAGE(Newmat REQUIRED)
FIND_PACKAGE(GoIgeslib REQUIRED)
FIND_PACKAGE(GoTopology REQUIRED)
FIND_PACKAGE(GoImplicitization REQUIRED)
FIND_PACKAGE(GoIntersections REQUIRED)

# Mimimum GoTools version
IF(GoTools_VERSION_MAJOR LESS 4 OR NOT GoTools_VERSION_MAJOR)
  MESSAGE(FATAL_ERROR "GoTools >= 4.0.0 required. bailing")
ENDIF(GoTools_VERSION_MAJOR LESS 4 OR NOT GoTools_VERSION_MAJOR)

# Order of GoTools libraries is not random!
SET(DEPLIBS ${PYTHON_LIBRARIES}
            ${GoCompositeModel_LIBRARIES}
            ${SISL_LIBRARIES}
            ${Parametrization_LIBRARIES}
            ${TTL_LIBRARIES}
            ${GoIgeslib_LIBRARIES}
            ${GoTopology_LIBRARIES}
            ${GoIntersections_LIBRARIES}
            ${GoImplicitization_LIBRARIES}
            ${Newmat_LIBRARIES}
            ${GoTools_LIBRARIES}
            ${GoTrivariate_LIBRARIES})


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GoTools_CXX_FLAGS} -fPIC")
SET(INCLUDE_DIRS ${GoTools_INCLUDE_DIRS}
                 ${PYTHON_INCLUDE_DIRS} 
                 ${PROJECT_SOURCE_DIR}/include)

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

# Additional compiler flags required by GoTools
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GoTools_CXX_FLAGS}")

FILE(GLOB_RECURSE SRCS src/*.C)

PYTHON_ADD_MODULE(GoTools ${SRCS})
TARGET_LINK_LIBRARIES(GoTools ${DEPLIBS})

ADD_EXECUTABLE(geoModeler ${SRCS} app/geoModeler.C)
TARGET_LINK_LIBRARIES(geoModeler ${DEPLIBS})

EXECUTE_PROCESS(COMMAND python -c "import sys; print sys.version[:3]" OUTPUT_VARIABLE PYTHON_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

ADD_CUSTOM_TARGET(doc PYTHONPATH=${CMAKE_BINARY_DIR} ${PYTHON_EXECUTABLE}
                      ${CMAKE_SOURCE_DIR}/doc/gendocumentation.py
                      DEPENDS GoTools
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Generating API documentation" VERBATIM)

ADD_CUSTOM_TARGET(pdfdoc pdflatex
                  api.tex
                  DEPENDS doc
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/latex
                  COMMENT "Compiling documentation to PDF (doc/latex/api.pdf)" VERBATIM)

SET(TAR_BALL GoTools-python-bindings-${GEOMODELLER_VERSION_MAJOR}.${GEOMODELLER_VERSION_MINOR}.${GEOMODELLER_VERSION_PATCH}.tar.bz2)
FILE(GLOB_RECURSE DIST_SRCS RELATIVE ${PROJECT_SOURCE_DIR} CMakeLists.txt src/*.C include/*.h include/*.in app/*.C doc/gendocumentation.py cmake/*.cmake)
STRING(REGEX REPLACE "include/geomodversion.h;" "" DIST_SRCS "${DIST_SRCS}")
ADD_CUSTOM_TARGET(dist tar jcvf ${TAR_BALL} ${DIST_SRCS} > /dev/null
                  DEPENDS GoTools
                  COMMENT "Building tar ball \"${TAR_BALL}\"" VERBATIM)
SET_DIRECTORY_PROPERTIES(${PROJECT_SOURCE_DIR} ADDITIONAL_MAKE_CLEAN_FILES ${TAR_BALL})
SET_DIRECTORY_PROPERTIES(${PROJECT_SOURCE_DIR}/doc ADDITIONAL_MAKE_CLEAN_DIRECTORIES html)


INSTALL(TARGETS GoTools DESTINATION lib/python${PYTHON_VERSION}/dist-packages/)
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/python-GoTools/html)
INSTALL(TARGETS geoModeler DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

ENABLE_TESTING()
FILE(GLOB_RECURSE GEOMOD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/tests ${PROJECT_SOURCE_DIR}/tests/*.reg)
FOREACH(TESTFILE ${GEOMOD_TESTFILES})
  STRING(REGEX REPLACE ".reg" ".py" TESTSCRIPT ${TESTFILE})
  ADD_TEST(${TESTFILE} ${PROJECT_SOURCE_DIR}/tests/dotest ${CMAKE_BINARY_DIR}/geoModeler ${PROJECT_SOURCE_DIR}/scripts/${TESTSCRIPT} ${PROJECT_SOURCE_DIR}/tests/${TESTFILE})
ENDFOREACH()
