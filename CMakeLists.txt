project(GoTools-python-bindings)

set(GEOMODELLER_VERSION_MAJOR 0)
set(GEOMODELLER_VERSION_MINOR 8)
set(GEOMODELLER_VERSION_PATCH 0)
set(GEOMODELLER_VERSION ${GEOMODELLER_VERSION_MAJOR}.${GEOMODELLER_VERSION_MINOR}.${GEOMODELLER_VERSION_PATCH})

cmake_minimum_required(VERSION 2.6)

enable_language(CXX)

option(USE_OPENNURBS "Use OpenNURBS?" ON)
option(USE_HDF5 "Use HDF5?" ON)
option(USE_GPM "Use GPM?" ON)
option(USE_VTF "Use VTF?" ON)
option(BUILD_PYTHON_MODULE "Build standalone python module?" ON)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_SOURCE_DIR}/cmake/Modules)

find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(GoTools REQUIRED)
find_package(GoTrivariate REQUIRED)
find_package(GoTrivariateModel REQUIRED)
find_package(GoCompositeModel REQUIRED)
find_package(Parametrization REQUIRED)
find_package(TTL REQUIRED)
find_package(Newmat REQUIRED)
find_package(GoIgeslib REQUIRED)
find_package(SISL REQUIRED)
find_package(GoTopology REQUIRED)
find_package(GoImplicitization REQUIRED)
find_package(GoIntersections REQUIRED)
find_package(LATEX)
find_package(EpyDoc)
find_package(LRSpline)

if(USE_OPENNURBS)
  find_package(OpenNURBS)
endif()
if(USE_HDF5)
  find_package(HDF5)
endif()

if(USE_GPM)
  find_package(GPM)
endif()

if(USE_VTF)
  find_package(VTFWriter)
endif()

# Generate header with version info
configure_file(${PROJECT_SOURCE_DIR}/include/geomodversion.h.in
               geomodversion.h @ONLY)

# Mimimum GoTools version
if(GoTools_VERSION_MAJOR LESS 4 OR NOT GoTools_VERSION_MAJOR)
  message(FATAL_ERROR "GoTools >= 4.0.0 required. bailing")
endif()

# Order of GoTools libraries is not random!
set(DEPLIBS ${PYTHON_LIBRARIES}
            ${GoTrivariateModel_LIBRARIES}
            ${GoCompositeModel_LIBRARIES}
            ${Parametrization_LIBRARIES}
            ${TTL_LIBRARIES}
            ${GoIgeslib_LIBRARIES}
            ${SISL_LIBRARIES}
            ${GoTopology_LIBRARIES}
            ${GoIntersections_LIBRARIES}
            ${GoImplicitization_LIBRARIES}
            ${Newmat_LIBRARIES}
            ${GoTools_LIBRARIES}
            ${GoTrivariate_LIBRARIES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GoTools_CXX_FLAGS} -fPIC")
SET(INCLUDE_DIRS ${GoTools_INCLUDE_DIRS}
                 ${PYTHON_INCLUDE_DIRS} 
                 ${PROJECT_SOURCE_DIR}/include)

if(OpenNURBS_LIBRARIES AND OpenNURBS_INCLUDE_DIRS)
  set(INCLUDE_DIRS ${INCLUDE_DIRS}
                   ${OpenNURBS_INCLUDE_DIRS})
  set(DEPLIBS ${DEPLIBS}
              ${OpenNURBS_LIBRARIES})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_OPENNURBS=1")
endif()

if(HDF5_LIBRARIES AND HDF5_INCLUDE_DIRS)
  set(DEPLIBS ${DEPLIBS}
              ${HDF5_LIBRARIES})
  set(INCLUDE_DIRS ${INCLUDE_DIRS}
                   ${HDF5_INCLUDE_DIRS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_HDF5=1")
endif()

if(GPM_LIBRARIES AND GPM_INCLUDES)
  if(GPM_VERSION_MINOR GREATER 3 OR GPM_VERSION_PATCH GREATER 0)
    list(APPEND DEPLIBS ${GPM_LIBRARIES})
    list(APPEND INCLUDE_DIRS ${GPM_INCLUDES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_GPM=1")
  endif()
endif()

if(LRSpline_FOUND)
  list(APPEND DEPLIBS ${LRSpline_LIBRARIES})
  list(APPEND INCLUDE_DIRS ${LRSpline_INCLUDE_DIRS})
endif()

if(VTFWRITER_LIBRARIES AND VTFWRITER_INCLUDES)
  list(APPEND DEPLIBS ${VTFWRITER_LIBRARIES})
  list(APPEND INCLUDE_DIRS ${VTFWRITER_INCLUDE_DIRS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_VTF=1")
endif()

include_directories(${INCLUDE_DIRS} ${CMAKE_BINARY_DIR})

# Additional compiler flags required by GoTools
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GoTools_CXX_FLAGS}")

file(GLOB_RECURSE SRCS src/*.C)

if(BUILD_PYTHON_MODULE)
  python_add_module(GoTools ${SRCS})
  target_link_libraries(GoTools ${DEPLIBS})
endif()

add_executable(geoModeler app/geoModeler.C ${SRCS})
target_link_libraries(geoModeler ${DEPLIBS})

execute_process(COMMAND python2 -c "import sys; print sys.version[:3]" OUTPUT_VARIABLE PYTHON_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

if(EPYDOC_EXECUTABLE)
  add_custom_target(doc PYTHONPATH=${CMAKE_BINARY_DIR}:${PROJECT_SOURCE_DIR} ${PYTHON_EXECUTABLE}
                        ${CMAKE_SOURCE_DIR}/doc/gendocumentation.py
                        DEPENDS GoTools
                        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                        COMMENT "Generating API documentation" VERBATIM)
endif()
if(PDFLATEX_COMPILER)
  add_custom_target(pdfdoc ${PDFLATEX_COMPILER} api.tex
                    DEPENDS doc
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/latex
                    COMMENT "Compiling documentation to PDF (doc/latex/api.pdf)" VERBATIM)
endif()

set(TAR_BALL GoTools-python-bindings-${GEOMODELLER_VERSION_MAJOR}.${GEOMODELLER_VERSION_MINOR}.${GEOMODELLER_VERSION_PATCH}.tar.bz2)
file(GLOB_RECURSE DIST_SRCS RELATIVE ${PROJECT_SOURCE_DIR} CMakeLists.txt src/*.C include/*.h include/*.in app/*.C doc/gendocumentation.py cmake/*.cmake)
string(REGEX REPLACE "include/geomodversion.h;" "" DIST_SRCS "${DIST_SRCS}")
add_custom_target(dist tar jcvf ${TAR_BALL} ${DIST_SRCS} > /dev/null
                  DEPENDS GoTools
                  COMMENT "Building tar ball \"${TAR_BALL}\"" VERBATIM)
set_directory_properties(${PROJECT_SOURCE_DIR} ADDITIONAL_MAKE_CLEAN_FILES ${TAR_BALL})
set_directory_properties(${PROJECT_SOURCE_DIR}/doc ADDITIONAL_MAKE_CLEAN_DIRECTORIES html)

# 'install' target
# Make install force doc build
if(EPYDOC_EXECUTABLE)
  install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_BUILD_TOOL} doc WORKING_DIRECTORY \"${CMAKE_CURRENT_BINARY_DIR}\")" COMPONENT doc)
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION share/doc/python-GoTools
          COMPONENT doc)
endif()
if(BUILD_PYTHON_MODULE)
  if(EXISTS "/etc/debian_version")
    install(TARGETS GoTools DESTINATION lib/python${PYTHON_VERSION}/dist-packages/ COMPONENT bin)
  else()
    install(TARGETS GoTools DESTINATION lib/python${PYTHON_VERSION}/site-packages/ COMPONENT bin)
  endif()
endif()

install(TARGETS geoModeler DESTINATION bin COMPONENT bin)

# Add correct path to the geomod binary for convenience
set(GEOMOD_BINARY ${CMAKE_BINARY_DIR}/geoModeler)
set(PROJECT_PATH ${PROJECT_SOURCE_DIR})
configure_file(${PROJECT_SOURCE_DIR}/tests/makereg.in
               makereg)
configure_file(${PROJECT_SOURCE_DIR}/tests/dotest.in
               dotest)

file(GLOB GEOUTILS_SRCS RELATIVE ${PROJECT_SOURCE_DIR} GeoUtils/*.py)
install(FILES ${GEOUTILS_SRCS} DESTINATION lib/python${PYTHON_VERSION}/dist-packages/GeoUtils
        COMPONENT bin)

enable_testing()
file(GLOB_RECURSE GEOMOD_TESTFILES RELATIVE ${PROJECT_SOURCE_DIR}/tests ${PROJECT_SOURCE_DIR}/tests/*.reg)
foreach(TESTFILE ${GEOMOD_TESTFILES})
  string(REGEX REPLACE ".reg" ".py" TESTSCRIPT ${TESTFILE})
  add_test(${TESTFILE} ${CMAKE_BINARY_DIR}/dotest ${CMAKE_BINARY_DIR}/geoModeler ${PROJECT_SOURCE_DIR}/scripts/${TESTSCRIPT} ${PROJECT_SOURCE_DIR}/tests/${TESTFILE})
endforeach()

add_custom_target(check-commits COMMAND ${CMAKE_COMMAND}
                                        -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
                                        -DCMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}
                                        -P ${PROJECT_SOURCE_DIR}/cmake/Scripts/CheckCommits.cmake)


# PIP/Distutils Package
configure_file(${CMAKE_SOURCE_DIR}/setup.py.in ${CMAKE_BINARY_DIR}/setup.py)
file(COPY ${CMAKE_SOURCE_DIR}/GeoUtils DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/README DESTINATION ${CMAKE_BINARY_DIR})
add_custom_target(pip DEPENDS GoTools COMMAND python2 setup.py sdist --dist-dir ${CMAKE_BINARY_DIR})
