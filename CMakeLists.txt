PROJECT(geomod)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

ENABLE_LANGUAGE(CXX)

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

INSTALL(TARGETS GoTools DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/dist-packages/)
INSTALL(TARGETS DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION}/dist-packages/)
FILE(GLOB html-docs ${PROJECT_SOURCE_DIR}/doc/html/*)
INSTALL(FILES ${html-docs} DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/python-GoTools/html)
