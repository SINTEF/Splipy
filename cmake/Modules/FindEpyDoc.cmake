include(FindPackageHandleStandardArgs)

if(NOT EPYDOC_EXECUTABLE)
  find_program(EPYDOC_EXECUTABLE NAMES epydoc)
  message(STATUS "Looking for Epydoc executable ")
else()
  if(NOT EXISTS ${EPYDOC_EXECUTABLE})
    set(EPYDOC_EXECUTABLE no)
  endif()
endif()

if(EPYDOC_EXECUTABLE)
  set(env_lc_all $ENV{LC_ALL})
  set(ENV{LC_ALL} "C")
  execute_process(COMMAND ${EPYDOC_EXECUTABLE} --version
    OUTPUT_VARIABLE stdout)
  string(REGEX REPLACE ".* ([0-9]+\\.[0-9]+(\\.[0-9]+)?).*" "\\1"
    EPYDOC_VERSION "${stdout}")
  set(ENV{LC_ALL} ${env_lc_all})
endif()

find_package_handle_standard_args(Epydoc
  REQUIRED_VARS EPYDOC_EXECUTABLE
  VERSION_VAR EPYDOC_VERSION
  ) 
