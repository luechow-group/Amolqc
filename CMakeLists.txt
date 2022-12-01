cmake_minimum_required(VERSION 3.10)

### options
option(CONFIG_FILE "ON - reading from cmake/config.cmake,   OFF - read no config file"               ON)
option(WARNINGS    "ON - display compiler warnings,         OFF - do not display compiler warnings"  OFF)
option(PFUNIT      "ON - using pFUnit for unit tests,       OFF - skip building unit tests"          ON)
option(OPENMP      "ON - enables adding of OpenMP targets,  OFF - disables adding of OpenMP targets" OFF)
option(MPI         "ON - compiling with MPI (-DMPI),        OFF - compiling the serial version"      OFF)
option(MKL         "ON - using LAPACK from MKL,             OFF - using compiled LAPACK"             OFF)
option(LINEINFO    "ON - compiling with line info (eg. -g), OFF - compiling without line info"       OFF)
option(CHKNANUP    "ON - check invalid updates in muldiDet, OFF - do not check for invalid updates"  OFF)
# ADDITIONAL_FLAGS sets additional compile flags                                          Default:
# WFMAXLINES sets maximum number of lines for the wf file                                 Default: 50000
# WFMAXCOLS sets maximum line length (# cols) for the wf file                             Default:   160
# MAXDETSPERCSF sets maximum number of dets per csf                                       Default:   120

# for parallel unit tests
set(NPROC 2)

message("")
message("---  STARTED AMOLQC BUILD PREPARATION  ---")

# apply parts of config.cmake before setting language
include(cmake/rules/preRules.cmake)

# define project name and used language
project(amolqc Fortran)

# apply rest of config.cmake
include(cmake/rules/postRules.cmake)

# set compile flags
include(cmake/utils/setFlags.cmake)

# allow for using add_test() and run tests with 'ctest'
enable_testing()

# add pFUnit support
include(cmake/utils/pFUnit.cmake)

# set path for executable amolqc
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin)

# add custom cmake functions
include(cmake/utils/fortranFunctions.cmake)

# add testsuite to tests
include(cmake/testsuite/testsuite.cmake)

# add subdirectories
add_subdirectory(utils)
add_subdirectory(src)

message("---  FINISHED AMOLQC BUILD PREPARATION ---")
message("")
