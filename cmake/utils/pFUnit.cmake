### for the use of pFUnit

if (${PFUNIT})
  enable_testing()
  include(ExternalProject)

  # python is needed for the preprocessing
  find_program(PYTHON python)

  if (DEFINED ENV{PFUNIT})
    set(PFUNIT_BUILD_PATH $ENV{PFUNIT})
    message("found PFUNIT in $ENV{PFUNIT}.")
    add_custom_target(pFUnit)
  else()
    set(PFUNIT_BUILD_PATH ${PROJECT_BINARY_DIR}/pFUnit)

    if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      # Avoiding global name length warnings
      set(PFUNIT_COMPILE_FLAGS "-diag-disable 5462")
    endif()

    # this way, pFUnit is configured during the cmake build step
    ExternalProject_Add(pFUnit
      PREFIX pFUnit
      SOURCE_DIR ${PROJECT_SOURCE_DIR}/pFUnit
      BINARY_DIR ${PFUNIT_BUILD_PATH}
      TMP_DIR ${PFUNIT_BUILD_PATH}/tmp
      STAMP_DIR ${PFUNIT_BUILD_PATH}/stamp
      CMAKE_ARGS
            "-DCMAKE_INSTALL_PREFIX=${PFUNIT_BUILD_PATH}"
            "-DINSTALL_PATH=${PFUNIT_BUILD_PATH}"
            "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
            "-DCMAKE_Fortran_FLAGS=${PFUNIT_COMPILE_FLAGS}"
            "-DMPI=${MPI}"
    )
  endif()

  function(add_pFUnit_test moduleName libraries)
  # Input:
  # - moduleName: name of a module, used for the name of the test
  # - libraries:  list of libraries, that are needed to build the *.pf files
  # The list of libraries has to be given as "${libraries}" when this function is used

    # testSuites.inc is read by driver.F90 to know, which tests to compile
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/testSuites.inc "")

    # storing *.pf files in testlist
    file(GLOB testlist RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*.pf)

    foreach(testfile ${testlist})
      # setting test to the basename of testfile
      string(REPLACE ".pf" "" test ${testfile})

      # the actual preprocessing with pFUnitParser.py
      add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${test}.f90
        COMMAND ${PYTHON} ${PFUNIT_BUILD_PATH}/bin/pFUnitParser.py ${CMAKE_CURRENT_LIST_DIR}/${test}.pf ${CMAKE_CURRENT_BINARY_DIR}/${test}.f90
        DEPENDS pFUnit ${CMAKE_CURRENT_LIST_DIR}/${test}.pf)

      # adding the test to testSuites.inc
      file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/testSuites.inc "ADD_TEST_SUITE(${test}_suite)\n")

      # appending to a list of preprocessed files
      list(APPEND ppSources ${CMAKE_CURRENT_BINARY_DIR}/${test}.f90)
    endforeach()

    add_executable(${moduleName}Tests
      ${PROJECT_SOURCE_DIR}/pFUnit/include/driver.F90
      ${ppSources})

    target_include_directories(${moduleName}Tests
                 # for the pFUnit modules
                 PUBLIC ${PFUNIT_BUILD_PATH}/mod
                 # for testSuite.inc
                 PUBLIC ${CMAKE_CURRENT_BINARY_DIR}
                 )

    target_link_libraries(${moduleName}Tests
      ${PFUNIT_BUILD_PATH}/lib/libpfunit.a
      ${libraries})

    add_test(NAME ${moduleName}Tests
            COMMAND ${EXECUTABLE_OUTPUT_PATH}/${moduleName}Tests --verbose)

    set_tests_properties(${moduleName}Tests
          PROPERTIES LABELS "${moduleName};UnitTest;${PROJECT_NAME}")
  endfunction(add_pFUnit_test)
else()
  function(add_pFUnit_test moduleName libraries)
  endfunction(add_pFUnit_test)
endif()
