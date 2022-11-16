# adding pFUnit support

if (${PFUNIT})
  if (DEFINED ENV{PFUNIT_DIR})
    message("using pFUnit installed at $ENV{PFUNIT_DIR}")
    find_package(PFUNIT)
  else ()
    set(PFUNIT_VERSION "v4.6.0")
    message("PFUNIT_DIR not set, downloading pFUnit ${PFUNIT_VERSION}")
    include(FetchContent)
    set(FETCHCONTENT_BASE_DIR ${CMAKE_BINARY_DIR}/pFUnit)
    FetchContent_Declare(
      PFUNIT
      GIT_REPOSITORY "https://github.com/Leonard-Reuter/pFUnit"
      GIT_TAG 286af89755d2c54558b7d43836da6d88593bb87c 
    )  
    FetchContent_MakeAvailable(PFUNIT)
  endif()
else()
  # dummy function replacing add_pfunit_ctest
  function(add_pfunit_ctest arg)
  endfunction()
endif()

