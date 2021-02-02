### sets linalg as a list of linear algebra libraries

# adds Amolqc/cmake dir to path that is searched for FindPackage.cmake files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# includes FindPackageHandleStandardArgs.cmake for custom FindPackage.cmake files
include(FindPackageHandleStandardArgs)

if (${MKL})

    if(NOT DEFINED ENV{MKLROOT})
        message( FATAL_ERROR "LAPACK mkl is given but $MKLROOT is not set.")
    endif()

    message("found MKLROOT $ENV{MKLROOT}")
    find_package(MKL REQUIRED)

    set(linalg
            ${MKL_INTERFACE_LIBRARY}
            ${MKL_SEQUENTIAL_LAYER_LIBRARY}
            ${MKL_CORE_LIBRARY}
            )
else()

    if(NOT DEFINED ENV{MATHLIBS})
        # using FindLAPACK.cmake file provided by cmake
        find_package(LAPACK REQUIRED)
        set(linalg
                ${LAPACK_LIBRARIES})
    else()
        message("found MATHLIBS $ENV{MATHLIBS}")
        # using self-written FindMATHLIBS.cmake file
        find_package(MATHLIBS REQUIRED)
        set(linalg
                ${BLAS_LIBRARY}
                ${LAPACK_LIBRARY}
                )
    endif()
endif()
