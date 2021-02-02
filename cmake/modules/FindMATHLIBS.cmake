### finds lapack, blas in $MATHLIBS

find_library(LAPACK_LIBRARY
        NAMES lapack
        PATHS
        $ENV{MATHLIBS}
        NO_DEFAULT_PATH
        )
find_library(BLAS_LIBRARY
        NAMES blas
        PATHS
        $ENV{MATHLIBS}
        NO_DEFAULT_PATH
        )

find_package_handle_standard_args(LAPACK DEFAULT_MSG
        LAPACK_LIBRARY
        BLAS_LIBRARY)
