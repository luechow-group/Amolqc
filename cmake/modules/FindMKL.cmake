### finds mkl in $MKLROOT
find_library(MKL_INTERFACE_LIBRARY
        NAMES mkl_intel_lp64
        PATHS
        $ENV{MKLROOT}/lib
        $ENV{MKLROOT}/lib/intel64
        NO_DEFAULT_PATH
        )

find_library(MKL_SEQUENTIAL_LAYER_LIBRARY
        NAMES mkl_sequential
        PATHS
        $ENV{MKLROOT}/lib
        $ENV{MKLROOT}/lib/intel64
        NO_DEFAULT_PATH)

find_library(MKL_CORE_LIBRARY
        NAMES mkl_core
        PATHS
        $ENV{MKLROOT}/lib
        $ENV{MKLROOT}/lib/intel64
        NO_DEFAULT_PATH)

find_package_handle_standard_args(MKL DEFAULT_MSG
        MKL_INTERFACE_LIBRARY
        MKL_SEQUENTIAL_LAYER_LIBRARY
        MKL_CORE_LIBRARY)
