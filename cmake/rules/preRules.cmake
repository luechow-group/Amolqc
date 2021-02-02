### applying rules from ../config.cmake

if(${CONFIG_FILE})
    include(cmake/config.cmake)

    # 1. MPI option is a boolean
    message("MPI is set to ${MPI} - from config.cmake.")

    # 2. OPTLEVEL only if CMAKE_BUILD_TYPE is not set.
    if(NOT CMAKE_BUILD_TYPE)
        include(cmake/rules/optlevel.cmake)
        message("using OPTLEVEL ${OPTLEVEL} (${CMAKE_BUILD_TYPE}) - from config.cmake.")
    else()
        message("using OPTLEVEL ${CMAKE_BUILD_TYPE} - from variable CMAKE_BUILD_TYPE.")
    endif()

    # 3. LAPACK checked if valid
    if (${MKL})
        message("using MKL - from config.cmake.")
    else()
        message("using LAPACK compiled instead of MKL - from config.cmake.")
    endif()

    # 4. COMPILER sets compiler
    set(CMAKE_Fortran_COMPILER ${COMPILER})
    message("using COMPILER ${COMPILER} - from config.cmake.")

    # warning/error, if COMPILER ist mpi, but MPI is off - or vice versa
    if ((${MPI}) AND (NOT (${CMAKE_Fortran_COMPILER} MATCHES "mpi")))
        message( FATAL_ERROR "MPI set to ${MPI}, but compiler ${CMAKE_Fortran_COMPILER} is non-mpi.")
    elseif((NOT (${MPI})) AND (${CMAKE_Fortran_COMPILER} MATCHES "mpi"))
        message("Warning: Compiler ${CMAKE_Fortran_COMPILER} contains mpi, but MPI is set to ${MPI}.")
    endif ()
else()
    if((NOT DEFINED ENV{FC}) AND (NOT CMAKE_Fortran_COMPILER))
        message( FATAL_ERROR "Could not determine Fortran compiler:
        config.cmake is not read,
        environment variable $FC is empty
        and CMAKE_Fortran_COMPILER is not set.")
    endif()
endif()
