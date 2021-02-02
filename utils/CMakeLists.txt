cmake_minimum_required(VERSION 3.10)

add_subdirectory(basics)
add_subdirectory(generic)
add_subdirectory(ddda)
add_subdirectory(optimizers)
add_subdirectory(random)
add_subdirectory(statistics)

add_fortran_library(utils STATIC
        cspline_m.f90
        parsing_m.f90
        atom_m.f90
        utils_m.f90
        sorting_m.f90
        gcube_m.f90
        bins3D_m.f90
        numderivs_m.f90
        string_utility_module.f90
        compilerStrings_m.f90
        linAlg_m.f90
        )

target_link_libraries(utils
        PUBLIC basics
        PUBLIC generic
        PUBLIC statistics
        PUBLIC ${linalg}
        INTERFACE ddda
        INTERFACE optimizers
        INTERFACE random
        )

add_subdirectory(tests)
