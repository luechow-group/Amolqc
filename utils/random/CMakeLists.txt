cmake_minimum_required(VERSION 3.10)

add_fortran_library(random STATIC
        mrg_m.f90
        mt_m.f90
        random_m.f90
        normal_m.f90
        )


target_link_libraries(random
        PRIVATE basics)

add_subdirectory(tests)
