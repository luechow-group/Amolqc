cmake_minimum_required(VERSION 3.10)

add_fortran_library(generic STATIC
        blockAllocator_m.f90
        genericFilter_m.f90
        plateau_m.f90)

target_link_libraries(generic
        PRIVATE basics
        )

add_subdirectory(tests)
