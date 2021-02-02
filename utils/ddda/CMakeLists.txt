cmake_minimum_required(VERSION 3.10)

add_fortran_library(ddda STATIC
        boltzmann_m.f90
        decorrelation_m.f90
        parallelDDDA_m.f90
        statistic_m.f90)

target_link_libraries(ddda
        PRIVATE basics
        PRIVATE optimizers
        PRIVATE generic
        )

add_subdirectory(tests)
