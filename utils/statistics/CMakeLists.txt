cmake_minimum_required(VERSION 3.10)

add_fortran_library(statistics STATIC
        intList_m.f90
        statistics_m.f90
        newStatistics_m.f90
        vectorStatistics_m.f90
        )

target_link_libraries(statistics
        PRIVATE basics)

add_subdirectory(tests)
