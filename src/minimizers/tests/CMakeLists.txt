cmake_minimum_required(VERSION 3.10)

add_fortran_library(minimizers_tests STATIC
        minimizer_tm.f90
        minimizer_ws_tm.f90
        myfctn_m.f90
        myfctn_tm.f90
        )

target_link_libraries(minimizers_tests
        PUBLIC minimizers
        PUBLIC utils
        )

add_executable(minimizersTests minimizersTests_p.f90)

target_link_libraries(minimizersTests
        minimizers_tests)

add_test(NAME MinimizersTests
         COMMAND minimizersTests)
set_tests_properties(MinimizersTests
        PROPERTIES LABELS "minimizers;moduletest;amolqc")
