cmake_minimum_required(VERSION 3.10)

add_fortran_library(utils_tests STATIC
        atom_tm.f90
        )

target_link_libraries(utils_tests
        PUBLIC utils
        )

add_executable(oldUtilsTests oldUtilsTests_p.f90)

target_link_libraries(oldUtilsTests
        utils_tests)

add_test(NAME OldUtilsTests
         COMMAND oldUtilsTests)
set_tests_properties(OldUtilsTests
        PROPERTIES LABELS "utils;moduletest;amolqc")

add_pFUnit_test("Utils" "utils")
