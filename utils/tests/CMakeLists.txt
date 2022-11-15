cmake_minimum_required(VERSION 3.10)

add_fortran_library(old_utils_tests STATIC
        atom_tm.f90
        )

target_link_libraries(old_utils_tests
        PUBLIC utils
        )

add_executable(oldUtilsTests oldUtilsTests_p.f90)

target_link_libraries(oldUtilsTests
        old_utils_tests)

add_test(NAME old_utils_tests
         COMMAND oldUtilsTests)
set_tests_properties(old_utils_tests
        PROPERTIES LABELS "utils;moduletest;Amolqc")

add_pFUnit_ctest(utils_tests
        TEST_SOURCES bins3D_tm.pf parsing_tm.pf sorting_tm.pf
        LINK_LIBRARIES utils
        LABELS utils moduletest Amolqc pFUnit)
