cmake_minimum_required(VERSION 3.10)

add_executable(randomTests randomTests_p.f90)

target_link_libraries(randomTests
        PUBLIC random
        PUBLIC basics)

add_test(NAME RandomTests
         COMMAND randomTests)

set_tests_properties(RandomTests
        PROPERTIES LABELS "utils;moduletest;Amolqc")

add_pFUnit_ctest(normal_test
        TEST_SOURCES normal_tm.pf
        LINK_LIBRARIES basics random
        LABELS utils moduletest Amolqc pFUnit)
