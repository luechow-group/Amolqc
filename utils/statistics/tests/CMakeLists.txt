cmake_minimum_required(VERSION 3.10)

add_fortran_library(statistics_tests STATIC
        intList_tm.f90
        newStatistics_tm.f90
        statistics_tm.f90
        )

target_link_libraries(statistics_tests
        PUBLIC statistics
        PUBLIC basics
        )

add_executable(statisticsTests statisticsTests_p.f90)

target_link_libraries(statisticsTests
        statistics_tests)

add_test(NAME StatisticsTests_serial
         COMMAND statisticsTests)
set_tests_properties(StatisticsTests_serial
        PROPERTIES LABELS "utils;moduletest;Amolqc;serial")

if(${MPI})
    add_test(NAME StatisticsTests_parallel
            COMMAND mpiexec -n ${NPROC} ${EXECUTABLE_OUTPUT_PATH}/statisticsTests)
    set_tests_properties(StatisticsTests_parallel
            PROPERTIES LABELS "utils;moduletest;Amolqc;parallel")
endif()
