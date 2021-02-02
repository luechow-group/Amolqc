file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/testsuite
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(BASH bash)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/utils/macros.cmake)

if (BASH)
    get_subdirectories(TESTLIST "${CMAKE_CURRENT_BINARY_DIR}/testsuite/run")
    foreach(TEST ${TESTLIST})
        add_test(NAME ${TEST}_serial
                COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/testsuite/testsuite.sh
                ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/bin/amolqc
                ${TEST}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(${TEST}_serial
                PROPERTIES LABELS "testsuite;serial;amolqc")
    endforeach(TEST)
    if (${MPI})
        foreach(TEST ${TESTLIST})
        add_test(NAME ${TEST}_parallel
                COMMAND ${BASH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/testsuite/testsuite.sh
                ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/bin/amolqc
                PARALLEL ${TEST}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
        set_tests_properties(${TEST}_parallel
                PROPERTIES LABELS "testsuite;parallel;amolqc")
        endforeach(TEST)
    endif()
endif()
