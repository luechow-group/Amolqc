cmake_minimum_required(VERSION 3.10)

add_fortran_library(optimizers STATIC
        lbfgsb3_m.f90
        nl2sol.f90
        nl2sol_i.f90
        )

target_link_libraries(optimizers
        PRIVATE basics)

if(${WARNINGS})
    set_file_compile_options(nl2sol.f90 ${NOWARN_FLAGS})
    set_file_compile_options(lbfgsb3_m.f90 ${NOWARN_FLAGS})
endif()