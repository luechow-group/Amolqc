cmake_minimum_required(VERSION 3.10)

add_fortran_library(minimizers STATIC
        fctn_m.f90
        minimizer_factory_m.f90
        minimizer_ws_bfgst_m.f90
        minimizer_fire_m.f90
        minimizer_bfgs_m.f90
        minimizer_ws_factory_m.f90
        minimizer_m.f90
        minimizer_ws_steep_desc_m.f90
        singularityCorrection_m.f90
        singularityParticles_m.f90
        minimizer_steep_desc_m.f90
        line_search_ws_simple_m.f90
        line_search_simple_m.f90
        minimizer_w_sing_m.f90
        minimizer_ws_fire_m.f90
        minimizer_ws_bfgs_m.f90
        minimizer_ws_newton_m.f90
        minimizer_ws_none_m.f90
        line_search_weak_wolfe_m.f90
        line_search_m.f90
        )

target_link_libraries(minimizers
        PRIVATE utils
        )

add_subdirectory(tests)
