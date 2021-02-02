cmake_minimum_required(VERSION 3.10)

add_fortran_library(basics STATIC
        error_m.f90
        globalUtils_m.f90
        kinds_m.f90
        verbosity_m.f90
        )
