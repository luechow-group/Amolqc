# https://developer.amd.com/wordpress/media/2013/12/AOCC-1.2-Flang-the%20Fortran%20Compiler.pdf

set(ALWAYS_FLAGS "-Mno-backslash")

set(WARN_FLAGS "-Wall")
set(NOWARN_FLAGS "-w")
set(PP_FLAGS "-cpp")
set(OPENMP_FLAGS "-mp")
set(LINEINFO_FLAGS "-g")
set(STANDARD_FLAGS "-std=f2008")

set(DEBUG_FLAGS "-fcheck=all")
set(RELEASE_FLAGS "-O3")
