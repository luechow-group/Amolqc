# https://software.intel.com/en-us/node/677967

set(WARN_FLAGS "-warn all")
set(NOWARN_FLAGS "-warn none")
set(PP_FLAGS "-fpp")
set(OPENMP_FLAGS "-fopenmp")
set(LINEINFO_FLAGS "-g -traceback")
set(VECTORINFO_FLAGS "-qopt-report -qopt-report-phase=vec -qopt-report-annotate=html")

set(DEBUG_FLAGS "-check all")
set(RELEASE_FLAGS "-O3")
