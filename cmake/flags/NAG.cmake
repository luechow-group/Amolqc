# https://www.nag.co.uk/nagware/np/r62_doc/manual/compiler_2_4.html

set(WARN_FLAGS "-w=x95 -w=x77")
set(NOWARN_FLAGS "-w=all")
set(PP_FLAGS "-fpp")
set(OPENMP_FLAGS "-openmp")
set(LINEINFO_FLAGS "-g -gline")

set(DEBUG_FLAGS "-C -C=alias -C=dangling")
set(RELEASE_FLAGS "-O3")
