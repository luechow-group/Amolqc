# https://www.pgroup.com/resources/docs/17.10/x86/pgi-ref-guide/index.htm#cmdln-options-ref

set(WARN_FLAGS "")
set(NOWARN_FLAGS "-w")
set(PP_FLAGS "-cpp")
set(OPENMP_FLAGS "-mp")
set(LINEINFO_FLAGS "-g")

set(DEBUG_FLAGS "-C -Minfo=all")
set(RELEASE_FLAGS "-O3")
