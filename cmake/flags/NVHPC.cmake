# https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html

set(WARN_FLAGS "")
set(NOWARN_FLAGS "-w")
set(PP_FLAGS "-cpp")
set(OPENMP_FLAGS "-mp")
set(LINEINFO_FLAGS "-g -traceback")

set(DEBUG_FLAGS "-C -Minfo=all")
set(RELEASE_FLAGS "-O2")
