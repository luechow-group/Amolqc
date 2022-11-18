# https://docs.nvidia.com/hpc-sdk/compilers/hpc-compilers-user-guide/index.html

set(WARN_FLAGS "")
set(NOWARN_FLAGS "-w")
set(PP_FLAGS "-cpp")
set(OPENMP_FLAGS "-mp")
set(LINEINFO_FLAGS "-g -traceback")
set(STANDARD_FLAGS "-Mstandard")

set(DEBUG_FLAGS "-Mbackslash -C -Minfo=all")
set(RELEASE_FLAGS "-Mbackslash -O2")
