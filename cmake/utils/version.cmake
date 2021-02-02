# creates the version file

message("getting version")
execute_process(
        COMMAND git describe --abbrev=6 --dirty --always --tags
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(VERSION_FLAG "-DVERSION=\\\"${VERSION}\\\"")
