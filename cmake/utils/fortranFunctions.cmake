function(add_fortran_library LIB)
    add_library(${LIB} ${ARGN})

    # set module path to libdir/mod
    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)

    # making LIB_DIR/mod available for libraries linking LIB
    target_include_directories(${LIB} PUBLIC ${LIB_DIR}/mod)
endfunction(add_fortran_library)

function(set_file_compile_options FILE)
    get_property(OLD_FLAGS
            SOURCE ${FILE}
            PROPERTY COMPILE_FLAGS)
    set_property(SOURCE ${FILE}
            APPEND_STRING PROPERTY
            COMPILE_FLAGS "${OLD_FLAGS} ${ARGN}"
            )
endfunction(set_file_compile_options)
