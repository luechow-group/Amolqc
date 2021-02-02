### sets OPTLEVEL

if(OPTLEVEL STREQUAL high)
    set(CMAKE_BUILD_TYPE Release)
elseif(OPTLEVEL STREQUAL debug)
    set(CMAKE_BUILD_TYPE Debug)
elseif(OPTLEVEL STREQUAL no)
    # CMAKE_BUILD_TYPE remains empty
else()
    message( FATAL_ERROR 'Optlevel not known')
endif()
