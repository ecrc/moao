################################ 
# load fits library (compulsory)
INCLUDE(FindPkgConfig)
pkg_check_modules(CFITSIO REQUIRED cfitsio)
if(CFITSIO_FOUND)
    set(CFITSIO_LIB ${CFITSIO_LIBRARY_DIRS})
    set(CFITSIO_INC ${CFITSIO_INCLUDE_DIRS})
else()

    set(CFITSIO_LIB $ENV{CFITSIO_LIB} CACHE PATH "Path to cfitsio libraries")
    set(CFITSIO_INC $ENV{CFITSIO_INC} CACHE PATH "Path to cfitsio header")
        if("${CFITSIO_LIB}" STREQUAL "")
    message( FATAL_ERROR "Environment variable CFITSIO_LIB is not set. Specify it with -DCFITSIO_LIB=path/to/cfitsio/libraries")
    endif()
    if("${CFITSIO_INC}" STREQUAL "")
        message( FATAL_ERROR "Environment variable CFITSIO_INC is not set. Specify it with -DCFITSIO_INC=path/to/cfitsio/header")
    endif()
    set(CFITSIO_LIBRARIES "cfitsio")
endif()
add_library(cfitsio UNKNOWN IMPORTED)
set_target_properties(cfitsio PROPERTIES
                        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                        IMPORTED_LOCATION "${CFITSIO_LIB}/libcfitsio.a"
                        INTERFACE_INCLUDE_DIRECTORIES "${CFITSIO_INC}"
                        )
