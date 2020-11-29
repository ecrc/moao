################################ 
# load fftw3 library
# no need to set this up if fftw is in /usr/lib
set(FFTW_LIB $ENV{FFTW_LIB} CACHE PATH "Path to fftw3 libraries")
set(FFTW_INC $ENV{FFTW_INC} CACHE PATH "Path to fftw3 headers")
if("${FFTW_LIB}" STREQUAL "")
    message( FATAL_ERROR "Environment variable FFTW_LIB  is not set. Specify it with -DFFTW_ROOT=path/to/fftw3/libraries")
endif()
if("${FFTW_INC}" STREQUAL "")
    message( FATAL_ERROR "Environment variable FFTW_INC  is not set. Specify it with -DFFTW_ROOT=path/to/fftw3/includes")
endif()
add_library(fftw3 UNKNOWN IMPORTED)
set_target_properties(fftw3 PROPERTIES
                        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                        IMPORTED_LOCATION "${FFTW_LIB}/libfftw3.so"
                        INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INC}"
                        )
