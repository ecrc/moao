    cmake_minimum_required(VERSION 3.2.3)


    #set common absolute path
    get_filename_component(COMMON "${CMAKE_SOURCE_DIR}/common" ABSOLUTE)

    include(${CMAKE_SOURCE_DIR}/${project}/CMake_project)
    message(STATUS "project ${PROJECT_NAME} with ${languages}")
    PROJECT(${PROJECT_NAME} ${languages})

    # directly make an error if in-source build
    if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
        message(FATAL_ERROR "In-source builds are not allowed.\n"
                            "Please create a build directory first and execute cmake configuration from "
                            "this directory. Example: mkdir build && cd build && cmake ..")
    endif()

    #########################
    ### VARIABLES ###########
    #########################
    option(MOAO_TESTING "Generates testing binaries" OFF)
    option(GPU "Enable GPU support" OFF)
    option(USE_INTERSAMPLE "Go through the intersample process and generates PSFs" OFF)

    # Version
    # Number 0.1.1 is temptative, may differ from actual release.
    SET( MOAO_VERSION_MAJOR 0 )
    SET( MOAO_VERSION_MINOR 1 )
    SET( MOAO_VERSION_PATCH 1 )
    set( MOAO_VERSION "${MOAO_VERSION_MAJOR}.${MOAO_VERSION_MINOR}.${MOAO_VERSION_PATCH}" )

    # Set the RPATH config
    # --------------------
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    # when building, use the install RPATH already
    # (automated test will need this)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    # the RPATH to be used when installing
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
     
    # Packaging (make package)
    SET(CPACK_PACKAGE_VERSION ${MOAO_VERSION})
    SET(CPACK_GENERATOR "TGZ")
    INCLUDE(CPack)

    # Define the common include directory of the project
    include_directories( ${CMAKE_SOURCE_DIR}/include )

    ## CMAKE MODULES :: ECRC
    ## REQUIRED FOR TESTS TO LINK LIBRARIES
    include(${PROJECT_SOURCE_DIR}/cmake_modules/CMake_initEcrc)


    list(APPEND link_lib "")
    ################################ 
    # load fits library (compulsory)
    set(CFITSIO_ROOT $ENV{CFITSIO_ROOT} CACHE PATH "Path to the root of cfitsio")
    if("${CFITSIO_ROOT}" STREQUAL "")
        message( FATAL_ERROR "Environment variable CFITSIO_ROOT is not set. Specify it with -DCFITSIO_ROOT=path/to/cfitsio")
    endif()
        message( WARNING "Environment variable CFITSIO_ROOT ${CFITSIO_ROOT}")
    include_directories(${CFITSIO_ROOT}/include)
    link_directories(${CFITSIO_ROOT})
    list(INSERT link_lib 0 "cfitsio" )

    
    ################################ 
    # load fftw3 library
    # no need to set this up if fftw is in /usr/lib
    set(FFTW_ROOT $ENV{FFTW_ROOT} CACHE PATH "Path to the root of fftw3")
    if("${FFTW_ROOT}" STREQUAL "")
        message( FATAL_ERROR "Environment variable FFTW_ROOT is not set. Specify it with -DFFTW_ROOT=path/to/fftw3")
    endif()
    include_directories(${FFTW_ROOT}/include)
    link_directories( ${FFTW_ROOT} )
    list(INSERT link_lib 0 "fftw3" )


    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ")
    set(extra_flags "")

    if(${USE_LAPACK})
        include (${CMAKE_SOURCE_DIR}/cmake_modules/CMake_loadLapack)
    	set(extra_flags "${extra_flags} -DUSE_LAPACK")     
	if(USE_OPENMP)
            find_package(OpenMP)
            set(extra_flags "${extra_flags} -DUSE_OPENMP")
        endif()
    endif()

    if(${USE_CHAMELEON})
        include (${CMAKE_SOURCE_DIR}/cmake_modules/CMake_loadChameleon)
    	set(extra_flags "${extra_flags} -DUSE_CHAMELEON")     
    endif()

    if(${USE_INTERSAMPLE})
        set(extra_flags "${extra_flags} -DUSE_INTERSAMPLE")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${extra_flags}")

    if(${GPU})
        set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${extra_flags} --shared --compiler-options '-fPIC' -gencode arch=compute_30,code=sm_30 --compiler-options -Wall")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_GPU ")
    endif()



    # MOAO code
    include(${CMAKE_SOURCE_DIR}/${project}/CMakeLists.txt)

    #building library
    #if(${GPU})
    #    CUDA_ADD_LIBRARY(moao_chameleon SHARED ${SRC} )
    #    set_target_properties(moao_chameleon PROPERTIES LINKER_LANGUAGE C )
    #else()
    #    ADD_LIBRARY(moao_chameleon SHARED ${SRC} )
    #endif()
    #target_link_libraries(moao_chameleon "${link_lib}" )

    #buiding executable
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_CHAMELEON -DUSE_MATCOV_TILED -DUSE_INTERSAMPLE -DCHECK_CHAMELEON")     
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_MATCOV_TILED ")     
    if(${GPU})
        CUDA_ADD_EXECUTABLE(${MAIN} ${SRC})
        set_target_properties(${MAIN} PROPERTIES LINKER_LANGUAGE C )
    else()
        ADD_EXECUTABLE(${MAIN} ${SRC})
    endif()
    target_link_libraries(${MAIN} "${link_lib}")
    #add_dependencies()


