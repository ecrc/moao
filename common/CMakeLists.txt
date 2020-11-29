    cmake_minimum_required(VERSION 3.5.1)

    PROJECT(common C CXX)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

    # directly make an error if in-source build
    if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
        message(FATAL_ERROR "In-source builds are not allowed.\n"
                            "Please create a build directory first and execute cmake configuration from "
                            "this directory. Example: mkdir build && cd build && cmake ..")
    endif()

    #########################
    ### VARIABLES ###########
    #########################
    option(COMMON_TESTING "Generates testing binaries" OFF)
    option(USE_GPU "Enable GPU support" OFF)
    set(CUDA_GENCODE "60 " CACHE STRING "cuda gencode")
    option(USE_OPENMP "Enable OpenMP support" OFF)
    option(USE_MATCOV_TILED "create covariance matrix using rodconan model (else random)" ON)
    option(USE_INTERSAMPLE "Go through the intersample process and generates PSFs" OFF)

if(USE_GPU)
    find_package(CUDA REQUIRED)
    set(GENCODE "")
    string(REPLACE " " ";" CUDA_GENCODE "${CUDA_GENCODE}")
    foreach(code ${CUDA_GENCODE} )
        set(GENCODE ${GENCODE} "-gencode arch=compute_${code},code=sm_${code} ")
    endforeach()
    message(WARNING ${GENCODE}) 
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} ${extra_flags} --shared --compiler-options '-fPIC' ${GENCODE} --compiler-options -Wall")
    message(STATUS "CUDA FLAGS: ${CUDA_NVCC_FLAGS}")
endif()


    # Version
    # Number 0.1.1 is temptative, may differ from actual release.
    SET( COMMON_VERSION_MAJOR 0 )
    SET( COMMON_VERSION_MINOR 1 )
    SET( COMMON_VERSION_PATCH 1 )
    set( COMMON_VERSION "${COMMON_VERSION_MAJOR}.${COMMON_VERSION_MINOR}.${COMMON_VERSION_PATCH}" )

    # Set the RPATH config
    # --------------------
    # use, i.e. don't skip the full RPATH for the build tree
    set(CMAKE_SKIP_BUILD_RPATH  FALSE)
    # when building, use the install RPATH already
    # (automated test will need this)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    # the RPATH to be used when installing
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    list(APPEND public_link_lib "")
    list(APPEND private_link_lib "m")
    ################################ 
    # load fits library (compulsory)
    include(cmake_modules/CMake_loadcfitsio.cmake)
    list(INSERT public_link_lib 0 "${CFITSIO_LIBRARIES}")

    ################################ 
    # load fftw3 library
    # no need to set this up if fftw is in /usr/lib
    include(cmake_modules/CMake_loadfftw.cmake)
    list(INSERT public_link_lib 0 "fftw3")

    if(USE_GPU)
        add_definitions(-DUSE_GPU)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_GPU ")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_GPU ")
    endif()

    set(extra_flags "")

    include(tomo_struct/CMakeLists.txt)
    include(io/CMakeLists.txt)
    include(kernels/CMakeLists.txt)
    include(compute/CMakeLists.txt)
    include(utils/CMakeLists.txt)


    set(SRC ${TOMO_STRUCT_SRC}
            ${IO_SRC}
            ${KERNELS_SRC}
            ${COMPUTE_SRC}
            ${UTILS_SRC}
        )

message(STATUS "source files: ${SRC}")
    #building library
    if(USE_GPU)
        CUDA_ADD_LIBRARY(${PROJECT_NAME} SHARED ${SRC})
        target_link_libraries(${PROJECT_NAME}               ${public_link_lib})
    else()
        ADD_LIBRARY(${PROJECT_NAME} SHARED ${SRC})
        target_link_libraries(${PROJECT_NAME}             PUBLIC  ${public_link_lib})
    endif()
    target_include_directories(${PROJECT_NAME}        PUBLIC ${CFITSIO_INC} ${FFTW_INC})
    target_include_directories(${PROJECT_NAME}        PUBLIC  tomo_struct compute kernels io ./ )

    

## Installing
# add installation for target common
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export
		RUNTIME DESTINATION bin 
		ARCHIVE DESTINATION lib
		LIBRARY DESTINATION lib
		INCLUDES DESTINATION include
)
# copy public include files
install(FILES 
	include/${PROJECT_NAME}/config.h 
	include/${PROJECT_NAME}/make_${PROJECT_NAME}.h
	DESTINATION include/${PROJECT_NAME})


    set(SRC_DEF default_sys/default_sys.cpp)
    PROJECT(default_sys CXX)
    ADD_EXECUTABLE(default_sys ${SRC_DEF})
    target_link_libraries(default_sys common )
