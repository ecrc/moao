    cmake_minimum_required(VERSION 3.5.1)


    #set common absolute path
    #get_filename_component(COMMON "${CMAKE_SOURCE_DIR}/common" ABSOLUTE)

    include(${CMAKE_SOURCE_DIR}/${project}/CMake_project)
    #message(STATUS "project ${PROJECT_NAME} with ${languages}")
    PROJECT(moao C CXX)
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
    option(MOAO_TESTING "Generates testing binaries" OFF)
    option(USE_GPU "Enable GPU support" OFF)
    set(CUDA_GENCODE "60" CACHE STRING "cuda gencode")
    option(USE_OPENMP "Enable OpenMP support" OFF)
    option(USE_INTERSAMPLE "Go through the intersample process and generates PSFs" OFF)
    option(USE_MATCOV_TILED "create covariance matrix using rodconan model (else random)" ON)
    option(DEBUG "use debug flags" OFF)

    set(GENCODE "")
    string(REPLACE " " ";" CUDA_GENCODE "${CUDA_GENCODE}")
    foreach(code ${CUDA_GENCODE} )
        set(GENCODE ${GENCODE} "-gencode arch=compute_${code},code=sm_${code} ")
    endforeach()
    message(WARNING ${GENCODE}) 

    if(${DEBUG})
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
    endif()
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
    #include_directories( ${CMAKE_SOURCE_DIR}/include )

    # CMAKE MODULES :: ECRC
    # REQUIRED FOR TESTS TO LINK LIBRARIES
    include(${PROJECT_SOURCE_DIR}/cmake_modules/CMake_initEcrc)



	#if(USE_OPENMP)#TODO
    #    find_package(OpenMP)
    #    add_definitions(-DUSE_OPENMP ${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS})
    #endif()

    if(USE_INTERSAMPLE)
        add_definitions(-DUSE_INTERSAMPLE)
    endif()
    if(USE_MATCOV_TILED)
        add_definitions(-DUSE_MATCOV_TILED)
    endif()
    if(${project} STREQUAL "chameleon")
        add_definitions(-DUSE_CHAMELEON)
    endif()
    if(LAPACKE_LIBRARIES MATCHES "mkl")
        add_definitions( -DUSE_MKL)
    endif()
    if(USE_GPU)
        add_definitions(-DUSE_GPU)
    endif()

    add_subdirectory(common)
    if(NOT ${project} STREQUAL "common")
        add_subdirectory(${project})
        set(SRC ${project}/moao_${project}.cpp)
    
        
    link_directories(${LINK_DIRS})
        if(USE_GPU)
            CUDA_ADD_EXECUTABLE(${PROJECT_NAME} ${SRC})
            set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C )
            ADD_EXECUTABLE(check_${PROJECT_NAME} ${project}/check_moao_${project}.cpp)
            set_target_properties(check_${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C )
        else()
            ADD_EXECUTABLE(${PROJECT_NAME} ${SRC})
            ADD_EXECUTABLE(check_${PROJECT_NAME} ${project}/check_moao_${project}.cpp)
        endif()
        target_include_directories(check_${PROJECT_NAME}  PRIVATE ${project} common)
        target_include_directories(      ${PROJECT_NAME}  PRIVATE ${project} common)
        target_link_libraries(check_${PROJECT_NAME} common moao_${project})
        target_link_libraries(      ${PROJECT_NAME} common moao_${project})

    endif()


    ## Installing
# add installation for target common
install(TARGETS ${PROJECT_NAME} EXPORT ${PROJECT_NAME}_export
		RUNTIME DESTINATION ${PREFIX}/bin 
		ARCHIVE DESTINATION ${PREFIX}/lib
		LIBRARY DESTINATION ${PREFIX}/lib
		INCLUDES DESTINATION ${PREFIX}/include
)
# copy public include files
message(WARNING "destination:include/${PROJECT_NAME}")
install(FILES 
	include/${PROJECT_NAME}/config.h 
	include/${PROJECT_NAME}/make_${PROJECT_NAME}.h
	DESTINATION include/${PROJECT_NAME})


# create export file which can be imported by other cmake projects
#install(EXPORT ${PROJECT_NAME}_export
#   NAMESPACE MOAO::
#   DESTINATION lib/cmake/${PROJECT_NAME}
#)
