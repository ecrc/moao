    cmake_minimum_required(VERSION 3.5.1)

    PROJECT(moao_chameleon C CXX)

    # directly make an error if in-source build
    if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
        message(FATAL_ERROR "In-source builds are not allowed.\n"
                            "Please create a build directory first and execute cmake configuration from "
                            "this directory. Example: mkdir build && cd build && cmake ..")
    endif()

    include (${CMAKE_SOURCE_DIR}/cmake_modules/CMake_loadChameleon)

    if(STARPU_FOUND)
        include(runtime/starpu/CMakeLists.txt)
        set(RUNTIME_INC runtime/starpu/codelet)
    endif()

    include(compute/CMakeLists.txt)
    set(SRC ${COMPUTE_SRC}
            ${CODELET_SRC}
    )

    link_directories(${CHAMELEON_LIBRARY_DIRS})
    set(LINK_DIRS "${CHAMELEON_LIBRARY_DIRS}" ) #PARENT_SCOPE)
    message(STATUS "LINK DIR: ${CHAMELEON_LIBRARY_DIRS}")
    set(LINK_DIRS ${LINK_DIRS} PARENT_SCOPE)
    message(STATUS "LINK DIR2: ${LINK_DIRS}")
    if(USE_GPU)
        find_package(CUDA REQUIRED)
        add_definitions(-DUSE_GPU)
        CUDA_ADD_LIBRARY(${PROJECT_NAME} SHARED ${SRC})
        set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE C )
        target_link_libraries(${PROJECT_NAME}        ${CHAMELEON_LIBRARIES_DEP} common)
    else()
        ADD_LIBRARY(${PROJECT_NAME} SHARED ${SRC})
        target_link_libraries(${PROJECT_NAME} PUBLIC ${CHAMELEON_LIBRARIES_DEP} common)
    endif()
        target_include_directories(${PROJECT_NAME} PUBLIC compute  ./ ${RUNTIME_INC} ${CHAMELEON_INCLUDE_DIRS_DEP} )

## Installing
# add installation for target lapack
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


## create export file which can be imported by other cmake projects
#install(EXPORT ${PROJECT_NAME}_export
#	NAMESPACE MOAO::
#	DESTINATION lib/cmake/${PROJECT_NAME}
#)
