    cmake_minimum_required(VERSION 3.5.1)

    PROJECT(moao_lapack C CXX )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

    # directly make an error if in-source build
    if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
        message(FATAL_ERROR "In-source builds are not allowed.\n"
                            "Please create a build directory first and execute cmake configuration from "
                            "this directory. Example: mkdir build && cd build && cmake ..")
    endif()


    include (${CMAKE_SOURCE_DIR}/cmake_modules/CMake_loadLapack)

    include(compute/CMakeLists.txt)
    set(SRC ${COMPUTE_SRC})

    ##static library
    #ADD_LIBRARY(${PROJECT_NAME}_static STATIC ${SRC})
    ##target_include_directories(${PROJECT_NAME}_static PUBLIC compute )
    #target_link_libraries(${PROJECT_NAME}_static PUBLIC common_static openblas)
    #dynamic library
    ADD_LIBRARY(${PROJECT_NAME} SHARED ${SRC})
    target_include_directories(${PROJECT_NAME} PUBLIC compute ${LAPACKE_INCLUDE_DIRS}  )
    target_link_libraries(${PROJECT_NAME} ${LAPACKE_LIBRARIES} common)



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
