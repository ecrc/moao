set(TOMO_STRUCT_PATH tomo_struct)

set(TOMO_STRUCT_SRC
    ${TOMO_STRUCT_PATH}/tomo_struct.cpp
    ${TOMO_STRUCT_PATH}/atmParams.cpp
    ${TOMO_STRUCT_PATH}/sysParams.cpp
    ${TOMO_STRUCT_PATH}/noise.cpp
    )

if(USE_INTERSAMPLE)
    set(TOMO_STRUCT_SRC ${TOMO_STRUCT_SRC}
        ${TOMO_STRUCT_PATH}/intersample.cpp)
endif()

