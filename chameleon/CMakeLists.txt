
set(CHAMELEON ${CMAKE_SOURCE_DIR}/chameleon)

set(MAIN moao)

if(LAPACKE_LIBRARIES MATCHES "mkl")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DUSE_MKL")
endif()

#TODO: include from common
set(SRC_COMMON ${COMMON}/arguments.c
               ${COMMON}/fits.c
               ${COMMON}/tomo_struct.c
               ${COMMON}/noise.c
               ${COMMON}/utils.c
    )
include_directories(${COMMON})

set(SRC ${CHAMELEON}/moao_chameleon.c
        ${COMMON}/reconstructor.c
        ${COMMON}/codelet_dscaldiag.c
        ${COMMON}/codelet_matcov.c
        ${COMMON}/dscaldiag.c
        ${COMMON}/matcov_kernels.c
        ${COMMON}/matcov.c
        ${COMMON}/cee_cvv.c
   ) 

if(${GPU})
set(SRC ${SRC} 
   ${COMMON}/matcov_kernels_gpu.cu
   )
endif()

include_directories(${CHAMELEON})

set(SRC_INTERSAMPLE ${COMMON}/intersample.c
   )

set(SRC ${SRC} ${SRC_COMMON})
if(USE_INTERSAMPLE)
set(SRC ${SRC} ${SRC_INTERSAMPLE})
endif()
