set(KERNELS_PATH kernels)

set(KERNELS_SRC
    ${KERNELS_PATH}/matcov_kernels.cpp
    )

if(USE_GPU)
    set(KERNELS_SRC ${KERNELS_SRC}
        ${KERNELS_PATH}/matcov_kernels_gpu.cu
        )
endif()
