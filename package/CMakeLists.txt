
PROJECT(PBIG_CUDA_release)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
FIND_PACKAGE(CUDA)


if(CUDA_FOUND)
	message("CUDA_TOOLKIT_INCLUDE=${CUDA_TOOLKIT_INCLUDE}")
	message("CUDA_CUDART_LIBRARY=${CUDA_CUDART_LIBRARY}")
	message("CUDA_TOOLKIT_ROOT_DIR=${CUDA_TOOLKIT_ROOT_DIR}")
	message("CUDA_SDK_ROOT_DIR=${CUDA_SDK_ROOT_DIR}")
endif(CUDA_FOUND)

INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/common/inc)
INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/../shared/inc)

find_library(CUDA_SDK_CUTIL_LIB cutil_x86_64 "${CUDA_SDK_ROOT_DIR}/lib")
find_library(CUDA_SDK_SHR_LIB shrutil_x86_64 "${CUDA_SDK_ROOT_DIR}/../shared/lib")
find_library(PBIGLIB NAMES pbig PATHS .)


set(CUDA_PROPAGATE_HOST_FLAGS off)
set(FLAGS "-O3 -Wall")
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
                OUTPUT_VARIABLE GCC_VERSION)
if (GCC_VERSION VERSION_GREATER 4.5 OR GCC_VERSION VERSION_EQUAL 4.5)
	set(FLAGS "${FLAGS} -flto")
endif()
set(ARCH  -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20)
		   
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -O3 ${ARCH})
set(CMAKE_CXX_FLAGS "${FLAGS}")


link_libraries( 
				${PBIGLIB}
				${CUDA_SDK_CUTIL_LIB} 
				${CUDA_SDK_SHR_LIB}
				${CUDA_CUDART_LIBRARY})
cuda_add_executable(pbig
						report.cc
						libpbig.a
						main.cc)


