PROJECT(Intersection_CUDA)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

INCLUDE(config.cmake)
FIND_PACKAGE(CUDA)

if(CUDA_FOUND)
	message("CUDA_TOOLKIT_INCLUDE=${CUDA_TOOLKIT_INCLUDE}")
	message("CUDA_CUDART_LIBRARY=${CUDA_CUDART_LIBRARY}")
	message("CUDA_TOOLKIT_ROOT_DIR=${CUDA_TOOLKIT_ROOT_DIR}")
	message("CUDA_SDK_ROOT_DIR=${CUDA_SDK_ROOT_DIR}")
endif(CUDA_FOUND)

if (ALGO STREQUAL "PBIG")
	add_definitions(-DUSE_PBIG)
endif(ALGO STREQUAL "PBIG")

if (ALGO STREQUAL "PRI")
	add_definitions(-DUSE_PBIG)
	set(NStreams 1)
	set(NOuts 1)
	set(BALANCE 0)
	set(COMPRESS 1)
endif(ALGO STREQUAL "PRI")


##### set alog's optimizations #####
add_definitions(-DNTPB=${NTPB})
add_definitions(-DNStreams=${NStreams})
add_definitions(-DBALANCE=${BALANCE})
add_definitions(-DMERGED=${MERGED})
add_definitions(-DNOuts=${NOuts})
add_definitions(-DDIMS=${DIMS})
add_definitions(-DVERIFY=${VERIFY})
add_definitions(-DPROFILING=${PROFILING})
if (COMPRESS EQUAL 2)
	add_definitions(-DCOMPRESSION=1)
	add_definitions(-DRICE_COMPRESSION=1)
elseif (COMPRESS EQUAL 1)
	add_definitions(-DCOMPRESSION=1)
endif(COMPRESS EQUAL 2)

#if (CMAKE_BUILD_TYPE STREQUAL "Release")
	add_definitions(-DRELEASE=1)
#else(CMAKE_BUILD_TYPE STREQUAL "Release")
#	message("No Release VERSION")
#	add_definitions(-DRELEASE=0)
#endif(CMAKE_BUILD_TYPE STREQUAL "Release")

INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/common/inc)
INCLUDE_DIRECTORIES(${CUDA_SDK_ROOT_DIR}/../shared/inc)
find_library(CUDA_SDK_CUTIL_LIB cutil_x86_64 "${CUDA_SDK_ROOT_DIR}/lib")
find_library(CUDA_SDK_SHR_LIB shrutil_x86_64 "${CUDA_SDK_ROOT_DIR}/../shared/lib")


set(CUDA_PROPAGATE_HOST_FLAGS off)
set(FLAGS "-O3 -Wall")
set(L1 -Xptxas -dlcm=cg)
set(REGS  -maxrregcount=30)
set(PTXSHOW  --ptxas-options=-v)
set(ARCH  
		-gencode=arch=compute_30,code=sm_30 
		-gencode=arch=compute_30,code=compute_30
		-gencode=arch=compute_20,code=sm_20 
		-gencode=arch=compute_20,code=compute_20
		-gencode=arch=compute_13,code=sm_13
		-gencode=arch=compute_13,code=compute_13 
)
		   
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} ${PTXSHOW} -O3 ${ARCH})
set(CMAKE_CXX_FLAGS "${FLAGS}")



######## for building PBIG and PRI ##################
if (${ALGO} STREQUAL "PBIG" OR ${ALGO} STREQUAL "PRI")
	cuda_add_library(pbig STATIC 
						pbig.cu
						decode.cc
						)
	link_libraries( pbig	
					${CUDA_SDK_CUTIL_LIB} 
					${CUDA_SDK_SHR_LIB}
					${CUDA_CUDART_LIBRARY})
	cuda_add_executable(check.${ALGO} 
						hybrid.cc
						report.cc
						#decode.cc
						main.cc)
endif()




#target_link_libraries(check.${ALGO}
#						hybrid
#						${CUDA_SDK_CUTIL_LIB} 
#						${CUDA_SDK_SHR_LIB}
#						${CUDA_CUDART_LIBRARY}
#						)




######  for installation ########################
install(FILES 
			pbig.h
			${CMAKE_CURRENT_BINARY_DIR}/libpbig.a 
		DESTINATION
			${CMAKE_CURRENT_SOURCE_DIR}/../../package)

