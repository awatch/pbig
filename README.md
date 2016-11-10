========================================================================================
Implementation of Pairwise Box Intersection checking on Graphics processing units (PBIG)
Date: Jan. 2013
Version: 0.2
Developer:
	Shih-Hsiang Lo (awatch@gmail.com)
========================================================================================

The implementation requires the packages:
	cmake ( >= 2.8.5 ) for cuda support
	g++ 
	nvcc 

Before compiliation, several nvironment variables below are required:
	export LD_LIBRARY_PATH='path to cuda runtime library'
		e.g., /usr/local/cuda/lib64
	export NVSDKCUDA_ROOT='path to cuda sdk'
		e.g., /usr/local/cuda/NVIDIA_GPU_Computing_SDK/C

The implementation needs several CUDA libraries:
	libcutil
	libshrutil

Project struture:
	src			: source files
	benchmarks	: three benchmarks
	package		: an example of using pbig library


Compilation:
	1. Proceed to src/pbig
	2. Please read src/pbig/README

Testing:
	1. A simple test by performing ./test.sh
	2. Benchmarks
		go to benchmarks/ and check the 'README' file inside
		
