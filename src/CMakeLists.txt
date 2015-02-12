cmake_minimum_required (VERSION 2.8)

# writegrid project (requires cuda)
project (writegrid)
find_package(CUDA)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -g -O3 -arch sm_20)
cuda_add_executable(writegrid writegrid.cu Framework.cc Forcefield.cc)

# henry project (requires cuda)
project (henry)
cuda_add_executable(henry henry.cu Framework.cc Forcefield.cc)

# gcmc project (requires C++11)
project (gcmc)
SET(CMAKE_CXX_FLAGS "-std=c++11")
add_executable(gcmc gcmc.cc Framework.cc Forcefield.cc pocketblocking.cc)

# put binaries one directory up
INSTALL(TARGETS writegrid henry gcmc DESTINATION ${writegrid_SOURCE_DIR}/..)
