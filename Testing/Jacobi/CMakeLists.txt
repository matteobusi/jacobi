cmake_minimum_required(VERSION 2.6)
project(Jacobi)

#SET(CMAKE_CXX_COMPILER icc)
set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp solvers/jacobisolver.cpp solvers/jacobisolver.h solvers/jacobireport.h solvers/jacobisequentialsolver.cpp solvers/jacobisequentialsolver.h solvers/jacobiffsolver.cpp solvers/jacobiffsolver.h solvers/jacobipthreadsolver.cpp solvers/jacobipthreadsolver.h)
add_executable(Jacobi ${SOURCE_FILES})

## Compiler flags
include_directories ("/home/caos/Master/SPM/SPM_Jacobi/Jacobi/fastflow")

#if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "-O3")        ## So vectorization is enabled
set(CMAKE_CXX_FLAGS "-DNO_DEFAULT_MAPPING")
TARGET_LINK_LIBRARIES(Jacobi pthread)
#endif()