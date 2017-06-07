cmake_minimum_required(VERSION 3.6)
project(JacobiFF)

include_directories ("${PROJECT_SOURCE_DIR}/fastflow/")

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(JacobiFF ${SOURCE_FILES})

## Compiler flags
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O3")        ## So vectorization is enabled
    TARGET_LINK_LIBRARIES(JacobiFF pthread)
endif()