cmake_minimum_required(VERSION 2.8.9)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(CMAKE_MACOSX_RPATH 1)
project (foo)

set(CMAKE_CXX_STANDARD 11) 

find_package( ZLIB REQUIRED )
find_package(HTSlib REQUIRED)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

file(GLOB lib_H
        "utils/*.h"
)

file(GLOB lib_SRC
        "align_stats/*.cpp"
        "utils/*.cpp"
)


message(${HTSlib_INCLUDE_DIR})
include_directories(${HTSlib_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

set (SOURCE AlignmentLibrary.cpp )
add_library(foo MODULE ${SOURCE} ${lib_SRC} ${lib_H})

TARGET_LINK_LIBRARIES(foo hts z)