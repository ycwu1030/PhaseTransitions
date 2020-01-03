CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

PROJECT(TraceMinimum)

find_package(GSL 2.1 REQUIRED)
include_directories(include ${GSL_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -stdlib=libc++")

# add_subdirectory(src)
add_subdirectory(unit_test)
