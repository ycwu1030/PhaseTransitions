CMAKE_MINIMUM_REQUIRED(VERSION 2.8 FATAL_ERROR)

PROJECT(PhaseTransitions CXX)


SET(PhaseTransition_VERSION "0.2")

find_package(GSL 2.1 REQUIRED)
include_directories(include ${GSL_INCLUDE_DIRS})
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -g -std=c++11 -stdlib=libc++")

add_subdirectory(src)
