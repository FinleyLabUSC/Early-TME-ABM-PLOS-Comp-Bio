cmake_minimum_required(VERSION 2.0)
project(eim)

set(CMAKE_CXX_COMPILER "g++-9")
set(CMAKE_CXX_FLAGS "-fopenmp -std=c++14 -O3")
include_directories(inc)

file(GLOB SOURCES "src/*.cpp")

add_executable(eim ${SOURCES})
