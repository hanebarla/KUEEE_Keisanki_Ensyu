cmake_minimum_required(VERSION 3.1)
project(Fisrt_Practice CXX)

find_package(OpenMP)

set (CMAKE_CXX_FLAGS "-std=c++17 -O3 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -U__STRICT_ANSI__ -mtune=native -march=native")

add_executable(main main.cpp)