cmake_minimum_required(VERSION 3.1)
project(Fisrt_Practice CXX)

find_package(OpenMP)

set (CMAKE_CXX_FLAGS "-std=c++17 -O3 ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -U__STRICT_ANSI__ -mtune=native -march=native")

add_library(utils SHARED libs/utils.cpp)

add_executable(main main.cpp)
add_executable(1_1_1_newton first_arxiv/1_1_1_newton.cpp)
add_executable(1_1_1_bisection first_arxiv/1_1_1_bisection.cpp)
add_executable(1_2_eigen first_arxiv/1_2_eigen.cpp)
add_executable(1_3_1_multinewton first_arxiv/1_3_1_multinewton.cpp)
add_executable(1_3_2_multinewton first_arxiv/1_3_2_multinewton.cpp)