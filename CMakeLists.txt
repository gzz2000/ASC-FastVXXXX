cmake_minimum_required(VERSION 3.10)

project(FastVXXXX LANGUAGES CXX)

add_subdirectory(pybind11)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

pybind11_add_module(hamming hamming_parallel.cpp)
pybind11_add_module(find_ps_sites find_ps_sites.cpp)
