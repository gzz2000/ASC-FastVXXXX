cmake_minimum_required(VERSION 3.10)

project(FastVXXXX LANGUAGES CXX)

add_subdirectory(pybind11)

pybind11_add_module(hamming hamming_parallel.cpp)
