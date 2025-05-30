cmake_minimum_required(VERSION 3.30)
project(cube_eq_decan)

set(CMAKE_CXX_STANDARD 20)

add_executable(cube_eq_decan main.cpp
        functions_for_modeling.h
        fluid.h
        "phase_diagram.h")
