cmake_minimum_required(VERSION 3.8)
project(untitled)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp
        functions.cpp
        functions.h
        matrix.cpp
        matrix.h
        matrix_fill.cpp
        matrix_fill.h
        residuals.cpp
        residuals.h)
add_executable(untitled ${SOURCE_FILES})