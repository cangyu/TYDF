cmake_minimum_required(VERSION 3.10)
project(PLOT3D)
set(CMAKE_CXX_STANDARD 17)
set(sources main.cc
        ../../src/common.cc
        ../../src/plot3d.cc)
add_executable(${PROJECT_NAME} ${sources})
