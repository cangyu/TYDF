cmake_minimum_required(VERSION 3.10)

project(Plot3DGridIO)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME}
		main.cc
		../../src/common.cc
		../../src/plot3d.cc)
