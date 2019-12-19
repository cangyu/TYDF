cmake_minimum_required(VERSION 3.10)

project(Grid-Glue)

set(CMAKE_CXX_STANDARD 17)

set(sources main.cc
	../../src/common.cc
	../../src/nmf.cc
	../../src/plot3d.cc
	../../src/xf.cc
	../../src/glue.cc )

add_executable(${PROJECT_NAME} ${sources})
