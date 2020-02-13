cmake_minimum_required(VERSION 3.10)

project(NeutralMapFileIO)

set(CMAKE_CXX_STANDARD 17)

add_executable(${PROJECT_NAME}
	main.cc 
	../../src/nmf.cc
	../../src/common.cc)
