cmake_minimum_required(VERSION 3.10)

project(NMF)

set(CMAKE_CXX_STANDARD 17)

set(sources main.cc 
	../../src/nmf.cc
	../../src/common.cc)

add_executable(${PROJECT_NAME} ${sources})
