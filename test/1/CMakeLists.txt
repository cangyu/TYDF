cmake_minimum_required(VERSION 3.10)

project(FluentMeshReader)

set(CMAKE_CXX_STANDARD 11)

include_directories(../../inc)

add_executable(${PROJECT_NAME} main.cc)