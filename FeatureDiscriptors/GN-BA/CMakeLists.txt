cmake_minimum_required(VERSION 3.10)
project(GN-BA)

set(CMAKE_CXX_FLAGS "-std=c++11")


include_directories( "/usr/local/include/eigen3" )
include_directories( "/usr/local/include/sophus" )


add_executable(GN-BA GN-BA.cpp)
