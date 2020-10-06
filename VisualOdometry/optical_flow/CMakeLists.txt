cmake_minimum_required(VERSION 3.10)
project(optical_flow)

set(CMAKE_CXX_FLAGS "-std=c++11")
#寻找opencv库
find_package( OpenCV 3.0 REQUIRED )
#添加头文件
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories( "/usr/local/include/eigen3" )

add_executable(optical_flow optical_flow.cpp)
target_link_libraries(optical_flow ${OpenCV_LIBS})