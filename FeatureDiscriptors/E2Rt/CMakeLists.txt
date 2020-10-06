cmake_minimum_required(VERSION 3.10)
project(E2Rt)

set(CMAKE_CXX_FLAGS "-std=c++11")
#寻找opencv库
find_package( OpenCV 3.0 REQUIRED )
#添加头文件
include_directories(${OpenCV_INCLUDE_DIRS})

include_directories( "/usr/local/include/eigen3" )
include_directories( "/usr/local/include/sophus" )


add_executable(E2Rt E2Rt.cpp)
target_link_libraries(E2Rt ${OpenCV_LIBS})