cmake_minimum_required(VERSION 3.16)
project(undistort_image)

set(CMAKE_CXX_FLAGS "-std=c++11")
#寻找opencv库
find_package( OpenCV 3.0 REQUIRED )
#添加头文件
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(undistort_image undistort_image.cpp)
target_link_libraries(undistort_image ${OpenCV_LIBS})

