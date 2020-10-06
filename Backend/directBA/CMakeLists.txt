cmake_minimum_required(VERSION 3.10)
project(directBA)

set(CMAKE_CXX_FLAGS "-std=c++11")
#寻找opencv库
find_package(OpenCV 3.0 REQUIRED )
find_Package(g2o REQUIRED)
find_package(Pangolin REQUIRED)
#添加头文件
include_directories(${OpenCV_INCLUDE_DIRS})
include_directories(${g2o_INCLUDE_DIRS})
include_directories(${Pangolin_INCLUDE_DIRS})

include_directories( "/usr/local/include/eigen3" )
include_directories( "/usr/local/include/sophus" )

if (POLICY CMP0072)
    set(OpenGL_GL_PREFERENCE LEGACY)
endif()

SET(g2o_LIBS g2o_cli g2o_ext_freeglut_minimal g2o_simulator g2o_solver_slam2d_linear
        g2o_types_icp g2o_types_slam2d g2o_types_sba g2o_types_slam3d g2o_core g2o_interface
        g2o_solver_csparse g2o_solver_structure_only g2o_csparse_extension g2o_opengl_helper g2o_solver_dense
        g2o_stuff g2o_types_sclam2d g2o_parser g2o_solver_pcg g2o_types_data g2o_types_sim3 cxsparse )


add_executable(${PROJECT_NAME} directBA.cpp)
target_link_libraries(${PROJECT_NAME} ${OpenCV_LIBS} ${g2o_LIBS} ${Pangolin_LIBS})