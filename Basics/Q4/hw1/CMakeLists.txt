cmake_minimum_required(VERSION 2.8)

project(SAYHELLO)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF()

include_directories(
${PROJECT_SOURCE_DIR}
${PROJECT_SOURCE_DIR}/include
)

ADD_LIBRARY(hello SHARED src/hello.c)

LINK_LIBRARIES(hello)

INSTALL(FILES include/hello.h DESTINATION /usr/local/include)

INSTALL(TARGETS hello LIBRARY DESTINATION /usr/local/lib)

SET(SRC_LIST useHello.c)
ADD_EXECUTABLE(sayhello useHello.c)
TARGET_LINK_LIBRARIES(sayhello hello)
