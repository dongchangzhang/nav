cmake_minimum_required(VERSION 3.10)
project(nav)

set(CMAKE_CXX_STANDARD 11)

add_executable(nav src/main.cpp 
    src/nav.h src/nav.cpp 
    src/tool.h src/tool.cpp 
    src/landmark.h src/landmark.cpp 
    src/orbit.h src/point.h src/point3.h src/constants.h)
