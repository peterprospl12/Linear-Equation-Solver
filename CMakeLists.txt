cmake_minimum_required(VERSION 3.28)
project(Solving_Linear_Equations)

set(CMAKE_CXX_STANDARD 17)

add_executable(Solving_Linear_Equations main.cpp
        Matrix.h
        Matrix.cpp)

find_package(Python COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})
target_link_libraries(Solving_Linear_Equations ${Python_LIBRARIES})