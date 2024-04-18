cmake_minimum_required(VERSION 3.28)
project(Solving_Linear_Equations)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")
add_executable(Solving_Linear_Equations main.cpp
        Matrix.h
        tasks.h
        plots.h
        tasks.cpp
        plots.cpp
)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
endif()

find_package(Python COMPONENTS Interpreter Development)
include_directories(${Python_INCLUDE_DIRS})
target_link_libraries(Solving_Linear_Equations ${Python_LIBRARIES})