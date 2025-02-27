cmake_minimum_required(VERSION 3.10)
project(project1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

# set(CMAKE_BUILD_TYPE Debug)
# set(CMAKE_CXX_FLAGS_DEBUG "-g")

add_executable(project1 main.cpp data_structures.cpp watch_list.cpp utilities.h decision.h vsids.cpp moms.cpp)