cmake_minimum_required(VERSION 3.6)
project(localwl)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O2 -Wall -Wextra -pedantic")

set(SOURCE_FILES main.cpp src/Graph.h src/Graph.cpp src/ColorRefinementKernel.h src/ColorRefinementKernel.cpp src/GraphletKernel.cpp src/ShortestPathKernel.cpp  src/AuxiliaryMethods.h src/AuxiliaryMethods.cpp src/GenerateTwo.cpp src/GenerateTwo.h src/GenerateThree.cpp src/GenerateThreeSampling.cpp src/GenerateThreeSampling.h src/GenerateThree.h)
add_executable(localwl ${SOURCE_FILES})