cmake_minimum_required(VERSION 3.16)
project(Graph_cpp)

set(CMAKE_CXX_STANDARD 20)

add_executable(Graph_cpp main.cpp Graph.h SimpleGraph.h WeightedGraph.h SnakesAndLadders.h BoardGame.h UnionFind.h)
