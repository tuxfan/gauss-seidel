cmake_minimum_required(VERSION 2.8)

project(solver)

option(ENABLE_RBGS "" OFF)

if(ENABLE_RBGS)
  add_definitions(-DRBGS)
else()
  add_definitions(-DGS)
endif()

add_executable(solver solver.cc)
