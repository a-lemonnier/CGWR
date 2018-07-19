cmake_minimum_required(VERSION 2.6)
project(wigner)

add_compile_options(-std=c++11 -o3 -Wall -pedantic)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(wigner main.cpp cgwr.cpp)

install(TARGETS wigner RUNTIME DESTINATION bin)
