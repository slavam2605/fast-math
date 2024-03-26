cmake_minimum_required(VERSION 3.28)
project(fast_math)

set(CMAKE_CXX_STANDARD 23)
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -march=native -flto")

add_subdirectory(big_int)
add_subdirectory(test)
