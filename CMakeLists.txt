cmake_minimum_required(VERSION 3.28)
project(fast_math)

set(CMAKE_CXX_STANDARD 23)
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -Ofast -march=native -flto")

add_executable(fast_math main.cpp
        big_int.h
        big_int.cpp
        big_int_ops.h
        big_int_ops.cpp
        karatsuba.cpp
        karatsuba.h)
