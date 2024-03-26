cmake_minimum_required(VERSION 3.28)
project(fast_math)

set(CMAKE_CXX_STANDARD 23)
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -Ofast -march=native -flto")

add_library(big_int_lib
        big_int/big_int.h
        big_int/big_int.cpp
        big_int/big_int_ops.h
        big_int/big_int_ops.cpp
        big_int/karatsuba.cpp
        big_int/karatsuba.h
        big_int/toom_cook.cpp
        big_int/toom_cook.h)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(big_int_tests test/bit_test_with_pow.cpp)
target_link_libraries(big_int_tests GTest::gtest_main big_int_lib)

include(GoogleTest)
gtest_discover_tests(big_int_tests)

