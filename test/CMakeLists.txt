include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
set(INSTALL_GTEST OFF)
set(INSTALL_GMOCK OFF)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_executable(big_int_tests big_test_with_pow.cpp
        simple_division_test.cpp)
target_link_libraries(big_int_tests GTest::gtest_main big_int_lib)

include(GoogleTest)
gtest_discover_tests(big_int_tests)