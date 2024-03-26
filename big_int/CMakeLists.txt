add_library(big_int_lib SHARED
        include/big_int/big_int.h
        src/big_int.cpp
        include/big_int/big_int_ops.h
        src/big_int_ops.cpp
        include/big_int/karatsuba.h
        src/karatsuba.cpp
        include/big_int/toom_cook.h
        src/toom_cook.cpp)

target_include_directories(big_int_lib PUBLIC include)