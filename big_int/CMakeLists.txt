add_library(big_int_lib SHARED
        include/big_int/big_int.h
        src/big_int.cpp
        include/big_int/big_int_ops.h
        src/big_int_ops.cpp
        include/big_int/karatsuba.h
        src/karatsuba.cpp
        include/big_int/toom_cook.h
        src/toom_cook.cpp
        src/burnikel_ziegler.cpp
        include/big_int/burnikel_ziegler.h
)

target_include_directories(big_int_lib PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

include(GNUInstallDirs)
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/BigIntLib)

set_target_properties(big_int_lib PROPERTIES EXPORT_NAME BigIntLib)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(TARGETS big_int_lib
    EXPORT big_int_lib_config
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(EXPORT big_int_lib_config
    FILE BigIntLibConfig.cmake
    NAMESPACE BigIntLib::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BigIntLib
)