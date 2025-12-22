cmake_minimum_required(VERSION 3.16)
project(SolarConstantCalculator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

option(USE_BUNDLED_QWT "Build Qwt from the bundled 3rdparty sources" ON)

if (USE_BUNDLED_QWT)
    include(ExternalProject)

    set(QWT_INSTALL_DIR "${CMAKE_BINARY_DIR}/qwt-install")
    set(QWT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/qwt")

    if (NOT EXISTS "${QWT_SOURCE_DIR}/qwt.pro")
        message(FATAL_ERROR "Bundled Qwt source not found at ${QWT_SOURCE_DIR}. Ensure 3rdparty/qwt contains Qwt sources.")
    endif()

    if (TARGET Qt6::qmake)
        get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
    else()
        find_program(QT_QMAKE_EXECUTABLE NAMES qmake6 qmake)
    endif()

    if (NOT QT_QMAKE_EXECUTABLE)
        message(FATAL_ERROR "Qt qmake executable not found. Install Qt development tools.")
    endif()

    set(QWT_INSTALL_PREFIX "${QWT_INSTALL_DIR}")
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qwtconfig.pri.in"
        "${CMAKE_BINARY_DIR}/qwtconfig.pri"
        @ONLY
    )

    ExternalProject_Add(qwt
        SOURCE_DIR "${QWT_SOURCE_DIR}"
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_BINARY_DIR}/qwtconfig.pri" <SOURCE_DIR>/qwtconfig.pri
        CONFIGURE_COMMAND "${QT_QMAKE_EXECUTABLE}" <SOURCE_DIR>/qwt.pro
        BUILD_COMMAND ${CMAKE_MAKE_PROGRAM}
        INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
        BUILD_IN_SOURCE 1
        UPDATE_COMMAND ""
    )

    add_library(Qwt::Qwt STATIC IMPORTED GLOBAL)
    set_target_properties(Qwt::Qwt PROPERTIES
        IMPORTED_LOCATION "${QWT_INSTALL_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}qwt${CMAKE_STATIC_LIBRARY_SUFFIX}"
        INTERFACE_INCLUDE_DIRECTORIES "${QWT_INSTALL_DIR}/include"
    )
    add_dependencies(Qwt::Qwt qwt)
else()
    find_package(Qwt REQUIRED)
endif()

add_executable(solar_constant
    src/planet_presets.cpp
    src/solar_calculator.cpp
    src/solar_display.cpp
    src/surface_temperature_calculator.cpp
    src/surface_temperature_plot.cpp
)

target_include_directories(solar_constant PRIVATE include)
target_link_libraries(solar_constant PRIVATE Qt6::Widgets Qwt::Qwt)

if (USE_BUNDLED_QWT)
    add_dependencies(solar_constant qwt)
endif()

if (MSVC)
    target_compile_options(solar_constant PRIVATE /W4 /permissive-)
else()
    target_compile_options(solar_constant PRIVATE -Wall -Wextra -pedantic)
endif()

install(TARGETS solar_constant)
