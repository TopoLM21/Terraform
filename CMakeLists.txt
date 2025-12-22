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

    find_program(QT_QMAKE_EXECUTABLE NAMES qmake6)

    if (NOT QT_QMAKE_EXECUTABLE)
        message(FATAL_ERROR "Qt qmake6 executable not found. Install Qt development tools.")
    endif()

    find_program(GNU_MAKE_EXECUTABLE NAMES gmake make)

    if (NOT GNU_MAKE_EXECUTABLE)
        message(FATAL_ERROR "GNU make executable not found. Install make/gmake to build Qwt.")
    endif()

    set(QWT_INSTALL_PREFIX "${QWT_INSTALL_DIR}")
    set(QWT_BINARY_DIR "${CMAKE_BINARY_DIR}/_deps/qwt-build")
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qwtconfig.pri.in"
        "${CMAKE_BINARY_DIR}/qwtconfig.pri"
        @ONLY
    )

    ExternalProject_Add(qwt
        SOURCE_DIR "${QWT_SOURCE_DIR}"
        BINARY_DIR "${QWT_BINARY_DIR}"
        PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_BINARY_DIR}/qwtconfig.pri" <SOURCE_DIR>/qwtconfig.pri
        CONFIGURE_COMMAND "${QT_QMAKE_EXECUTABLE}" -o <BINARY_DIR>/Makefile <SOURCE_DIR>/qwt.pro
        BUILD_COMMAND ${GNU_MAKE_EXECUTABLE}
        INSTALL_COMMAND ${GNU_MAKE_EXECUTABLE} install
        UPDATE_COMMAND ""
    )
    ExternalProject_Add_Step(qwt create_moc_dir
        COMMAND ${CMAKE_COMMAND} -E make_directory <BINARY_DIR>/src/moc
        DEPENDEES configure
        DEPENDERS build
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
