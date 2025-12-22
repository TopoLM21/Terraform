cmake_minimum_required(VERSION 3.16)
project(SolarConstantCalculator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)
find_package(Qwt REQUIRED)

add_executable(solar_constant
    src/planet_presets.cpp
    src/solar_calculator.cpp
    src/solar_display.cpp
    src/surface_temperature_calculator.cpp
    src/surface_temperature_plot.cpp
)

target_include_directories(solar_constant PRIVATE include)
target_link_libraries(solar_constant PRIVATE Qt6::Widgets Qwt::Qwt)

if (MSVC)
    target_compile_options(solar_constant PRIVATE /W4 /permissive-)
else()
    target_compile_options(solar_constant PRIVATE -Wall -Wextra -pedantic)
endif()

install(TARGETS solar_constant)
