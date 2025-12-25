cmake_minimum_required(VERSION 3.16)
project(SolarConstantCalculator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets Concurrent Charts)
# find_package( REQUIRED)
# target_link_libraries()
add_executable(solar_constant
    include/atmosphere_chart_widget.h
    include/atmosphere_model.h
    include/atmospheric_radiation_model.h
    include/atmospheric_circulation_model.h
    include/atmosphere_widget.h
    include/segment_selector_widget.h
    src/atmosphere_chart_widget.cpp
    src/atmosphere_model.cpp
    src/atmospheric_radiation_model.cpp
    src/atmospheric_circulation_model.cpp
    src/atmosphere_widget.cpp
    src/mode_illustration_widget.cpp
    src/mode_illustration_widget.h
    src/orbit_segment_calculator.cpp
    src/planet_presets.cpp
    src/segment_selector_widget.cpp
    src/solar_calculator.cpp
    src/solar_display.cpp
    src/surface_temperature_calculator.cpp
    src/surface_temperature_plot.cpp
    src/temperature_plot_tracker.cpp
)

target_include_directories(solar_constant PRIVATE include src)
target_link_libraries(solar_constant PRIVATE Qt5::Widgets Qt5::Concurrent Qt5::Charts qwt-qt5)

if (MSVC)
    target_compile_options(solar_constant PRIVATE /W4 /permissive-)
else()
    target_compile_options(solar_constant PRIVATE -Wall -Wextra -pedantic)
endif()

install(TARGETS solar_constant)
