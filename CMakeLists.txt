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
    include/atmospheric_lapse_rate_model.h
    include/atmosphere_widget.h
    include/equal_area_projection.h
    include/mollweide_projection.h
    include/planet_surface_grid.h
    include/planet_surface_point.h
    include/segment_selector_widget.h
    include/surface_point_state.h
    include/surface_temperature_state.h
    include/surface_temperature_scale_widget.h
    include/surface_map_widget.h
    include/surface_map_interpolator.h
    include/surface_map_cache.h
    include/surface_cell.h
    include/surface_height_model.h
    include/surface_heightmap.h
    include/surface_globe_widget.h
    include/surface_point_status_dialog.h
    include/temperature_color_scale.h
    src/atmosphere_chart_widget.cpp
    src/atmosphere_model.cpp
    src/atmospheric_radiation_model.cpp
    src/atmospheric_circulation_model.cpp
    src/atmospheric_lapse_rate_model.cpp
    src/atmosphere_widget.cpp
    src/mollweide_projection.cpp
    src/mode_illustration_widget.cpp
    src/mode_illustration_widget.h
    src/orbit_segment_calculator.cpp
    src/planet_surface_grid.cpp
    src/planet_presets.cpp
    src/segment_selector_widget.cpp
    src/solar_calculator.cpp
    src/solar_display.cpp
    src/surface_map_widget.cpp
    src/surface_map_interpolator.cpp
    src/surface_map_cache.cpp
    src/surface_globe_widget.cpp
    src/surface_height_model.cpp
    src/surface_heightmap.cpp
    src/surface_point_status_dialog.cpp
    src/surface_point_state.cpp
    src/surface_temperature_calculator.cpp
    src/surface_temperature_scale_widget.cpp
    src/surface_temperature_state.cpp
    src/surface_temperature_plot.cpp
    src/temperature_plot_tracker.cpp
    src/temperature_color_scale.cpp
)

target_include_directories(solar_constant PRIVATE include src)
target_link_libraries(solar_constant PRIVATE Qt5::Widgets Qt5::Concurrent Qt5::Charts qwt-qt5)

if (MSVC)
    target_compile_options(solar_constant PRIVATE /W4 /permissive-)
else()
    target_compile_options(solar_constant PRIVATE -Wall -Wextra -pedantic)
endif()

install(TARGETS solar_constant)
