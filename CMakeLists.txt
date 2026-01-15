cmake_minimum_required(VERSION 3.16)
project(SolarConstantCalculator VERSION 0.1.77 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets Concurrent Charts)
# find_package( REQUIRED)
# target_link_libraries()
add_executable(solar_constant
    include/atmosphere_chart_widget.h
    include/atmosphere_model.h
    include/atmospheric_cell_state.h
    include/atmospheric_radiation_model.h
    include/layered_radiation_model.h
    include/radiation_model.h
    include/radiation_model_utils.h
    include/atmospheric_lapse_rate_model.h
    include/atmospheric_thermodynamics.h
    include/radiative_convective_profile.h
    include/atmosphere_widget.h
    include/co2_convective_adjustment.h
    include/emission_layer_model.h
    include/equal_area_projection.h
    include/geothermal_flux_model.h
    include/mollweide_projection.h
    include/planet_surface_grid.h
    include/planet_surface_point.h
    include/segment_selector_widget.h
    include/height_color_scale.h
    include/surface_point_state.h
    include/surface_temperature_state.h
    include/spectral_band_model.h
    include/subsurface_grid.h
    include/subsurface_temperature_solver.h
    include/thermal_conductivity_model.h
    include/surface_temperature_scale_widget.h
    include/surface_height_scale_widget.h
    include/surface_wind_scale_widget.h
    include/surface_pressure_scale_widget.h
    include/surface_realistic_color.h
    include/surface_realistic_scale_widget.h
    include/surface_tile_temperature_calculator.h
    include/surface_map_mode.h
    include/surface_map_widget.h
    include/surface_map_interpolator.h
    include/surface_map_cache.h
    include/star_color.h
    include/StarSystemTopView.h
    include/surface_cell.h
    include/surface_advection_model.h
    include/surface_pressure_transport_model.h
    include/surface_height_model.h
    include/surface_heightmap.h
    include/surface_globe_widget.h
    include/surface_point_status_dialog.h
    include/surface_atmosphere_coupler.h
    include/temperature_color_scale.h
    src/atmospheric_pressure_model.h
    src/atmospheric_cell_state.cpp
    src/atmosphere_chart_widget.cpp
    src/atmosphere_model.cpp
    src/atmospheric_pressure_model.cpp
    src/atmospheric_radiation_model.cpp
    src/geothermal_flux_model.cpp
    src/layered_radiation_model.cpp
    src/radiation_model_factory.cpp
    src/radiation_model_utils.cpp
    src/radiative_convective_profile.cpp
    src/spectral_band_model.cpp
    src/atmospheric_lapse_rate_model.cpp
    src/atmospheric_thermodynamics.cpp
    src/atmosphere_widget.cpp
    src/co2_convective_adjustment.cpp
    src/emission_layer_model.cpp
    src/mollweide_projection.cpp
    src/mode_illustration_widget.cpp
    src/mode_illustration_widget.h
    src/orbit_animation_model.cpp
    src/orbit_animation_model.h
    src/orbit_segment_calculator.cpp
    src/planet_surface_grid.cpp
    src/planet_presets.cpp
    src/segment_selector_widget.cpp
    src/solar_calculator.cpp
    src/solar_display.cpp
    src/height_color_scale.cpp
    src/surface_map_widget.cpp
    src/surface_map_interpolator.cpp
    src/surface_map_cache.cpp
    src/surface_globe_widget.cpp
    src/star_color.cpp
    src/StarSystemTopView.cpp
    src/surface_advection_model.cpp
    src/surface_pressure_transport_model.cpp
    src/surface_height_model.cpp
    src/surface_heightmap.cpp
    src/surface_point_status_dialog.cpp
    src/surface_point_state.cpp
    src/surface_height_scale_widget.cpp
    src/surface_wind_scale_widget.cpp
    src/surface_pressure_scale_widget.cpp
    src/surface_atmosphere_coupler.cpp
    src/surface_realistic_color.cpp
    src/surface_realistic_scale_widget.cpp
    src/surface_tile_temperature_calculator.cpp
    src/surface_temperature_scale_widget.cpp
    src/surface_temperature_state.cpp
    src/subsurface_grid.cpp
    src/subsurface_temperature_solver.cpp
    src/thermal_conductivity_model.cpp
    src/surface_temperature_plot.cpp
    src/temperature_plot_tracker.cpp
    src/temperature_color_scale.cpp
    src/wind_field_model.cpp
    src/wind_field_model.h
)

target_include_directories(solar_constant PRIVATE include src)
target_link_libraries(solar_constant PRIVATE Qt5::Widgets Qt5::Concurrent Qt5::Charts qwt-qt5)

if (MSVC)
    target_compile_options(solar_constant PRIVATE /W4 /permissive-)
else()
    target_compile_options(solar_constant PRIVATE -Wall -Wextra -pedantic)
endif()

install(TARGETS solar_constant)
