cmake_minimum_required(VERSION 3.15)
project(leo_orbit_propagator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_compile_options(/O2 /W4)
else()
    add_compile_options(-O3 -Wall -Wextra -march=native)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(third_party/pybind11)

add_library(leo_propagator_core STATIC
    src/tle_parser.cpp
    src/coordinates.cpp
    src/gravity.cpp
    src/integrator.cpp
    src/propagator.cpp
)

pybind11_add_module(_leo_propagator bindings/wrap_propagator.cpp)
target_link_libraries(_leo_propagator PRIVATE leo_propagator_core)

set_target_properties(_leo_propagator PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/python/leo_propagator"
)
