cmake_minimum_required(VERSION 3.14)
project(IKDH CXX)

set(CMAKE_CXX_STANDARD 11)

# Default build type (can be overridden by the user)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ── HuPf IK solver (Husty–Pfurner, integrated) ───────────────────────────────
add_library(hupf STATIC
    src/hupf/Calculate.cpp
    src/hupf/ik.cpp
)
target_include_directories(hupf PUBLIC src)
target_compile_definitions(hupf PRIVATE LIBHUPF_EXPORTS)

# Silence warnings from third-party code (GCC/Clang only)
target_compile_options(hupf PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
        -Wno-unknown-pragmas
        -Wno-deprecated-declarations
        -Wno-parentheses
        -Wno-unused-function
    >
)

# ── IKDH library ──────────────────────────────────────────────────────────────
add_library(ikdh STATIC src/ikdh.cpp)
target_include_directories(ikdh
    PUBLIC  include
    PRIVATE src
)
target_link_libraries(ikdh PRIVATE hupf)

# ── Demo executable ───────────────────────────────────────────────────────────
# ROBOTS_DIR defaults to "robots" (relative to the working directory).
# Run executables from the repository root: ./build/demo
add_executable(demo examples/demo.cpp)
target_include_directories(demo PRIVATE include)
target_link_libraries(demo PRIVATE ikdh)

add_executable(poly examples/poly.cpp)
target_include_directories(poly PRIVATE include)
target_link_libraries(poly PRIVATE ikdh)

add_executable(benchmark examples/benchmark.cpp)
target_include_directories(benchmark PRIVATE include)
target_link_libraries(benchmark PRIVATE ikdh)

