cmake_minimum_required(VERSION 3.15)
project(IKDH VERSION 1.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# Required so static libs can be linked into shared Python extension modules
# (critical on aarch64/ARM where non-PIC code cannot be used in .so files)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(Threads REQUIRED)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Default build type (can be overridden by the user)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# MSVC does not expose M_PI by default — this enables it across all targets
if(MSVC)
    add_compile_definitions(_USE_MATH_DEFINES)
endif()

# ── HuPf IK solver (Husty–Pfurner, integrated) ───────────────────────────────
add_library(hupf STATIC
    src/hupf/Calculate.cpp
    src/hupf/ik.cpp
)
target_include_directories(hupf PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_definitions(hupf PUBLIC LIBHUPF_EXPORTS)
target_link_libraries(hupf PRIVATE Threads::Threads)

# Silence warnings from third-party code and enable native SIMD (GCC/Clang only)
target_compile_options(hupf PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
        -march=native
        -Wno-unknown-pragmas
        -Wno-deprecated-declarations
        -Wno-parentheses
        -Wno-unused-function
    >
)

# ── IKDH library ──────────────────────────────────────────────────────────────
add_library(ikdh STATIC src/ikdh.cpp)
target_include_directories(ikdh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE src
)
target_compile_options(ikdh PRIVATE
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-march=native>
)
target_link_libraries(ikdh PRIVATE hupf)

# ── Demo executable ───────────────────────────────────────────────────────────
# Run all executables from the repository root (e.g. ./build/demo),
# so that relative paths like "robots/gofa5.yaml" resolve correctly.
add_executable(demo examples/cpp/demo.cpp)
target_include_directories(demo PRIVATE include examples/cpp)
target_link_libraries(demo PRIVATE ikdh)

add_executable(run_from_DH examples/cpp/run_from_DH.cpp)
target_include_directories(run_from_DH PRIVATE include examples/cpp)
target_link_libraries(run_from_DH PRIVATE ikdh)

add_executable(run_from_yaml examples/cpp/run_from_yaml.cpp)
target_include_directories(run_from_yaml PRIVATE include examples/cpp)
target_link_libraries(run_from_yaml PRIVATE ikdh)

add_executable(derive_poly tools/derive_poly.cpp)
target_include_directories(derive_poly PRIVATE include)
target_link_libraries(derive_poly PRIVATE ikdh)

add_executable(run_benchmark tools/run_benchmark.cpp)
target_include_directories(run_benchmark PRIVATE include)
target_link_libraries(run_benchmark PRIVATE ikdh)

# ── Python bindings (optional, requires pybind11) ─────────────────────────────
option(IKDH_BUILD_PYTHON "Build Python bindings via pybind11" OFF)

if(IKDH_BUILD_PYTHON)
    set(PYBIND11_FINDPYTHON ON)  # use CMake FindPython instead of legacy FindPythonLibsNew
    include(FetchContent)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.13.6
    )
    FetchContent_MakeAvailable(pybind11)

    pybind11_add_module(ikdh_py src/ikdh_bindings.cpp)
    target_include_directories(ikdh_py PRIVATE include)
    target_link_libraries(ikdh_py PRIVATE ikdh)
    set_target_properties(ikdh_py PROPERTIES OUTPUT_NAME ikdh)
    # Suppress the MSVC import .lib — it collides with the static ikdh.lib
    # Python loads .pyd files directly so no import library is needed
    if(MSVC)
        target_link_options(ikdh_py PRIVATE /NOIMPLIB)
    endif()

    # Installed by `pip install .` (see pyproject.toml) directly into site-packages root.
    install(TARGETS ikdh_py DESTINATION .)
endif()

# ── Install rules ─────────────────────────────────────────────────────────────
# Only add install rules when this is the top-level project (not a subdirectory),
# and not when building the Python wheel (SKBUILD, set by scikit-build-core) —
# the wheel should contain only the compiled extension, not C++ dev artifacts.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND NOT SKBUILD)

    # Install libraries
    install(TARGETS ikdh hupf
        EXPORT IKDHTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

    # Install public header
    install(FILES include/ikdh.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    # Export targets file (IKDHTargets.cmake)
    install(EXPORT IKDHTargets
        FILE IKDHTargets.cmake
        NAMESPACE IKDH::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IKDH
    )

    # Generate IKDHConfig.cmake from template
    configure_package_config_file(
        cmake/IKDHConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/IKDHConfig.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IKDH
    )

    # Generate IKDHConfigVersion.cmake
    write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/IKDHConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/IKDHConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/IKDHConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/IKDH
    )

endif()
