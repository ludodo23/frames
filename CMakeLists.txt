cmake_minimum_required(VERSION 3.15)
project(frames CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -------------------------
# LIB HEADER-ONLY
# -------------------------
add_library(frames INTERFACE)
target_include_directories(frames INTERFACE include/)

if(FRAMES_WITH_EIGEN)
    find_package(Eigen3 REQUIRED)
    target_link_libraries(frames INTERFACE Eigen3::Eigen)
    target_compile_definitions(frames INTERFACE FRAMES_WITH_EIGEN)
endif()

# -------------------------
# TESTS
# -------------------------
enable_testing()
add_subdirectory(tests)