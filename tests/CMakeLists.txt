cmake_minimum_required(VERSION 3.15)

project(frames_tests)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(Catch2 REQUIRED)

add_executable(test_frames test_frames.cpp test_benchmark.cpp)

target_include_directories(test_frames PRIVATE ../include)

target_link_libraries(test_frames
    Eigen3::Eigen
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)

catch_discover_tests(test_frames)
