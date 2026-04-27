find_package(Catch2 REQUIRED)
find_package(interpolation REQUIRED)

add_executable(test_frames test_frames.cpp)

target_link_libraries(test_frames PRIVATE
    frames                          # hérite de Eigen + FRAMES_WITH_EIGEN via INTERFACE
    interpolation::interpolation
    Catch2::Catch2WithMain
)

include(CTest)
include(Catch)
catch_discover_tests(test_frames)