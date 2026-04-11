#include <catch2/catch_all.hpp>

#include <iostream>
#include <chrono>
#include <random>
#include <vector>

#include "frames.hpp"

using namespace frames;


void bench(int n_frames, int n_samples, int n_eval) {
std::cout << n_frames << ", " << n_samples << ", " << n_eval << ", ";

    // ============================
    // Génération des données
    // ============================
    std::vector<double> t(n_samples);
    std::vector<Vector3> pos(n_samples);
    std::vector<Quaternion> rot(n_samples);

    for (int i = 0; i < n_samples; ++i) {
        double ti = i * 0.01;
        t[i] = ti;

        pos[i] = Vector3(
            std::sin(ti),
            std::cos(ti),
            std::sin(2 * ti)
        );

        rot[i] = Quaternion(
            Eigen::AngleAxisd(ti, Vector3::UnitZ()) *
            Eigen::AngleAxisd(0.5 * ti, Vector3::UnitY())
        );
    }

    // ============================
    // Construction FrameGraph
    // ============================
    FrameGraph fg;

    int parent = 0;

    for (int i = 0; i < n_frames; ++i) {
        SampledTranslation trans(t, pos);
        SampledRotation rota(t, rot);

        parent = fg.add_frame(parent, rota, trans);
    }

    // ============================
    // Benchmark
    // ============================
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n_eval; ++i) {
        double time = (i % n_samples) * 0.01;

        fg.update(time);

        // accès pour éviter optimisation
        volatile auto T = fg.to_world_from(n_frames - 1);
        // (void)T;
    }

    auto end = std::chrono::high_resolution_clock::now();

    double elapsed = std::chrono::duration<double, std::milli>(end - start).count();

    std::cout  << elapsed << ", " << elapsed / n_eval << std::endl;
}

TEST_CASE("Benchmark sampled number of frame", "[benchmark]")
{
    constexpr int N_SAMPLES = 1000;
    constexpr int N_EVAL = 1000;
    std::vector<int> N {1000, 10000, 50000, 100000};

    for (int i : N) {
            bench(i, N_SAMPLES, N_EVAL);
    }
}

TEST_CASE("Benchmark sampled number of sample", "[benchmark]")
{
    constexpr int N_FRAMES = 1000;
    constexpr int N_EVAL = 1000;
    std::vector<int> N {1000, 10000, 50000, 100000};

    for (int i : N) {
            bench(N_FRAMES, i, N_EVAL);
    }
}

TEST_CASE("Benchmark sampled number of eval", "[benchmark]")
{
    constexpr int N_FRAMES = 1000;
    constexpr int N_SAMPLES = 1000;
    std::vector<int> N {1000, 10000, 50000, 100000};

    for (int i : N) {
            bench(N_FRAMES, N_SAMPLES, i);
    }
}