#include <iostream>
#include <chrono>
#include <random>
#include "frame_graph.hpp"

using Clock = std::chrono::high_resolution_clock;

fg::Sampled makeSampled(int n)
{
    fg::Sampled s;
    s.t.resize(n);
    s.q.resize(n);
    s.p.resize(n);

    for (int i = 0; i < n; ++i)
    {
        s.t[i] = i * 10.0;

        s.q[i] = Eigen::Quaterniond::UnitRandom();
        s.p[i] = Eigen::Vector3d::Random() * 1e6;
    }
    return s;
}

int main()
{
    fg::FrameGraph g;

    const int N = 1000;      // nombre de repères
    const int STEPS = 1000;  // frames simulées

    int root = g.addRoot();

    std::vector<int> frames;
    frames.push_back(root);

    // ============================
    // Construction du graphe
    // ============================

    for (int i = 1; i < N; ++i)
    {
        int parent = frames[rand() % frames.size()];

        if (i % 3 == 0)
        {
            auto s = makeSampled(50);
            frames.push_back(g.addSampled(parent, s));
        }
        else if (i % 3 == 1)
        {
            Eigen::Isometry3d T = Eigen::Isometry3d::Identity();
            T.translation() = Eigen::Vector3d::Random();

            frames.push_back(g.addConstant(parent, T));
        }
        else
        {
            int dep = parent;
            int deps[] = {dep};

            frames.push_back(
                g.addAnalytic(parent, deps, 1, fg::lvlh_func)
            );
        }
    }

    std::cout << "Frames: " << g.size() << std::endl;

    // ============================
    // Benchmark update
    // ============================

    auto t0 = Clock::now();

    double time = 0.0;

    for (int i = 0; i < STEPS; ++i)
    {
        g.update(time);

        // lecture (important pour éviter l'optimisation)
        volatile double acc = 0.0;

        for (int j = 0; j < N; ++j)
        {
            acc += g.get(j).translation().x();
        }

        time += 1.0;
    }

    auto t1 = Clock::now();

    double dt = std::chrono::duration<double>(t1 - t0).count();

    std::cout << "Total time: " << dt << " s\n";
    std::cout << "Per frame: " << dt / STEPS * 1e6 << " us\n";
    std::cout << "Per frame per node: "
              << (dt / STEPS / N) * 1e9 << " ns\n";

    return 0;
}