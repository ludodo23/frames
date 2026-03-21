/*
 * License: CeCILL-C
 * 
 * Copyright (c) 2026 Ludovic Andrieux
 * contributor(s): Ludovic Andrieux (2026)
 * 
 * ludovic.andrieux23@gmail.com
 * 
 * This software is a header-only C++ frames conversion library provided as a
 * single header file. [TODO].
 * 
 * This software is governed by the CeCILL-C license under French law and
 * abiding by the rules of distribution of free software. You can use,
 * modify and/or redistribute the software under the terms of the CeCILL-C
 * license as circulated by CEA, CNRS and INRIA at the following URL:
 * https://www.cecill.info
 * 
 * As a counterpart to the access to the source code and rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty and the software's author, the holder of the
 * economic rights, and the successive licensors have only limited
 * liability.
 * 
 * In this respect, the user's attention is drawn to the risks associated
 * with loading, using, modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean that it is complicated to manipulate, and that also
 * therefore means that it is reserved for developers and experienced
 * professionals having in-depth computer knowledge.
 * 
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL-C license and that you accept its terms.
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <cassert>
#include <limits>
#include <cmath>

namespace frames
{

using Transform = Eigen::Isometry3d;
using Quaternion = Eigen::Quaterniond;
using Vector3 = Eigen::Vector3d;

// ============================================================
// Utils
// ============================================================

inline Transform makeTransform(const Quaternion& q,
                               const Vector3& t)
{
    Transform T = Transform::Identity();
    T.linear() = q.toRotationMatrix();
    T.translation() = t;
    return T;
}

// ============================================================
// Sampled
// ============================================================

struct Sampled
{
    std::vector<double> t;
    std::vector<Quaternion> q;
    std::vector<Vector3> p;

    int find(double time) const
    {
        int lo = 0;
        int hi = (int)t.size() - 1;

        while (hi - lo > 1)
        {
            int mid = (lo + hi) / 2;
            if (t[mid] <= time) lo = mid;
            else hi = mid;
        }
        return lo;
    }

    Transform eval(double time) const
    {
        if (t.empty())
            return Transform::Identity();

        if (time <= t.front())
            return makeTransform(q.front(), p.front());

        if (time >= t.back())
            return makeTransform(q.back(), p.back());

        int i = find(time);

        double t0 = t[i];
        double t1 = t[i+1];
        double a = (time - t0) / (t1 - t0);

        Quaternion qr = q[i].slerp(a, q[i+1]);
        Vector3 pr = (1.0 - a) * p[i] + a * p[i+1];

        return makeTransform(qr, pr);
    }
};

// ============================================================
// Sources
// ============================================================

enum class SourceType : uint8_t
{
    Identity,
    Constant,
    Sampled,
    External,
    Analytic
};

// TODO revoir External (appel à spice par ex) et Analytic (définition user)

struct External
{
    Transform (*func)(double t, void* user);
    void* user;
};

struct Analytic
{
    const int* deps;
    int ndeps;

    Transform (*func)(
        double t,
        const Transform* world,
        const int* deps,
        int ndeps,
        void* user
    );

    void* user;
};

// ============================================================
// FrameGraph
// ============================================================

class FrameGraph
{
public:

    // TODO constructeur ?
    int addRoot()
    {
        int id = size();

        parent.push_back(-1);
        type.push_back(SourceType::Identity);
        index.push_back(0);

        return id;
    }

    int addConstant(int p, const Transform& T)
    {
        int id = size();
        assert(p < id);

        parent.push_back(p);
        type.push_back(SourceType::Constant);
        index.push_back(constants.size());

        constants.push_back(T);

        return id;
    }

    int addSampled(int p, const Sampled& s)
    {
        int id = size();
        assert(p < id);

        parent.push_back(p);
        type.push_back(SourceType::Sampled);
        index.push_back(sampled.size());

        sampled.push_back(s);

        return id;
    }

    // int addExternal(int p, External e)
    // {
    //     int id = size();
    //     assert(p < id);

    //     parent.push_back(p);
    //     type.push_back(SourceType::External);
    //     index.push_back(external.size());

    //     external.push_back(e);

    //     return id;
    // }

    // int addAnalytic(int p,
    //                 const int* deps,
    //                 int ndeps,
    //                 Transform (*func)(double,
    //                                   const Transform*,
    //                                   const int*,
    //                                   int,
    //                                   void*),
    //                 void* user = nullptr)
    // {
    //     int id = size();
    //     assert(p < id);

    //     for (int i = 0; i < ndeps; ++i)
    //         assert(deps[i] < id);

    //     parent.push_back(p);
    //     type.push_back(SourceType::Analytic);
    //     index.push_back(analytic.size());

    //     analytic.push_back({deps, ndeps, func, user});

    //     return id;
    }

    // ========================================================
    // Update
    // ========================================================

    void update(double t)
    {
        if (!dirty && t == last_time)
            return;

        last_time = t;
        dirty = false;

        size_t N = parent.size();
        world.resize(N);

        for (size_t i = 0; i < N; ++i)
        {
            Transform local;

            switch (type[i])
            {
                case SourceType::Identity:
                    local = Transform::Identity();
                    break;

                case SourceType::Constant:
                    local = constants[index[i]];
                    break;

                case SourceType::Sampled:
                    local = sampled[index[i]].eval(t);
                    break;

                // case SourceType::External:
                // {
                //     auto& e = external[index[i]];
                //     local = e.func(t, e.user);
                //     break;
                // }

                // case SourceType::Analytic:
                // {
                //     auto& a = analytic[index[i]];
                //     local = a.func(t, world.data(),
                //                    a.deps, a.ndeps, a.user);
                //     break;
                // }
            }

            int p = parent[i];

            if (p < 0) {
                world[i] = local;
            } else {
                // T_world_i = T_world_parent ∘ T_parent_i
                world[i] = world[p] * local;
            }
        }
    }

    const Transform& to_world_from(int id) const {
        return world[id];
    }

    Transform transform(int to, int from) {
        // T_to_from = T_to_world ∘ T_world_from
        return world[to].inverse() * world[from]
    }

    Vector3 position(int frame, int reference_frame, int expression_frame=-1) {
        // TODO
    }

    Quaternon attitude(int frame, int reference_frame) {
        // TODO
    }

    int size() const { return (int)parent.size(); }

private:

    std::vector<int> parent;
    std::vector<SourceType> type;
    std::vector<uint32_t> index;

    std::vector<Transform> world;

    std::vector<Transform> constants;
    std::vector<Sampled> sampled;
    std::vector<External> external;
    std::vector<Analytic> analytic;

    double last_time = std::numeric_limits<double>::quiet_NaN();
    bool dirty = true;
};

}