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
using Matrix3 = Eigen::Matrix3d;



// Axes
enum class Axis { X, Y, Z };

// Tags
struct Intrinsic {};
struct Extrinsic {};

// ============================================================
// Rotation élémentaire (compile-time)
// ============================================================

template<Axis A>
inline Matrix3 rot(double a);

template<>
inline Matrix3 rot<Axis::X>(double a)
{
    double c = std::cos(a), s = std::sin(a);
    Matrix3 R;
    R << 1,0,0,
         0,c,-s,
         0,s,c;
    return R;
}

template<>
inline Matrix3 rot<Axis::Y>(double a)
{
    double c = std::cos(a), s = std::sin(a);
    Matrix3 R;
    R <<  c,0,s,
          0,1,0,
         -s,0,c;
    return R;
}

template<>
inline Matrix3 rot<Axis::Z>(double a)
{
    double c = std::cos(a), s = std::sin(a);
    Matrix3 R;
    R << c,-s,0,
         s, c,0,
         0, 0,1;
    return R;
}

// ============================================================
// Euler → Quaternion (24 conventions)
// ============================================================

template<Axis A1, Axis A2, Axis A3, typename Mode>
inline Quaternion eulerToQuaternion(double a1, double a2, double a3)
{
    if constexpr (std::is_same_v<Mode, Intrinsic>)
    {
        Matrix3 R =
            rot<A1>(a1) *
            rot<A2>(a2) *
            rot<A3>(a3);

        return Quaternion(R);
    }
    else // Extrinsic
    {
        Matrix3 R =
            rot<A3>(a3) *
            rot<A2>(a2) *
            rot<A1>(a1);

        return Quaternion(R);
    }
}

// ============================================================
// makeTransform overloads
// ============================================================

// Quaternion
inline Transform makeTransform(const Quaternion& q,
                               const Vector3& t)
{
    Transform T = Transform::Identity();
    T.linear() = q.toRotationMatrix();
    T.translation() = t;
    return T;
}

// Rotation matrix (i -> world)
inline Transform makeTransform(const Matrix3& R,
                               const Vector3& t)
{
    Transform T = Transform::Identity();
    T.linear() = R;
    T.translation() = t;
    return T;
}

// Euler direct
template<Axis A1, Axis A2, Axis A3, typename Mode>
inline Transform makeTransform(double a1, double a2, double a3,
                               const Vector3& t)
{
    Quaternion q = eulerToQuaternion<A1,A2,A3,Mode>(a1,a2,a3);
    return makeTransform(q, t);
}

// Helpers lisibles
template<Axis A1, Axis A2, Axis A3>
inline Transform makeIntrinsic(double a1, double a2, double a3,
                               const Vector3& t)
{
    return makeTransform<A1,A2,A3,Intrinsic>(a1,a2,a3,t);
}

template<Axis A1, Axis A2, Axis A3>
inline Transform makeExtrinsic(double a1, double a2, double a3,
                               const Vector3& t)
{
    return makeTransform<A1,A2,A3,Extrinsic>(a1,a2,a3,t);
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
        if (t.empty()) {
            return Transform::Identity();
        }

        if (time <= t.front()) {
            return makeTransform(q.front(), p.front());
        }

        if (time >= t.back()) {
            return makeTransform(q.back(), p.back());
        }

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
    FixedAtEpoch  // TODO
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

    /// @brief get position of frame a w.r.t frame b, eventually projected in frame c.
    /// @param a frame id to get the position.
    /// @param b frame id w.r.t. get the position
    /// @param c frame id for projection. with -1, b is used instead.
    /// @return position of frame a in frame b (projected on frame c).
    Vector3 position(int a, int b, int c=-1) {
        // Tba = T_b_world T_world_a
        const Transform Twa & world[a];
        const Transform Twb & world[b];
        Vector3 BAb = (Twb.inverse() * Twa]).translation();
        if (c == -1) {
            return BAb;
        } else {
            const Transform & Twc (world[c]);
            return (Twc.linear().transpose() * Twb.linear()) * BAb; // BAc
        }
    }

    /// @brief get attitude of frame a w.r.t. frame b.
    /// @param a frame to get the attitude.
    /// @param b frame w.r.t. get the attitute.
    /// @return attitude quaterniin of a w.r.t. b.
    Quaternion attitude(int a, int b) {
        // Tba = T_b_world T_world_a
        const Transform Twa & world[a];
        const Transform Twb & world[b];
        return Quaternion((Twb.inverse() * Twa).linear())
    }

    int size() const { return (int)parent.size(); }

private:

    std::vector<int> parent;
    std::vector<SourceType> type;
    std::vector<uint32_t> index;

    std::vector<Transform> world;

    std::vector<Transform> constants;
    std::vector<Sampled> sampled;

    double last_time = std::numeric_limits<double>::quiet_NaN();
    bool dirty = true;
};

}