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
#include <memory>
#include <stdexcept>

#define FRAMES_VERSION "0.1.0"

namespace frames
{

using Transform = Eigen::Isometry3d;
using Quaternion = Eigen::Quaterniond;
using Vector3 = Eigen::Vector3d;
using Matrix3 = Eigen::Matrix3d;

// Axes
enum class Axis
{
    X,
    Y,
    Z
};

// Tags
struct Intrinsic
{
};
struct Extrinsic
{
};

// ============================================================
//  elementary rotation
// ============================================================

template <Axis A>
inline Matrix3 rot(double a);

template <>
inline Matrix3 rot<Axis::X>(double a)
{
    double c = std::cos(a), s = std::sin(a);
    Matrix3 R;
    R << 1, 0, 0,
        0, c, -s,
        0, s, c;
    return R;
}

template <>
inline Matrix3 rot<Axis::Y>(double a)
{
    double c = std::cos(a), s = std::sin(a);
    Matrix3 R;
    R << c, 0, s,
        0, 1, 0,
        -s, 0, c;
    return R;
}

template <>
inline Matrix3 rot<Axis::Z>(double a)
{
    double c = std::cos(a), s = std::sin(a);
    Matrix3 R;
    R << c, -s, 0,
        s, c, 0,
        0, 0, 1;
    return R;
}

// ============================================================
// Euler → Quaternion (24 conventions)
// ============================================================

template <Axis A1, Axis A2, Axis A3, typename Mode>
inline Quaternion eulerToQuaternion(double a1, double a2, double a3)
{
    if constexpr (std::is_same_v<Mode, Intrinsic>) {
        Matrix3 R =
            rot<A1>(a1) *
            rot<A2>(a2) *
            rot<A3>(a3);

        return Quaternion(R);
    } else {
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
inline Transform makeTransform(const Quaternion &q,
                                const Vector3 &t) {
    Transform T = Transform::Identity();
    T.linear() = q.toRotationMatrix();
    T.translation() = t;
    return T;
}

// Rotation matrix (i -> world)
inline Transform makeTransform(const Matrix3 &R,
                                const Vector3 &t) {
    Transform T = Transform::Identity();
    T.linear() = R;
    T.translation() = t;
    return T;
}

// Euler direct
template <Axis A1, Axis A2, Axis A3, typename Mode>
inline Transform makeTransform(double a1, double a2, double a3,
                                const Vector3 &t) {
    Quaternion q = eulerToQuaternion<A1, A2, A3, Mode>(a1, a2, a3);
    return makeTransform(q, t);
}

// Helpers lisibles
template <Axis A1, Axis A2, Axis A3>
inline Transform makeIntrinsic(double a1, double a2, double a3,
                                const Vector3 &t) {
    return makeTransform<A1, A2, A3, Intrinsic>(a1, a2, a3, t);
}

template <Axis A1, Axis A2, Axis A3>
inline Transform makeExtrinsic(double a1, double a2, double a3,
                                const Vector3 &t)
{
    return makeTransform<A1, A2, A3, Extrinsic>(a1, a2, a3, t);
}

// ============================================================
// CRTP for strategies definitions.
// ============================================================

// Forward
class FrameGraph;

template <typename T>
T interp(double a, const T &v0, const T &v1) {
    return (1.0 - a) * v0 + a * v1;
}

template <>
Quaternion interp(double a, const Quaternion &q0, const Quaternion &q1) {
    return q0.slerp(a, q1);
}

template <typename T>
struct Constant {
    T value;
    T operator()(double t, const FrameGraph &fg) const {
        return value;
    }
};

template <typename T>
struct FixedAtEpoch {
    double epoch;
    T operator()(double t, const FrameGraph &fg) const {
        throw std::runtime_error("Not implemented !");
    }
};

template <typename T>
struct Sampled {
    std::vector<double> t;
    std::vector<T> value;

    int find(double time) const
    {
        int lo = 0;
        int hi = (int)t.size() - 1;

        while (hi - lo > 1)
        {
            int mid = (lo + hi) / 2;
            if (t[mid] <= time)
                lo = mid;
            else
                hi = mid;
        }
        return lo;
    }

    T operator()(double time, const FrameGraph &fg) const {
        if (t.empty()) {
            return T{};
        }

        if (time <= t.front()) {
            return value.front();
        }

        if (time >= t.back()) {
            return value.back();
        }

        int i = find(time);

        double t0 = t[i];
        double t1 = t[i + 1];
        double a = (time - t0) / (t1 - t0);
        T value0 = value[i];
        T value1 = value[i + 1];

        return interp(a, value0, value1);
    }
};

typedef Constant<Quaternion> ConstantRotation;
typedef FixedAtEpoch<Quaternion> FixedAtEpochRotation;
typedef Sampled<Quaternion> SampledRotation;

typedef Constant<Vector3> ConstantTranslation;
typedef FixedAtEpoch<Vector3> FixedAtEpochTranslation;
typedef Sampled<Vector3> SampledTranslation;



// ============================================================
// wrapper for custom types.
// ============================================================
template <typename T>
class Interface {
    struct Concept {
        virtual ~Concept() = default;
        virtual T eval(double t, const FrameGraph& g) const = 0;
    };

    template <typename Derived>
    struct Model : Concept {
        Derived impl;

        Model(Derived v) : impl(std::move(v)) {}

        T eval(double t, const FrameGraph& g) const override {
            return impl(t, g);
        }
    };
    
    std::unique_ptr<Concept> self;

public:
    template <typename Derived>
    Interface(Derived v) : self(std::make_unique<Model<Derived>>(std::move(v))) {}

    T eval(double t, const FrameGraph& g) const {
        return self->eval(t, g);
    }
    
};

// ============================================================
// FrameGraph
// ============================================================

class FrameGraph {
public:

    FrameGraph() :
        _clean(true) {
        _add_root(); 
    }

    virtual ~FrameGraph() {
    }

    template <typename RotationType, typename TranslationType>
    int add_frame(int p, const RotationType& rotation, const TranslationType& translation) {
        assert(p < size());
        _add_rotation<RotationType>(p, rotation);
        _add_translation<TranslationType>(p, translation);
        return size() - 1;
    }

    void update(double t) {
        if (t == _last_time && _clean) {
            return;
        }

        _last_time = t;

        _update(t);

        _clean = true;
    }

    Transform to_world_from(int id) const {
        return makeTransform(_world_rotation[id], _world_position[id]);
    }

    Transform transform(int to, int from) {
        // T_to_from = T_to_world ∘ T_world_from
        return to_world_from(to).inverse() * to_world_from(from);
    }

    Vector3 position(int a, int b) { 
        // Tba = T_b_world T_world_a
        const Vector3 & OAw (_world_position[a]);
        const Vector3 & OBw (_world_position[b]);
        const Quaternion & Qwb (_world_rotation[b]);
        Vector3 BAb = Qwb.toRotationMatrix().transpose() * (OAw - OBw);
        return BAb;
    }

    /// @brief get position of frame a w.r.t frame b, eventually projected in frame c.
    /// @param a frame id to get the position.
    /// @param b frame id w.r.t. get the position
    /// @param c frame id for projection. with -1, b is used instead.
    /// @return position of frame a in frame b (projected on frame c).
    Vector3 position(int a, int b, int c) {  // TODO two overloads 
        const Quaternion & Qwb (_world_rotation[b]);
        const Quaternion & Qwc (_world_rotation[c]);

        Vector3 BAb = position(a, b);
        
        return (Qwc.conjugate() * Qwb).toRotationMatrix() * BAb; // BAc
    }

    /// @brief get attitude of frame a w.r.t. frame b.
    /// @param a frame to get the attitude.
    /// @param b frame w.r.t. get the attitute.
    /// @return attitude quaternion of a w.r.t. b.
    Quaternion attitude(int a, int b) {
        // Tba = T_b_world T_world_a
        const Quaternion & Qwa (_world_rotation[a]);
        const Quaternion & Qwb (_world_rotation[b]);
        return Qwb.conjugate() * Qwa;
    }

    int size() const {
        return (int)_world_rotation.size();
    }

private:

    void _update(double t) {
        for (size_t i = 1; i < size(); ++i) { // 0 -> world = Id
            Quaternion q = _rot_fn[i].eval(t, *this);
            Vector3 p    = _pos_fn[i].eval(t, *this);

            int pr_id = _rot_parent[i];
            int pp_id = _pos_parent[i];

            // Q_world_i = Q_world_parent * Q_parent_i
            _world_rotation[i] = _world_rotation[pr_id] * q;
            // OI_world = OP_world + R_world_parent * PI_parent
            _world_position[i] = _world_position[pp_id] + _world_rotation[pp_id].toRotationMatrix() * p;
        }

    }

    void _update_rotation(double t) {
        _clean = false;
        for (size_t i = 1; i < _rot_fn.size(); ++i) { // 0 -> world = Id
            Quaternion q = _rot_fn[i].eval(t, *this);

            int pr_id = _rot_parent[i];

            // Q_world_i = Q_world_parent * Q_parent_i
            _world_rotation[i] = _world_rotation[pr_id] * q;
    }
}

    void _update_translation(double t) {
        _clean = false;
        _update(t);
    }
    
    int _add_root() {
        int id = size();

        _add_rotation(-1, ConstantRotation{Quaternion::Identity()});
        _add_translation(-1, ConstantTranslation{Vector3::Zero()});

        return id;
    }

    template <typename RotationType>
    void _add_rotation(int parent, const RotationType& rotation) {
        if constexpr (std::is_same_v<RotationType, FixedAtEpochRotation>) {
            int parent_of_parent = _rot_parent[parent];
            _update_rotation(rotation.epoch);
            _add_rotation(parent_of_parent, ConstantRotation{attitude(parent, parent_of_parent)});
        } else {
            _rot_parent.push_back(parent);
            _rot_fn.emplace_back(rotation);
             _world_rotation.emplace_back(Quaternion::Identity());
        }
    }

    template <typename TranslationType>
   void _add_translation(int parent, const TranslationType& translation) {
        if constexpr (std::is_same_v<TranslationType, FixedAtEpochTranslation>) {
            int parent_of_parent = _pos_parent[parent];
            _update_translation(translation.epoch);
        _add_translation(parent_of_parent, ConstantTranslation{position(parent, parent_of_parent)});
        } else {
            _pos_parent.push_back(parent);
            _pos_fn.emplace_back(translation);
            _world_position.emplace_back(Vector3::Zero());
        }
    }

    std::vector<Quaternion> _world_rotation;
    std::vector<Vector3> _world_position;

    std::vector<int> _rot_parent;
    std::vector<int> _pos_parent;

    std::vector<Interface<Quaternion>> _rot_fn;
    std::vector<Interface<Vector3>> _pos_fn;

    double _last_time = std::numeric_limits<double>::quiet_NaN();
    bool _clean;
};

}