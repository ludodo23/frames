
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
  
#ifdef FRAMES_WITH_EIGEN  
#include <Eigen/Dense>  
#include <Eigen/Geometry>  
#endif  
#include <vector>  
#include <cassert>  
#include <limits>  
#include <cmath>  
#include <memory>  
#include <stdexcept>
#include <algorithm>
#include <type_traits>
  
#include <interpolation.hpp>  
  
#define FRAMES_VERSION "0.1.0"  
  
namespace frames  
{  
  
// ============================================================  
// strategies definitions.  
// ============================================================  
  
// Forward  
template <typename Backend>  
class FrameGraph;  
  
template <typename Backend>  
struct BConstantRotation { 
using Quat = typename Backend::Quaternion; 
    Quat value;  
    Quat operator()(double t, const FrameGraph<Backend> &fg) const {  
        return value;  
    }  
};  

template <typename Backend>  
struct BConstantTranslation { 
using Vec3 = typename Backend::Vector3; 
    Vec3 value;  
    Vec3 operator()(double t, const FrameGraph<Backend> &fg) const {  
        return value;  
    }  
};  

template <typename Backend>
class BFixedAtEpochRotation {
private:
    using Quat = typename Backend::Quaternion;
    double epoch;
    Quat Q_pp_n;
    int parent;
    int parent_of_parent;

    Quat get_rotation(const FrameGraph<Backend>& fg, double t, int a, int b) const {
        if (fg.is_updated(t)) {
            return fg.attitude(b, a);
        } else {
            Quat Qwa = fg.eval_rotation(epoch, a);
            Quat Qwb = fg.eval_rotation(epoch, b);
            return Backend::compose_rotation(Backend::inverse_rotation(Qwa), Qwb);
        }
    }

public:
    BFixedAtEpochRotation() = delete;
    BFixedAtEpochRotation(double epoch_, int parent_, const FrameGraph<Backend>& fg) : 
        epoch(epoch_),
        parent(parent_) {
        parent_of_parent = fg.get_parent(parent);
        Q_pp_n = get_rotation(fg, epoch, parent_of_parent, parent);
    }

    Quat operator()(double t, const FrameGraph<Backend>& fg) const {
        Quat Q_p_pp = get_rotation(fg, t, parent, parent_of_parent);
        return Backend::compose_rotation(Q_p_pp, Q_pp_n);
    }
};


template <typename Backend>
class BFixedAtEpochTranslation {
private:
    using Vec3 = typename Backend::Vector3;
    using Quat = typename Backend::Quaternion;

    double epoch;
    int parent;
    int parent_of_parent;
    Vec3 P_pp_n;

    Quat get_rotation(const FrameGraph<Backend>& fg, double t, int a, int b) const {
        if (fg.is_updated(t)) {
            return fg.attitude(b, a);
        } else {
            Quat Qwa = fg.eval_rotation(epoch, a);
            Quat Qwb = fg.eval_rotation(epoch, b);
            return Backend::compose_rotation(Backend::inverse_rotation(Qwa), Qwb);
        }
    }

    Vec3 get_position(const FrameGraph<Backend>& fg, double t, int a, int b) const {
        if (fg.is_updated(t)) {
            return fg.position(a, b);
        } else {
            Vec3 BAw = fg.eval_translation(t, a) - fg.eval_translation(t, b);
            Quat Qwb = fg.eval_rotation(t, b);
            return Backend::rotate(Qwb, BAw);
        }
    }

public:
    BFixedAtEpochTranslation() = delete;
    BFixedAtEpochTranslation(double epoch_, int parent_, const FrameGraph<Backend>& fg) : 
        epoch(epoch_),
        parent(parent_) {
        parent_of_parent = fg.get_parent(parent);
        P_pp_n = get_position(fg, epoch, parent, parent_of_parent);
    }

    Vec3 operator()(double t, const FrameGraph<Backend>& fg) const {

        Quat Q_p_pp = get_rotation(fg, t, parent, parent_of_parent);
        Vec3 P_p_pp = get_position(fg, t, parent, parent_of_parent);

        return P_p_pp + Backend::rotate(Q_p_pp, P_pp_n);
    }

};

// ============================================================
// Sampled strategies
// ============================================================
template <typename T>  
struct SampledData {  
    std::vector<double> t;  
    std::vector<T> value;  
    std::vector<T> derivative;  
  
private:  
    SampledData() = delete;  
  
public:  
    SampledData(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_  
    ) : SampledData(t_, value_, {}) {}  
  
    SampledData(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_,  
        const std::vector<T> & derivative_  
    ) : t(t_), value(value_), derivative(derivative_) {}  
  
};  
  
template <typename Backend>  
struct BSampledRotation {  
private:  
    using T = typename Backend::Quaternion;  
    SampledData<T> _data;  
    std::unique_ptr<interpolation::IntervalSearch> _search;  
public:  
    BSampledRotation(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_  
    ) :   
        _data(SampledData(t_, value_)) {  
        _search = std::make_unique<interpolation::CachedInterval>(  
            std::make_shared<const std::vector<double>>(_data.t)  
        );  
    }  
  
    T operator()(double time, const FrameGraph<Backend> &fg) const {  
        int i = _search->find(time);  
        double alpha = (time - _data.t[i]) / (_data.t[i+1] - _data.t[i]);  
        return Backend::slerp(_data.value[i], _data.value[i + 1], alpha);  
    }  
};  
  
template <typename Backend>  
struct BSampledTranslation {  
private:  
    using T = typename Backend::Vector3;  
    SampledData<T> _data;  
    std::unique_ptr<interpolation::Interpolator<T>> _interp;  
  
public:  
    BSampledTranslation(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_  
    ) : _data(SampledData(t_, value_)) {  
        _interp = std::make_unique<interpolation::CatmullRomInterpolator<T>>(  
            std::make_shared<const std::vector<double>>(_data.t),  
            std::make_shared<const std::vector<T>>(_data.value)  
        );  
    }  
  
    BSampledTranslation(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_,  
        const std::vector<T> & derivative_  
    ) : _data(SampledData(t_, value_, derivative_)) {  
       _interp = std::make_unique<interpolation::CubicHermiteInterpolator<T>>(  
            std::make_shared<const std::vector<double>>(_data.t),  
            std::make_shared<const std::vector<T>>(_data.value),  
            std::make_shared<const std::vector<T>>(_data.derivative)  
        );  
    }  
    T operator()(double time, const FrameGraph<Backend> &fg) const {  
        return _interp->eval(time);  
    }  
};  
  
  
// ============================================================  
// FrameGraph  
// ============================================================  
  
template <typename Backend>  
class FrameGraph {  
    using B = Backend;  
  
    // ============================================================  
    // wrapper for custom types.  
    // ============================================================  
    template <typename T>
    class Interface {
        struct Concept {
            virtual ~Concept() = default;
            virtual T eval(double t, const FrameGraph<Backend>& g) const = 0;
        };

        template <typename Derived>
        struct Model final : Concept {
            Derived impl;

            template <typename U>
            Model(U&& v)
                : impl(std::forward<U>(v)) {}

            T eval(double t, const FrameGraph<Backend>& g) const override {
                return impl(t, g);
            }
        };

        std::unique_ptr<Concept> self;

    public:
        Interface() = default;

        // PERFECT FORWARDING CONSTRUCTOR
        template <typename Derived>
        Interface(Derived&& v)
            : self(std::make_unique<Model<std::decay_t<Derived>>>(
                std::forward<Derived>(v))) {}

        // move only
        Interface(Interface&&) noexcept = default;
        Interface& operator=(Interface&&) noexcept = default;

        // copie deleted
        Interface(const Interface&) = delete;
        Interface& operator=(const Interface&) = delete;

        T eval(double t, const FrameGraph<Backend>& g) const {
            return self->eval(t, g);
        }
    };
  
public:  
    using Vector3 = typename B::Vector3;  
    using Quaternion = typename B::Quaternion;  
  
    FrameGraph() :  
        _clean(true) {  
        _add_root();   
    }  
  
    virtual ~FrameGraph() = default;  
  
    template <typename RotationType, typename TranslationType>  
    int add_frame(int p, RotationType&& rotation, TranslationType&& translation) {  
        int id;
        if (!_free_list.empty() && _free_list. back() > p) {  
            id = _free_list.back();
            _free_list.pop_back();  
            _set_rotation<RotationType>(id, std::move(rotation));  
            _set_translation<TranslationType>(id, std::move(translation));  
            _parent[id] = p;  
            _children[id].clear();  
            _alive[id] = 1;  
        } else { 
            assert(p < size());
            id = size();  
 
            _add_rotation<RotationType>(std::move(rotation));  
            _add_translation<TranslationType>(std::move(translation));  
            _id.push_back(id);
            _parent.push_back(p);  
            _children.emplace_back();    
            _alive.push_back(1);  
        } 
        _children[p].push_back(id);
        _clean = false;
        return id;
    }  

    int get_parent(int id) const {
        return _parent[id];
    }

    void remove_frame(int id) {
        _remove_subtree(id);
    }
    // ----------------------------------------------------  
    // ANCESTORS (O(h))  
    // includes root (0)  
    // ----------------------------------------------------  
    std::vector<int> get_ancestors(int id) const {  
        std::vector<int> ancestors;  

        while (id != -1) {  
            ancestors.push_back(id);  
            id = _parent[id];  
        }  

        std::reverse(ancestors.begin(), ancestors.end());  
        return ancestors;  
    }  
  
    // ----------------------------------------------------  
    // CHILDREN (O(1))  
    // ----------------------------------------------------  
    const std::vector<int>& get_children(int id) const {  
        return _children[id];  
    }  

    void update(double t) {  
        if (is_updated(_last_time) && _clean) {  
            return;  
        }  
        _last_time = t;  
        _update(t, 1, _id);  
        _clean = true;  
    }  

    void update(double t, int id) {
        if (is_updated(t)) {  
            return;  
        }  
        _update(t, 1, get_ancestors(id)); 
        _clean = false;
    }

    Quaternion eval_rotation(double time, int id) const {
        if (id == 0) {
            return B::quat_identity();
        } else {
            return B::compose_rotation(this->eval_rotation(time, _parent[id]), _rot_fn[id].eval(time, *this));
        }
    }

    Vector3 eval_translation(double time, int id) const {
        if (id == 0) {
            return B::vec_zero();
        } else {
            Quaternion Qwp = this->eval_rotation(time, _parent[id]); 
            Vector3 PIp =  _pos_fn[id].eval(time, *this);
            return this->eval_translation(time, _parent[id]) + B::rotate(Qwp, PIp);
        }
    }
  
    Vector3 position(int a, int b) const {   
        // Tba = T_b_world T_world_a  
        const Vector3 & OAw (_world_position[a]);  
        const Vector3 & OBw (_world_position[b]);  
        const Quaternion & Qwb (_world_rotation[b]);  
        return B::rotate(  
            B::inverse_rotation(Qwb),  
            OAw - OBw  
        );  
    }  
  
    /// @brief get position of frame a w.r.t frame b, eventually projected in frame c.  
    /// @param a frame id to get the position.  
    /// @param b frame id w.r.t. get the position  
    /// @param c frame id for projection. with -1, b is used instead.  
    /// @return position of frame a in frame b (projected on frame c).  
    Vector3 position(int a, int b, int c) const {    
        const Quaternion & Qwb (_world_rotation[b]);  
        const Quaternion & Qwc (_world_rotation[c]);  
  
        Vector3 BAb = position(a, b);  
  
        return B::rotate(  
            B::compose_rotation(  
                B::inverse_rotation(Qwb),  
                Qwc  
            ),  
            BAb  
        );  
    }  
  
    /// @brief get attitude of frame a w.r.t. frame b.  
    /// @param a frame to get the attitude.  
    /// @param b frame w.r.t. get the attitute.  
    /// @return attitude quaternion of a w.r.t. b.  
    Quaternion attitude(int a, int b) const {  
        // Tba = T_b_world T_world_a  
        const Quaternion & Qwa (_world_rotation[a]);  
        const Quaternion & Qwb (_world_rotation[b]);  
        return B::compose_rotation(  
            B::inverse_rotation(Qwb),  
            Qwa  
        );  
    }  
  
    int size() const {  
        return (int)_world_rotation.size();  
    }  

    bool is_alive(int id) const {  
        return _alive[id];  
    }

    bool is_updated(double time) const {
        return time == _last_time;
    }
  
private:  
  
    void _update(double t, int begin, const std::vector<int> & list_to_update) {  
        for (size_t i = begin; i < list_to_update.size(); ++i) { // 0 -> world = Id
            int id = list_to_update[i];
            if (_alive[id]) {
                Quaternion q = _rot_fn[id].eval(t, *this);  
                Vector3 p = _pos_fn[id].eval(t, *this);  
    
                int p_id = _parent[id];  
    
                // Q_world_i = Q_world_parent * Q_parent_i  
                _world_rotation[id] = B::compose_rotation(_world_rotation[p_id], q);  
                // OI_world = OP_world + R_world_parent * PI_parent  
                _world_position[id] = _world_position[p_id] + B::rotate(_world_rotation[p_id], p);  
            }
        }
    }  
      
    int _add_root() {  
        int id = size();  
        _add_rotation(BConstantRotation<Backend>{B::quat_identity()});  
        _add_translation(BConstantTranslation<Backend>{B::vec_zero()});  
        return id;  
    }  
  
    template <typename RotationType>  
    void _add_rotation( RotationType&& rotation) {  
            _rot_fn.emplace_back(std::move(rotation)); 
           _world_rotation.emplace_back(B::quat_identity());  
    }  
  
    template <typename TranslationType>  
    void _add_translation( TranslationType&& translation) {
            _pos_fn.emplace_back(std::move(translation));  
            _world_position.emplace_back(B::vec_zero());  
    }  

    template <typename RotationType>
    void _set_rotation(int id, RotationType&& rotation) {
        _rot_fn[id] = Interface<Quaternion>(std::move(rotation));
    }

    template <typename TranslationType>
    void _set_translation(int id, TranslationType&& translation) {
        _pos_fn[id] = Interface<Vector3>(std::move(translation));
    }

    // ----------------------------------------------------  
    // REMOVE SUBTREE (iterative DFS, O(k))  
    // ----------------------------------------------------  
    void _remove_subtree(int root) {  
        std::vector<int> stack;  
        stack.push_back(root);  

        while (!stack.empty()) {  
            int node = stack.back();  
            stack.pop_back();  

            if (!_alive[node]) continue;  

            _alive[node] = 0;  
            _free_list.push_back(node);  

            for (int child : _children[node]) {  
                stack.push_back(child);  
            }  
            _children[node].clear();  
        }  
    }  
  
    std::vector<Quaternion> _world_rotation; 
    std::vector<Vector3> _world_position;
    std::vector<int> _id;
    std::vector<int> _parent;
    std::vector<std::vector<int>> _children;
    std::vector<Interface<Quaternion>> _rot_fn;  
    std::vector<Interface<Vector3>> _pos_fn;
    std::vector<uint8_t> _alive; // 1 = alive, 0 = dead  
    std::vector<int> _free_list; // recycled ids    
  
    double _last_time = std::numeric_limits<double>::quiet_NaN();  
    bool _clean;  
            
};  
  
// ============================================================  
// Eigen Backend  
// ============================================================  
  
#ifdef FRAMES_WITH_EIGEN  
  
struct EigenBackend {  
    using Vector3 = Eigen::Vector3d;  
    using Quaternion = Eigen::Quaterniond;  
  
    static Vector3 vec_zero() { return Vector3::Zero(); }  
    static Quaternion quat_identity() { return Quaternion::Identity(); }  
  
    static Quaternion compose_rotation(const Quaternion& a, const Quaternion& b) {  
        return a * b;  
    }  
  
    static Quaternion inverse_rotation(const Quaternion& q) {  
        return q.conjugate();  
    }  
  
    static Vector3 rotate(const Quaternion& q, const Vector3& v) {  
        return q * v;  
    }  
  
    static Quaternion slerp(const Quaternion& q0,  
                            const Quaternion& q1,  
                            double t) {  
        return q0.slerp(t, q1);  
    }  
};  
  
using EigenFrameGraph = FrameGraph<EigenBackend>;  
  
typedef BConstantRotation<EigenBackend> ConstantRotation;  
typedef BFixedAtEpochRotation<EigenBackend> FixedAtEpochRotation;  
typedef BSampledRotation<EigenBackend> SampledRotation;  
  
typedef BConstantTranslation<EigenBackend> ConstantTranslation;  
typedef BFixedAtEpochTranslation<EigenBackend> FixedAtEpochTranslation;  
typedef BSampledTranslation<EigenBackend> SampledTranslation;  
  
using Matrix3 = Eigen::Matrix3d;  
using Vector3 = Eigen::Vector3d;  
using Quaternion = Eigen::Quaterniond;  
  
// Axes  
enum class Axis {X, Y, Z};  
  
// Tags  
struct Intrinsic {};  
struct Extrinsic {};  
  
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
  
#endif  
  
}
