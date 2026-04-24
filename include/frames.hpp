
/*  
 * License: CeCILL-C  
 *  
 * Copyright (c) 2026 Ludovic Andrieux  
 * contributor(s): Ludovic Andrieux (2026)  
 *  
 * ludovic.andrieux23@gmail.com  
 *  
 * This software is a header-only C++ frames conversion library provided as a  
 * single header file. It provides a frame graph data structure to represent a
 * hierarchy of frames, and supports various strategies for defining the 
 * rotation and translation of each frame. The library is designed to be
 * flexible and extensible, allowing users to define custom strategies for
 * frame transformations. It also includes an Eigen backend for efficient
 * mathematical operations.
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

/**  
 * @brief The main namespace for the frames library.  
 */
namespace frames  
{  
  
// ============================================================  
// strategies definitions.  
// ============================================================  
  
// Forward  
template <typename Backend>  
class FrameGraph;  

/**  
 * @brief A strategy for constant rotation.
 * 
 * @tparam Backend The backend type that defines the Quaternion type and related operations.
 */
template <typename Backend>  
struct BConstantRotation { 
using Quat = typename Backend::Quaternion; 
    /**
     * @brief The constant rotation value.
     */
    Quat value; 
    /**
     * @brief Evaluates the rotation at a given time.
     * @param t The time at which to evaluate the rotation.
     * @param fg The frame graph.
     * @return The constant rotation value.
     */
    Quat operator()(double t, const FrameGraph<Backend> &fg) const {  
        return value;  
    }  
};  

/**
 * @brief A strategy for constant translation.
 * 
 * @tparam Backend The backend type that defines the Vector3 type and related operations.
 */
template <typename Backend>  
struct BConstantTranslation { 
using Vec3 = typename Backend::Vector3;
    /**
     * @brief The constant translation value.
     */
    Vec3 value;
    /**
     * @brief Evaluates the translation at a given time.
     * @param t The time at which to evaluate the translation.
     * @param fg The frame graph.
     * @return The constant translation value.
     */
    Vec3 operator()(double t, const FrameGraph<Backend> &fg) const {  
        return value;  
    }  
};  

/**
 * @brief A strategy for fixed rotation at a specific epoch.
 * @tparam Backend The backend type that defines the Quaternion type and related operations.
 */
template <typename Backend>
class BFixedAtEpochRotation {
public:
    using Quat = typename Backend::Quaternion;
    /**
     * @brief The epoch at which the rotation is fixed.
     */
    double epoch;
    /**
     * @brief The fixed rotation at the specified epoch.
     */
    Quat Qwi;
    /**
     * @brief The parent frame index.
     */
    int parent;

public:
    BFixedAtEpochRotation() = delete;
    /**
     * @brief Constructs a BFixedAtEpochRotation strategy.
     * @param epoch_ The epoch at which the rotation is fixed.
     * @param parent_ The parent frame index.
     * @param fg The frame graph to evaluate the rotation at the specified epoch.
     */
    BFixedAtEpochRotation(double epoch_, int parent_, const FrameGraph<Backend>& fg) : 
        epoch(epoch_),
        parent(parent_) {
        Qwi = fg.eval_rotation(epoch, parent);
    }

    /**
     * @brief Evaluates the rotation at a given time.
     * @param t The time at which to evaluate the rotation.
     * @param fg The frame graph.
     * @return The interpolated rotation value.
     */
    Quat operator()(double t, const FrameGraph<Backend>& fg) const {
        Quat Qpw = Backend::inverse_rotation(fg.get_rotation(t, parent));
        return Backend::compose_rotation(Qpw, Qwi);
    }
};

/**
 * @brief A strategy for fixed translation at a specific epoch.
 * @tparam Backend The backend type that defines the Vector3 and Quaternion types and related operations
 */
template <typename Backend>
class BFixedAtEpochTranslation {
private:
    using Vec3 = typename Backend::Vector3;
    using Quat = typename Backend::Quaternion;

    /**
     * @brief The epoch at which the translation is fixed.
     */
    double epoch;
    /**
     * @brief The fixed translation at the specified epoch.
     */
    Vec3 Pos; // OIw
    /**
     * @brief The parent frame index.
     */
    int parent;

public:
    BFixedAtEpochTranslation() = delete;
    /**
     * @brief Constructs a BFixedAtEpochTranslation strategy.
     * @param epoch_ The epoch at which the translation is fixed.
     * @param parent_ The parent frame index.
     * @param fg The frame graph to evaluate the translation at the specified epoch.
     */
    BFixedAtEpochTranslation(double epoch_, int parent_, const FrameGraph<Backend>& fg) : 
        epoch(epoch_),
        parent(parent_) {
        Pos = fg.eval_translation(epoch, parent); // OIw = OPw(epoch)
    }

    /**
     * @brief Evaluates the translation at a given time.
     * @param t The time at which to evaluate the translation.
     * @param fg The frame graph.
     * @return The interpolated translation value.
     */
    Vec3 operator()(double t, const FrameGraph<Backend>& fg) const {

        Quat Qwp = fg.get_rotation(t, parent);
        Vec3 OPw = fg.get_position(t, parent);
        // PIp = Qwp * (OIw - OPw)
        return Backend::rotate(Qwp, Pos - OPw);
    }

};

// ============================================================
// Sampled strategies
// ============================================================

/**
 * @brief A structure to hold sampled data.
 * @tparam T The type of the sampled data.
 */
template <typename T>  
struct SampledData { 
    /** @brief The time points. */
    std::vector<double> t;  
    /** @brief The sampled values. */
    std::vector<T> value;  
    /** @brief The derivatives at the time points. */
    std::vector<T> derivative;  
  
private:  
    SampledData() = delete;  
  
public:
    /** 
     * @brief Constructs a SampledData structure with time points and values.
     * @param t_ The time points.
     * @param value_ The sampled values at the time points.
     */
    SampledData(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_  
    ) : SampledData(t_, value_, {}) {}  
  
    /** 
     * @brief Constructs a SampledData structure with time points, values, and derivatives.
     * @param t_ The time points.
     * @param value_ The sampled values at the time points.
     * @param derivative_ The derivatives at the time points.
     */
    SampledData(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_,  
        const std::vector<T> & derivative_  
    ) : t(t_), value(value_), derivative(derivative_) {}  
  
};  

/** 
 * @brief A strategy for sampled rotation.
 * @tparam Backend The backend type that defines the Quaternion type and related operations.
 */
template <typename Backend>  
struct BSampledRotation {  
private:  
    using T = typename Backend::Quaternion;  
    /** @brief The sampled data for rotation. */
    SampledData<T> _data;  
    /** @brief The interval search object for finding the appropriate time interval. */
    std::unique_ptr<interpolation::IntervalSearch> _search;  
public:
    /** 
     * @brief Constructs a BSampledRotation strategy with time points and values.
     * @param t_ The time points.
     * @param value_ The sampled rotation values at the time points.
     */
    BSampledRotation(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_  
    ) :   
        _data(SampledData(t_, value_)) {  
        _search = std::make_unique<interpolation::LinearCachedIntervalSearch>(  
            std::make_shared<const std::vector<double>>(_data.t)  
        );  
    }  
  
    /** 
     * @brief Evaluates the rotation at a given time.
     * @param time The time at which to evaluate the rotation.
     * @param fg The frame graph.
     * @return The interpolated rotation value.
     */
    T operator()(double time, const FrameGraph<Backend> &fg) const {  
        int i = _search->find(time);  
        double alpha = (time - _data.t[i]) / (_data.t[i+1] - _data.t[i]);  
        return Backend::slerp(_data.value[i], _data.value[i + 1], alpha);  
    }  
};  

/** 
 * @brief A strategy for sampled translation.
 * @tparam Backend The backend type that defines the Vector3 type and related operations.
 */
template <typename Backend>  
struct BSampledTranslation {  
private:  
    using T = typename Backend::Vector3;
    /** @brief The sampled data for translation. */
    SampledData<T> _data;  
    /** @brief The interpolator for evaluating the translation at a given time. */
    std::unique_ptr<interpolation::Interpolator<T>> _interp;  
  
public:  
    /** 
    * @brief Constructs a BSampledTranslation strategy with time points and values.
    * 
    * A CatmullRomInterpolator is used, which provides a smooth interpolation.
    * 
    * @param t_ The time points.
    * @param value_ The sampled translation values at the time points.
    */
    BSampledTranslation(  
        const std::vector<double> & t_,  
        const std::vector<T> & value_  
    ) : _data(SampledData(t_, value_)) {  
        _interp = std::make_unique<interpolation::CatmullRomInterpolator<T>>(  
            std::make_shared<const std::vector<double>>(_data.t),  
            std::make_shared<const std::vector<T>>(_data.value)  
        );  
    }  
  
    /** 
     * @brief Constructs a BSampledTranslation strategy with time points, values, and derivatives.
     * 
     * A CubicHermiteInterpolator is used, allowing for smoother interpolation
     * that takes into account the rate of change at the sampled points.
     * 
     * @param t_ The time points.
     * @param value_ The sampled translation values at the time points.
     * @param derivative_ The derivatives at the time points.
     */
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
    /** 
     * @brief Evaluates the translation at a given time.
     * @param time The time at which to evaluate the translation.
     * @param fg The frame graph.
     * @return The interpolated translation value.
     */
    T operator()(double time, const FrameGraph<Backend> &fg) const {  
        return _interp->eval(time);  
    }  
};  
  
  
// ============================================================  
// FrameGraph  
// ============================================================  
/** 
 * @brief A graph of frames with time-dependent transformations.
 * 
 * The FrameGraph class represents a hierarchy of frames, where each frame can have a parent and multiple children. Each frame's rotation and translation can be defined by various strategies, allowing for flexible and dynamic transformations over time. The class provides methods to add and remove frames, retrieve ancestors and children, and evaluate the position and attitude of frames at specific time points.
 * 
 * @tparam Backend The backend type that defines the Vector3 and Quaternion types and related operations.
 */
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
  
    /** 
     * @brief Constructs a new frame graph.
     */
    FrameGraph() :  
        _clean(true) {  
        _add_root();   
    }  
  
    /** 
     * @brief Destroys the frame graph.
     */
    virtual ~FrameGraph() = default;  
  
    /** 
     * @brief Adds a new frame to the graph.
     * @param p The parent frame id.
     * @param rotation The rotation strategy for the frame.
     * @param translation The translation strategy for the frame.
     * @return The id of the newly added frame.
     */
    template <typename RotationType, typename TranslationType>  
    int add_frame(int p, RotationType&& rotation, TranslationType&& translation) {  
        int id;
        bool reuse = false;
        if (!_free_list.empty()) {
            if(_free_list.back() > p) {
                reuse = true;
            }
        }
        // std::cout << "reuse : " << reuse << std::endl;
        if (reuse) {  
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

            // std::cout << "id = " << id << "\n";
 
            _add_rotation<RotationType>(std::move(rotation)); 
            // std::cout << "rotation added." << std::endl; 
            _add_translation<TranslationType>(std::move(translation));
            // std::cout << "translation added." << std::endl;
            _id.push_back(id);
            // std::cout << "id added = " << _id.back() << std::endl;
            _parent.push_back(p);
            // std::cout << "parent added = " << _parent.back() << std::endl;
            _children.emplace_back();  
            // std::cout << "children added" << std::endl;  
            _alive.push_back(1); 
            // std::cout << "alive added = " << _alive.back() << std::endl;

        } 
        if (p >= 0) {
            _children[p].push_back(id);
        }
        _clean = false;
        return id;
    }  

    /** 
     * @brief Gets the parent frame id.
     * @param id The frame id.
     * @return The parent frame id.
     */
    int get_parent(int id) const {
        return _parent[id];
    }

    /** 
     * @brief Removes a frame from the graph.
     * @param id The frame id of the frame to remove.
     */
    void remove_frame(int id) {
        _remove_subtree(id);
    }
    /** 
     * @brief Gets the ancestors of a frame.
     * 
     * This method returns a vector containing the ids of the ancestors of
     * the specified frame, starting from the root (id 0) down to the frame
     * itself. The ancestors are ordered from the root to the frame,
     * allowing for easy traversal of the hierarchy. The method operates 
     * in O(h) time complexity, where h is the height of the tree from
     * the specified frame to the root.
     * 
     * @param id The frame id.
     * @return A vector containing the ids of the ancestors.
     */
    std::vector<int> get_ancestors(int id) const {  
        std::vector<int> ancestors;  

        while (id != -1) {  
            ancestors.push_back(id);  
            id = _parent[id];  
        }  

        std::reverse(ancestors.begin(), ancestors.end());  
        return ancestors;  
    }  
  
    /**
     * @brief Gets the children of a frame.
     * @param id The frame id.
     * @return A reference to the vector of child frame ids.
     */
    const std::vector<int>& get_children(int id) const {  
        return _children[id];  
    }  

    /** 
     * @brief Updates the frame graph.
     * 
     * This method updates the frame graph at a given time. If the graph is 
     * already updated at the specified time and is clean, it will return 
     * immediately without performing any updates. Otherwise, it will compute
     * the new world rotations and positions for all alive frames based on
     * their respective rotation and translation strategies. After updating,
     * it marks the graph as clean, indicating that it is entirely up-to-date
     * with the given time. The method operates in O(n) time complexity, where 
     * n is the number of frames in the graph.
     * 
     * @param t The time at which to update.
     */
    void update(double t) {  
        if (is_updated(t) && _clean) {  
            return;  
        }  
        _last_time = t;  
        _update(t);
        _clean = true;  
    }

    /** 
     * @brief Evaluates the rotation of a frame at a given time.
     * 
     * This method computes the rotation of a frame at a specified time by
     * recursively composing the rotations of its ancestors.
     * It is important to note that this method does not perform any caching,
     * so it may be inefficient if called multiple times for the same time
     * and frame. For optimal performance, it is recommended to call the
     * `update` method before evaluating rotations, which will cache the
     * results and allow for O(1) retrieval of rotations using the
     * `get_rotation` or `attitude`methods.
     * 
     * @param time The time at which to evaluate.
     * @param id The frame id of the frame for which to evaluate the rotation.
     * @return The rotation quaternion.
     */
    Quaternion eval_rotation(double time, int id) const {
        if (id == 0) {
            return B::quat_identity();
        } else {
            return B::compose_rotation(this->eval_rotation(time, _parent[id]), _rot_fn[id].eval(time, *this));
        }
    }

    /** 
     * @brief Evaluates the translation of a frame at a given time.
     * 
     * This method computes the translation of a frame at a specified time by
     * recursively composing the translations of its ancestors.
     * It is important to note that this method does not perform any caching,
     * so it may be inefficient if called multiple times for the same time
     * and frame. For optimal performance, it is recommended to call the
     * `update` method before evaluating translations, which will cache the
     * results and allow for O(1) retrieval of translations using the
     * `get_position` or `position` methods.
     * 
     * @param time The time at which to evaluate.
     * @param id The frame id of the frame for which to evaluate the translation.
     * @return The translation vector.
     */
    Vector3 eval_translation(double time, int id) const {
        if (id == 0) {
            return B::vec_zero();
        } else {
            Quaternion Qwp = this->eval_rotation(time, _parent[id]); 
            Vector3 PIp =  _pos_fn[id].eval(time, *this);
            return this->eval_translation(time, _parent[id]) + B::rotate(Qwp, PIp);
        }
    }
  
    /** 
     * @brief Gets the rotation of a frame at a given time.
     * 
     * This method retrieves the rotation of a frame at a specified time from the cache
     * if the graph is updated at that time. Otherwise, it evaluates the rotation
     * using the `eval_rotation` method.
     * 
     * @param t The time at which to get the rotation.
     * @param a The frame id of the frame for which to get the rotation.
     * @return The rotation quaternion.
     */
    Quaternion get_rotation(double t, int a) const {
        // return Qwa
        if (is_updated(t)) {
            return attitude(a, 0);
        } else {
            return eval_rotation(t, a);
        }
    }

    /** 
     * @brief Gets the translation of a frame at a given time.
     * 
     * This method retrieves the translation of a frame at a specified time from the cache
     * if the graph is updated at that time. Otherwise, it evaluates the translation
     * using the `eval_translation` method.
     * 
     * @param t The time at which to get the translation.
     * @param a The frame id of the frame for which to get the translation.
     * @return The translation vector.
     */
    Vector3 get_position(double t, int a) const {
        // return OAw
        if (is_updated(t)) {
            return position(a, 0);
        } else {
            return eval_translation(t, a);
        }
    }

    /**
     * @brief get position of frame a w.r.t frame b.
     * 
     * The FrameGraph should be updated at the time of interest before calling
     * this method to ensure that the world positions and rotations are up-to-date.
     * 
     * @param a frame id to get the position.
     * @param b frame id w.r.t. get the position.
     * @return position vector of frame a in frame b.
     */
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
  
    /**
     * @brief get position of frame a w.r.t frame b projected in frame c.
     * 
     * The FrameGraph should be updated at the time of interest before calling
     * this method to ensure that the world positions and rotations are up-to-date.
     *  
     * @param a frame id to get the position.  
     * @param b frame id w.r.t. get the position  
     * @param c frame id for projection.  
     * @return position of frame a in frame b projected on frame c.  
     */
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
  
    /**
     * @brief get attitude of frame a w.r.t. frame b. 
     * 
     * The FrameGraph should be updated at the time of interest before calling
     * this method to ensure that the world rotations are up-to-date.
     * 
     * @param a frame to get the attitude.  
     * @param b frame w.r.t. get the attitute.  
     * @return attitude quaternion of a w.r.t. b.  
     */
    Quaternion attitude(int a, int b) const {  
        // Tba = T_b_world T_world_a  
        const Quaternion & Qwa (_world_rotation[a]);  
        const Quaternion & Qwb (_world_rotation[b]);  
        return B::compose_rotation(  
            B::inverse_rotation(Qwb),  
            Qwa  
        );  
    }  
  
    /**
     * @brief get the size of the frame graph.
     * @return the number of frames in the graph.
     */
    int size() const {  
        return (int)_world_rotation.size();  
    }  

    /**
     * @brief check if a frame is alive (not removed).
     * @param id frame id to check.
     * @return true if the frame is alive, false otherwise.
     */
    bool is_alive(int id) const {  
        return _alive[id];  
    }

    /**
     * @brief check if the frame graph is updated at a given time.
     * @param time the time at which to check.
     * @return true if the frame graph is updated at the given time, false otherwise.
     */
    bool is_updated(double time) const {
        return time == _last_time;
    }
  
private:  
  
    void _update(double t) {
        for (size_t id = 1; id < size(); ++id) {
            if (_alive[id]) {
                Quaternion q = _rot_fn[id].eval(t, *this);  
                Vector3 p = _pos_fn[id].eval(t, *this);

                int p_id = _parent[id];
                // std::cout << " frame " << id << " (parent = " << p_id << ") --> alive" << std::endl;

    
                // Q_world_i = Q_world_parent * Q_parent_i  
                _world_rotation[id] = B::compose_rotation(_world_rotation[p_id], q);  
                // OI_world = OP_world + R_world_parent * PI_parent  
                _world_position[id] = _world_position[p_id] + B::rotate(_world_rotation[p_id], p);
                // std::cout << "PIp\n";
                // std::cout << p;
                // std::cout << "\nOPw\n";
                // std::cout << _world_position[p_id];
                // std::cout << "\nQwp\n";
                // std::cout << _world_rotation[p_id];
                // std::cout << " --> OIw\n";
                // std::cout << _world_position[id];
                // std::cout << "\n---------------------------------\n";
            }
            else {
                // std::cout << id << "/" << size() << " --> not alive" << std::endl;

            }
        }
    }  
      
    int _add_root() {  
        int id = size();
        // std::cout << "root id = " << id << std::endl;
        _children.emplace_back();
        add_frame(-1, BConstantRotation<Backend>{B::quat_identity()}, BConstantTranslation<Backend>{B::vec_zero()});
        // std::cout << "root added." << std::endl;
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
  
    /** @brief World rotations of the frames. */
    std::vector<Quaternion> _world_rotation; 
    /** @brief World positions of the frames. */
    std::vector<Vector3> _world_position;
    /** @brief IDs of the frames. */
    std::vector<int> _id;
    /** @brief Parent IDs of the frames. */
    std::vector<int> _parent;
    /** @brief Child IDs of the frames. */
    std::vector<std::vector<int>> _children;
    /** @brief Rotation functions for the frames. */
    std::vector<Interface<Quaternion>> _rot_fn;  
    /** @brief Translation functions for the frames. */
    std::vector<Interface<Vector3>> _pos_fn;
    /** @brief Alive status of the frames. */
    std::vector<int> _alive; // 1 = alive, 0 = dead  
    /** @brief List of free IDs for recycling. */
    std::vector<int> _free_list; // recycled ids    
  
    /** @brief Last update time. */
    double _last_time = std::numeric_limits<double>::quiet_NaN();  
    /**
     * @brief Clean flag.
     * @details a frame add results in a clean flag being set to false.
     */
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
