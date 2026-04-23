# Getting Started — `frames.hpp`

> A header-only C++ reference frame management library by [Ludovic Andrieux](https://github.com/ludodo23).  
> Version: `0.1.0` — License: CeCILL-C  
> Depends on: [`interpolation.hpp`](https://github.com/ludodo23/interpolation)

---

## What is it?

`frames.hpp` lets you build a **tree of reference frames**, where each frame has a position and orientation relative to its parent. Frames can be:

- **static** — fixed rotation/translation
- **frozen at an epoch** — locked to a parent's state at a given time
- **sampled** — driven by time-series data, interpolated automatically

The library is backend-agnostic (you bring your own math types), with a ready-to-use **Eigen backend** included.

---

## Installation

Both libraries are header-only. Place them so they are both accessible:

```
your_project/
├── include/
│   ├── interpolation.hpp   ← required dependency
│   └── frames.hpp
```

Then include in your source:

```cpp
// Optional: enable the built-in Eigen backend
#define FRAMES_WITH_EIGEN
#include "frames.hpp"
```

> Requires **C++17** and **Eigen 3** (only if using `FRAMES_WITH_EIGEN`).

---

## Core Concepts

| Concept | Description |
|---|---|
| `FrameGraph<Backend>` | The tree of frames. Root frame (id=0) is always the world frame. |
| Frame id | An `int` returned by `add_frame()`. Use it to query positions/attitudes. |
| Rotation strategy | A callable that returns `Quaternion(t, graph)` |
| Translation strategy | A callable that returns `Vector3(t, graph)` |
| `update(t)` | Propagates the whole tree at time `t` (cached — only recomputes when needed) |

---

## 1. Setup with the Eigen Backend

```cpp
#define FRAMES_WITH_EIGEN
#include "frames.hpp"
using namespace frames;

EigenFrameGraph fg;
// Frame 0 is the world frame — already created automatically
```

---

## 2. Adding a Static Frame

```cpp
// A frame with constant orientation and position
int body = fg.add_frame(
    0,                                        // parent: world
    ConstantRotation{ Quaternion::Identity() },
    ConstantTranslation{ Vector3(1.0, 0.0, 0.0) }
);
```

---

## 3. Adding a Sampled Frame (time-series)

For a frame whose position evolves over time, provide sampled data. `frames.hpp` uses `interpolation.hpp` internally:

```cpp
std::vector<double> times = {0.0, 1.0, 2.0, 3.0};

// Sampled translation — uses CatmullRom interpolation automatically
std::vector<Vector3> positions = {
    Vector3(0, 0, 0),
    Vector3(1, 0, 0),
    Vector3(2, 1, 0),
    Vector3(3, 1, 0)
};

// Sampled rotation — uses SLERP automatically
std::vector<Quaternion> rotations = {
    Quaternion::Identity(),
    Quaternion(Eigen::AngleAxisd(0.5, Vector3::UnitZ())),
    Quaternion(Eigen::AngleAxisd(1.0, Vector3::UnitZ())),
    Quaternion(Eigen::AngleAxisd(1.5, Vector3::UnitZ()))
};

int satellite = fg.add_frame(
    0,                                      // parent: world
    SampledRotation(times, rotations),
    SampledTranslation(times, positions)
);
```

If you have known derivatives (velocities), pass them for Hermite interpolation:

```cpp
std::vector<Vector3> velocities = { ... }; // dy/dt at each node

int frame = fg.add_frame(
    0,
    SampledRotation(times, rotations),
    SampledTranslation(times, positions, velocities)  // uses CubicHermite
);
```

---

## 4. Freezing a Frame at an Epoch

A frame that follows its parent up to a given time, then stays fixed in the world:

```cpp
// Freeze the satellite's attitude at t=1.5
int frozen = fg.add_frame(
    satellite,
    FixedAtEpochRotation(1.5, satellite, fg),
    FixedAtEpochTranslation(1.5, satellite, fg)
);
```

---

## 5. Querying the Tree

Always call `update(t)` first for efficient cached evaluation:

```cpp
double t = 1.5;
fg.update(t);

// Position of `satellite` relative to world (frame 0)
Vector3 pos = fg.position(satellite, 0);

// Attitude of `satellite` relative to world
Quaternion att = fg.attitude(satellite, 0);

// Position of `satellite` w.r.t. another frame, projected in a third frame
Vector3 pos_in_c = fg.position(satellite, /*w.r.t.=*/0, /*projected in=*/frozen);
```

For one-off evaluations at arbitrary times (no caching):

```cpp
Vector3  pos = fg.eval_translation(t, satellite);
Quaternion q = fg.eval_rotation(t, satellite);
```

---

## 6. Tree Navigation

```cpp
int parent = fg.get_parent(satellite);

const std::vector<int>& children = fg.get_children(satellite);

std::vector<int> ancestors = fg.get_ancestors(satellite); // root → satellite
```

---

## 7. Removing a Frame

Removing a frame also recursively removes its entire subtree:

```cpp
fg.remove_frame(satellite); // satellite and all its children are removed
```

Removed ids are recycled for future `add_frame()` calls.

---

## 8. Euler Angles to Quaternion (Eigen backend)

The library provides a utility for all 24 Euler conventions:

```cpp
using namespace frames;

// ZYX intrinsic (aerospace convention)
Quaternion q = eulerToQuaternion<Axis::Z, Axis::Y, Axis::X, Intrinsic>(
    yaw, pitch, roll
);

// XYZ extrinsic
Quaternion q2 = eulerToQuaternion<Axis::X, Axis::Y, Axis::Z, Extrinsic>(
    a1, a2, a3
);
```

---

## 9. Custom Backend

Implement your own backend to use any math library (glm, custom types, etc.):

```cpp
struct MyBackend {
    using Vector3   = MyVec3;
    using Quaternion = MyQuat;

    static Vector3    vec_zero()                                   { return {0,0,0}; }
    static Quaternion quat_identity()                              { return MyQuat::identity(); }
    static Quaternion compose_rotation(const Quaternion& a, const Quaternion& b) { return a * b; }
    static Quaternion inverse_rotation(const Quaternion& q)        { return q.inverse(); }
    static Vector3    rotate(const Quaternion& q, const Vector3& v){ return q.rotate(v); }
    static Quaternion slerp(const Quaternion& q0, const Quaternion& q1, double t) {
        return MyQuat::slerp(q0, q1, t);
    }
};

FrameGraph<MyBackend> fg;
```

---

## Quick Reference

| Strategy | Rotation | Translation |
|---|---|---|
| Static | `ConstantRotation{ q }` | `ConstantTranslation{ v }` |
| Frozen at epoch | `FixedAtEpochRotation(t, parent, fg)` | `FixedAtEpochTranslation(t, parent, fg)` |
| Time-series | `SampledRotation(times, quats)` | `SampledTranslation(times, vecs)` |
| Time-series + derivatives | — | `SampledTranslation(times, vecs, derivs)` |

| Method | Description |
|---|---|
| `fg.update(t)` | Propagate the full tree at time `t` (cached) |
| `fg.position(a, b)` | Position of frame `a` w.r.t. frame `b` |
| `fg.position(a, b, c)` | Position of `a` w.r.t. `b`, projected in `c` |
| `fg.attitude(a, b)` | Attitude quaternion of `a` w.r.t. `b` |
| `fg.eval_translation(t, id)` | One-off position evaluation (no cache) |
| `fg.eval_rotation(t, id)` | One-off rotation evaluation (no cache) |
| `fg.get_parent(id)` | Parent frame id |
| `fg.get_children(id)` | Direct children ids |
| `fg.get_ancestors(id)` | All ancestors from root to frame |
| `fg.remove_frame(id)` | Remove frame and its subtree |
