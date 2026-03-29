#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "frames.hpp"

using namespace frames;

constexpr double PI = 3.14159265358979323846;

// ============================================================
// Helpers
// ============================================================

static bool isApprox(const Vector3& a, const Vector3& b, double eps = 1e-9)
{
    return (a - b).norm() < eps;
}

static bool isApprox(const Quaternion& q1, const Quaternion& q2, double eps = 1e-9)
{
    return std::abs(q1.dot(q2)) > (1.0 - eps);
}

// ============================================================
// Tests de base
// ============================================================

TEST_CASE("Root frame is identity", "[basic]")
{
    FrameGraph g;

    REQUIRE(g.size() == 1);

    g.update(0.0);

    const Transform& T = g.to_world_from(0);

    REQUIRE(isApprox(T.translation(), Vector3::Zero()));
    REQUIRE(isApprox(Quaternion(T.linear()), Quaternion::Identity()));
}

// ============================================================
// Constant
// ============================================================

TEST_CASE("Constant frame", "[constant]")
{
    FrameGraph g;

    Quaternion q = Quaternion::Identity();
    Vector3 p(1.0, 2.0, 3.0);

    int f = g.add_frame(0,
        ConstantRotation{q},
        ConstantTranslation{p});

    g.update(0.0);

    const Transform& T = g.to_world_from(f);

    REQUIRE(isApprox(T.translation(), p));
    REQUIRE(isApprox(Quaternion(T.linear()), q));
}

// ============================================================
// Chaîne de transformations
// ============================================================

TEST_CASE("Chain composition", "[composition]")
{
    FrameGraph g;

    int f1 = g.add_frame(0,
        ConstantRotation{Quaternion::Identity()},
        ConstantTranslation{Vector3(1,0,0)});

    int f2 = g.add_frame(f1,
        ConstantRotation{Quaternion::Identity()},
        ConstantTranslation{Vector3(0,1,0)});

    g.update(0.0);

    Vector3 p = g.position(f2, 0);

    REQUIRE(isApprox(p, Vector3(1,1,0)));
}

// ============================================================
// Sampled interpolation
// ============================================================

TEST_CASE("Sampled interpolation", "[sampled]")
{
    FrameGraph g;

    SampledTranslation st;
    st.t = {0.0, 10.0};
    st.value = {Vector3(0,0,0), Vector3(10,0,0)};

    int f = g.add_frame(0,
        ConstantRotation{Quaternion::Identity()},
        st);

    g.update(5.0);

    Vector3 p = g.position(f, 0);

    REQUIRE(isApprox(p, Vector3(5,0,0)));
}

// ============================================================
// Rotation interpolation
// ============================================================

TEST_CASE("Quaternion slerp", "[rotation]")
{
    FrameGraph g;

    SampledRotation sr;
    sr.t = {0.0, 1.0};
    sr.value = {
        Quaternion::Identity(),
        Quaternion(Eigen::AngleAxisd(PI, Vector3::UnitZ()))
    };

    int f = g.add_frame(0, sr, ConstantTranslation{Vector3::Zero()});

    g.update(0.5);

    Quaternion q = g.attitude(f, 0);

    REQUIRE(std::abs(q.angularDistance(
        Quaternion(Eigen::AngleAxisd(PI/2, Vector3::UnitZ()))
    )) < 1e-6);
}

// ============================================================
// Snapshot (FixedAtEpoch via node)
// ============================================================

TEST_CASE("FixedAtEpoch via snapshot", "[snapshot]")
{
    FrameGraph g;

    // R1 bouge
    SampledTranslation st;
    st.t = {0.0, 10.0};
    st.value = {Vector3(0,0,0), Vector3(10,0,0)};

    int R1 = g.add_frame(0,
        ConstantRotation{Quaternion::Identity()},
        st);

    // snapshot à t0 = 5
    double t0 = 5.0;
    g.update(t0);

    Transform Twp = g.to_world_from(R1);

    int frozen = g.add_frame(
        0,
        ConstantRotation{Quaternion(Twp.linear())},
        ConstantTranslation{Twp.translation()}
    );

    // R2 fixé au snapshot
    int R2 = g.add_frame(frozen,
        ConstantRotation{Quaternion::Identity()},
        ConstantTranslation{Vector3(1,0,0)});

    // maintenant on bouge le temps
    g.update(8.0);

    // R2 doit rester constant dans world
    Vector3 p = g.position(R2, 0);

    REQUIRE(isApprox(p, Vector3(5 + 1, 0, 0)));
}

// ============================================================
// Transform cohérence
// ============================================================

TEST_CASE("Transform consistency", "[transform]")
{
    FrameGraph g;

    int f1 = g.add_frame(0,
        ConstantRotation{Quaternion::Identity()},
        ConstantTranslation{Vector3(1,0,0)});

    int f2 = g.add_frame(0,
        ConstantRotation{Quaternion::Identity()},
        ConstantTranslation{Vector3(0,1,0)});

    g.update(0.0);

    Transform T = g.transform(f1, f2);

    Vector3 p = T.translation();

    REQUIRE(isApprox(p, Vector3(-1,1,0)));
}