#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

#include "frames.hpp"

#include <iostream>

using namespace frames;
using frames::EigenBackend;

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
    EigenFrameGraph g;

    REQUIRE(g.size() == 1);

    g.update(0.0);

    Quaternion Q = g.attitude(0, 0);

    CAPTURE(Q);
    REQUIRE(isApprox(Q, EigenBackend::quat_identity()));
}

// ============================================================
// Constant
// ============================================================

TEST_CASE("Constant frame", "[constant]")
{
    EigenFrameGraph g;

    Quaternion q = EigenBackend::quat_identity();
    Vector3 p(1.0, 2.0, 3.0);

    int f = g.add_frame(0,
        ConstantRotation{q},
        ConstantTranslation{p});

    g.update(5.0);

    Quaternion Q = g.attitude(f, 0);
    Vector3 P = g.position(f, 0);

    REQUIRE(isApprox(Q, EigenBackend::quat_identity()));
    REQUIRE(isApprox(P, p));
}

// ============================================================
// Chaîne de transformations
// ============================================================

TEST_CASE("Chain composition", "[composition]")
{
    EigenFrameGraph g;

    int f1 = g.add_frame(0,
        ConstantRotation{EigenBackend::quat_identity()},
        ConstantTranslation{Vector3(1,0,0)});

    int f2 = g.add_frame(f1,
        ConstantRotation{EigenBackend::quat_identity()},
        ConstantTranslation{Vector3(0,1,0)});
    
    g.update(0.0);

    Vector3 p = g.position(f2, 0);

    CAPTURE(p);
    REQUIRE(isApprox(p, Vector3(1,1,0)));

    g.update(7.0);

    REQUIRE(isApprox(p, Vector3(1,1,0)));
}

// ============================================================
// Sampled interpolation
// ============================================================

TEST_CASE("Sampled translation", "[sampled]")
{
    EigenFrameGraph g;

    SampledTranslation st({0.0, 10.0}, {Vector3(0,0,0), Vector3(10,0,0)});

    int f = g.add_frame(0,
        ConstantRotation{EigenBackend::quat_identity()},
        st);
    
    for (int i = 0; i <= 10; ++i) {
        Vector3 p = g.eval_translation(i, f);
        CAPTURE(i);
        CAPTURE(p);
        REQUIRE(isApprox(Vector3(i, 0., 0.), p));
    }

    for (int i = 0; i <= 10; ++i) {
        g.update(i);
        CAPTURE(i);
        Vector3 p = g.position(f, 0);
        REQUIRE(isApprox(Vector3(i, 0., 0.), p));
    }
}

// ============================================================
// Rotation interpolation
// ============================================================

TEST_CASE("Sampled rotation", "[sampled]")
{
    EigenFrameGraph g;

    SampledRotation sr({0.0, 1.0}, {EigenBackend::quat_identity(), Quaternion(Eigen::AngleAxisd(PI, Vector3::UnitZ()))});

    int f = g.add_frame(0, sr, ConstantTranslation{EigenBackend::vec_zero()});

    g.update(0.5);

    Quaternion q = g.attitude(f, 0);

    REQUIRE(std::abs(q.angularDistance(
        Quaternion(Eigen::AngleAxisd(PI/2, Vector3::UnitZ()))
    )) < 1e-6);
}

// ============================================================
// Snapshot (FixedAtEpoch via node)
// ============================================================

TEST_CASE("FixedAtEpochTranslation", "[fixed_at_epoch]")
{
    EigenFrameGraph g;

    // R1 bouge
    SampledTranslation st({0.0, 10.0}, {Vector3(0,0,0), Vector3(10,0,0)});

    int R1 = g.add_frame(0,
        ConstantRotation{EigenBackend::quat_identity()},
        st);

    // snapshot à t0 = 5
    double t0 = 5.0;

    g.update(2.0);

    FixedAtEpochTranslation ft(t0, 1, g);

    int R2 = g.add_frame(
        1,
        ConstantRotation{EigenBackend::quat_identity()},
        ft
    );

    for (int i = 0; i < 11; ++i) {
        g.update(i);
        Vector3 p = g.position(R2, 0);
        CAPTURE(p);
        CAPTURE(i);
        REQUIRE(isApprox(p, Vector3(5, 0, 0)));
        REQUIRE(isApprox(g.position(R2, 0), g.position(R1, 0)) == (i == 5));
    }
}

TEST_CASE("FixedAtEpochRotation", "[fixed_at_epoch]")
{
    EigenFrameGraph g;

    // R1 tourne autour de Z de 0 → PI
    SampledRotation sr(
        {0.0, 10.0},
        {
            EigenBackend::quat_identity(),
            Quaternion(Eigen::AngleAxisd(PI, Vector3::UnitZ()))
        }
    );

    int R1 = g.add_frame(
        0,
        sr,
        ConstantTranslation{EigenBackend::vec_zero()}
    );

    // snapshot à t0 = 5 → rotation = PI/2
    double t0 = 5.0;

    g.update(2.0);

    FixedAtEpochRotation fr(t0, R1, g);

    int R2 = g.add_frame(
        R1,
        fr,
        ConstantTranslation{EigenBackend::vec_zero()}
    );

    for (int i = 0; i < 11; ++i) {
        g.update(i);

        Quaternion q = g.attitude(R2, 0);

        CAPTURE(i);
        // attendu : rotation constante = PI/2 autour de Z
        Quaternion expected(
            Eigen::AngleAxisd(PI/2, Vector3::UnitZ())
        );

        REQUIRE(isApprox(q, expected));

        // égal à R1 uniquement à t0
        bool same = isApprox(g.attitude(R2, 0), g.attitude(R1, 0));

        REQUIRE(same == (i == 5));
    }
}

// ============================================================
// Transform cohérence
// ============================================================

TEST_CASE("Complex chain composition", "[compositon]")
{
    EigenFrameGraph g;

    int f1 = g.add_frame(0,
        ConstantRotation{EigenBackend::quat_identity()},
        ConstantTranslation{Vector3(1,0,0)});

    int f2 = g.add_frame(0,
        ConstantRotation{EigenBackend::quat_identity()},
        ConstantTranslation{Vector3(0,1,0)});

    g.update(9.0);

    Vector3 p = g.position(f1, f2);

    CAPTURE(p);
    REQUIRE(isApprox(p, Vector3(1,-1,0)));
}