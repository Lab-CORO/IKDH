#include <ikdh.h>
#include <cmath>
#include <cstdio>

int main()
{
    const double pi = M_PI;

    // ── ABB GoFa CRB 15000-5 — DH parameters ─────────────────────────────────
    //        Joint:        1          2          3          4          5          6
    IKDH::DHTable dh;
    dh.a    [0] = 0.000; dh.a    [1] = 0.444; dh.a    [2] = 0.110;
    dh.a    [3] = 0.000; dh.a    [4] = 0.080; dh.a    [5] = 0.000; // metres
    dh.d    [0] = 0.265; dh.d    [1] = 0.000; dh.d    [2] = 0.000;
    dh.d    [3] = 0.470; dh.d    [4] = 0.000; dh.d    [5] = 0.101; // metres
    dh.alpha[0] = -pi/2; dh.alpha[1] = 0;     dh.alpha[2] = -pi/2;
    dh.alpha[3] =  pi/2; dh.alpha[4] = -pi/2; dh.alpha[5] = 0;    // radians
    dh.theta[0] = 0;     dh.theta[1] = -pi/2; dh.theta[2] = 0;
    dh.theta[3] = 0;     dh.theta[4] = 0;     dh.theta[5] = pi;   // radians (offsets)
    for (int i = 0; i < 6; ++i) dh.revolute[i] = true;

    // ── Joint limits (degrees) ────────────────────────────────────────────────
    IKDH::JointLimits limits({
        {-180,  180},   // J1
        {-180,  180},   // J2
        {-225,   85},   // J3
        {-180,  180},   // J4
        {-180,  180},   // J5
        {-180,  180},   // J6
    });

    IKDH::Solver solver(dh, limits);

    // ── Solve IK for two end-effector poses ───────────────────────────────────
    // poseFromXYZRPW: x y z in mm, rx ry rz in degrees (RoboDK convention Rz*Ry*Rx)
    IKDH::Transform poses[] = {
        IKDH::poseFromXYZRPW(200.0, 0.0, 600.0,   0.0, 90.0, 0.0),
        IKDH::poseFromXYZRPW(400.0, 0.0, 300.0, 180.0,  0.0, 0.0),
    };

    for (const auto& ee : poses) {
        auto sols = solver.solve(ee);

        printf("%zu solution(s) found\n", sols.size());
        for (size_t i = 0; i < sols.size(); ++i) {
            printf("  [%2zu]", i);
            for (double v : sols[i]) printf("  %7.3f", v);
            printf("   FK err = %.1e\n", IKDH::fkError(ee, IKDH::forwardKin(dh, sols[i])));
        }
        printf("\n");
    }
}
