#pragma once

#include <ikdh.h>
#include <robots.h>

#include <cstdio>

// Solve IK for `ee` and print all solutions. Compatible with robodk_verify.py
// (solution lines contain "FK err = <val>").
static void printIKResults(const Robots::Robot& robot, const IKDH::Transform& ee)
{
    IKDH::Solver solver(robot.dh, robot.limits);
    auto sols = solver.solve(ee);

    printf("%zu solution(s)\n", sols.size());
    for (size_t i = 0; i < sols.size(); ++i) {
        IKDH::Transform check = IKDH::forwardKin(robot.dh, sols[i]);
        double err = IKDH::fkError(ee, check);
        printf("  [%2zu]", i);
        for (double v : sols[i]) printf("  %7.3f", v);
        printf("   FK err = %.1e\n", err);
    }
    printf("\n");
}

// Solve IK from a joint configuration.
// Header: "-> <label|robot_name>" + reference joints subtitle.
static void runFromJoints(const Robots::Robot& robot,
                          const IKDH::JointConfig& q_ref,
                          const char* label = nullptr)
{
    printf("-> %s\n", label ? label : robot.name.c_str());
    printf("  Joints:");
    for (double v : q_ref) printf("  %.2f", v);
    printf(" deg\n");
    printIKResults(robot, IKDH::forwardKin(robot.dh, q_ref));
}

// Solve IK from a Cartesian pose (mm, degrees, RoboDK convention).
// Header: "-> <label|robot_name>" + pose subtitle.
static void runFromPose(const Robots::Robot& robot,
                        double x_mm, double y_mm, double z_mm,
                        double rx_deg, double ry_deg, double rz_deg,
                        const char* label = nullptr)
{
    printf("-> %s\n", label ? label : robot.name.c_str());
    printf("  Pose: (%.3f, %.3f, %.3f) mm  Rx=%.1f  Ry=%.1f  Rz=%.1f deg\n",
           x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg);
    printIKResults(robot, IKDH::poseFromXYZRPW(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg));
}
