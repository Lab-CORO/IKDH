#include <ikdh.h>
#include <robots.h>
#include "example_common.h"

#include <cstdio>
#include <string>

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

static std::string robotPath(const char* filename)
{
    return std::string(ROBOTS_DIR) + "/" + filename;
}

// Solve IK for `ee` and print all solutions. Compatible with robodk_verify.py
// (solution lines contain "FK err = <val>").
static void printIKResults(const Robots::Robot& robot, const IKDH::Transform& ee,
                           bool expand_wraps = false)
{
    IKDH::Solver solver(robot.dh, robot.limits);
    printSolutions(solver, robot.dh, ee, expand_wraps);
}

// Solve IK from a joint configuration.
// Header: "-> <label|robot_name>" + reference joints subtitle.
static void runFromJoints(const Robots::Robot& robot,
                          const IKDH::JointConfig& q_ref,
                          const char* label = nullptr,
                          bool expand_wraps = false)
{
    printf("-> %s\n", label ? label : robot.name.c_str());
    printf("  Joints:");
    for (double v : q_ref) printf("  %.2f", v);
    printf(" deg\n");
    printIKResults(robot, IKDH::forwardKin(robot.dh, q_ref), expand_wraps);
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

int main()
{
    // UR5e
    {
        auto robot = Robots::loadRobot(robotPath("ur5e.yaml"));
        runFromJoints(robot, {5.73, -68.75, 85.94, -103.13, -85.94, 28.65});
    }

    // ABB GoFa CRB 15000-5
    {
        auto robot = Robots::loadRobot(robotPath("gofa5.yaml"));
        runFromJoints(robot, {0, 45.0, 0, 0, 45.0, 0});
        runFromPose(robot, 571.0, 0.0, 899.0, 0.0, 90.0, 0.0);
    }

    // ABB CRB 15000-10
    {
        auto robot = Robots::loadRobot(robotPath("crb15000_10.yaml"));
        runFromJoints(robot, {0, 45.0, 0, 0, 45.0, 180.0});
    }

    // ABB GoFa CRB 15000-12
    {
        auto robot = Robots::loadRobot(robotPath("gofa12.yaml"));
        runFromJoints(robot, {0, 45.0, 0, 0, 45.0, 180.0});
        runFromPose(robot, 700.0, 0.0, 400.0, 0.0, 90.0, 0.0);
    }

    // FANUC CRX-10iA
    {
        auto robot = Robots::loadRobot(robotPath("fanuc_crx_10ia.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // FANUC CRX-5iA
    {
        auto robot = Robots::loadRobot(robotPath("fanuc_crx_5ia.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // Doosan A0509
    {
        auto robot = Robots::loadRobot(robotPath("doosan_a0509.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // Auctech ACR-12
    {
        auto robot = Robots::loadRobot(robotPath("auctech_acr_12.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // expand_wraps comparison (robots with limits wider than ±180°)
    // For each solution, allWraps() generates every k*360° shift that stays
    // within the joint limits. These are distinct robot configurations that
    // reach the same end-effector pose  -  useful for motion planners.
    printf("=== expand_wraps=true comparison ===\n\n");

    {
        auto robot = Robots::loadRobot(robotPath("doosan_a0509.yaml"));
        printf("--- Doosan A0509 (limits ±360° on all joints) ---\n");
        printf("Without wrap expansion:\n");
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0}, nullptr, false);
        printf("With wrap expansion:\n");
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0}, nullptr, true);
    }

    {
        auto robot = Robots::loadRobot(robotPath("gofa12.yaml"));
        printf("--- ABB GoFa 12 (J3 limit -225 deg) ---\n");
        printf("Without wrap expansion:\n");
        runFromJoints(robot, {0.0, 45.0, 0.0, 0.0, 45.0, 180.0}, nullptr, false);
        printf("With wrap expansion:\n");
        runFromJoints(robot, {0.0, 45.0, 0.0, 0.0, 45.0, 180.0}, nullptr, true);
    }

    return 0;
}
