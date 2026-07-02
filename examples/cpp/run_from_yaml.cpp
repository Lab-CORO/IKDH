#include <ikdh.h>
#include <robots.h>
#include "example_common.h"

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

int main()
{
    // ── Load robot from YAML ──────────────────────────────────────────────────
    auto robot = Robots::loadRobot(ROBOTS_DIR "/gofa5.yaml");
    IKDH::Solver solver(robot.dh, robot.limits);

    // ── Solve IK for two end-effector poses ───────────────────────────────────
    // poseFromXYZRPW: x y z in mm, rx ry rz in degrees (RoboDK convention Rz*Ry*Rx)
    IKDH::Transform poses[] = {
        IKDH::poseFromXYZRPW(200.0, 0.0, 600.0,   0.0, 90.0, 0.0),
        IKDH::poseFromXYZRPW(400.0, 0.0, 300.0, 180.0,  0.0, 0.0),
    };

    for (const auto& ee : poses) printSolutions(solver, robot.dh, ee);
}
