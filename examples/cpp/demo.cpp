#include <ikdh.h>
#include <robots.h>
#include "example_common.h"

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

int main()
{
    // Forward kinematics: joint angles (deg) -> end-effector pose.
    auto robot = Robots::loadRobot(ROBOTS_DIR "/fanuc_crx_10ia.yaml");
    IKDH::Solver solver(robot.dh, robot.limits);

    IKDH::JointConfig q = {30.0, 45.0, 60.0, -20.0, 30.0, 10.0};
    IKDH::Transform ee = IKDH::forwardKin(robot.dh, q);

    printf("-> %s (from joint angles)\n", robot.name.c_str());
    printSolutions(solver, robot.dh, ee);

    // Inverse kinematics: end-effector pose (mm, deg) -> all joint solutions.
    auto robot2 = Robots::loadRobot(ROBOTS_DIR "/gofa5.yaml");
    IKDH::Solver solver2(robot2.dh, robot2.limits);

    IKDH::Transform target = IKDH::poseFromXYZRPW(571.0, 0.0, 899.0, 0.0, 90.0, 0.0);

    printf("-> %s (from a Cartesian pose)\n", robot2.name.c_str());
    printSolutions(solver2, robot2.dh, target);

    return 0;
}
