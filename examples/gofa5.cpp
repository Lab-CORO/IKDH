#include "ik_print.h"

#include <string>

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

int main()
{
    auto robot = Robots::loadRobot(std::string(ROBOTS_DIR) + "/gofa5.yaml");

    // find the 14 solutions here
    runFromPose(robot, 200.0, 0, 600.0, 0.0, 90.0, 0.0);
    // find the 10 solutions here
    runFromPose(robot, 400.0, 0.0, 300.0, 180.0, 0.0, 0.0);

    return 0;
}
