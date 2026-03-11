#include "ik_print.h"

#include <string>

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

static std::string robotPath(const char* filename)
{
    return std::string(ROBOTS_DIR) + "/" + filename;
}

int main()
{
    // ── UR5e ─────────────────────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("ur5e.yaml"));
        runFromJoints(robot, {5.73, -68.75, 85.94, -103.13, -85.94, 28.65});
    }

    // ── ABB GoFa CRB 15000-5 ─────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("gofa5.yaml"));
        runFromJoints(robot, {0, 45.0, 0, 0, 45.0, 0});
        runFromPose(robot, 571.0, 0.0, 899.0, 0.0, 90.0, 0.0);
    }

    // ── ABB CRB 15000-10 ─────────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("crb15000_10.yaml"));
        runFromJoints(robot, {0, 45.0, 0, 0, 45.0, 180.0});
    }

    // ── ABB GoFa CRB 15000-12 ────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("gofa12.yaml"));
        runFromJoints(robot, {0, 45.0, 0, 0, 45.0, 180.0});
        runFromPose(robot, 700.0, 0.0, 400.0, 0.0, 90.0, 0.0);
    }

    // ── FANUC CRX-10iA ───────────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("fanuc_crx_10ia.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── FANUC CRX-5iA ────────────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("fanuc_crx_5ia.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── Doosan A0509 ─────────────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("doosan_a0509.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── Auctech ACR-12 ───────────────────────────────────────────────────────
    {
        auto robot =Robots::loadRobot(robotPath("auctech_acr_12.yaml"));
        runFromJoints(robot, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    return 0;
}
