#include "ik_print.h"

#include <cmath>
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
        auto r = Robots::loadRobot(robotPath("ur5e.yaml"));
        runFromJoints(r, {5.73, -68.75, 85.94, -103.13, -85.94, 28.65});
    }

    // ── ABB GoFa CRB 15000-5 ─────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("gofa5.yaml"));
        runFromJoints(r, {0, 45.0, 0, 0, 45.0, 0});
        runFromPose(r, 571.0, 0.0, 899.0, 0.0, 90.0, 0.0);
    }

    // ── ABB CRB 15000-10 ─────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("crb15000_10.yaml"));
        runFromJoints(r, {0, 45.0, 0, 0, 45.0, M_PI});
    }

    // ── ABB GoFa CRB 15000-12 ────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("gofa12.yaml"));
        runFromJoints(r, {0, 45.0, 0, 0, 45.0, M_PI});
        runFromPose(r, 700.0, 0.0, 400.0, 0.0, 90.0, 0.0);
    }

    // ── FANUC CRX-10iA ───────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("fanuc_crx_10ia.yaml"));
        runFromJoints(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── FANUC CRX-5iA ────────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("fanuc_crx_5ia.yaml"));
        runFromJoints(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── Doosan A0509 ─────────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("doosan_a0509.yaml"));
        runFromJoints(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── Auctech ACR-12 ───────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("auctech_acr_12.yaml"));
        runFromJoints(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    return 0;
}
