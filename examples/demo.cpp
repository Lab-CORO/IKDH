#include <ikdh.h>
#include <robots.h>

#include <cmath>
#include <cstdio>
#include <string>

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

static double fkError(const IKDH::Transform& A, const IKDH::Transform& B)
{
    double err = 0;
    for (int k = 0; k < 16; ++k) { double d = A[k] - B[k]; err += d*d; }
    return std::sqrt(err);
}


static void run(const Robots::Robot& robot, const IKDH::JointConfig& q_ref)
{
    printf("-> %s\n", robot.name.c_str());

    IKDH::Solver    solver(robot.dh, robot.limits);
    IKDH::Transform ee   = IKDH::forwardKin(robot.dh, q_ref);
    auto            sols = solver.solve(ee);

    printf("%zu solution(s)\n", sols.size());
    for (size_t i = 0; i < sols.size(); ++i) {
        IKDH::Transform check = IKDH::forwardKin(robot.dh, sols[i]);
        double err = fkError(ee, check);
        printf("  [%2zu]", i);
        for (double v : sols[i]) printf("  %8.3f", v);
        printf("   FK err = %.1e\n", err);
    }
    printf("\n");
}

static void run_from_pose(const Robots::Robot& robot,
                           double x_mm, double y_mm, double z_mm,
                           double rx_deg, double ry_deg, double rz_deg)
{
    printf("-> %s\n", robot.name.c_str());

    IKDH::Solver    solver(robot.dh, robot.limits);
    IKDH::Transform ee   = IKDH::poseFromXYZRPW(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg);
    auto            sols = solver.solve(ee);

    printf("%zu solution(s)\n", sols.size());
    for (size_t i = 0; i < sols.size(); ++i) {
        IKDH::Transform check = IKDH::forwardKin(robot.dh, sols[i]);
        double err = fkError(ee, check);
        printf("  [%2zu]", i);
        for (double v : sols[i]) printf("  %8.3f", v);
        printf("   FK err = %.1e\n", err);
    }
    printf("\n");
}

static std::string robotPath(const char* filename)
{
    return std::string(ROBOTS_DIR) + "/" + filename;
}

int main()
{
    // ── UR5e ─────────────────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("ur5e.yaml"));
        run(r, {5.73, -68.75, 85.94, -103.13, -85.94, 28.65});
    }

    // ── ABB GoFa CRB 15000-5 ─────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("gofa5.yaml"));
        run(r, {0, 45.0, 0, 0, 45.0, 0});
        run_from_pose(r, 571.0, 0.0, 899.0, 0.0, 90.0, 0.0);
    }

    // ── ABB CRB 15000-10 ─────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("crb15000_10.yaml"));
        run(r, {0, 45.0, 0, 0, 45.0, M_PI});
    }

    // ── ABB GoFa CRB 15000-12 ────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("gofa12.yaml"));
        run(r, {0, 45.0, 0, 0, 45.0, M_PI});
        run_from_pose(r, 700.0, 0.0, 400.0, 0.0, 90.0, 0.0);
    }

    // ── FANUC CRX-10iA ───────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("fanuc_crx_10ia.yaml"));
        run(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── FANUC CRX-5iA ────────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("fanuc_crx_5ia.yaml"));
        run(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    // ── Doosan A0509 ─────────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("doosan_a0509.yaml"));
        run(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }
        // ── Doosan A0509 ─────────────────────────────────────────────────────────
    {
        auto r = Robots::loadRobot(robotPath("auctech_acr_12.yaml"));
        run(r, {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    }

    return 0;
}
