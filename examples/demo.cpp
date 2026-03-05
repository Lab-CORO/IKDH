#include <ikdh.h>
#include "robots.h"

#include <cmath>
#include <cstdio>

static double fkError(const IKDH::Transform& A, const IKDH::Transform& B)
{
    double err = 0;
    for (int k = 0; k < 16; ++k) { double d = A[k] - B[k]; err += d*d; }
    return std::sqrt(err);
}

// Build a 4x4 transform from a RoboDK-style pose:
//   x, y, z  in millimetres  →  converted to metres internally
//   rx, ry, rz  in degrees   →  intrinsic Rz*Ry*Rx rotation (RoboDK convention)
static IKDH::Transform poseFromXYZRPW(double x_mm, double y_mm, double z_mm,
                                       double rx_deg, double ry_deg, double rz_deg)
{
    const double deg = M_PI / 180.0;
    double rx = rx_deg * deg, ry = ry_deg * deg, rz = rz_deg * deg;

    double cx = std::cos(rx), sx = std::sin(rx);
    double cy = std::cos(ry), sy = std::sin(ry);
    double cz = std::cos(rz), sz = std::sin(rz);

    // R = Rz * Ry * Rx  (intrinsic ZYX — RoboDK convention)
    IKDH::Transform T;
    T[ 0] = cy*cz;              T[ 1] = cz*sx*sy - cx*sz;   T[ 2] = cx*cz*sy + sx*sz;   T[ 3] = x_mm * 1e-3;
    T[ 4] = cy*sz;              T[ 5] = cx*cz + sx*sy*sz;   T[ 6] = cx*sy*sz - cz*sx;   T[ 7] = y_mm * 1e-3;
    T[ 8] = -sy;                T[ 9] = cy*sx;               T[10] = cx*cy;               T[11] = z_mm * 1e-3;
    T[12] = 0;                  T[13] = 0;                   T[14] = 0;                   T[15] = 1;
    return T;
}

static void run(const char* name,
                const IKDH::DHTable&     dh,
                const IKDH::JointLimits& limits,
                const IKDH::JointConfig& q_ref)
{
    printf("-> %s\n", name);

    IKDH::Solver    solver(dh, limits);
    IKDH::Transform ee   = IKDH::forwardKin(dh, q_ref);
    auto            sols = solver.solve(ee);

    printf("%zu solution(s)\n", sols.size());
    for (size_t i = 0; i < sols.size(); ++i) {
        IKDH::Transform check = IKDH::forwardKin(dh, sols[i]);
        double err = fkError(ee, check);
        printf("  [%2zu]", i);
        for (double v : sols[i]) printf("  %8.3f", v);
        printf("   FK err = %.1e\n", err);
    }
    printf("\n");
}

static void run_from_pose(const char* name,
                           const IKDH::DHTable&     dh,
                           const IKDH::JointLimits& limits,
                           double x_mm, double y_mm, double z_mm,
                           double rx_deg, double ry_deg, double rz_deg)
{
    printf("-> %s\n", name);

    IKDH::Solver    solver(dh, limits);
    IKDH::Transform ee   = poseFromXYZRPW(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg);
    auto            sols = solver.solve(ee);

    printf("%zu solution(s)\n", sols.size());
    for (size_t i = 0; i < sols.size(); ++i) {
        IKDH::Transform check = IKDH::forwardKin(dh, sols[i]);
        double err = fkError(ee, check);
        printf("  [%2zu]", i);
        for (double v : sols[i]) printf("  %8.3f", v);
        printf("   FK err = %.1e\n", err);
    }
    printf("\n");
}

int main()
{
    // ── UR5e ─────────────────────────────────────────────────────────────────
    run("UR5e",
        Robots::ur5e_dh(),
        Robots::ur5e_limits(),
        {5.73, -68.75, 85.94, -103.13, -85.94, 28.65});

    // ── ABB GoFa CRB 15000-5 ─────────────────────────────────────────────────
    run("ABB GoFa CRB 15000-5",
        Robots::gofa5_dh(),
        Robots::gofa5_limits(),
        {0, 45.0, 0, 0, 45.0, 0});
    
    // ── ABB GoFa CRB 15000-5 — from TCP pose (RoboDK format) ─────────────────
    // Pose: x=571 mm, y=0, z=899 mm, rx=0°, ry=90°, rz=0°
    run_from_pose("ABB GoFa CRB 15000-5 (pose)",
        Robots::gofa5_dh(),
        Robots::gofa5_limits(),
        571.0, 0.0, 899.0,   // mm
        0.0, 90.0, 0.0);     // degrees (rx, ry, rz)

    // ── ABB GoFa CRB 15000 10 ────────────────────────────────────────────────
    run("ABB CRB 15000 10",
        Robots::abb_crb_15000_10_dh(),
        Robots::abb_crb_15000_10_limits(),
        {0, 45.0, 0, 0, 45.0, M_PI});

    // ── ABB GoFa CRB 15000 12 ────────────────────────────────────────────────
    run("ABB CRB 15000 12",
        Robots::abb_crb_15000_gofa_12_dh(),
        Robots::abb_crb_15000_gofa_12_limits(),
        {0, 45.0, 0, 0, 45.0, M_PI});
    
    run_from_pose("ABB GoFa CRB 15000-12 (pose)",
        Robots::abb_crb_15000_gofa_12_dh(),
        Robots::abb_crb_15000_gofa_12_limits(),
        700.0, 0.0, 400.0,   // mm
        0.0, 90.0, 0.0);     // degrees (rx, ry, rz)

    // ── FANUC CRX-10iA ───────────────────────────────────────────────────────
    run("FANUC CRX-10iA",
        Robots::fanuc_crx_10ia_dh(),
        Robots::fanuc_crx_10ia_limits(),
        {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});

    // ── FANUC CRX-5iA ────────────────────────────────────────────────────────
    run("FANUC CRX-5iA",
        Robots::fanuc_crx_5ia_dh(),
        Robots::fanuc_crx_5ia_limits(),
        {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});
    
    run("Doosan A0509",
        Robots::doosan_robotics_a0509_white_dh(),
        Robots::doosan_robotics_a0509_white_limits(),
        {30.0, 45.0, 60.0, -20.0, 30.0, 10.0});

    return 0;
}
