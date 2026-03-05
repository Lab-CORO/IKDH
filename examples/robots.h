#pragma once

#include <ikdh.h>
#include <cmath>

// Predefined DH tables and joint limits for common 6-DOF robots.
// Units: a, d in meters; alpha, theta in radians; joint limits in degrees.

namespace Robots {

// ── Universal Robots UR5e ─────────────────────────────────────────────────────
inline IKDH::DHTable ur5e_dh()
{
    IKDH::DHTable dh;
    const double a[]     = { 0.0,    -0.425,  -0.3922,  0.0,     0.0,     0.0    };
    const double d[]     = { 0.1625,  0.0,     0.0,     0.1333,  0.0997,  0.0996 };
    const double alpha[] = { M_PI/2,  0.0,     0.0,     M_PI/2, -M_PI/2,  0.0   };
    const double theta[] = { 0.0,     0.0,     0.0,     0.0,     0.0,     0.0   };
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits ur5e_limits()
{
    return IKDH::JointLimits({
        {-360.0, 360.0},   // J1
        {-360.0, 360.0},   // J2
        {-360.0, 360.0},   // J3
        {-360.0, 360.0},   // J4
        {-360.0, 360.0},   // J5
        {-360.0, 360.0},   // J6
    });
}

// ── ABB GoFa CRB 15000-5 ──────────────────────────────────────────────────────
inline IKDH::DHTable gofa5_dh()
{
    IKDH::DHTable dh;
    const double a[]     = { 0.0,    0.444,  0.110,  0.0,    0.080,  0.0  };
    const double d[]     = { 0.265,  0.0,    0.0,    0.470,  0.0,    0.101};
    const double alpha[] = {-M_PI/2, 0.0,   -M_PI/2, M_PI/2,-M_PI/2, 0.0 };
    const double theta[] = { 0.0,   -M_PI/2, 0.0,    0.0,    0.0,    M_PI};
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits gofa5_limits()
{
    return IKDH::JointLimits({
        {-180.0,  180.0},  // J1
        {-180.0,  180.0},  // J2
        {-225.0,   85.0},  // J3
        {-180.0,  180.0},  // J4
        {-180.0,  180.0},  // J5
        {-180.0,  180.0},  // J6
    });
}

// ── ABB CRB 15000 10 ──────────────────────────────────────────────────────────
inline IKDH::DHTable abb_crb_15000_10_dh()
{
    IKDH::DHTable dh;
    const double a[]     = { 0.15 , 0.707, 0.11 , 0.0  , 0.08 , 0.0   };
    const double d[]     = { 0.4  , 0.0  , 0.0  , 0.637, 0.0  , 0.101 };
    const double alpha[] = { -M_PI/2, 0.0    , -M_PI/2, M_PI/2 , -M_PI/2, 0.0     };
    const double theta[] = { 0.0    , -M_PI/2, 0.0    , 0.0    , 0.0    , M_PI    };
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits abb_crb_15000_10_limits()
{
    return IKDH::JointLimits({
        { -270.0,   270.0},  // J1
        { -180.0,   180.0},  // J2
        { -225.0,    85.0},  // J3
        { -180.0,   180.0},  // J4
        { -180.0,   180.0},  // J5
        { -270.0,   270.0}   // J6
    });
}

inline IKDH::DHTable abb_crb_15000_gofa_12_dh()
{
    IKDH::DHTable dh;
    const double a[]     = { 0.0  , 0.707, 0.11 , 0.0  , 0.08 , 0.0   };
    const double d[]     = { 0.338, 0.0  , -0.0 , 0.534, 0.0  , 0.101 };
    const double alpha[] = { -M_PI/2, 0.0    , -M_PI/2, M_PI/2 , -M_PI/2, 0.0     };
    const double theta[] = { 0.0    , -M_PI/2, 0.0    , 0.0    , 0.0    , M_PI    };
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits abb_crb_15000_gofa_12_limits()
{
    return IKDH::JointLimits({
        { -270.0,   270.0},  // J1
        { -180.0,   180.0},  // J2
        { -225.0,    85.0},  // J3
        { -180.0,   180.0},  // J4
        { -180.0,   180.0},  // J5
        { -270.0,   270.0}   // J6
    });
}

// ── Fanuc CRX-5iA ─────────────────────────────────────────────────────────────
inline IKDH::DHTable fanuc_crx_5ia_dh()
{
    IKDH::DHTable dh;
    const double a[]     = { 0.0 , 0.41, 0.0 , 0.0 , 0.0 , 0.0  };
    const double d[]     = { 0.185, 0.0  , 0.0  , 0.43 , -0.13, 0.145 };
    const double alpha[] = { -M_PI/2, 0.0    , -M_PI/2, M_PI/2 , -M_PI/2, 0.0     };
    const double theta[] = { 0.0    , -M_PI/2, 0.0    , 0.0    , 0.0    , 0.0     };
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits fanuc_crx_5ia_limits()
{
    return IKDH::JointLimits({
        { -200.0,   200.0},  // J1
        { -180.0,   180.0},  // J2
        { -310.0,   310.0},  // J3
        { -190.0,   190.0},  // J4
        { -180.0,   180.0},  // J5
        { -225.0,   225.0}   // J6
    });
}

// ── Fanuc CRX-10iA ────────────────────────────────────────────────────────────
inline IKDH::DHTable fanuc_crx_10ia_dh()
{
    IKDH::DHTable dh;
    // Converted from Modified DH (RoboDK) to Standard DH used by this solver.
    // Rule: alpha^s[i] = alpha^m[i+1], a^s[i] = a^m[i+1], theta/d unchanged.
    const double a[]     = { 0.0,      0.540,    0.0,      0.0,      0.0,      0.0     };
    const double d[]     = { 0.245,    0.0,      0.0,      0.540,   -0.150,    0.160   };
    const double alpha[] = {-M_PI/2,   0.0,     -M_PI/2,  M_PI/2, -M_PI/2,   0.0     };
    const double theta[] = { 0.0,     -M_PI/2,   0.0,      0.0,      0.0,      0.0     };
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits fanuc_crx_10ia_limits()
{
    return IKDH::JointLimits({
        { -180.0,   180.0},  // J1
        { -180.0,   180.0},  // J2
        { -360.0,   430.0},  // J3
        { -190.0,   190.0},  // J4
        { -180.0,   180.0},  // J5
        { -190.0,   190.0}   // J6
    });
}

// ── Doosan Robotics A0509 White ───────────────────────────────────────────────
inline IKDH::DHTable doosan_robotics_a0509_white_dh()
{
    IKDH::DHTable dh;
    const double a[]     = { 0.0  , 0.409, 0.0  , 0.0  , 0.0  , 0.0   };
    const double d[]     = { 0.155, 0.0  , -0.0 , 0.367, 0.0  , 0.124 };
    const double alpha[] = { -M_PI/2, 0.0    , -M_PI/2, M_PI/2 , -M_PI/2, 0.0     };
    const double theta[] = { 0.0    , -M_PI/2, -M_PI/2, 0.0    , 0.0    , M_PI/2  };
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i];
        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    return dh;
}

inline IKDH::JointLimits doosan_robotics_a0509_white_limits()
{
    return IKDH::JointLimits({
        { -360.0,   360.0},  // J1
        { -360.0,   360.0},  // J2
        { -150.0,   150.0},  // J3
        { -360.0,   360.0},  // J4
        { -360.0,   360.0},  // J5
        { -360.0,   360.0}   // J6
    });
}

} // namespace Robots
