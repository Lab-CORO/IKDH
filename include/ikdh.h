#pragma once

#include <array>
#include <cmath>
#include <initializer_list>
#include <utility>
#include <vector>

namespace IKDH {

// DH parameters for a 6-DOF serial robot.
// Convention: T_i = Rz(θ_i + θ_off[i]) · Tz(d_i) · Tx(a_i) · Rx(α_i)
// Units: a, d in metres; alpha, theta in radians.
struct DHTable {
    double a[6];        // link lengths (m)
    double d[6];        // link offsets (m)
    double alpha[6];    // link twists (rad)
    double theta[6];    // joint angle offsets (rad)
    bool   revolute[6]; // true = revolute, false = prismatic
};

// Per-joint motion limits in degrees (revolute) or metres (prismatic).
// Solutions outside [lo, hi] are wrapped by multiples of 360°; others are discarded.
struct JointLimits {
    double lo[6];
    double hi[6];

    JointLimits() { for (int i = 0; i < 6; ++i) { lo[i] = -180.0; hi[i] = 180.0; } }
    JointLimits(std::initializer_list<std::pair<double,double>> pairs) {
        int i = 0;
        for (auto& p : pairs) { lo[i] = p.first; hi[i++] = p.second; }
    }
};

// Row-major 4×4 homogeneous transform. Index: [row*4 + col].
using Transform   = std::array<double, 16>;

// Joint configuration: degrees (revolute) or metres (prismatic).
using JointConfig = std::array<double, 6>;

// IK solver — construct once per robot; thread-safe per instance.
class Solver {
public:
    explicit Solver(const DHTable& dh, const JointLimits& limits = JointLimits{});
    ~Solver();

    // Return all IK solutions within joint limits (up to 16 for a general 6R robot).
    // If expand_wraps is true, also includes ±360° equivalents within limits.
    std::vector<JointConfig> solve(const Transform& ee,
                                   bool expand_wraps = false) const;

    // Warm-start IK: refine seed to ee via damped Newton-Raphson, bypassing the
    // algebraic solver. Returns one solution if convergence is reached within
    // limits, otherwise empty. Recommended for path planning (O(1) per waypoint).
    std::vector<JointConfig> solveFromSeed(const Transform& ee,
                                           const JointConfig& seed,
                                           int max_iter = 100) const;

    // Coefficients (ascending) of the Sylvester resultant polynomial from the
    // last solve() call on this thread. Empty if solve() has not been called.
    std::vector<double> lastPolynomial() const;

private:
    void*       _impl;
    double      _theta_off[6];
    JointLimits _limits;
    DHTable     _dh;
};

// Forward kinematics. q in degrees (revolute) or metres (prismatic).
Transform forwardKin(const DHTable& dh, const JointConfig& q);

// Build a Transform from a RoboDK pose: x, y, z in mm; rx, ry, rz in degrees.
// Rotation convention: intrinsic Rz·Ry·Rx.
inline Transform poseFromXYZRPW(double x_mm, double y_mm, double z_mm,
                                 double rx_deg, double ry_deg, double rz_deg)
{
    const double deg = M_PI / 180.0;
    double rx = rx_deg*deg, ry = ry_deg*deg, rz = rz_deg*deg;
    double cx = std::cos(rx), sx = std::sin(rx);
    double cy = std::cos(ry), sy = std::sin(ry);
    double cz = std::cos(rz), sz = std::sin(rz);
    Transform T;
    T[ 0]=cy*cz;  T[ 1]=cz*sx*sy-cx*sz;  T[ 2]=cx*cz*sy+sx*sz;  T[ 3]=x_mm*1e-3;
    T[ 4]=cy*sz;  T[ 5]=cx*cz+sx*sy*sz;  T[ 6]=cx*sy*sz-cz*sx;  T[ 7]=y_mm*1e-3;
    T[ 8]=-sy;    T[ 9]=cy*sx;            T[10]=cx*cy;            T[11]=z_mm*1e-3;
    T[12]=0;      T[13]=0;                T[14]=0;                T[15]=1;
    return T;
}

// Σ(A_ij − B_ij)² over all 16 elements.
inline double fkError(const Transform& A, const Transform& B)
{
    double e = 0;
    for (int k = 0; k < 16; ++k) { double d = A[k]-B[k]; e += d*d; }
    return e;
}

} // namespace IKDH
