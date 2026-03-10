#pragma once

#include <array>
#include <cmath>
#include <initializer_list>
#include <utility>
#include <vector>

namespace IKDH {

// Denavit-Hartenberg parameters for a 6-DOF serial robot.
// Convention: T_i = Rz(theta_i + theta_offset[i]) * Tz(d[i]) * Tx(a[i]) * Rx(alpha[i])
// Units: a, d in meters; alpha, theta in radians.
struct DHTable {
    double a[6];        // link lengths
    double d[6];        // link offsets
    double alpha[6];    // link twists
    double theta[6];    // joint angle offsets
    bool   revolute[6]; // true = revolute, false = prismatic
};

// Joint limits in degrees (revolute) or meters (prismatic).
// Angles outside [lo, hi] are wrapped by multiples of 360°.
// Solutions that cannot be mapped into the range are discarded.
struct JointLimits {
    double lo[6];
    double hi[6];

    // Default: ±180° for all joints.
    JointLimits() {
        for (int i = 0; i < 6; ++i) { lo[i] = -180.0; hi[i] = 180.0; }
    }
    JointLimits(std::initializer_list<std::pair<double,double>> pairs) {
        int i = 0;
        for (auto& p : pairs) { lo[i] = p.first; hi[i++] = p.second; }
    }
};

// 4x4 homogeneous transformation matrix, row-major.
// [ R  t ]   index: [row*4 + col]
// [ 0  1 ]
using Transform = std::array<double, 16>;

// Joint configuration: 6 values (degrees for revolute, meters for prismatic).
using JointConfig = std::array<double, 6>;

// IK solver — construct once per robot, call solve() for each end-effector pose.
class Solver {
public:
    // Construct with DH parameters and optional joint limits.
    // Returned solutions are guaranteed to lie within the specified limits.
    explicit Solver(const DHTable& dh, const JointLimits& limits = JointLimits{});
    ~Solver();

    // Solve IK for the given end-effector transform.
    // Returns all real solutions within joint limits.
    std::vector<JointConfig> solve(const Transform& ee) const;

    // Returns the coefficients of the univariate resultant polynomial computed
    // during the most recent call to solve() on this thread.
    // coefficient[k] is the coefficient of t^k (ascending order).
    //
    // This is the Sylvester resultant of the two bivariate Study-quadric
    // constraints; its degree is robot-dependent (typically 56 for a general
    // 6R robot).  The actual IK solutions correspond to the real roots of this
    // polynomial (up to 16 for a 6R robot) after removal of spurious roots.
    //
    // Returns an empty vector if solve() has not yet been called on this thread.
    std::vector<double> lastPolynomial() const;

private:
    void*       _impl;
    double      _theta_off[6];
    JointLimits _limits;
    DHTable     _dh;   // stored for Newton-Raphson IK refinement
};

// Compute forward kinematics.
// q[i] is in degrees (revolute) or meters (prismatic).
Transform forwardKin(const DHTable& dh, const JointConfig& q);

// Build a Transform from a RoboDK-style pose:
//   x, y, z   in millimetres  →  stored in metres
//   rx, ry, rz in degrees     →  intrinsic Rz*Ry*Rx rotation (RoboDK convention)
inline Transform poseFromXYZRPW(double x_mm, double y_mm, double z_mm,
                                 double rx_deg, double ry_deg, double rz_deg)
{
    const double deg = M_PI / 180.0;
    double rx = rx_deg * deg, ry = ry_deg * deg, rz = rz_deg * deg;
    double cx = std::cos(rx), sx = std::sin(rx);
    double cy = std::cos(ry), sy = std::sin(ry);
    double cz = std::cos(rz), sz = std::sin(rz);
    Transform T;
    T[ 0] = cy*cz;            T[ 1] = cz*sx*sy - cx*sz; T[ 2] = cx*cz*sy + sx*sz; T[ 3] = x_mm * 1e-3;
    T[ 4] = cy*sz;            T[ 5] = cx*cz + sx*sy*sz; T[ 6] = cx*sy*sz - cz*sx; T[ 7] = y_mm * 1e-3;
    T[ 8] = -sy;              T[ 9] = cy*sx;             T[10] = cx*cy;             T[11] = z_mm * 1e-3;
    T[12] = 0; T[13] = 0; T[14] = 0; T[15] = 1;
    return T;
}

// Sum of squared element-wise differences between two 4x4 transforms.
inline double fkError(const Transform& A, const Transform& B)
{
    double e = 0;
    for (int k = 0; k < 16; ++k) { double d = A[k] - B[k]; e += d*d; }
    return e;
}

} // namespace IKDH
