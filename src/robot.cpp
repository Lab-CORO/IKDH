#include <robot.h>
#include <robots.h>

#include <algorithm>
#include <cmath>

namespace IKDH {

namespace {

// Forward kinematics helpers, kept local to this file (mirrors ikdh.cpp's
// internal dhMatrix/mul4, but returns Transform directly for convenience).

Transform dhMatrix(double theta, double d, double a, double alpha)
{
    double ct = std::cos(theta), st = std::sin(theta);
    double ca = std::cos(alpha), sa = std::sin(alpha);
    Transform T;
    T[0]=ct;  T[1]=-st*ca;  T[2]= st*sa;  T[3]=a*ct;
    T[4]=st;  T[5]= ct*ca;  T[6]=-ct*sa;  T[7]=a*st;
    T[8]=0;   T[9]=sa;      T[10]=ca;     T[11]=d;
    T[12]=0;  T[13]=0;      T[14]=0;      T[15]=1;
    return T;
}

Transform mul4(const Transform& A, const Transform& B)
{
    Transform C;
    for (int i = 0; i < 3; ++i) {
        C[i*4+0] = A[i*4+0]*B[0] + A[i*4+1]*B[4] + A[i*4+2]*B[8];
        C[i*4+1] = A[i*4+0]*B[1] + A[i*4+1]*B[5] + A[i*4+2]*B[9];
        C[i*4+2] = A[i*4+0]*B[2] + A[i*4+1]*B[6] + A[i*4+2]*B[10];
        C[i*4+3] = A[i*4+0]*B[3] + A[i*4+1]*B[7] + A[i*4+2]*B[11] + A[i*4+3];
    }
    C[12]=0; C[13]=0; C[14]=0; C[15]=1;
    return C;
}

// Cumulative transforms T_accum[0..6]: T_accum[j] is the pose of frame j-1
// in base coordinates (T_accum[0] = identity, T_accum[6] = end-effector).
std::array<Transform, 7> computeFrames(const DHTable& dh, const JointConfig& q)
{
    std::array<Transform, 7> T_accum;
    T_accum[0] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};
    for (int j = 0; j < 6; ++j) {
        double val = dh.revolute[j] ? q[j] * (M_PI / 180.0) : q[j];
        double jt  = dh.theta[j] + (dh.revolute[j] ? val : 0.0);
        double jd  = dh.d[j]     + (dh.revolute[j] ? 0.0 : val);
        T_accum[j+1] = mul4(T_accum[j], dhMatrix(jt, jd, dh.a[j], dh.alpha[j]));
    }
    return T_accum;
}

// Geometric Jacobian: linear rows z_j x (p_e - p_j), angular rows z_j.
Matrix6 geometricJacobian(const DHTable& dh, const JointConfig& q)
{
    auto frames = computeFrames(dh, q);
    const Transform& Tcur = frames[6];
    double p_e[3] = { Tcur[3], Tcur[7], Tcur[11] };

    Matrix6 J{};
    for (int j = 0; j < 6; ++j) {
        const Transform& Tj = frames[j];
        double zj[3] = { Tj[2], Tj[6], Tj[10] };

        if (!dh.revolute[j]) {
            for (int r = 0; r < 3; ++r) J[r][j] = zj[r];
            continue;
        }

        double pj[3] = { Tj[3], Tj[7], Tj[11] };
        double rv[3] = { p_e[0]-pj[0], p_e[1]-pj[1], p_e[2]-pj[2] };
        J[0][j] = zj[1]*rv[2] - zj[2]*rv[1];
        J[1][j] = zj[2]*rv[0] - zj[0]*rv[2];
        J[2][j] = zj[0]*rv[1] - zj[1]*rv[0];
        J[3][j] = zj[0]; J[4][j] = zj[1]; J[5][j] = zj[2];
    }
    return J;
}

// Small dense 6x6 linear algebra, no external dependency  -  same spirit as
// the QR solve in ikdh.cpp, kept separate since this file's operations
// (determinant, inverse) are different from the damped least-squares solve
// used by the IK refinement loop.

double det6(const Matrix6& M)
{
    Matrix6 A = M;
    double s = 1.0;
    for (int c = 0; c < 6; ++c) {
        int p = c;
        for (int r = c+1; r < 6; ++r) if (std::fabs(A[r][c]) > std::fabs(A[p][c])) p = r;
        if (p != c) { std::swap(A[c], A[p]); s = -s; }
        if (std::fabs(A[c][c]) < 1e-14) return 0.0;
        for (int r = c+1; r < 6; ++r) {
            double f = A[r][c] / A[c][c];
            for (int k = c; k < 6; ++k) A[r][k] -= f * A[c][k];
        }
    }
    double d = s;
    for (int i = 0; i < 6; ++i) d *= A[i][i];
    return d;
}

Matrix6 transpose6(const Matrix6& M)
{
    Matrix6 T{};
    for (int r = 0; r < 6; ++r) for (int c = 0; c < 6; ++c) T[c][r] = M[r][c];
    return T;
}

Matrix6 matmul6(const Matrix6& A, const Matrix6& B)
{
    Matrix6 C{};
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c) {
            double s = 0;
            for (int k = 0; k < 6; ++k) s += A[r][k] * B[k][c];
            C[r][c] = s;
        }
    return C;
}

// Gauss-Jordan inverse with partial pivoting. Returns false (leaving out
// untouched) if M is singular to working precision.
bool invert6(const Matrix6& M, Matrix6& out)
{
    Matrix6 A = M;
    Matrix6 I{};
    for (int i = 0; i < 6; ++i) I[i][i] = 1.0;

    for (int c = 0; c < 6; ++c) {
        int p = c;
        for (int r = c+1; r < 6; ++r) if (std::fabs(A[r][c]) > std::fabs(A[p][c])) p = r;
        if (std::fabs(A[p][c]) < 1e-12) return false;
        if (p != c) { std::swap(A[c], A[p]); std::swap(I[c], I[p]); }

        double piv = A[c][c];
        for (int k = 0; k < 6; ++k) { A[c][k] /= piv; I[c][k] /= piv; }

        for (int r = 0; r < 6; ++r) {
            if (r == c) continue;
            double f = A[r][c];
            if (f == 0.0) continue;
            for (int k = 0; k < 6; ++k) { A[r][k] -= f*A[c][k]; I[r][k] -= f*I[c][k]; }
        }
    }
    out = I;
    return true;
}

// Axis-angle orientation error (radians): the rotation that takes `current`
// onto `target`, i.e. the log map of R_target * R_current^T.
std::array<double, 3> orientationError(const Transform& target, const Transform& current)
{
    double Rt[9] = { target[0],target[1],target[2], target[4],target[5],target[6], target[8],target[9],target[10] };
    double Rc[9] = { current[0],current[1],current[2], current[4],current[5],current[6], current[8],current[9],current[10] };

    double Rerr[9] = {};
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) {
            double s = 0;
            for (int k = 0; k < 3; ++k) s += Rt[r*3+k] * Rc[c*3+k];  // Rt * Rc^T
            Rerr[r*3+c] = s;
        }

    double tr = Rerr[0] + Rerr[4] + Rerr[8];
    double angle = std::acos(std::max(-1.0, std::min(1.0, (tr - 1.0) / 2.0)));

    if (angle > 1e-6) {
        double s2 = 2.0 * std::sin(angle);
        return {
            (Rerr[7] - Rerr[5]) / s2 * angle,
            (Rerr[2] - Rerr[6]) / s2 * angle,
            (Rerr[3] - Rerr[1]) / s2 * angle,
        };
    }
    return {0.0, 0.0, 0.0};
}

double norm3(const std::array<double, 3>& v)
{
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Entrywise lerp between two transforms  -  exact at both ends, only used to
// guide intermediate moveL steps (not a slerp; the position/rotation error
// each step is computed against, and Jacobian-corrected toward, is what
// actually matters, not a perfectly geodesic waypoint).
std::vector<Transform> lerpTransforms(const Transform& A, const Transform& B, int n_points)
{
    std::vector<Transform> out;
    if (n_points < 1) return out;
    for (int i = 0; i < n_points; ++i) {
        double t = (n_points == 1) ? 0.0 : double(i) / (n_points - 1);
        Transform T;
        for (int k = 0; k < 16; ++k) T[k] = A[k] + t * (B[k] - A[k]);
        out.push_back(T);
    }
    return out;
}

} // namespace

Robot::Robot() : name(), dh{}, limits{}, _solver(dh, limits)
{
}

void Robot::defineDH(const std::array<double, 6>& a, const std::array<double, 6>& d,
                      const std::array<double, 6>& alpha, const std::array<double, 6>& theta,
                      const JointLimits& lim)
{
    for (int i = 0; i < 6; ++i) {
        dh.a[i] = a[i]; dh.d[i] = d[i]; dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];
        dh.revolute[i] = true;
    }
    limits = lim;
    _solver = Solver(dh, limits);
}

Robot Robot::loadYAML(const std::string& path)
{
    Robots::Robot r = Robots::loadRobot(path);
    Robot robot;
    robot.name = r.name;
    robot.dh = r.dh;
    robot.limits = r.limits;
    robot._solver = Solver(robot.dh, robot.limits);
    return robot;
}

Transform Robot::forwardKinematics(const JointConfig& q) const
{
    return forwardKin(dh, q);
}

std::vector<JointConfig> Robot::inverseKinematics(const Transform& pose, int n_seeds, bool expand_wraps) const
{
    return _solver.solve(pose, expand_wraps, n_seeds);
}

Matrix6 Robot::jacobian(const JointConfig& q) const
{
    return geometricJacobian(dh, q);
}

double Robot::jacobianDeterminant(const JointConfig& q) const
{
    return det6(jacobian(q));
}

double Robot::manipulability(const JointConfig& q) const
{
    Matrix6 J = jacobian(q);
    return std::sqrt(std::fabs(det6(matmul6(J, transpose6(J)))));
}

std::vector<JointConfig> Robot::moveL(const Transform& poseA, const Transform& poseB,
                                       const JointConfig& q0, int n_points) const
{
    std::vector<Transform> waypoints = lerpTransforms(poseA, poseB, n_points);
    std::vector<JointConfig> result;

    JointConfig q = q0;
    Transform curr = forwardKinematics(q);

    for (std::size_t i = 0; i + 1 < waypoints.size(); ++i) {
        const Transform& target = waypoints[i+1];

        std::array<double, 3> dp = { target[3]-curr[3], target[7]-curr[7], target[11]-curr[11] };
        std::array<double, 3> dw = orientationError(target, curr);
        double dx[6] = { dp[0],dp[1],dp[2], dw[0],dw[1],dw[2] };

        Matrix6 J = jacobian(q);
        Matrix6 Jinv;
        if (!invert6(J, Jinv)) break;

        double dq_rad[6] = {};
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c)
                dq_rad[r] += Jinv[r][c] * dx[c];

        JointConfig q_try = q;
        for (int j = 0; j < 6; ++j) q_try[j] += dq_rad[j] * (180.0 / M_PI);
        Transform curr_try = forwardKinematics(q_try);

        double trans_error = norm3({ target[3]-curr_try[3], target[7]-curr_try[7], target[11]-curr_try[11] });
        double orien_error_deg = norm3(orientationError(target, curr_try)) * (180.0 / M_PI);

        if (trans_error > 0.001) break;
        if (orien_error_deg > 1.0) break;

        bool in_limits = true;
        for (int j = 0; j < 6; ++j)
            if (q_try[j] < limits.lo[j] || q_try[j] > limits.hi[j]) { in_limits = false; break; }
        if (!in_limits) break;

        q = q_try;
        curr = curr_try;
        result.push_back(q);
    }

    return result;
}

std::vector<JointConfig> Robot::moveJ(const JointConfig& qA, const JointConfig& qB, int n_points) const
{
    std::vector<JointConfig> result;
    if (n_points < 1) return result;
    for (int i = 0; i < n_points; ++i) {
        double t = (n_points == 1) ? 0.0 : double(i) / (n_points - 1);
        JointConfig q;
        for (int j = 0; j < 6; ++j) q[j] = qA[j] + t * (qB[j] - qA[j]);
        result.push_back(q);
    }
    return result;
}

} // namespace IKDH
