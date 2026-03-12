//  path_planning.cpp — Cartesian path planning with IK branch tracking.
//
//  Interpolates a straight-line Cartesian path (linear translation + SLERP
//  rotation) between two poses, solves IK at every step, and tracks all
//  solution branches simultaneously.  A branch is killed (entry set to empty)
//  as soon as the closest available IK solution requires a joint-space jump
//  larger than JUMP_THRESHOLD_DEG degrees (Euclidean norm over the 6 joints).
//
//  Build (from the repository root, after cmake --build build):
//    g++ -std=c++17 -Iinclude examples/cpp/path_planning.cpp \
//        build/libikdh.a -o path_planning
//  or add it to CMakeLists.txt alongside the other examples.

#include <ikdh.h>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <optional>
#include <vector>

// ── Quaternion helpers ────────────────────────────────────────────────────────

struct Quat { double w, x, y, z; };

static Quat mat_to_quat(const IKDH::Transform& T)
{
    // Extract rotation from row-major 4x4 transform
    double m00=T[0], m01=T[1], m02=T[2];
    double m10=T[4], m11=T[5], m12=T[6];
    double m20=T[8], m21=T[9], m22=T[10];

    double trace = m00 + m11 + m22;
    Quat q;
    if (trace > 0) {
        double s = 0.5 / std::sqrt(trace + 1.0);
        q.w = 0.25 / s;
        q.x = (m21 - m12) * s;
        q.y = (m02 - m20) * s;
        q.z = (m10 - m01) * s;
    } else if (m00 > m11 && m00 > m22) {
        double s = 2.0 * std::sqrt(1.0 + m00 - m11 - m22);
        q.w = (m21 - m12) / s;
        q.x = 0.25 * s;
        q.y = (m01 + m10) / s;
        q.z = (m02 + m20) / s;
    } else if (m11 > m22) {
        double s = 2.0 * std::sqrt(1.0 + m11 - m00 - m22);
        q.w = (m02 - m20) / s;
        q.x = (m01 + m10) / s;
        q.y = 0.25 * s;
        q.z = (m12 + m21) / s;
    } else {
        double s = 2.0 * std::sqrt(1.0 + m22 - m00 - m11);
        q.w = (m10 - m01) / s;
        q.x = (m02 + m20) / s;
        q.y = (m12 + m21) / s;
        q.z = 0.25 * s;
    }
    return q;
}

static Quat slerp(Quat a, Quat b, double t)
{
    double dot = a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z;
    if (dot < 0.0) { b.w=-b.w; b.x=-b.x; b.y=-b.y; b.z=-b.z; dot=-dot; }
    if (dot > 0.9995) {
        // Linear interpolation for near-identical quaternions
        return { a.w + t*(b.w-a.w), a.x + t*(b.x-a.x),
                 a.y + t*(b.y-a.y), a.z + t*(b.z-a.z) };
    }
    double theta_0 = std::acos(dot);
    double theta   = theta_0 * t;
    double s0 = std::cos(theta) - dot * std::sin(theta) / std::sin(theta_0);
    double s1 = std::sin(theta) / std::sin(theta_0);
    return { s0*a.w + s1*b.w, s0*a.x + s1*b.x,
             s0*a.y + s1*b.y, s0*a.z + s1*b.z };
}

static IKDH::Transform quat_to_transform(const Quat& q,
                                          double tx, double ty, double tz)
{
    double w=q.w, x=q.x, y=q.y, z=q.z;
    IKDH::Transform T{};
    T[ 0]=1-2*(y*y+z*z); T[ 1]=2*(x*y-w*z);   T[ 2]=2*(x*z+w*y);   T[ 3]=tx;
    T[ 4]=2*(x*y+w*z);   T[ 5]=1-2*(x*x+z*z); T[ 6]=2*(y*z-w*x);   T[ 7]=ty;
    T[ 8]=2*(x*z-w*y);   T[ 9]=2*(y*z+w*x);   T[10]=1-2*(x*x+y*y); T[11]=tz;
    T[12]=0; T[13]=0; T[14]=0; T[15]=1;
    return T;
}

// ── Cartesian interpolation ───────────────────────────────────────────────────

static std::vector<IKDH::Transform>
interpolate_poses(const IKDH::Transform& A, const IKDH::Transform& B, int n)
{
    Quat qa = mat_to_quat(A), qb = mat_to_quat(B);
    double ta[3] = { A[3], A[7], A[11] };
    double tb[3] = { B[3], B[7], B[11] };

    std::vector<IKDH::Transform> poses;
    poses.reserve(n);
    for (int i = 0; i < n; ++i) {
        double s = static_cast<double>(i) / (n - 1);
        Quat q   = slerp(qa, qb, s);
        double tx = (1.0-s)*ta[0] + s*tb[0];
        double ty = (1.0-s)*ta[1] + s*tb[1];
        double tz = (1.0-s)*ta[2] + s*tb[2];
        poses.push_back(quat_to_transform(q, tx, ty, tz));
    }
    return poses;
}

// ── Branch tracking ───────────────────────────────────────────────────────────

// A branch entry: either a valid joint config, or nullopt (dead).
using Entry  = std::optional<IKDH::JointConfig>;
using Branch = std::vector<Entry>;

static double joint_dist(const IKDH::JointConfig& a, const IKDH::JointConfig& b)
{
    double sum = 0;
    for (int i = 0; i < 6; ++i) { double d = a[i]-b[i]; sum += d*d; }
    return std::sqrt(sum);
}

static std::vector<Branch>
plan_paths(const IKDH::Solver& solver,
           const std::vector<IKDH::Transform>& poses,
           double jump_threshold_deg = 30.0)
{
    auto sols0 = solver.solve(poses[0]);
    if (sols0.empty()) { printf("[ERROR] No IK solution at pose A\n"); return {}; }

    // One branch per solution at pose A
    std::vector<Branch> branches;
    for (const auto& q : sols0)
        branches.push_back({ Entry{q} });

    // Precompute IK for all remaining poses
    std::vector<std::vector<IKDH::JointConfig>> all_solutions;
    all_solutions.reserve(poses.size() - 1);
    for (size_t i = 1; i < poses.size(); ++i)
        all_solutions.push_back(solver.solve(poses[i]));

    // Extend each branch independently
    for (auto& branch : branches) {
        for (const auto& solutions : all_solutions) {
            if (!branch.back().has_value() || solutions.empty()) {
                branch.push_back(std::nullopt);
                continue;
            }

            const auto& prev_q = branch.back().value();
            double best_dist = std::numeric_limits<double>::max();
            size_t best_idx  = 0;
            for (size_t i = 0; i < solutions.size(); ++i) {
                double d = joint_dist(solutions[i], prev_q);
                if (d < best_dist) { best_dist = d; best_idx = i; }
            }

            if (best_dist > jump_threshold_deg)
                branch.push_back(std::nullopt);
            else
                branch.push_back(Entry{solutions[best_idx]});
        }
    }
    return branches;
}

// ── Main ──────────────────────────────────────────────────────────────────────

int main()
{
    const int    N_POINTS           = 50;
    const double JUMP_THRESHOLD_DEG = 30.0;

    IKDH::Transform pose_A = IKDH::poseFromXYZRPW(400.0,  100.0, 600.0,  0.0, 90.0,  0.0);
    IKDH::Transform pose_B = IKDH::poseFromXYZRPW(400.0, -100.0, 400.0,  0.0, 90.0, 45.0);

    // ── Robot (DH + limits) ───────────────────────────────────────────────────
    const double pi = M_PI;
    IKDH::DHTable dh;
    dh.a    [0]=0.000; dh.a    [1]=0.444; dh.a    [2]=0.110;
    dh.a    [3]=0.000; dh.a    [4]=0.080; dh.a    [5]=0.000;
    dh.d    [0]=0.265; dh.d    [1]=0.000; dh.d    [2]=0.000;
    dh.d    [3]=0.470; dh.d    [4]=0.000; dh.d    [5]=0.101;
    dh.alpha[0]=-pi/2; dh.alpha[1]=0;     dh.alpha[2]=-pi/2;
    dh.alpha[3]= pi/2; dh.alpha[4]=-pi/2; dh.alpha[5]=0;
    dh.theta[0]=0;     dh.theta[1]=-pi/2; dh.theta[2]=0;
    dh.theta[3]=0;     dh.theta[4]=0;     dh.theta[5]=pi;
    for (int i = 0; i < 6; ++i) dh.revolute[i] = true;

    IKDH::JointLimits limits({ {-180,180},{-180,180},{-225,85},
                                {-180,180},{-180,180},{-180,180} });
    IKDH::Solver solver(dh, limits);

    // ── Plan ──────────────────────────────────────────────────────────────────
    auto poses    = interpolate_poses(pose_A, pose_B, N_POINTS);
    auto branches = plan_paths(solver, poses, JUMP_THRESHOLD_DEG);

    printf("%zu branch(es) from A to B  (%d poses, threshold=%.0f deg)\n\n",
           branches.size(), N_POINTS, JUMP_THRESHOLD_DEG);

    int n_full = 0;
    for (size_t b = 0; b < branches.size(); ++b) {
        const auto& branch = branches[b];
        int n_valid = 0;
        int died_at = -1;
        for (size_t i = 0; i < branch.size(); ++i) {
            if (branch[i].has_value()) ++n_valid;
            else if (died_at < 0)      died_at = static_cast<int>(i);
        }
        if (died_at < 0) {
            printf("  Branch %2zu: %2d/%d valid — reaches B ✓\n",
                   b, n_valid, N_POINTS);
            ++n_full;
        } else {
            printf("  Branch %2zu: %2d/%d valid — dies at step %d\n",
                   b, n_valid, N_POINTS, died_at);
        }
    }

    printf("\n%d complete path(s):\n", n_full);
    for (size_t b = 0; b < branches.size(); ++b) {
        const auto& branch = branches[b];
        if (!branch.back().has_value()) continue;
        const auto& q0 = branch.front().value();
        const auto& qN = branch.back().value();
        printf("  Branch %2zu  start=[", b);
        for (int j = 0; j < 6; ++j) printf("%7.2f%s", q0[j], j<5?", ":"");
        printf("]  end=[");
        for (int j = 0; j < 6; ++j) printf("%7.2f%s", qN[j], j<5?", ":"");
        printf("]\n");
    }
}