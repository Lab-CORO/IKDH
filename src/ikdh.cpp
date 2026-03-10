#include <ikdh.h>
#include <hupf/ik.h>
#include <hupf/poly_capture.h>

#include <cmath>

namespace IKDH {

// ── Internal helpers ──────────────────────────────────────────────────────────

// Wrap angle (degrees) into [lo, hi] by adding/subtracting multiples of 360°.
// Returns false if the angle cannot be mapped into the interval.
static bool wrapAngle(double angle, double lo, double hi, double& out)
{
    // Normalise to (-180, 180]
    angle = std::fmod(angle, 360.0);
    if (angle >  180.0) angle -= 360.0;
    if (angle <= -180.0) angle += 360.0;

    // Try k=0 first (keep the normalised angle if already in range), then k=±1.
    // This prefers the most compact representation and avoids unnecessary wrapping
    // (e.g. 36° with limits [-360,360] stays 36°, not -324°).
    for (int k : {0, -1, 1}) {
        double a = angle + k * 360.0;
        if (a >= lo && a <= hi) { out = a; return true; }
    }
    return false;
}

static void mul4(const double A[16], const double B[16], double C[16])
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            C[i*4+j] = 0;
            for (int k = 0; k < 4; ++k)
                C[i*4+j] += A[i*4+k] * B[k*4+j];
        }
}

// T_i = Rz(theta) * Tz(d) * Tx(a) * Rx(alpha)
static void dhMatrix(double theta, double d, double a, double alpha, double T[16])
{
    double ct = std::cos(theta), st = std::sin(theta);
    double ca = std::cos(alpha), sa = std::sin(alpha);

    T[0]  = ct;  T[1]  = -st*ca;  T[2]  =  st*sa;  T[3]  = a*ct;
    T[4]  = st;  T[5]  =  ct*ca;  T[6]  = -ct*sa;  T[7]  = a*st;
    T[8]  = 0;   T[9]  =  sa;     T[10] =  ca;      T[11] = d;
    T[12] = 0;   T[13] = 0;       T[14] = 0;        T[15] = 1;
}

// Solve the 6×6 linear system A x = b via Gaussian elimination
// with partial pivoting.  Returns false if A is (near-)singular.
static bool solveLinear6(double A[6][6], double b[6], double x[6])
{
    double M[6][7];
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) M[i][j] = A[i][j];
        M[i][6] = b[i];
    }
    for (int col = 0; col < 6; ++col) {
        int piv = col;
        for (int row = col + 1; row < 6; ++row)
            if (std::fabs(M[row][col]) > std::fabs(M[piv][col])) piv = row;
        if (std::fabs(M[piv][col]) < 1e-12) return false;
        for (int j = 0; j <= 6; ++j) std::swap(M[col][j], M[piv][j]);
        double pivot = M[col][col];
        for (int row = 0; row < 6; ++row) {
            if (row == col) continue;
            double f = M[row][col] / pivot;
            for (int j = col; j <= 6; ++j) M[row][j] -= f * M[col][j];
        }
    }
    for (int i = 0; i < 6; ++i) x[i] = M[i][6] / M[i][i];
    return true;
}

// Newton-Raphson refinement: adjust q (degrees) until forwardKin(dh, q) ≈ T_target.
//
// Uses the analytical geometric Jacobian (computed from intermediate DH
// transforms accumulated during the FK pass) instead of finite differences.
// This replaces 6 extra FK calls per iteration with pure matrix operations,
// giving ~6× speedup on the inner loop at identical convergence.
//
// Damped least-squares (Tikhonov) step keeps wrist singularities stable.
// Returns true when the Frobenius FK error drops below tol.
static bool refineIK(const DHTable& dh, JointConfig& q, const Transform& T_target,
                     int maxIter = 40, double tol = 1e-9)
{
    const double lambda = 1e-6;  // Tikhonov damping

    for (int iter = 0; iter < maxIter; ++iter) {

        // ── FK with intermediate transforms ───────────────────────────────────
        // T_accum[j] = cumulative transform up to (but not including) joint j,
        // i.e. the pose of frame j-1 in base coordinates.
        //   T_accum[0] = I        (base frame, provides rotation axis for joint 0)
        //   T_accum[1] = DH_0     (frame 0,    provides rotation axis for joint 1)
        //   ...
        //   T_accum[6] = DH_0*…*DH_5  = end-effector pose
        double T_accum[7][16];

        // Initialise T_accum[0] = Identity
        for (int k = 0; k < 16; ++k) T_accum[0][k] = 0.0;
        T_accum[0][0] = T_accum[0][5] = T_accum[0][10] = T_accum[0][15] = 1.0;

        for (int j = 0; j < 6; ++j) {
            double val = dh.revolute[j] ? q[j] * (M_PI / 180.0) : q[j];
            double jt  = dh.theta[j] + (dh.revolute[j] ? val : 0.0);
            double jd  = dh.d[j]     + (dh.revolute[j] ? 0.0 : val);
            double Ti[16];
            dhMatrix(jt, jd, dh.a[j], dh.alpha[j], Ti);
            mul4(T_accum[j], Ti, T_accum[j + 1]);
        }

        // End-effector pose = T_accum[6]
        const double* Tcur = T_accum[6];

        // 12-component error: first 3 rows of (T_target - T_current)
        double err[12];
        double fe = 0.0;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 4; ++c) {
                double e = T_target[r*4+c] - Tcur[r*4+c];
                err[r*4+c] = e;
                fe += e * e;
            }
        if (fe < tol * tol) return true;

        // ── Analytical geometric Jacobian J (12×6) ────────────────────────────
        // For revolute joint j rotating about z_{j-1} (column 2 of T_accum[j]):
        //
        //   z_j  = T_accum[j] col-2 = [T2, T6, T10]
        //   p_j  = T_accum[j] col-3 = [T3, T7, T11]
        //   p_e  = T_accum[6] col-3
        //   R_e  = T_accum[6] 3×3 rotation block
        //
        //   ∂p_e / ∂q_j  =  z_j × (p_e − p_j)               (position rows, c=3)
        //   ∂R_e[:,c] / ∂q_j  =  z_j × R_e[:,c]   c=0,1,2   (rotation rows)
        double J[12][6] = {};

        const double p_e[3] = { Tcur[3], Tcur[7], Tcur[11] };

        for (int j = 0; j < 6; ++j) {
            const double* Tj = T_accum[j];
            const double zj[3] = { Tj[2], Tj[6], Tj[10] };

            if (!dh.revolute[j]) {
                // Prismatic: only translates along z_j, no rotation change.
                for (int r = 0; r < 3; ++r) J[r*4+3][j] = zj[r];
                continue;
            }

            const double pj[3] = { Tj[3], Tj[7], Tj[11] };

            // Position part: z_j × (p_e − p_j)
            const double rv[3] = { p_e[0]-pj[0], p_e[1]-pj[1], p_e[2]-pj[2] };
            J[0*4+3][j] = zj[1]*rv[2] - zj[2]*rv[1];
            J[1*4+3][j] = zj[2]*rv[0] - zj[0]*rv[2];
            J[2*4+3][j] = zj[0]*rv[1] - zj[1]*rv[0];

            // Rotation part: z_j × R_e[:,c]  for each column c
            for (int c = 0; c < 3; ++c) {
                const double rc[3] = { Tcur[c], Tcur[4+c], Tcur[8+c] };
                J[0*4+c][j] = zj[1]*rc[2] - zj[2]*rc[1];
                J[1*4+c][j] = zj[2]*rc[0] - zj[0]*rc[2];
                J[2*4+c][j] = zj[0]*rc[1] - zj[1]*rc[0];
            }
        }

        // Damped normal equations:  (J^T J + λ I) dq_rad = J^T err
        double JTJ[6][6] = {};
        double JTe[6]    = {};
        for (int a = 0; a < 6; ++a) {
            for (int b = 0; b < 6; ++b)
                for (int k = 0; k < 12; ++k)
                    JTJ[a][b] += J[k][a] * J[k][b];
            JTJ[a][a] += lambda;
            for (int k = 0; k < 12; ++k)
                JTe[a] += J[k][a] * err[k];
        }

        double dq_rad[6];
        if (!solveLinear6(JTJ, JTe, dq_rad)) return false;

        for (int j = 0; j < 6; ++j) {
            q[j] += dq_rad[j] * (180.0 / M_PI);
            // Normalise revolute joints to (-180, 180] after each step.
            // This prevents "winding" (accumulating thousands of degrees) when
            // the starting configuration is far from the target in joint space.
            if (dh.revolute[j]) {
                q[j] = std::fmod(q[j], 360.0);
                if (q[j] >  180.0) q[j] -= 360.0;
                if (q[j] <= -180.0) q[j] += 360.0;
            }
        }
    }
    return false;
}

// ── Solver ────────────────────────────────────────────────────────────────────

Solver::Solver(const DHTable& dh, const JointLimits& limits)
    : _limits(limits), _dh(dh)
{
    double a[6], d[6], theta[6], alpha[6];
    bool   rots[6];
    for (int i = 0; i < 6; ++i) {
        a[i]          = dh.a[i];
        d[i]          = dh.d[i];
        theta[i]      = dh.theta[i];
        alpha[i]      = dh.alpha[i];
        rots[i]       = dh.revolute[i];
        _theta_off[i] = dh.theta[i];
    }
    _impl = new LibHUPF::ik_solver(a, d, theta, alpha, rots);
}

Solver::~Solver()
{
    delete static_cast<LibHUPF::ik_solver*>(_impl);
}

std::vector<JointConfig> Solver::solve(const Transform& ee) const
{
    auto* iks = static_cast<LibHUPF::ik_solver*>(_impl);
    const double rad2deg = 180.0 / M_PI;

    // Solve IK for a flat 4×4 matrix, subtract theta offsets, wrap, filter limits.
    auto processRaw = [&](double* mat) {
        auto raw = iks->solve(mat);
        std::vector<JointConfig> res;
        res.reserve(raw.size());
        for (const auto& sol : raw) {
            JointConfig q;
            bool valid = true;
            for (int i = 0; i < 6; ++i) {
                // libhupf returns total DH angles (offset + user); subtract offset.
                double angle = sol[i] - _theta_off[i] * rad2deg;
                if (!wrapAngle(angle, _limits.lo[i], _limits.hi[i], q[i])) {
                    valid = false; break;
                }
            }
            if (valid) res.push_back(q);
        }
        return res;
    };

    double flat[16];
    for (int i = 0; i < 16; ++i) flat[i] = ee[i];
    auto result = processRaw(flat);

    // Snapshot the resultant polynomial from the primary solve before the
    // fallback overwrites poly_capture_ref() with perturbed-pose polynomials.
    LibHUPF::primary_poly_ref() = LibHUPF::poly_capture_ref();

    // Post-refine primary HuPf solutions.  Near degenerate poses the algebraic
    // roots can have poor precision (FK error ~1e-5 to 1e-7).  A few Newton
    // iterations bring them to < 1e-9 with negligible extra cost.
    //
    // After Newton, angles are normalised to (-180, 180] by the per-step
    // wrap inside refineIK.  For joints whose limit range extends below -180°
    // (e.g. J3 ∈ [-225, 85]) this can flip a valid -203° into +157°, which
    // is outside the limit.  Re-wrap to fix the representation; drop any
    // solution that genuinely cannot be mapped into limits.
    {
        std::vector<JointConfig> refined;
        refined.reserve(result.size());
        for (auto q : result) {
            refineIK(_dh, q, ee, 40);
            JointConfig qw;
            bool valid = true;
            for (int i = 0; i < 6; ++i)
                if (!wrapAngle(q[i], _limits.lo[i], _limits.hi[i], qw[i]))
                    { valid = false; break; }
            if (valid) refined.push_back(qw);
        }
        result = std::move(refined);
    }

    // ── Degeneracy / infinity fallback ────────────────────────────────────────
    // The HuPf solver uses the right-chain substitution u4 = 1/v4 where
    // v4 = tan(θ4/2).  This means:
    //   • J4 = ±180°  →  v4 = ±∞  →  u4 = 0   →  found by the polynomial
    //   • J4 =   0°   →  v4 =  0  →  u4 = ∞   →  NOT found (root at infinity)
    // Additionally, rank-deficient KinematicSurface matrices (e.g. pure-Ry
    // poses on a robot with α₅=0) further reduce the number of roots found.
    //
    // Strategy: right-multiply the target by small rotation perturbations to
    // break both issues, solve the non-degenerate perturbed problem, then
    // refine each new candidate back onto the original target via damped
    // Newton-Raphson.  Six independent directions × three magnitudes give good
    // coverage of all solution basins.
    //
    // We target up to 16 solutions (theoretical algebraic maximum for a 6R robot).
    // Stop early once we reach that count.
    {
        auto tryRefineAndAdd = [&](JointConfig q) {
            if (!refineIK(_dh, q, ee, 100)) return;
            JointConfig qw;
            bool valid = true;
            for (int i = 0; i < 6; ++i) {
                if (!wrapAngle(q[i], _limits.lo[i], _limits.hi[i], qw[i])) {
                    valid = false; break;
                }
            }
            if (!valid) return;
            // Deduplication via circular distance (mod 360°).
            // Threshold: sum-of-squares < 1 deg² (≈ 0.4° RMS per joint).
            for (const auto& r : result) {
                double d2 = 0;
                for (int i = 0; i < 6; ++i) {
                    double dd = std::fmod(std::fabs(r[i] - qw[i]), 360.0);
                    if (dd > 180.0) dd = 360.0 - dd;
                    d2 += dd * dd;
                }
                if (d2 < 1.0) return;
            }
            result.push_back(qw);
        };

        // Try a seed and its shoulder/wrist flip variants.
        // Flipping J1 by ±180° (shoulder flip), J4 by ±180° (wrist flip), or
        // J6 by ±180° (tool flip) jumps to solution basins unreachable by Newton
        // alone.  All 8 combinations of the three binary flips are tried.
        auto tryWithFlips = [&](JointConfig q) {
            for (int mask = 0; mask < 8; ++mask) {
                JointConfig qf = q;
                if (mask & 1) qf[0] += 180.0;  // J1 flip
                if (mask & 2) qf[3] += 180.0;  // J4 flip
                if (mask & 4) qf[5] += 180.0;  // J6 flip
                tryRefineAndAdd(qf);
            }
        };

        // Helper: run one perturbation magnitude (rotation only).
        auto runRotPerturbations = [&](double eps) {
            double ce = std::cos(eps), se = std::sin(eps);

            double Rs[6][16] = {
                { ce,-se, 0, 0,  se, ce, 0, 0,  0,  0, 1, 0,  0, 0, 0, 1 }, // Rz+
                { ce, se, 0, 0, -se, ce, 0, 0,  0,  0, 1, 0,  0, 0, 0, 1 }, // Rz−
                {  1,  0, 0, 0,   0, ce,-se, 0,  0, se,ce, 0,  0, 0, 0, 1 }, // Rx+
                {  1,  0, 0, 0,   0, ce, se, 0,  0,-se,ce, 0,  0, 0, 0, 1 }, // Rx−
                { ce,  0,se, 0,   0,  1, 0, 0, -se,  0,ce, 0,  0, 0, 0, 1 }, // Ry+
                { ce,  0,-se,0,   0,  1, 0, 0,  se,  0,ce, 0,  0, 0, 0, 1 }, // Ry−
            };
            for (auto& R : Rs) {
                if (result.size() >= 16) break;
                double perturbed[16];
                mul4(flat, R, perturbed);
                for (auto q : processRaw(perturbed)) tryWithFlips(q);
            }

            double Rls[6][9] = {
                { ce,-se, 0,  se, ce, 0,  0, 0, 1 }, // Left-Rz+
                { ce, se, 0, -se, ce, 0,  0, 0, 1 }, // Left-Rz−
                {  1,  0, 0,   0, ce,-se, 0,se,ce }, // Left-Rx+
                {  1,  0, 0,   0, ce, se, 0,-se,ce}, // Left-Rx−
                { ce,  0,se,   0,  1,  0,-se, 0,ce}, // Left-Ry+
                { ce,  0,-se,  0,  1,  0, se, 0,ce}, // Left-Ry−
            };
            for (auto& Rl : Rls) {
                if (result.size() >= 16) break;
                double perturbed[16];
                for (int k = 0; k < 16; ++k) perturbed[k] = flat[k];
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j) {
                        perturbed[i*4+j] = 0;
                        for (int k = 0; k < 3; ++k)
                            perturbed[i*4+j] += Rl[i*3+k] * flat[k*4+j];
                    }
                for (auto q : processRaw(perturbed)) tryWithFlips(q);
            }
        };

        // ── Phase 1 (always): standard magnitudes 1°, 5°, 10° ────────────────
        for (double eps : { 1.0 * M_PI/180.0, 5.0 * M_PI/180.0, 10.0 * M_PI/180.0 }) {
            if (result.size() >= 16) break;
            runRotPerturbations(eps);
        }

        // ── Phase 2 (last resort): joint-space Halton sampling ───────────────
        // Independent of libhupf entirely.  Activated when Phase 1+2 still leave
        // fewer than 4 solutions (e.g. robots where the HuPf polynomial returns
        // zero real roots for all poses, such as the Doosan A0509).
        // Halton sequence (bases 2,3,5,7,11,13) gives quasi-uniform coverage of
        // the 6D joint-space box [lo, hi]^6.
        if (result.size() < 4) {
            auto halton = [](int idx, int base) -> double {
                double f = 1.0, r = 0.0;
                for (; idx > 0; idx /= base) { f /= base; r += f * (idx % base); }
                return r;
            };
            static const int bases[6] = { 2, 3, 5, 7, 11, 13 };
            for (int s = 1; s <= 100 && result.size() < 16; ++s) {
                JointConfig q_seed;
                for (int j = 0; j < 6; ++j)
                    q_seed[j] = _limits.lo[j]
                              + halton(s, bases[j])
                              * (_limits.hi[j] - _limits.lo[j]);
                tryWithFlips(q_seed);
            }
        }
    }

    return result;
}

std::vector<double> Solver::lastPolynomial() const
{
    return LibHUPF::primary_poly_ref();
}

// ── Forward kinematics ────────────────────────────────────────────────────────

Transform forwardKin(const DHTable& dh, const JointConfig& q)
{
    double T[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};

    for (int i = 0; i < 6; ++i) {
        double val = dh.revolute[i] ? q[i] * M_PI / 180.0 : q[i];

        double joint_theta = dh.theta[i] + (dh.revolute[i] ? val : 0.0);
        double joint_d     = dh.d[i]     + (dh.revolute[i] ? 0.0 : val);

        double Ti[16], Tnew[16];
        dhMatrix(joint_theta, joint_d, dh.a[i], dh.alpha[i], Ti);
        mul4(T, Ti, Tnew);
        for (int k = 0; k < 16; ++k) T[k] = Tnew[k];
    }

    Transform result;
    for (int k = 0; k < 16; ++k) result[k] = T[k];
    return result;
}

} // namespace IKDH
