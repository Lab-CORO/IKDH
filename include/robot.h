#pragma once

#include <ikdh.h>

#include <array>
#include <string>
#include <vector>

namespace IKDH {

// 6x6 matrix, row-major access as J[row][col].
using Matrix6 = std::array<std::array<double, 6>, 6>;

// Central robot object: DH table + joint limits + everything derived from them
// (FK, IK, Jacobian/manipulability, Cartesian and joint-space path planning).
class Robot {
public:
    Robot();

    // Set the DH parameters and joint limits directly.
    void defineDH(const std::array<double, 6>& a, const std::array<double, 6>& d,
                  const std::array<double, 6>& alpha, const std::array<double, 6>& theta,
                  const JointLimits& limits = JointLimits{});

    // Load DH parameters and joint limits from a robot YAML file.
    static Robot loadYAML(const std::string& path);

    Transform forwardKinematics(const JointConfig& q) const;

    // n_seeds sets how many Halton seeds cover the joint-space cube before
    // Newton refinement (see Solver::solve).
    std::vector<JointConfig> inverseKinematics(const Transform& pose,
                                                int n_seeds = 64,
                                                bool expand_wraps = false) const;

    // Geometric Jacobian at q (degrees). Rows: vx,vy,vz,wx,wy,wz. Columns: joints.
    // Derivatives are with respect to radians, the standard robotics convention.
    Matrix6 jacobian(const JointConfig& q) const;
    double  jacobianDeterminant(const JointConfig& q) const;
    double  manipulability(const JointConfig& q) const;  // sqrt(det(J * J^T))

    // Cartesian path from poseA to poseB via inverse-Jacobian steps, n_points
    // waypoints including both ends. Stops early (returning what was found so
    // far) if a step's translation error exceeds 1mm, its rotation error
    // exceeds 1 degree, or it leaves the joint limits.
    std::vector<JointConfig> moveL(const Transform& poseA, const Transform& poseB,
                                    const JointConfig& q0, int n_points) const;

    // Joint-space path from qA to qB via straight linear interpolation,
    // n_points waypoints including both ends.
    std::vector<JointConfig> moveJ(const JointConfig& qA, const JointConfig& qB,
                                    int n_points) const;

    std::string name;
    DHTable     dh;
    JointLimits limits;

private:
    Solver _solver;
};

} // namespace IKDH
