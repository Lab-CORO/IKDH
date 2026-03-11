#include <ikdh.h>
#include <robots.h>
#include <cstdio>

int main()
{
    // ── Load robot from YAML ──────────────────────────────────────────────────
    auto robot = Robots::loadRobot("robots/gofa5.yaml");
    IKDH::Solver solver(robot.dh, robot.limits);

    // ── Solve IK for two end-effector poses ───────────────────────────────────
    // poseFromXYZRPW: x y z in mm, rx ry rz in degrees (RoboDK convention Rz*Ry*Rx)
    IKDH::Transform poses[] = {
        IKDH::poseFromXYZRPW(200.0, 0.0, 600.0,   0.0, 90.0, 0.0),
        IKDH::poseFromXYZRPW(400.0, 0.0, 300.0, 180.0,  0.0, 0.0),
    };

    for (const auto& ee : poses) {
        auto sols = solver.solve(ee);

        printf("%zu solution(s) found\n", sols.size());
        for (size_t i = 0; i < sols.size(); ++i) {
            printf("  [%2zu]", i);
            for (double v : sols[i]) printf("  %7.3f", v);
            printf("   FK err = %.1e\n", IKDH::fkError(ee, IKDH::forwardKin(robot.dh, sols[i])));
        }
        printf("\n");
    }
}
