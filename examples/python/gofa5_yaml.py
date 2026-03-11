import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build'))
import ikdh

# ── Load robot from YAML ──────────────────────────────────────────────────────
robot  = ikdh.load_robot("robots/gofa5.yaml")
solver = ikdh.Solver(robot.dh, robot.limits)

# ── Solve IK for two end-effector poses ───────────────────────────────────────
# pose_from_xyzrpw: x y z in mm, rx ry rz in degrees (RoboDK convention Rz*Ry*Rx)
poses = [
    ikdh.pose_from_xyzrpw(200.0, 0.0, 600.0,   0.0, 90.0, 0.0),
    ikdh.pose_from_xyzrpw(400.0, 0.0, 300.0, 180.0,  0.0, 0.0),
]

for ee in poses:
    sols = solver.solve(ee)
    print(f"{len(sols)} solution(s) found")
    for i, q in enumerate(sols):
        err = ikdh.fk_error(ee, ikdh.forward_kin(robot.dh, q))
        print(f"  [{i:2d}]  {' '.join(f'{v:7.3f}' for v in q)}   FK err = {err:.1e}")
    print()
