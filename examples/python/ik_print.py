"""
ik_print.py — Python equivalent of examples/cpp/ik_print.h
Shared helpers for printing IK results.
"""

import ikdh


def run_from_joints(robot, q_ref, label=None):
    name = label or robot.name
    print(f"-> {name}")
    print(f"  Joints: {' '.join(f'{v:7.2f}' for v in q_ref)} deg")

    ee = ikdh.forward_kin(robot.dh, q_ref)
    solver = ikdh.Solver(robot.dh, robot.limits)
    sols = solver.solve(ee)

    print(f"{len(sols)} solution(s) found")
    for i, q in enumerate(sols):
        err = ikdh.fk_error(ee, ikdh.forward_kin(robot.dh, q))
        print(f"  [{i:2d}]  {' '.join(f'{v:7.3f}' for v in q)}   FK err = {err:.1e}")
    print()


def run_from_pose(robot, x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg, label=None):
    name = label or robot.name
    print(f"-> {name}")
    print(f"  Pose: ({x_mm:.3f}, {y_mm:.3f}, {z_mm:.3f}) mm  "
          f"Rx={rx_deg:.1f}  Ry={ry_deg:.1f}  Rz={rz_deg:.1f} deg")

    ee = ikdh.pose_from_xyzrpw(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg)
    solver = ikdh.Solver(robot.dh, robot.limits)
    sols = solver.solve(ee)

    print(f"{len(sols)} solution(s) found")
    for i, q in enumerate(sols):
        err = ikdh.fk_error(ee, ikdh.forward_kin(robot.dh, q))
        print(f"  [{i:2d}]  {' '.join(f'{v:7.3f}' for v in q)}   FK err = {err:.1e}")
    print()
