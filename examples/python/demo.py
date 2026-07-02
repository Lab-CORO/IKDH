"""
demo.py — Mirror of examples/cpp/demo.cpp
Runs IK for all robots from a known joint configuration and prints solutions.

Usage (from the repository root):
    python3 examples/python/demo.py
"""

from _common import print_solutions
import ikdh

def run_from_joints(robot, q_ref, label=None):
    name = label or robot.name
    print(f"-> {name}")
    print(f"  Joints: {' '.join(f'{v:7.2f}' for v in q_ref)} deg")

    ee = ikdh.forward_kin(robot.dh, q_ref)
    solver = ikdh.Solver(robot.dh, robot.limits)
    print_solutions(solver, robot.dh, ee)


def run_from_pose(robot, x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg, label=None):
    name = label or robot.name
    print(f"-> {name}")
    print(f"  Pose: ({x_mm:.3f}, {y_mm:.3f}, {z_mm:.3f}) mm  "
          f"Rx={rx_deg:.1f}  Ry={ry_deg:.1f}  Rz={rz_deg:.1f} deg")

    ee = ikdh.pose_from_xyzrpw(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg)
    solver = ikdh.Solver(robot.dh, robot.limits)
    print_solutions(solver, robot.dh, ee)


def main():
    # ── UR5e ──────────────────────────────────────────────────────────────────
    r = ikdh.load_robot("robots/ur5e.yaml")
    run_from_joints(r, [5.73, -68.75, 85.94, -103.13, -85.94, 28.65])

    # ── ABB GoFa CRB 15000-5 ──────────────────────────────────────────────────
    r = ikdh.load_robot("robots/gofa5.yaml")
    run_from_joints(r, [0, 45.0, 0, 0, 45.0, 0])
    run_from_pose(r, 571.0, 0.0, 899.0, 0.0, 90.0, 0.0)

    # ── ABB CRB 15000-10 ──────────────────────────────────────────────────────
    r = ikdh.load_robot("robots/crb15000_10.yaml")
    run_from_joints(r, [0, 45.0, 0, 0, 45.0, 180.0])

    # ── ABB GoFa CRB 15000-12 ─────────────────────────────────────────────────
    r = ikdh.load_robot("robots/gofa12.yaml")
    run_from_joints(r, [0, 45.0, 0, 0, 45.0, 180.0])
    run_from_pose(r, 700.0, 0.0, 400.0, 0.0, 90.0, 0.0)

    # ── FANUC CRX-10iA ────────────────────────────────────────────────────────
    r = ikdh.load_robot("robots/fanuc_crx_10ia.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])

    # ── FANUC CRX-5iA ─────────────────────────────────────────────────────────
    r = ikdh.load_robot("robots/fanuc_crx_5ia.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])

    # ── Doosan A0509 ──────────────────────────────────────────────────────────
    r = ikdh.load_robot("robots/doosan_a0509.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])

    # ── Auctech ACR-12 ────────────────────────────────────────────────────────
    r = ikdh.load_robot("robots/auctech_acr_12.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, 0.0, 30.0, 10.0])


if __name__ == "__main__":
    main()
