"""
run_from_yaml.py — Minimal example: load a robot from YAML and solve IK.

Usage (from the repository root):
    python3 examples/python/run_from_yaml.py
"""

from _common import print_solutions
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
    print_solutions(solver, robot.dh, ee)
