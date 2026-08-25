"""
demo.py  -  Forward and inverse kinematics for a couple of robots.

Usage (from the repository root):
    python3 examples/python/demo.py
"""

from _common import print_solutions
import ikdh

# Forward kinematics: joint angles (deg) -> end-effector pose.
robot = ikdh.load_robot("robots/fanuc_crx_10ia.yaml")
q = [30.0, 45.0, 60.0, -20.0, 30.0, 10.0]
ee = ikdh.forward_kin(robot.dh, q)
solver = ikdh.Solver(robot.dh, robot.limits)

print(f"-> {robot.name} (from joint angles)")
print_solutions(solver, robot.dh, ee)

# Inverse kinematics: end-effector pose (mm, deg) -> all joint solutions.
robot2 = ikdh.load_robot("robots/gofa5.yaml")
target = ikdh.pose_from_xyzrpw(571.0, 0.0, 899.0, 0.0, 90.0, 0.0)
solver2 = ikdh.Solver(robot2.dh, robot2.limits)

print(f"-> {robot2.name} (from a Cartesian pose)")
print_solutions(solver2, robot2.dh, target)
