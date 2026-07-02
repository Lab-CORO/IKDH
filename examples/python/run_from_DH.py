"""
run_from_DH.py — Minimal example: define a DH table directly and solve IK.

Usage (from the repository root):
    python3 examples/python/run_from_DH.py
"""

import math

from _common import print_solutions
import ikdh

# ABB GoFa CRB 15000-5 — DH table (a/d in metres, alpha/theta in radians)
dh = ikdh.DHTable(
    a     = [0.0,        0.444,      0.110,      0.0,       0.080,      0.0  ],
    d     = [0.265,      0.0,        0.0,        0.470,     0.0,        0.101],
    alpha = [-math.pi/2, 0,         -math.pi/2,  math.pi/2,-math.pi/2,  0   ],
    theta = [0,         -math.pi/2,  0,           0,         0,          math.pi],
)

# Joint limits (degrees)
limits = ikdh.JointLimits([(-180,180), (-180,180), (-225,85),
                           (-180,180), (-180,180), (-180,180)])

solver = ikdh.Solver(dh, limits)

# Solve for two poses (x, y, z in mm — rx, ry, rz in degrees)
poses = [
    ikdh.pose_from_xyzrpw(200.0, 0.0, 600.0,   0.0, 90.0, 0.0),
    ikdh.pose_from_xyzrpw(400.0, 0.0, 300.0, 180.0,  0.0, 0.0),
]

for ee in poses:
    print_solutions(solver, dh, ee)
