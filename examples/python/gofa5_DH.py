"""
gofa5_DH.py — Minimal example: define a DH table directly and solve IK.

Usage (from the repository root):
    python3 examples/python/gofa5_DH.py
"""

import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build'))
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
    sols = solver.solve(ee)
    print(f"{len(sols)} solution(s) found")
    for i, q in enumerate(sols):
        err = ikdh.fk_error(ee, ikdh.forward_kin(dh, q))
        print(f"  [{i:2d}]  {' '.join(f'{v:7.3f}' for v in q)}   FK err = {err:.1e}")
    print()
