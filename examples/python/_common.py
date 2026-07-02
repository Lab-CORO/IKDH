"""
_common.py — Shared setup and printing helpers for the examples in this directory.

Not a public API; imported directly by sibling scripts run as
`python3 examples/python/<script>.py`.
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build'))
import ikdh


def print_solutions(solver, dh, ee):
    sols = solver.solve(ee)
    print(f"{len(sols)} solution(s) found")
    for i, q in enumerate(sols):
        err = ikdh.fk_error(ee, ikdh.forward_kin(dh, q))
        print(f"  [{i:2d}]  {' '.join(f'{v:7.3f}' for v in q)}   FK err = {err:.1e}")
    print()
