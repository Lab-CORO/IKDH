"""
demo.py — Mirror of examples/cpp/demo.cpp
Runs IK for all robots from a known joint configuration and prints solutions.

Usage (from the repository root):
    python3 examples/python/demo.py
"""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build'))
import ikdh
from ik_print import run_from_joints, run_from_pose

ROBOTS_DIR = "robots"


def main():
    # ── UR5e ──────────────────────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/ur5e.yaml")
    run_from_joints(r, [5.73, -68.75, 85.94, -103.13, -85.94, 28.65])

    # ── ABB GoFa CRB 15000-5 ──────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/gofa5.yaml")
    run_from_joints(r, [0, 45.0, 0, 0, 45.0, 0])
    run_from_pose(r, 571.0, 0.0, 899.0, 0.0, 90.0, 0.0)

    # ── ABB CRB 15000-10 ──────────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/crb15000_10.yaml")
    run_from_joints(r, [0, 45.0, 0, 0, 45.0, 180.0])

    # ── ABB GoFa CRB 15000-12 ─────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/gofa12.yaml")
    run_from_joints(r, [0, 45.0, 0, 0, 45.0, 180.0])
    run_from_pose(r, 700.0, 0.0, 400.0, 0.0, 90.0, 0.0)

    # ── FANUC CRX-10iA ────────────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/fanuc_crx_10ia.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])

    # ── FANUC CRX-5iA ─────────────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/fanuc_crx_5ia.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])

    # ── Doosan A0509 ──────────────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/doosan_a0509.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])

    # ── Auctech ACR-12 ────────────────────────────────────────────────────────
    r = ikdh.load_robot(f"{ROBOTS_DIR}/auctech_acr_12.yaml")
    run_from_joints(r, [30.0, 45.0, 60.0, -20.0, 30.0, 10.0])


if __name__ == "__main__":
    main()
