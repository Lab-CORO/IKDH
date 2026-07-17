#!/usr/bin/env python3
"""
regenerate_jacobians.py — Re-run derive_jacobian.py's Jacobian computation over
every robots/*.yaml (or a given list), overwriting each file's 'jacobian' block
in place. The DH table and limits are left untouched.

Usage (from the repository root):
    python3 tools/regenerate_jacobians.py                  # every robots/*.yaml
    python3 tools/regenerate_jacobians.py robots/ur5e.yaml robots/gofa5.yaml
"""

import glob
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import derive_jacobian as dj


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    paths = sys.argv[1:] or sorted(glob.glob(os.path.join(root, "robots", "*.yaml")))

    ok, failed = 0, []
    t0 = time.time()
    for i, path in enumerate(paths):
        print(f"[{i+1}/{len(paths)}] {path}")
        try:
            a, d, alpha, theta = dj._load_dh(path)
            J = dj.compute_jacobian(a, d, alpha, theta)
            dj.write_jacobian(path, J)
            ok += 1
        except Exception as e:
            print(f"  FAILED: {e}")
            failed.append((path, str(e)))

    print(f"\n{ok} regenerated, {len(failed)} failed, in {time.time()-t0:.0f}s.")
    for p, err in failed:
        print(f"  FAILED: {p}: {err}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
