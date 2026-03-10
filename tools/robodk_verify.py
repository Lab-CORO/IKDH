#!/usr/bin/env python3
"""Verify IK solutions against RoboDK.

Runs the demo binary, parses its output for a given robot, streams each
IK solution to a live RoboDK session, and reports the Cartesian position
deviation measured by RoboDK's forward kinematics.

Usage:
    python3 tools/robodk_verify.py [path/to/robot.yaml] [robodk_robot_name]

    Default robot file: robots/gofa5.yaml
    Default RoboDK name: taken from the YAML (can differ from the scene name)

Prerequisites:
    pip install robodk
    cmake --build build   # build/demo must exist
"""

import math, re, sys, time, subprocess

try:
    from robodk.robolink import Robolink, ITEM_TYPE_ROBOT
except ImportError:
    print("robodk package not found. Install it with:  pip install robodk")
    sys.exit(1)

# ── Load robot YAML ───────────────────────────────────────────────────────────

def parse_pi_expr(s):
    s = s.strip()
    neg = s.startswith('-')
    t = s.lstrip('-').strip()
    if t == 'pi':
        val = math.pi
    elif t.startswith('pi/'):
        val = math.pi / float(t[3:])
    elif t.startswith('pi*'):
        val = math.pi * float(t[3:])
    else:
        val = float(t)
    return -val if neg else val

def parse_array(line):
    lb, rb = line.find('['), line.find(']')
    return [parse_pi_expr(x) for x in line[lb+1:rb].split(',')]

def load_robot_yaml(path):
    robot = {'name': '', 'dh': {}, 'limits': []}
    section = None
    with open(path) as f:
        for line in f:
            s = line.split('#')[0].rstrip()
            if not s.strip():
                continue
            if 'name:' in s:
                robot['name'] = s.split(':', 1)[1].strip()
            elif s.strip() == 'dh:':
                section = 'dh'
            elif s.strip() == 'limits:':
                section = 'limits'
            elif section == 'dh' and '[' in s:
                if 'a:' in s:         robot['dh']['a']     = parse_array(s)
                elif 'd:' in s:       robot['dh']['d']     = parse_array(s)
                elif 'alpha:' in s:   robot['dh']['alpha'] = parse_array(s)
                elif 'theta:' in s:   robot['dh']['theta'] = parse_array(s)
            elif section == 'limits' and '[' in s:
                row = parse_array(s)
                robot['limits'].append((row[0], row[1]))
    return robot

# ── Configuration ──────────────────────────────────────────────────────────────

YAML_FILE = sys.argv[1] if len(sys.argv) > 1 else "robots/gofa5.yaml"
robot_cfg = load_robot_yaml(YAML_FILE)

ROBOT_LABEL  = robot_cfg['name']        # must match "-> <name>" header in demo output
ROBOT_RDK    = sys.argv[2] if len(sys.argv) > 2 else robot_cfg['name']  # robot name in RoboDK scene
JOINT_LIMITS = robot_cfg['limits']      # [(lo, hi), ...] in degrees
INTERVAL     = 2                      # seconds between solutions
DEMO_BIN     = 'build/gofa5'

# ── Helpers ───────────────────────────────────────────────────────────────────

def wrap_angle(angle, lo, hi):
    """Map angle into [lo, hi] by ±360 steps. Returns None if impossible."""
    angle = angle % 360.0
    if angle > 180.0:
        angle -= 360.0
    for offset in (0, 360, -360):
        v = angle + offset
        if lo <= v <= hi:
            return v
    return None

def pos_xyz(mat):
    """Extract (x, y, z) mm from a RoboDK Mat (row-major .rows)."""
    return mat.rows[0][3], mat.rows[1][3], mat.rows[2][3]

# ── 1. Parse solutions from demo binary ───────────────────────────────────────

proc   = subprocess.run([DEMO_BIN], capture_output=True, text=True)
output = proc.stdout

solutions_raw = []
in_block = False
for line in output.splitlines():
    if line.strip() == f'-> {ROBOT_LABEL}':
        in_block = True
        continue
    if in_block and line.strip().startswith('->'):
        break
    if in_block:
        m = re.search(r'\[\s*(\d+)\s*\]((?:\s+[-\d.]+){6})', line)
        if m:
            angles = [float(x) for x in m.group(2).split()]
            err_m  = re.search(r'FK err\s*=\s*([0-9e.+\-]+)', line)
            err    = float(err_m.group(1)) if err_m else None
            solutions_raw.append((angles, err))

if not solutions_raw:
    print(f"No solutions found for '{ROBOT_LABEL}' in demo output.")
    print("Raw output:\n", output)
    sys.exit(1)

# Apply joint limits
solutions = []
for angles, err in solutions_raw:
    wrapped = [wrap_angle(a, lo, hi) for a, (lo, hi) in zip(angles, JOINT_LIMITS)]
    if all(w is not None for w in wrapped):
        solutions.append((wrapped, err))

print(f"{len(solutions_raw)} raw solution(s), {len(solutions)} within joint limits.\n")

# ── 2. Connect to RoboDK ──────────────────────────────────────────────────────

RDK   = Robolink()
robot = RDK.Item(ROBOT_RDK, ITEM_TYPE_ROBOT)
if not robot.Valid():
    print(f"Robot '{ROBOT_RDK}' not found in RoboDK.")
    sys.exit(1)

# ── 3. Reference FK from RoboDK ───────────────────────────────────────────────

ref_joints  = solutions[0][0] if solutions else None
pose_ref    = robot.SolveFK(ref_joints) if ref_joints else None
xr, yr, zr  = pos_xyz(pose_ref) if pose_ref else (None, None, None)

print(f"{'Sol':>4}  {'joint angles (deg)':^54}  {'FK err (internal)':>18}  {'delta pos RoboDK':>16}")
print("-" * 104)

# ── 4. Display each solution ──────────────────────────────────────────────────

for i, (angles, err) in enumerate(solutions):
    robot.setJoints(angles)

    pose_sol = robot.SolveFK(angles)
    x, y, z  = pos_xyz(pose_sol)
    delta    = ((x-xr)**2 + (y-yr)**2 + (z-zr)**2)**0.5 if xr is not None else float('nan')

    angles_str = '  '.join(f'{a:8.3f}' for a in angles)
    print(f"  {i:>2}  [{angles_str}]  {err:.1e}  delta={delta:.2f} mm")

    time.sleep(INTERVAL)

print("\nDone.")
