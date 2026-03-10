#!/usr/bin/env python3
"""Extract DH parameters from a RoboDK robot and save as a YAML robot file.

Usage:
    python tools/robodk_dh.py [robot_name]

Output:
    robots/<slug>.yaml   — created/overwritten automatically
    Also prints the YAML to stdout.
"""

import math, os, re, sys
from robodk.robolink import Robolink, ITEM_TYPE_ROBOT

# ── Matrix helpers ─────────────────────────────────────────────────────────────

def mm4(A, B):
    return [[sum(A[i][k]*B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def inv4(M):
    R  = [M[i][:3] for i in range(3)]
    t  = [M[i][3]  for i in range(3)]
    RT = [[R[j][i] for j in range(3)] for i in range(3)]
    p  = [-sum(RT[i][k]*t[k] for k in range(3)) for i in range(3)]
    return [RT[i] + [p[i]] for i in range(3)] + [[0,0,0,1]]

# ── Select robot ───────────────────────────────────────────────────────────────

RDK    = Robolink()
robots = RDK.ItemList(ITEM_TYPE_ROBOT)
arg    = " ".join(sys.argv[1:]).strip()

if arg:
    robot = next(r for r in robots if arg.lower() in r.Name().lower())
elif len(robots) == 1:
    robot = robots[0]
else:
    print("Robots:", [r.Name() for r in robots]); sys.exit(1)

robot_name = robot.Name()

# ── Extract Modified DH from JointPoses ───────────────────────────────────────

def mat_to_m4(m):
    return [[float(m[i, j]) for j in range(4)] for i in range(4)]

q0 = [0.0] * 6

all_poses_raw = robot.JointPoses(q0)
all_poses = [mat_to_m4(p) for p in all_poses_raw]
tcp_frame = mat_to_m4(robot.SolveFK(q0))

T_prev = all_poses
T_cum  = all_poses[1:] + [tcp_frame]

dh_mod = []
for i in range(6):
    T     = mm4(inv4(T_prev[i]), T_cum[i])
    alpha = math.atan2(-T[1][2],  T[2][2])
    theta = math.atan2(-T[0][1],  T[0][0])
    a_mm  = T[0][3]
    ca, sa = math.cos(alpha), math.sin(alpha)
    d_mm  = T[2][3] * ca - T[1][3] * sa
    dh_mod.append((theta, d_mm, a_mm, alpha))

# Joint limits (degrees)
lim   = robot.JointLimits()
lower = [float(v) for row in lim[0].rows for v in row]
upper = [float(v) for row in lim[1].rows for v in row]

# ── Convert Modified DH → Standard DH ────────────────────────────────────────

dh_params = []
for i in range(6):
    theta, d_mm, _, _ = dh_mod[i]
    a_s     = dh_mod[i + 1][2] * 1e-3 if i < 5 else 0.0
    alpha_s = dh_mod[i + 1][3]        if i < 5 else 0.0
    dh_params.append((theta, d_mm * 1e-3, a_s, alpha_s))

# ── Formatting ─────────────────────────────────────────────────────────────────

PI_FRACS = [
    (math.pi,     "pi"),    (-math.pi,     "-pi"),
    (math.pi/2,   "pi/2"),  (-math.pi/2,   "-pi/2"),
    (math.pi/3,   "pi/3"),  (-math.pi/3,   "-pi/3"),
    (math.pi/4,   "pi/4"),  (-math.pi/4,   "-pi/4"),
    (0.0, "0"),
]

def fmt_rad(v):
    for val, s in PI_FRACS:
        if abs(v - val) < 1e-9:
            return s
    return f"{v:.6f}"

def fmt_m(v):
    s = f"{v:.4f}".rstrip('0').rstrip('.')
    return s if '.' in s else s + '.0'

# ── Build YAML ─────────────────────────────────────────────────────────────────

a_vals     = [fmt_m  (p[2]) for p in dh_params]
d_vals     = [fmt_m  (p[1]) for p in dh_params]
alpha_vals = [fmt_rad(p[3]) for p in dh_params]
theta_vals = [fmt_rad(p[0]) for p in dh_params]

def yaml_list(vals):
    return "[" + ", ".join(vals) + "]"

yaml_lines = [
    f"name: {robot_name}",
    "dh:",
    f"  a:     {yaml_list(a_vals)}",
    f"  d:     {yaml_list(d_vals)}",
    f"  alpha: {yaml_list(alpha_vals)}",
    f"  theta: {yaml_list(theta_vals)}",
    "limits:",
]
for i in range(6):
    yaml_lines.append(f"  J{i+1}: [{lower[i]:7.1f}, {upper[i]:7.1f}]")

yaml_content = "\n".join(yaml_lines) + "\n"

# ── Print + save ──────────────────────────────────────────────────────────────

print(yaml_content)

slug = re.sub(r'_+', '_', re.sub(r'[^a-z0-9_]', '', re.sub(r'[\s\-\./]+', '_', robot_name.lower()))).strip('_')
out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "robots")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, slug + ".yaml")

with open(out_path, "w") as f:
    f.write(yaml_content)

print(f"Saved → {out_path}")
