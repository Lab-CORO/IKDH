#!/usr/bin/env python3
"""Print DH table and joint limits of a RoboDK robot as C++ (robots.h format).

Usage: python tools/robodk_dh.py [robot_name]
"""

import math, re, sys
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

# ── Extract DH modified from JointPoses ──────────────────────────────────────
#
# robot.JointPoses(q) returns the CUMULATIVE transform (base → frame_i) for
# each joint i.  To get the Modified DH parameters we need the INTERMEDIATE
# (relative) transform between consecutive frames:
#
#   T_rel[i] = inv(T_cum[i-1]) · T_cum[i],   T_cum[-1] = identity (base)
#
# Each relative transform satisfies the Modified DH equation:
#
#   T = Rx(α) · Tx(a) · Rz(θ) · Tz(d)
#
#   ┌  cθ    -sθ     0    a   ┐
#   │  cα·sθ  cα·cθ -sα  -sα·d│
#   │  sα·sθ  sα·cθ  cα   cα·d│
#   └  0      0      0    1   ┘
#
# Reading off the four DH parameters from this matrix:
#   α = atan2(-T[1][2],  T[2][2])
#   θ = atan2(-T[0][1],  T[0][0])
#   a = T[0][3]                          (mm, RoboDK units)
#   d = T[2][3]·cos(α) − T[1][3]·sin(α) (dot product recovers d)

def mat_to_m4(m):
    """Convert a robodk Mat to a plain 4×4 list-of-lists."""
    return [[float(m[i, j]) for j in range(4)] for i in range(4)]

q0 = [0.0] * 6

# robot.JointPoses(q) returns n+1 cumulative frames:
#   [base_frame, joint_1_frame, ..., joint_(n-1)_frame]
# For a 6-DOF robot this gives 6 elements (base + joints 1..5).
# The last joint frame (J6/TCP) is absent; we fetch it via SolveFK.
#
# Relative transform for joint i:
#   T_rel[i] = inv(all_poses[i]) · all_poses[i+1]   for i = 0..4
#   T_rel[5] = inv(all_poses[5]) · tcp_frame
all_poses_raw = robot.JointPoses(q0)
all_poses = [mat_to_m4(p) for p in all_poses_raw]   # 6 frames: base + joints 1..5
tcp_frame = mat_to_m4(robot.SolveFK(q0))            # joint 6 = TCP

T_prev = all_poses                       # [base, j1, j2, j3, j4, j5]
T_cum  = all_poses[1:] + [tcp_frame]    # [j1,   j2, j3, j4, j5, tcp]

dh_mod = []   # [(theta_rad, d_mm, a_mm, alpha_rad), ...]  — Modified DH

for i in range(6):
    T     = mm4(inv4(T_prev[i]), T_cum[i])
    alpha = math.atan2(-T[1][2],  T[2][2])
    theta = math.atan2(-T[0][1],  T[0][0])
    a_mm  = T[0][3]
    ca, sa = math.cos(alpha), math.sin(alpha)
    d_mm  = T[2][3] * ca - T[1][3] * sa
    dh_mod.append((theta, d_mm, a_mm, alpha))

# Joint limits (degrees) — JointLimits() returns two 1×n Mat objects
lim   = robot.JointLimits()
lower = [float(v) for row in lim[0].rows for v in row]
upper = [float(v) for row in lim[1].rows for v in row]

# ── Convert DH standard ───────────────────────────────────────────────────────
#
# Modified DH → Standard DH (convention used by robots.h / libhupf):
#
#   Standard DH:  T = Rz(θ) · Tz(d) · Tx(a) · Rx(α)
#
# In Modified DH, the (α, a) of row i describe the link between axis i-1 and
# axis i.  In Standard DH, (α, a) of row i describe the link between axis i
# and axis i+1.  Therefore:
#
#   α_s[i] = α_m[i+1],   a_s[i] = a_m[i+1]   ← shift one row forward
#   θ_s[i] = θ_m[i],     d_s[i] = d_m[i]      ← unchanged
#
# The last joint gets α_s[5] = a_s[5] = 0 (end-effector frame is aligned).

dh_params = []   # [(theta_rad, d_m, a_m, alpha_rad), ...]  — Standard DH
for i in range(6):
    theta, d_mm, _, _ = dh_mod[i]
    a_s     = dh_mod[i + 1][2] * 1e-3 if i < 5 else 0.0   # mm → m
    alpha_s = dh_mod[i + 1][3]        if i < 5 else 0.0
    dh_params.append((theta, d_mm * 1e-3, a_s, alpha_s))

# ── Formatting ─────────────────────────────────────────────────────────────────

PI_FRACS = [(math.pi, "M_PI"), (-math.pi, "-M_PI"),
            (math.pi/2, "M_PI/2"), (-math.pi/2, "-M_PI/2"),
            (math.pi/3, "M_PI/3"), (-math.pi/3, "-M_PI/3"),
            (math.pi/4, "M_PI/4"), (-math.pi/4, "-M_PI/4"), (0.0, "0.0")]

def fmt_rad(v):
    for val, s in PI_FRACS:
        if abs(v - val) < 1e-9: return s
    return f"{v:.6f}"

def fmt_m(v):
    s = f"{v:.3f}".rstrip('0').rstrip('.')
    return s if '.' in s else s + '.0'

def arr(strs):
    w = max(len(s) for s in strs)
    return "{ " + ", ".join(s.ljust(w) for s in strs) + " }"

func = re.sub(r'_+', '_', re.sub(r'[^a-z0-9_]', '', re.sub(r'[\s\-\./]+', '_', robot_name.lower()))).strip('_')

# ── Output ─────────────────────────────────────────────────────────────────────

dashes = "─" * max(0, 80 - 6 - len(robot_name))
print(f"// ── {robot_name} {dashes}")
print(f"inline IKDH::DHTable {func}_dh()")
print("{")
print("    IKDH::DHTable dh;")
print(f"    const double a[]     = {arr([fmt_m  (p[2]) for p in dh_params])};")
print(f"    const double d[]     = {arr([fmt_m  (p[1]) for p in dh_params])};")
print(f"    const double alpha[] = {arr([fmt_rad(p[3]) for p in dh_params])};")
print(f"    const double theta[] = {arr([fmt_rad(p[0]) for p in dh_params])};")
print("    for (int i = 0; i < 6; ++i) {")
print("        dh.a[i] = a[i]; dh.d[i] = d[i];")
print("        dh.alpha[i] = alpha[i]; dh.theta[i] = theta[i];")
print("        dh.revolute[i] = true;")
print("    }")
print("    return dh;")
print("}\n")
print(f"inline IKDH::JointLimits {func}_limits()")
print("{")
print("    return IKDH::JointLimits({")
for i in range(6):
    print(f"        {{{lower[i]:7.1f}, {upper[i]:7.1f}}}{',' if i<5 else ' '}  // J{i+1}")
print("    });")
print("}")