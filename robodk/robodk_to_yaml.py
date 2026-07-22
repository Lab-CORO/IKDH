"""
robodk_to_yaml.py  -  Converts a RoboDK .robot binary file to an IKDH-compatible .yaml file.

Binary layout (reverse-engineered from UR5e.robot):
  Header  : 4 bytes, uint32 big-endian = decompressed size
  Payload : zlib-compressed blob, little-endian float64 values
  DH      : 6 joints × 160-byte blocks starting at offset 1624
                +0   → d      [m]
                +136 → alpha  [rad]
                +144 → a      [m]
            theta = 0 at home position
  Limits  : 6 × float64 at offset 4800 (lower) and 4960 (upper), in degrees

Usage (from the repository root):
    python3 robodk/robodk_to_yaml.py robodk/UR5e.robot
    python3 robodk/robodk_to_yaml.py robodk/UR5e.robot output_dir/
"""

import os
import struct
import sys
import zlib

from _format import fmt_angle


# Binary layout constants
HEADER_SIZE = 4
N_JOINTS    = 6
JOINT_BLOCK = 160
FIRST_JOINT = 1624
FIELD_D     = 0
FIELD_ALPHA = 136
FIELD_A     = 144
LIMITS_LO   = 4800   # 6 × float64, degrees
LIMITS_HI   = 4960   # 6 × float64, degrees
#


def load_robot_file(path):
    with open(path, "rb") as f:
        raw = f.read()
    expected = struct.unpack(">I", raw[:HEADER_SIZE])[0]
    payload  = zlib.decompress(raw[HEADER_SIZE:])
    if len(payload) != expected:
        raise ValueError(f"Size mismatch: expected {expected} bytes, got {len(payload)}")
    return payload


def _f64(data, offset):
    return struct.unpack_from("<d", data, offset)[0]


def extract_dh(data):
    rows = []
    for j in range(N_JOINTS):
        base = FIRST_JOINT + j * JOINT_BLOCK
        rows.append({
            "d":     _f64(data, base + FIELD_D),
            "a":     _f64(data, base + FIELD_A),
            "alpha": _f64(data, base + FIELD_ALPHA),
            "theta": 0.0,
        })
    return rows


def extract_limits(data):
    lo = [_f64(data, LIMITS_LO + j * 8) for j in range(N_JOINTS)]
    hi = [_f64(data, LIMITS_HI + j * 8) for j in range(N_JOINTS)]
    return lo, hi


def _fmt_float(v):
    """Compact float representation (strip trailing zeros)."""
    return f"{v:.6g}"


def save_yaml(path, name, dh, lo, hi):
    # d and a are in mm in the binary  -  convert to metres for IKDH.
    # Note: sign conventions may differ per manufacturer; verify after conversion.
    lines = []
    lines.append(f"name: {name}")
    lines.append("dh:")

    def _row(key, values, fmt):
        inner = ", ".join(fmt(v) for v in values)
        lines.append(f"  {key}: [{inner}]")

    _row("a",     [r["a"] / 1000.0 for r in dh], _fmt_float)
    _row("d",     [r["d"] / 1000.0 for r in dh], _fmt_float)
    _row("alpha", [r["alpha"] for r in dh], fmt_angle)
    _row("theta", [r["theta"] for r in dh], fmt_angle)

    lines.append("limits:")
    for j in range(N_JOINTS):
        lines.append(f"  J{j+1}: [{lo[j]:.1f}, {hi[j]:.1f}]")

    lines.append("")   # trailing newline
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def convert(robot_path, out_dir=None):
    """Convert a .robot file to .yaml. Returns the output path."""
    data = load_robot_file(robot_path)
    dh   = extract_dh(data)
    lo, hi = extract_limits(data)

    name     = os.path.splitext(os.path.basename(robot_path))[0]
    out_dir  = out_dir or os.path.dirname(os.path.abspath(robot_path))
    out_path = os.path.join(out_dir, name + ".yaml")

    save_yaml(out_path, name, dh, lo, hi)
    return out_path


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    robot_path = sys.argv[1]
    out_dir    = sys.argv[2] if len(sys.argv) > 2 else None

    print(f"Loading: {robot_path}")
    data     = load_robot_file(robot_path)
    dh       = extract_dh(data)
    lo, hi   = extract_limits(data)

    print(f"Decompressed payload: {len(data):,} bytes\n")

    print(f"{'Joint':^7} | {'d (m)':^10} | {'a (m)':^10} | {'alpha':^10} | {'theta':^10}")
    print("-" * 58)
    for j, r in enumerate(dh):
        print(f"  J{j+1}   | {r['d']/1000:10.5f} | {r['a']/1000:10.5f} | {fmt_angle(r['alpha']):^10} | {fmt_angle(r['theta']):^10}")

    print()
    print(f"{'Joint':^7} | {'lower (°)':^12} | {'upper (°)':^12}")
    print("-" * 38)
    for j in range(N_JOINTS):
        print(f"  J{j+1}   | {lo[j]:12.2f} | {hi[j]:12.2f}")

    out_path = convert(robot_path, out_dir)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
