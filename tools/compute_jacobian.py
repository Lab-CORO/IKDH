#!/usr/bin/env python3
"""
compute_jacobian.py — Compute the symbolic geometric Jacobian for a 6R robot
from a DH YAML file and write the result back to the same file.

Usage (from the repository root):
    python3 tools/compute_jacobian.py robots/gofa5.yaml

Requires:
    pip install sympy

The Jacobian is appended to the YAML file as a 'jacobian' block.
Each of the 6×6 elements is stored as a Python expression string using:
    s1=sin(q1), c1=cos(q1), ..., s6=sin(q6), c6=cos(q6)  (radians)

To evaluate in Python:
    import math, yaml

    with open('robots/gofa5.yaml') as f:
        data = yaml.safe_load(f)

    q = [0.0, 0.785, 0.0, 0.0, 0.785, 0.0]  # joint angles in radians
    env = {f's{i+1}': math.sin(q[i]) for i in range(6)}
    env.update({f'c{i+1}': math.cos(q[i]) for i in range(6)})
    J = [[eval(e, {}, env) for e in row] for row in data['jacobian']['rows']]
"""

import sys
import os
import fractions

try:
    import sympy as sp
    from sympy.printing.str import StrPrinter
except ImportError:
    print("sympy is required:  pip install sympy")
    sys.exit(1)


# ── Custom printer: rationals and floats as compact decimals ──────────────────

class _CompactPrinter(StrPrinter):
    def _print_Rational(self, expr):
        if expr.q == 1:
            return str(expr.p)
        return f"{float(expr):.6g}"

    def _print_Float(self, expr):
        return f"{float(expr):.6g}"


_printer = _CompactPrinter()


# ── DH YAML parser (mirrors robots.h, handles pi expressions) ─────────────────

def _parse_pi_expr(s):
    s = s.strip()
    neg = s.startswith('-')
    if neg:
        s = s[1:].strip()
    if s == 'pi':
        val = sp.pi
    elif s.startswith('pi/'):
        val = sp.pi / sp.Integer(s[3:])
    elif s.startswith('pi*'):
        val = sp.pi * sp.Integer(s[3:])
    else:
        # Use exact Rational to avoid float contamination of pi
        frac = fractions.Fraction(s)
        val = sp.Rational(frac.numerator, frac.denominator)
    return -val if neg else val


def _parse_array(s):
    s = s.strip()
    if s.startswith('[') and s.endswith(']'):
        s = s[1:-1]
    return [_parse_pi_expr(t) for t in s.split(',')]


def _load_dh(yaml_path):
    a = d = alpha = theta = None
    with open(yaml_path) as f:
        for line in f:
            t = line.strip()
            if t.startswith('a:')     and '[' in t: a     = _parse_array(t[2:])
            if t.startswith('d:')     and '[' in t: d     = _parse_array(t[2:])
            if t.startswith('alpha:') and '[' in t: alpha = _parse_array(t[6:])
            if t.startswith('theta:') and '[' in t: theta = _parse_array(t[6:])
    if None in (a, d, alpha, theta):
        raise ValueError(f"Incomplete DH table in {yaml_path}")
    return a, d, alpha, theta


# ── Symbolic DH matrix ────────────────────────────────────────────────────────

def _dh_matrix(a, d, alpha, q_sym, theta_off):
    """4×4 symbolic DH transform: Rz(q+offset) * Tz(d) * Tx(a) * Rx(alpha)."""
    t = q_sym + theta_off
    ct, st = sp.cos(t), sp.sin(t)
    ca, sa = sp.cos(alpha), sp.sin(alpha)
    return sp.Matrix([
        [ct, -st*ca,  st*sa,  a*ct],
        [st,  ct*ca, -ct*sa,  a*st],
        [ 0,     sa,     ca,     d],
        [ 0,      0,      0,     1],
    ])


# ── Geometric Jacobian computation ────────────────────────────────────────────

def compute_jacobian(a, d, alpha, theta_off):
    """
    Returns the 6×6 symbolic geometric Jacobian as a sympy Matrix.
    Columns: one per joint (revolute).
    Rows 0-2: linear velocity  (z_{i-1} × (p_e - p_{i-1}))
    Rows 3-5: angular velocity (z_{i-1})
    """
    n = 6
    q = sp.symbols(f'q1:{n+1}')   # q1 … q6

    # Accumulate transforms:  T[i] = base → frame i
    T = [sp.eye(4)]
    for i in range(n):
        T.append(T[-1] * _dh_matrix(a[i], d[i], alpha[i], q[i], theta_off[i]))

    p_e = T[n][:3, 3]             # end-effector position

    J = sp.zeros(6, n)
    for i in range(n):
        z = T[i][:3, 2]           # z-axis of frame i-1
        p = T[i][:3, 3]           # origin of frame i-1
        J[:3, i] = z.cross(p_e - p)
        J[3:, i] = z

    print("  Applying trigsimp...")
    J = J.applyfunc(sp.trigsimp)
    return J


# ── Expression → compact string  (s1,c1 … s6,c6 notation) ───────────────────

def _to_str(expr):
    s = _printer.doprint(expr)
    for i in range(6, 0, -1):
        s = s.replace(f'sin(q{i})', f's{i}')
        s = s.replace(f'cos(q{i})', f'c{i}')
    return s


# ── YAML update ───────────────────────────────────────────────────────────────

def _remove_existing_jacobian(content):
    """Strip any existing 'jacobian:' block from the file content."""
    lines = content.splitlines()
    result, skip = [], False
    for line in lines:
        if line.startswith('jacobian:'):
            skip = True
            continue
        if skip and line and not line.startswith(' '):
            skip = False
        if not skip:
            result.append(line)
    # Remove trailing blank lines
    while result and not result[-1].strip():
        result.pop()
    return '\n'.join(result)


def write_jacobian(yaml_path, J):
    with open(yaml_path) as f:
        content = f.read()

    content = _remove_existing_jacobian(content)

    rows_yaml = []
    for row in range(6):
        exprs = [f'"{_to_str(J[row, col])}"' for col in range(6)]
        rows_yaml.append(f'    - [{", ".join(exprs)}]')

    block = [
        '',
        'jacobian:',
        '  rows:',
    ] + rows_yaml

    with open(yaml_path, 'w') as f:
        f.write(content + '\n' + '\n'.join(block) + '\n')

    print(f"  Written to {yaml_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: python3 {sys.argv[0]} robots/<robot>.yaml")
        sys.exit(1)

    yaml_path = sys.argv[1]
    if not os.path.exists(yaml_path):
        print(f"File not found: {yaml_path}")
        sys.exit(1)

    print(f"Loading DH from {yaml_path}...")
    a, d, alpha, theta = _load_dh(yaml_path)

    print("Computing symbolic Jacobian...")
    J = compute_jacobian(a, d, alpha, theta)

    print("Writing result...")
    write_jacobian(yaml_path, J)

    print("Done.")
