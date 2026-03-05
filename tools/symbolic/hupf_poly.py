#!/usr/bin/env python3
"""
tools/symbolic/hupf_poly.py

Symbolically reduce the HuPf degree-56 IK polynomial to the genuine
degree-16 polynomial for the ABB GoFa CRB 15000-5 (6R robot).

Method:
  1. Port the HuPf hyperplane equations from C++ to SymPy
  2. Build the kinematic surface r[0..7] via 7x7 symbolic determinants
  3. Form e1 (Study quadric) and e2 (8th hyperplane constraint)
  4. Compute resultant(e1, e2, u) → univariate polynomial in v
     with coefficients that are polynomials in t0..t7 (Study parameters of TCP)
  5. Factor over Q(t0..t7) → extract the degree-16 factor
  6. Save coefficient expressions (for fast runtime evaluation)

Runtime: potentially several hours. Checkpoints are saved automatically.

Usage:
  pip install sympy
  python hupf_poly.py [--checkpoint <file>] [--robot gofa5]
"""

import sympy as sp
from sympy import Rational, symbols, Matrix, Poly, factor_list, resultant
import pickle, time, os, sys, argparse

###############################################################################
# Checkpoint helpers
###############################################################################

CHECKPOINT_DIR = os.path.join(os.path.dirname(__file__), "checkpoints")
os.makedirs(CHECKPOINT_DIR, exist_ok=True)

def save(name, obj):
    path = os.path.join(CHECKPOINT_DIR, name + ".pkl")
    with open(path, "wb") as f:
        pickle.dump(obj, f)
    print(f"  [saved {path}]")

def load(name):
    path = os.path.join(CHECKPOINT_DIR, name + ".pkl")
    if os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    return None

def timed(label, func):
    print(f"\n=== {label} ===")
    t0 = time.time()
    result = func()
    print(f"  done in {time.time()-t0:.1f}s")
    return result

###############################################################################
# GoFa5 DH parameters (exact rationals)
# Convention: al_i = tan(alpha_i/2), alpha = [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]
###############################################################################

al1 = Rational(-1)       # tan(-pi/4)
al2 = Rational(0)
al3 = Rational(-1)       # tan(-pi/4)
al4 = Rational(1)        # tan(pi/4)
al5 = Rational(-1)       # tan(-pi/4)
al6 = Rational(0)

a1 = Rational(0)
a2 = Rational(444, 1000)
a3 = Rational(110, 1000)
a4 = Rational(0)
a5 = Rational(80,  1000)
a6 = Rational(0)

d2 = Rational(0)
d3 = Rational(0)
d4 = Rational(470, 1000)
d5 = Rational(0)
d6 = Rational(101, 1000)

###############################################################################
# Symbolic variables
###############################################################################

# Study parameters of TCP pose (8D projective dual-quaternion representation)
t0, t1, t2, t3, t4, t5, t6, t7 = symbols('t0:8')

# Half-angle substitution variables
u = symbols('u')   # left chain:  v2 = tan(theta2/2)
v = symbols('v')   # right chain: v4 = tan(theta4/2)

###############################################################################
# LEFT CHAIN hyperplanes: h_tc_v2
# (GoFa5 left chain uses Tv2 parameterization — see Hyperplane.h logic)
# Each function returns a list of 8 SymPy expressions in u.
# C++ Polynomial(c0, c1) → c0 + c1*u
###############################################################################

def h1_tc_v2():
    return [
        (-2*al3*al1*d2 + 2*al3*al2*d2 - 2*al3*d3*al2 - 2*al3*d3*al1)
        + (2*al3*a2 + 2*al3*al1*al2*a2 - 2*al3*al1*a1*al2 - 2*al3*a1 + 2*a3*al1 - 2*al2*a3)*u,

        (2*al1*d2 - 2*al2*d2 + 2*d3*al2 + 2*d3*al1)
        + (-2*a2 - 2*al1*al2*a2 + 2*al1*a1*al2 + 2*a1 + 2*al3*a3*al1 - 2*a3*al3*al2)*u,

        (-2*a1 - 2*a2 + 2*a2*al1*al2 + 2*a1*al1*al2 - 2*a3*al3*al2 - 2*a3*al3*al1)
        + (2*d2*al1 + 2*d2*al2 + 2*d3*al1 - 2*d3*al2)*u,

        (2*al3*a1 + 2*al3*a2 - 2*al3*al1*al2*a2 - 2*al3*al1*al2*a1 - 2*al2*a3 - 2*al1*a3)
        + (-2*al1*d2*al3 - 2*d2*al3*al2 - 2*d3*al1*al3 + 2*al3*d3*al2)*u,

        sp.Integer(0)  + (-4*al3*(al1 - al2))*u,
        sp.Integer(0)  + ( 4*al1 - 4*al2)*u,
        (-4*al2 - 4*al1),
        ( 4*al3*(al1 + al2)),
    ]

def h2_tc_v2():
    return [
        (-2*d2 - 2*al1*al2*d2 + 2*d3*al1*al2 - 2*d3)
        + (2*al2*a2 + 2*al1*a1 - 2*al1*a2 - 2*al2*a1 - 2*al3*a3*al1*al2 - 2*al3*a3)*u,

        (-2*al3*d2 - 2*al3*al1*al2*d2 + 2*al3*d3*al1*al2 - 2*al3*d3)
        + (2*a2*al3*al2 + 2*al3*al1*a1 - 2*al3*al1*a2 - 2*a1*al3*al2 + 2*a3*al1*al2 + 2*a3)*u,

        (2*a2*al3*al1 + 2*a1*al3*al1 + 2*a1*al3*al2 + 2*a2*al3*al2 - 2*a3*al1*al2 + 2*a3)
        + (2*d2*al3 - 2*al1*d2*al3*al2 + 2*al3*d3*al1*al2 + 2*d3*al3)*u,

        (2*al1*a2 + 2*al1*a1 + 2*al2*a1 + 2*al2*a2 + 2*al3*al1*al2*a3 - 2*al3*a3)
        + (2*d2 - 2*al1*d2*al2 + 2*d3*al1*al2 + 2*d3)*u,

        sp.Integer(0) + (-4*al1*al2 - 4)*u,
        sp.Integer(0) + (-4*al3*(1 + al1*al2))*u,
        4*al3*(-1 + al1*al2),
        4*al1*al2 - 4,
    ]

def h3_tc_v2():
    return [
        (2*al3*a1 + 2*al3*a2 - 2*al3*al1*al2*a2 - 2*al3*al1*al2*a1 - 2*al2*a3 - 2*al1*a3)
        + (-2*al1*d2*al3 - 2*d2*al3*al2 - 2*d3*al1*al3 + 2*al3*d3*al2)*u,

        (-2*a1 - 2*a2 + 2*a2*al1*al2 + 2*a1*al1*al2 - 2*a3*al3*al2 - 2*a3*al3*al1)
        + (2*d2*al1 + 2*d2*al2 + 2*d3*al1 - 2*d3*al2)*u,

        (2*al2*d2 - 2*al1*d2 - 2*d3*al2 - 2*d3*al1)
        + (2*a2 - 2*a1 + 2*al1*al2*a2 - 2*al1*a1*al2 + 2*a3*al3*al2 - 2*al3*a3*al1)*u,

        (-2*al3*al2*d2 + 2*al3*al1*d2 + 2*al3*d3*al2 + 2*al3*d3*al1)
        + (-2*al3*a2 + 2*al3*a1 - 2*al3*al1*al2*a2 + 2*al3*al1*a1*al2 + 2*al2*a3 - 2*a3*al1)*u,

        4*al3*(al1 + al2),
        (-4*al2 - 4*al1),
        sp.Integer(0) + (4*al2 - 4*al1)*u,
        sp.Integer(0) + (4*al3*(al1 - al2))*u,
    ]

def h4_tc_v2():
    return [
        (-2*al2*a1 - 2*al1*a1 - 2*al2*a2 - 2*al1*a2 + 2*al3*a3 - 2*al3*al1*al2*a3)
        + (-2*d2 + 2*al1*d2*al2 - 2*d3*al1*al2 - 2*d3)*u,

        (-2*a1*al3*al2 - 2*a1*al3*al1 - 2*a2*al3*al2 - 2*a2*al3*al1 - 2*a3 + 2*a3*al1*al2)
        + (-2*d2*al3 + 2*al1*d2*al3*al2 - 2*al3*d3*al1*al2 - 2*d3*al3)*u,

        (-2*al3*d2 - 2*al3*al1*al2*d2 + 2*al3*d3*al1*al2 - 2*al3*d3)
        + (2*a2*al3*al2 + 2*al3*al1*a1 - 2*al3*al1*a2 - 2*a1*al3*al2 + 2*a3*al1*al2 + 2*a3)*u,

        (-2*d2 - 2*al1*al2*d2 + 2*d3*al1*al2 - 2*d3)
        + (2*al2*a2 + 2*al1*a1 - 2*al1*a2 - 2*al2*a1 - 2*al3*a3*al1*al2 - 2*al3*a3)*u,

        4 - 4*al1*al2,
        -(4*al3*(-1 + al1*al2)),
        sp.Integer(0) + (-(4*al3*(1 + al1*al2)))*u,
        sp.Integer(0) + (-4*al1*al2 - 4)*u,
    ]

###############################################################################
# RIGHT CHAIN hyperplanes: h_v4q
# (GoFa5 right chain uses Tv4 parameterization — see Hyperplane.h logic)
# Each function returns a list of 8 SymPy expressions in v,
# with coefficients linear in t0..t7.
# C++ Polynomial(c0_expr, c1_expr) → c0_expr + c1_expr*v
###############################################################################

def h1_v4q():
    return [
        (-t0*al6*a6 - 2*t4 + t2*al4*d5 - 2*t4*al4*al6 + t3*d4 + t3*d6 + t0*al4*a6
         - t2*al4*d4 - t0*al4*a4 + t0*a5*al5 + t0*al6*a4 + t2*al6*d5 + t2*al6*d4
         + t2*al6*d6 + t2*d6*al4 - t1*a5*al5*al4 - t1*al6*al4*a4 + t1*al6*a5*al5
         + t1*al6*a6*al4 - t3*al6*al4*d5 + t3*al6*al4*d4 - t3*al6*al4*d6
         - t1*a4 + t1*a6 + t3*d5 + t0*al6*a5*al5*al4 + 2*t5*al4 - 2*t5*al6)
        + (t1*al6*d4 - t3*a5*al5 - 2*t6*al4 + 2*t7 + t0*d5 - t3*al6*a5*al5*al4
           - t2*al6*a6*al4 + t2*a4 - t0*al6*al4*d5 + t1*al4*d5 + t1*al6*d5
           + t1*al6*d6 + t3*al4*a4 - t3*al6*a4 + t3*al6*a6 - t3*a6*al4
           + 2*t7*al6*al4 + t2*a5*al5*al4 - t2*al6*a5*al5 + t2*al6*al4*a4
           - t1*al4*d4 + t1*al4*d6 + t0*al6*al4*d4 - t0*al6*al4*d6
           + 2*t6*al6 + t0*d6 + t0*d4 - t2*a6)*v,

        (t1*a5*al5 + t1*al4*a6 - 2*t5 - 2*t4*al4 - t2*d6 - t1*al6*a6
         - t0*al6*a5*al5 + t3*al6*d5 - 2*t5*al4*al6 - t3*al4*d4 + t3*al6*d4
         + t3*al4*d5 + t0*a5*al5*al4 + t0*al6*al4*a4 - t0*al6*a6*al4
         + t2*al6*al4*d5 - t2*al6*al4*d4 - t1*al4*a4 + t2*al6*al4*d6
         + t1*al6*a4 + t3*al6*d6 + t3*d6*al4 + t0*a4 - t0*a6 - t2*d5
         - t2*d4 + 2*t4*al6 + t1*al6*a5*al5*al4)
        + (-t3*al6*a5*al5 - 2*t6 - t2*al6*a6 - t3*a6 + 2*t7*al6 + t2*a5*al5
           + t3*a4 + t1*d5 + t2*al6*a5*al5*al4 - t0*al4*d5 - t0*al6*d5
           - t0*al6*d6 - t2*al4*a4 + t2*al6*a4 + t2*a6*al4 - 2*t6*al6*al4
           - 2*t7*al4 - t1*al6*al4*d5 + t3*a5*al5*al4 + t3*al6*al4*a4
           - t3*al6*a6*al4 + t1*d6 + t1*d4 + t1*al6*al4*d4 - t1*al6*al4*d6
           + t0*al4*d4 - t0*al6*d4 - t0*al4*d6)*v,

        (-t3*al6*a5*al5 - 2*t6 - t2*al6*a6 - t3*a6 + 2*t7*al6 + t2*a5*al5
         + t3*a4 + t1*d5 + t2*al6*a5*al5*al4 - t0*al4*d5 - t0*al6*d5
         - t0*al6*d6 - t2*al4*a4 + t2*al6*a4 + t2*a6*al4 - 2*t6*al6*al4
         - 2*t7*al4 - t1*al6*al4*d5 + t3*a5*al5*al4 + t3*al6*al4*a4
         - t3*al6*a6*al4 + t1*d6 + t1*d4 + t1*al6*al4*d4 - t1*al6*al4*d6
         + t0*al4*d4 - t0*al6*d4 - t0*al4*d6)
        + (-t1*a5*al5 - t1*al4*a6 + 2*t5 + 2*t4*al4 + t2*d6 + t1*al6*a6
           + t0*al6*a5*al5 - t3*al6*d5 + 2*t5*al4*al6 + t3*al4*d4 - t3*al6*d4
           - t3*al4*d5 - t0*a5*al5*al4 - t0*al6*al4*a4 + t0*al6*a6*al4
           - t2*al6*al4*d5 + t2*al6*al4*d4 + t1*al4*a4 - t2*al6*al4*d6
           - t1*al6*a4 - t3*al6*d6 - t3*d6*al4 - t0*a4 + t0*a6 + t2*d5
           + t2*d4 - 2*t4*al6 - t1*al6*a5*al5*al4)*v,

        (-t1*al6*d4 + t3*a5*al5 + 2*t6*al4 - 2*t7 - t0*d5 + t3*al6*a5*al5*al4
         + t2*al6*a6*al4 - t2*a4 + t0*al6*al4*d5 - t1*al4*d5 - t1*al6*d5
         - t1*al6*d6 - t3*al4*a4 + t3*al6*a4 - t3*al6*a6 + t3*a6*al4
         - 2*t7*al6*al4 - t2*a5*al5*al4 + t2*al6*a5*al5 - t2*al6*al4*a4
         + t1*al4*d4 - t1*al4*d6 - t0*al6*al4*d4 + t0*al6*al4*d6
         - 2*t6*al6 - t0*d6 - t0*d4 + t2*a6)
        + (-t0*al6*a6 - 2*t4 + t2*al4*d5 - 2*t4*al4*al6 + t3*d4 + t3*d6
           + t0*al4*a6 - t2*al4*d4 - t0*al4*a4 + t0*a5*al5 + t0*al6*a4
           + t2*al6*d5 + t2*al6*d4 + t2*al6*d6 + t2*d6*al4 - t1*a5*al5*al4
           - t1*al6*al4*a4 + t1*al6*a5*al5 + t1*al6*a6*al4 - t3*al6*al4*d5
           + t3*al6*al4*d4 - t3*al6*al4*d6 - t1*a4 + t1*a6 + t3*d5
           + t0*al6*a5*al5*al4 + 2*t5*al4 - 2*t5*al6)*v,

        (-2*t0 - 2*t0*al6*al4 + 2*t1*al4 - 2*t1*al6)
        + (-2*t2*al4 + 2*t2*al6 + 2*t3 + 2*t3*al6*al4)*v,

        (-2*t1 - 2*t1*al6*al4 - 2*t0*al4 + 2*t0*al6)
        + (-2*t3*al4 + 2*t3*al6 - 2*t2 - 2*t2*al6*al4)*v,

        (-2*t3*al4 + 2*t3*al6 - 2*t2 - 2*t2*al6*al4)
        + (2*t0*al4 - 2*t0*al6 + 2*t1 + 2*t1*al6*al4)*v,

        (-2*t3 - 2*t3*al6*al4 + 2*t2*al4 - 2*t2*al6)
        + (-2*t0 - 2*t0*al6*al4 + 2*t1*al4 - 2*t1*al6)*v,
    ]

def h2_v4q():
    return [
        (2*t4*al5*al6 - t1*a5 - 2*t5*al5 - t3*al5*al4*d5 + t1*al6*al5*a4
         - t1*al6*a6*al5 - t1*al5*al4*a4 + t3*al5*al4*d4 - t2*al6*al5*al4*d5
         + t2*al6*al5*al4*d4 - t0*al6*al5*al4*a6 + t2*al5*d5
         + t2*al6*al5*al4*d6 - t3*al6*al5*d5 + t1*a6*al5*al4
         + t0*al6*al5*al4*a4 - 2*t5*al6*al5*al4 - t3*al6*al5*d4
         + t3*d6*al5*al4 + t0*al6*a5 - t1*al6*a5*al4 + t3*al6*al5*d6
         - 2*t4*al5*al4 - t0*al5*a6 + t2*al5*d4 - t2*al5*d6
         - t0*a5*al4 + t0*al5*a4)
        + (2*t6*al5 + t3*a6*al5 + 2*t7*al5*al4 - t3*al5*a4 - 2*t7*al6*al5
           - t1*al5*d6 + t3*a5*al4 - t0*al6*al5*d5 + t1*al6*d6*al5*al4
           - t3*al6*a5 + t2*al5*al4*a4 - t0*al6*al5*d4 + t0*al5*al4*d6
           - t2*al5*al6*a6 + t2*al5*al4*a6 + 2*t6*al5*al4*al6  # Note: sign in C++: +2*t6*al5*al4*al6
           # Actually from C++: +2*t6*al5*al4*al6 -- I'll trust the C++ here
           + t0*al5*al4*d4 + t2*a5 + t0*al6*al5*d6 - t0*al5*al4*d5
           + t3*al6*a6*al5*al4 - t1*al6*al5*al4*d5 + t1*al5*d5
           + t1*al5*d4 + t2*al6*a5*al4 + t1*al6*al5*al4*d4
           - t3*al6*al5*al4*a4 - t2*al6*al5*a4)*v,

        (2*t4*al5 + t3*al5*d4 + t1*al6*a5 - t1*a6*al5 + t0*a5
         + 2*t5*al6*al5 - t2*al6*al5*d6 - t2*al5*d6*al4 - t3*d6*al5
         + t2*al6*al5*d5 - t3*al6*al5*al4*d5 + t0*al5*al4*a4
         - t0*al6*al5*a4 - t1*a5*al4 + t3*al6*al5*al4*d4 + t0*al6*a5*al4
         + t0*al6*a6*al5 + t1*al5*a4 + t2*al6*al5*d4
         + 2*t4*al5*al6*al4 + t3*al6*al5*al4*d6 - t0*al5*a6*al4
         - 2*t5*al5*al4 + t3*al5*d5 + t1*al6*al5*al4*a4
         + t2*al5*al4*d5 - t2*al5*al4*d4 - t1*al6*al5*al4*a6)
        + (t3*a5 + t1*al6*d6*al5 + t2*al6*a5 + 2*t7*al5 - 2*t6*al5*al4
           + t3*al6*a5*al4 + t3*al5*al4*a4 + 2*t6*al6*al5
           - t2*al5*al6*a6*al4 - t0*al5*d4 + t0*al5*d6
           + t2*al6*al5*al4*a4 - t1*al5*al4*d5 + t0*al6*al5*al4*d5
           - t3*al6*al5*a4 - t0*al6*al5*d6*al4 - t1*al6*al5*d5
           - t2*a6*al5 - t0*al6*al5*al4*d4 + t3*al5*al6*a6
           - t3*al5*al4*a6 + t1*al5*al4*d6 + 2*t7*al5*al4*al6
           - t1*al5*al6*d4 + t1*al5*al4*d4 - t0*al5*d5
           - t2*a5*al4 + t2*al5*a4)*v,

        (t3*a5 + t1*al6*d6*al5 + t2*al6*a5 + 2*t7*al5 - 2*t6*al5*al4
         + t3*al6*a5*al4 + t3*al5*al4*a4 + 2*t6*al6*al5
         - t2*al5*al6*a6*al4 - t0*al5*d4 + t0*al5*d6
         + t2*al6*al5*al4*a4 - t1*al5*al4*d5 + t0*al6*al5*al4*d5
         - t3*al6*al5*a4 - t0*al6*al5*d6*al4 - t1*al6*al5*d5
         - t2*a6*al5 - t0*al6*al5*al4*d4 + t3*al5*al6*a6
         - t3*al5*al4*a6 + t1*al5*al4*d6 + 2*t7*al5*al4*al6
         - t1*al5*al6*d4 + t1*al5*al4*d4 - t0*al5*d5
         - t2*a5*al4 + t2*al5*a4)
        + (-2*t4*al5 - t3*al5*d4 - t1*al6*a5 + t1*a6*al5 - t0*a5
           - 2*t5*al6*al5 + t2*al6*al5*d6 + t2*al5*d6*al4 + t3*d6*al5
           - t2*al6*al5*d5 + t3*al6*al5*al4*d5 - t0*al5*al4*a4
           + t0*al6*al5*a4 + t1*a5*al4 - t3*al6*al5*al4*d4 - t0*al6*a5*al4
           - t0*al6*a6*al5 - t1*al5*a4 - t2*al6*al5*d4
           - 2*t4*al5*al6*al4 - t3*al6*al5*al4*d6 + t0*al5*a6*al4
           + 2*t5*al5*al4 - t3*al5*d5 - t1*al6*al5*al4*a4
           - t2*al5*al4*d5 + t2*al5*al4*d4 + t1*al6*al5*al4*a6)*v,

        (-2*t6*al5 - t3*a6*al5 - 2*t7*al5*al4 + t3*al5*a4 + 2*t7*al6*al5
         + t1*al5*d6 - t3*a5*al4 + t0*al6*al5*d5 - t1*al6*d6*al5*al4
         + t3*al6*a5 - t2*al5*al4*a4 + t0*al6*al5*d4 - t0*al5*al4*d6
         - t2*al5*al6*a6 + t2*al5*al4*a6 - 2*t6*al5*al4*al6
         - t0*al5*al4*d4 - t2*a5 - t0*al6*al5*d6 + t0*al5*al4*d5
         - t3*al6*a6*al5*al4 + t1*al6*al5*al4*d5 - t1*al5*d5
         - t1*al5*d4 - t2*al6*a5*al4 - t1*al6*al5*al4*d4
         + t3*al6*al5*al4*a4 + t2*al6*al5*a4)
        + (2*t4*al5*al6 - t1*a5 - 2*t5*al5 - t3*al5*al4*d5 + t1*al6*al5*a4
           - t1*al6*a6*al5 - t1*al5*al4*a4 + t3*al5*al4*d4
           - t2*al6*al5*al4*d5 + t2*al6*al5*al4*d4 - t0*al6*al5*al4*a6
           + t2*al5*d5 + t2*al6*al5*al4*d6 - t3*al6*al5*d5
           + t1*a6*al5*al4 + t0*al6*al5*al4*a4 - 2*t5*al6*al5*al4
           - t3*al6*al5*d4 + t3*d6*al5*al4 + t0*al6*a5 - t1*al6*a5*al4
           + t3*al6*al5*d6 - 2*t4*al5*al4 - t0*al5*a6 + t2*al5*d4
           - t2*al5*d6 - t0*a5*al4 + t0*al5*a4)*v,

        (-2*t0*al5*al4 + 2*t0*al6*al5 - 2*t1*al5 - 2*t1*al6*al5*al4)
        + (2*t2*al5 + 2*t2*al6*al5*al4 + 2*t3*al5*al4 - 2*t3*al6*al5)*v,

        (-2*t1*al5*al4 + 2*t1*al6*al5 + 2*t0*al5 + 2*t0*al6*al5*al4)
        + (2*t3*al5 + 2*t3*al6*al5*al4 - 2*t2*al5*al4 + 2*t2*al6*al5)*v,

        (2*t3*al5 + 2*t3*al6*al5*al4 - 2*t2*al5*al4 + 2*t2*al6*al5)
        + (-2*t0*al5 - 2*t0*al6*al5*al4 + 2*t1*al5*al4 - 2*t1*al6*al5)*v,

        (-2*t3*al5*al4 + 2*t3*al6*al5 - 2*t2*al5 - 2*t2*al6*al5*al4)
        + (-2*t0*al5*al4 + 2*t0*al6*al5 - 2*t1*al5 - 2*t1*al6*al5*al4)*v,
    ]

def h3_v4q():
    return [
        (-2*t6*al5 - t3*a6*al5 + 2*t7*al5*al4 - t3*al5*a4 + 2*t7*al6*al5
         + t1*al5*d6 + t3*a5*al4 + t0*al6*al5*d5 + t1*al6*d6*al5*al4
         + t3*al6*a5 - t2*al5*al4*a4 + t0*al6*al5*d4 + t0*al5*al4*d6
         - t2*al5*al6*a6 - t2*al5*al4*a6 + 2*t6*al5*al4*al6
         + t0*al5*al4*d4 - t2*a5 - t0*al6*al5*d6 - t0*al5*al4*d5
         + t3*al6*a6*al5*al4 - t1*al6*al5*al4*d5 - t1*al5*d5
         - t1*al5*d4 + t2*al6*a5*al4 + t1*al6*al5*al4*d4
         + t3*al6*al5*al4*a4 - t2*al6*al5*a4)
        + (2*t4*al5*al6 - t1*a5 - 2*t5*al5 + t3*al5*al4*d5 - t1*al6*al5*a4
           - t1*al6*a6*al5 - t1*al5*al4*a4 - t3*al5*al4*d4
           + t2*al6*al5*al4*d5 - t2*al6*al5*al4*d4 + t0*al6*al5*al4*a6
           + t2*al5*d5 - t2*al6*al5*al4*d6 - t3*al6*al5*d5
           - t1*a6*al5*al4 + t0*al6*al5*al4*a4 + 2*t5*al6*al5*al4
           - t3*al6*al5*d4 - t3*d6*al5*al4 + t0*al6*a5 + t1*al6*a5*al4
           + t3*al6*al5*d6 + 2*t4*al5*al4 - t0*al5*a6 + t2*al5*d4
           - t2*al5*d6 + t0*a5*al4 - t0*al5*a4)*v,

        (-t3*a5 - t1*al6*d6*al5 - t2*al6*a5 - 2*t7*al5 - 2*t6*al5*al4
         + t3*al6*a5*al4 - t3*al5*al4*a4 - 2*t6*al6*al5
         - t2*al5*al6*a6*al4 + t0*al5*d4 - t0*al5*d6
         - t2*al6*al5*al4*a4 - t1*al5*al4*d5 + t0*al6*al5*al4*d5
         - t3*al6*al5*a4 - t0*al6*al5*d6*al4 + t1*al6*al5*d5
         + t2*a6*al5 - t0*al6*al5*al4*d4 - t3*al5*al6*a6
         - t3*al5*al4*a6 + t1*al5*al4*d6 + 2*t7*al5*al4*al6
         + t1*al5*al6*d4 + t1*al5*al4*d4 + t0*al5*d5 - t2*a5*al4 + t2*al5*a4)
        + (2*t4*al5 + t3*al5*d4 + t1*al6*a5 - t1*a6*al5 + t0*a5
           + 2*t5*al6*al5 - t2*al6*al5*d6 + t2*al5*d6*al4 - t3*d6*al5
           + t2*al6*al5*d5 + t3*al6*al5*al4*d5 + t0*al5*al4*a4
           + t0*al6*al5*a4 + t1*a5*al4 - t3*al6*al5*al4*d4 - t0*al6*a5*al4
           + t0*al6*a6*al5 - t1*al5*a4 + t2*al6*al5*d4 - 2*t4*al5*al6*al4
           - t3*al6*al5*al4*d6 + t0*al5*a6*al4 + 2*t5*al5*al4 + t3*al5*d5
           + t1*al6*al5*al4*a4 - t2*al5*al4*d5 + t2*al5*al4*d4
           + t1*al6*al5*al4*a6)*v,

        (t3*al5*d4 + 2*t4*al5 + t2*al5*al4*d4 + t1*a5*al4 + t1*al6*al5*al4*a4
         - t1*al5*a6 + t3*al6*al5*al4*d5 + t1*al6*al5*al4*a6 - t1*al5*a4
         + t2*d6*al5*al4 + t1*al6*a5 - t2*al5*al4*d5 + t2*al6*al5*d5
         - t3*al6*al5*al4*d6 + 2*t5*al6*al5 + t0*a5 - t3*d6*al5
         + t0*al5*al4*a4 + t0*al6*al5*a4 - t3*al6*al5*al4*d4
         - t0*al6*a5*al4 + t0*al6*a6*al5 + t2*al6*al5*d4
         - t2*al6*al5*d6 - 2*t4*al6*al5*al4 + t0*al5*a6*al4
         + 2*t5*al5*al4 + t3*al5*d5)
        + (t3*a5 + t1*al6*d6*al5 + t2*al6*a5 + 2*t7*al5 + 2*t6*al5*al4
           - t3*al6*a5*al4 + t3*al5*al4*a4 + 2*t6*al6*al5
           + t2*al5*al6*a6*al4 - t0*al5*d4 + t0*al5*d6
           + t2*al6*al5*al4*a4 + t1*al5*al4*d5 - t0*al6*al5*al4*d5
           + t3*al6*al5*a4 + t0*al6*al5*d6*al4 - t1*al6*al5*d5
           - t2*a6*al5 + t0*al6*al5*al4*d4 + t3*al5*al6*a6 + t3*al5*al4*a6
           - t1*al5*al4*d6 - 2*t7*al5*al4*al6 - t1*al5*al6*d4
           - t1*al5*al4*d4 - t0*al5*d5 + t2*a5*al4 - t2*al5*a4)*v,

        (-t0*al6*a5 - t0*al6*al5*al4*a4 + t1*a5 + t2*d6*al5 + t1*al5*al4*a4
         - t2*al5*d5 - 2*t5*al6*al5*al4 + 2*t5*al5 + t2*al6*al5*al4*d6
         + t2*al6*al5*al4*d4 + t3*al6*al5*d4 + t1*al5*a6*al4
         + t0*al5*a4 - t3*al5*al4*d5 - t2*al6*al5*al4*d5
         - t1*al6*a5*al4 + t3*al5*al4*d4 + t1*al6*al5*a4
         + t1*al6*a6*al5 - t2*al5*d4 + t3*d6*al5*al4 + t3*al6*al5*d5
         - 2*t4*al5*al4 - t0*a5*al4 + t0*al5*a6 - t3*al6*al5*d6
         - 2*t4*al6*al5 - t0*al6*al5*al4*a6)
        + (2*t7*al6*al5 + 2*t7*al5*al4 - t2*a5 - t2*al6*a6*al5
           + t1*al6*al5*al4*d4 + t3*al6*a5 - t0*al6*d6*al5 - 2*t6*al5
           - t1*al6*al5*al4*d5 + t3*a5*al4 - t1*al5*d5 + t1*al6*d6*al5*al4
           + t0*al5*al4*d4 + t0*al6*al5*d4 - t2*al5*al4*a6
           + t0*al5*al4*d6 + 2*t6*al5*al4*al6 - t2*al6*al5*a4
           + t2*al6*a5*al4 + t3*al5*al6*a6*al4 - t2*al5*al4*a4
           - t3*a6*al5 + t3*al6*al5*al4*a4 - t3*al5*a4
           - t0*al5*al4*d5 - t1*al5*d4 + t1*al5*d6 + t0*al6*al5*d5)*v,

        (-2*t2*al5 + 2*t2*al6*al5*al4 + 2*t3*al5*al4 + 2*t3*al6*al5)
        + (2*t0*al5*al4 + 2*t0*al6*al5 - 2*t1*al5 + 2*t1*al6*al5*al4)*v,

        (-2*t3*al5 + 2*t3*al6*al5*al4 - 2*t2*al5*al4 - 2*t2*al6*al5)
        + (2*t1*al5*al4 + 2*t1*al6*al5 + 2*t0*al5 - 2*t0*al6*al5*al4)*v,

        (2*t1*al5*al4 + 2*t1*al6*al5 + 2*t0*al5 - 2*t0*al6*al5*al4)
        + (2*t2*al5*al4 + 2*t2*al6*al5 + 2*t3*al5 - 2*t3*al6*al5*al4)*v,

        (2*t1*al5 - 2*t1*al6*al5*al4 - 2*t0*al5*al4 - 2*t0*al6*al5)
        + (-2*t2*al5 + 2*t2*al6*al5*al4 + 2*t3*al5*al4 + 2*t3*al6*al5)*v,
    ]

def h4_v4q():
    return [
        (2*t6*al6 + 2*t7 - t2*a6 + t1*al6*d6 + t0*d4 - t1*al4*d5 + t3*al6*a6
         - t2*a4 + t3*al6*a4 + t1*al4*d4 + 2*t6*al4 + t0*d6 + t3*al4*a6
         + t0*d5 + t0*al6*al4*d6 - t1*d6*al4 + t1*al6*d5 - t3*a5*al5
         + t0*al6*al4*d5 - 2*t7*al4*al6 + t1*al6*d4 - t0*al6*al4*d4
         + t2*al6*a6*al4 + t2*al6*al4*a4 - t2*al6*a5*al5 + t3*al4*a4
         + t3*al6*a5*al5*al4 - t2*a5*al5*al4)
        + (t2*al4*d5 - 2*t4*al6*al4 + t0*al4*a4 - t1*a6 + 2*t4 + t0*al6*a4
           - t0*a5*al5 - t3*d4 - t3*d5 - t3*d6 + t3*al6*al4*d4
           - t3*al6*al4*d6 - t2*al6*d5 + 2*t5*al6 + t1*al6*a6*al4
           + t2*al4*d6 - t2*al6*d4 - t2*al4*d4 + t1*al6*al4*a4
           - t1*a5*al5*al4 - t1*a4 + 2*t5*al4 + t0*al6*a5*al5*al4
           + t0*al6*a6 - t3*al6*al4*d5 + t0*a6*al4
           - t1*al6*a5*al5 - t2*al6*d6)*v,

        (-t2*al6*a6 + 2*t7*al6 - t0*al6*d6 - 2*t6 - t2*al6*a4 - t2*al4*a6
         - t0*al6*d4 - t2*al4*a4 + t1*al6*al4*d5 - t3*a6 + t2*a5*al5
         + t1*d5 - t1*al6*al4*d4 + t0*d6*al4 + t1*d4 - t0*al4*d4
         + 2*t6*al4*al6 - t0*al6*d5 - t3*a4 + t1*d6
         - t2*al6*a5*al5*al4 + t1*al6*al4*d6 - t3*a5*al5*al4
         - 2*t7*al4 - t3*al6*a5*al5 + t3*al6*al4*a4
         + t0*al4*d5 + t3*al6*a6*al4)
        + (2*t5 + t2*d6 - 2*t4*al6 + t1*al6*a5*al5*al4 - t3*al6*d4
           - t3*al4*d4 + t3*al4*d6 + t2*d5 + t0*a6 + t0*al6*a5*al5
           + t0*a5*al5*al4 - t3*al6*d5 - t3*al6*d6 + t1*al4*a4
           + t1*al6*a4 + t2*d4 - 2*t5*al6*al4 + t3*al4*d5
           - t0*al6*al4*a4 - t0*al6*a6*al4 - 2*t4*al4 + t0*a4
           + t1*al6*a6 + t2*al6*al4*d5 + t1*a6*al4 + t2*al6*al4*d6
           - t2*al6*al4*d4 - t1*a5*al5)*v,

        (2*t5 + t2*d6 - 2*t4*al6 + t1*al6*a5*al5*al4 - t3*al6*d4
         - t3*al4*d4 + t3*al4*d6 + t2*d5 + t0*a6 + t0*al6*a5*al5
         + t0*a5*al5*al4 - t3*al6*d5 - t3*al6*d6 + t1*al4*a4
         + t1*al6*a4 + t2*d4 - 2*t5*al6*al4 + t3*al4*d5
         - t0*al6*al4*a4 - t0*al6*a6*al4 - 2*t4*al4 + t0*a4
         + t1*al6*a6 + t2*al6*al4*d5 + t1*a6*al4 + t2*al6*al4*d6
         - t2*al6*al4*d4 - t1*a5*al5)
        + (t2*al6*a6 - 2*t7*al6 + t0*al6*d6 + 2*t6 + t2*al6*a4 + t2*al4*a6
           + t0*al6*d4 + t2*al4*a4 - t1*al6*al4*d5 + t3*a6 - t2*a5*al5
           - t1*d5 + t1*al6*al4*d4 - t0*d6*al4 - t1*d4 + t0*al4*d4
           - 2*t6*al4*al6 + t0*al6*d5 + t3*a4 - t1*d6 + t2*al6*a5*al5*al4
           - t1*al6*al4*d6 + t3*a5*al5*al4 + 2*t7*al4 + t3*al6*a5*al5
           - t3*al6*al4*a4 - t0*al4*d5 - t3*al6*a6*al4)*v,

        (-t2*al4*d5 + 2*t4*al6*al4 - t0*al4*a4 + t1*a6 - 2*t4 - t0*al6*a4
         + t0*a5*al5 + t3*d4 + t3*d5 + t3*d6 - t3*al6*al4*d4
         + t3*al6*al4*d6 + t2*al6*d5 - 2*t5*al6 - t1*al6*a6*al4
         - t2*al4*d6 + t2*al6*d4 + t2*al4*d4 - t1*al6*al4*a4
         + t1*a5*al5*al4 + t1*a4 - 2*t5*al4 - t0*al6*a5*al5*al4
         - t0*al6*a6 + t3*al6*al4*d5 - t0*a6*al4
         + t1*al6*a5*al5 + t2*al6*d6)
        + (2*t6*al6 + 2*t7 - t2*a6 + t1*al6*d6 + t0*d4 - t1*al4*d5 + t3*al6*a6
           - t2*a4 + t3*al6*a4 + t1*al4*d4 + 2*t6*al4 + t0*d6 + t3*al4*a6
           + t0*d5 + t0*al6*al4*d6 - t1*d6*al4 + t1*al6*d5 - t3*a5*al5
           + t0*al6*al4*d5 - 2*t7*al4*al6 + t1*al6*d4 - t0*al6*al4*d4
           + t2*al6*a6*al4 + t2*al6*al4*a4 - t2*al6*a5*al5 + t3*al4*a4
           + t3*al6*a5*al5*al4 - t2*a5*al5*al4)*v,

        (2*t2*al4 + 2*t2*al6 + 2*t3 - 2*t3*al6*al4)
        + (2*t0 - 2*t0*al6*al4 + 2*t1*al4 + 2*t1*al6)*v,

        (2*t3*al4 + 2*t3*al6 - 2*t2 + 2*t2*al6*al4)
        + (2*t1 - 2*t1*al6*al4 - 2*t0*al4 - 2*t0*al6)*v,

        (2*t1 - 2*t1*al6*al4 - 2*t0*al4 - 2*t0*al6)
        + (2*t2 - 2*t2*al6*al4 - 2*t3*al4 - 2*t3*al6)*v,

        (-2*t1*al4 - 2*t1*al6 - 2*t0 + 2*t0*al6*al4)
        + (2*t2*al4 + 2*t2*al6 + 2*t3 - 2*t3*al6*al4)*v,
    ]

###############################################################################
# Build the 8x8 hyperplane matrix H[i][j]
# Row 0-3: left chain (h1_tc_v2 .. h4_tc_v2)
# Row 4-7: right chain (h1_v4q .. h4_v4q)
# Column j=0: free term (x0 coefficient)
# Column j=1..7: coefficients of x1..x7
###############################################################################

def build_hyperplane_matrix():
    left  = [h1_tc_v2(), h2_tc_v2(), h3_tc_v2(), h4_tc_v2()]
    right = [h1_v4q(),   h2_v4q(),   h3_v4q(),   h4_v4q()  ]
    H = left + right  # 8 rows, each with 8 entries
    return H

###############################################################################
# Kinematic surface r[0..7] via Cramer's rule
# A[i][j] = H[i][j+1]  (7x7 matrix, x0 = 1)
# b[i]    = -H[i][0]   (RHS vector)
# r[k] = det(A with column k-1 replaced by b)  for k = 1..7
# r[0] = det(A)  (common denominator)
###############################################################################

def compute_kinematic_surface(H):
    """
    Returns r[0..7] as SymPy expressions in u, v, t0..t7.
    r[0] = det(A)
    r[1..7] = Cramer numerators for x1..x7
    """
    # Build 7x7 matrix A and RHS b
    n = 7
    A = [[sp.expand(H[i][j+1]) for j in range(n)] for i in range(n)]
    b = [sp.expand(-H[i][0]) for i in range(n)]

    r = [None] * 8

    print("  Computing det(A) [r0] ...")
    Amat = Matrix(A)
    r[0] = sp.expand(Amat.det())
    print(f"    r[0]: degree_u={sp.degree(r[0], u)}, degree_v={sp.degree(r[0], v)}")

    for k in range(1, 8):
        print(f"  Computing Cramer det for r[{k}] ...")
        Ak = [row[:] for row in A]      # copy
        for i in range(n):
            Ak[i][k-1] = b[i]          # replace column k-1 with b
        r[k] = sp.expand(Matrix(Ak).det())
        print(f"    r[{k}]: degree_u={sp.degree(r[k], u)}, degree_v={sp.degree(r[k], v)}")

    return r

###############################################################################
# Build e1 and e2, then compute the resultant
###############################################################################

def build_e1(r):
    """Study quadric: r[0]*r[4] + r[1]*r[5] + r[2]*r[6] + r[3]*r[7]"""
    e = (r[0]*r[4] + r[1]*r[5] + r[2]*r[6] + r[3]*r[7])
    return sp.expand(e)

def build_e2(H, r):
    """8th hyperplane h[7] dotted with r[0..7]"""
    e = sum(H[7][j] * r[j] for j in range(8))
    return sp.expand(e)

###############################################################################
# Main pipeline
###############################################################################

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--skip-to', default='start',
                        choices=['start', 'r', 'e1e2', 'resultant', 'factor'],
                        help='Skip to a checkpoint stage')
    args = parser.parse_args()

    print("=" * 60)
    print("IKDH Symbolic Polynomial Reduction — GoFa5")
    print("=" * 60)

    # ── Stage 1: Hyperplane matrix ────────────────────────────
    print("\n[1] Building hyperplane matrix ...")
    H = build_hyperplane_matrix()
    # Evaluate numeric DH values (substitution is already done by Python eval)
    # Simplify each entry
    H = [[sp.nsimplify(sp.expand(e), rational=True) for e in row] for row in H]
    save("H", H)
    print("  Done.")

    # ── Stage 2: Kinematic surface r[0..7] ───────────────────
    r = load("r")
    if r is None or args.skip_to == 'start':
        r = timed("Kinematic surface r[0..7]", lambda: compute_kinematic_surface(H))
        save("r", r)
    else:
        print("\n[2] Loaded r[] from checkpoint.")

    # ── Stage 3: e1 and e2 ───────────────────────────────────
    e1e2 = load("e1e2")
    if e1e2 is None or args.skip_to in ('start', 'r'):
        e1 = timed("e1 (Study quadric)", lambda: build_e1(r))
        e2 = timed("e2 (8th hyperplane)", lambda: build_e2(H, r))
        print(f"  e1: deg_u={sp.degree(e1, u)}, deg_v={sp.degree(e1, v)}")
        print(f"  e2: deg_u={sp.degree(e2, u)}, deg_v={sp.degree(e2, v)}")
        e1e2 = (e1, e2)
        save("e1e2", e1e2)
    else:
        print("\n[3] Loaded e1, e2 from checkpoint.")
        e1, e2 = e1e2

    # ── Stage 4: Resultant in u ───────────────────────────────
    real_pol = load("real_pol")
    if real_pol is None or args.skip_to in ('start', 'r', 'e1e2'):
        print("\n[4] Computing resultant(e1, e2, u) ...")
        print("    (this is the expensive step — may take hours)")
        t_start = time.time()
        # Use Poly for efficiency
        gens = (u, v, t0, t1, t2, t3, t4, t5, t6, t7)
        p1 = Poly(e1, *gens, domain='QQ')
        p2 = Poly(e2, *gens, domain='QQ')
        real_pol = p1.resultant(p2)   # resultant in u (first generator)
        print(f"    Done in {time.time()-t_start:.1f}s")
        print(f"    Degree in v: {real_pol.degree(v)}")
        save("real_pol", real_pol)
    else:
        print("\n[4] Loaded resultant from checkpoint.")

    # ── Stage 5: Factor ───────────────────────────────────────
    print("\n[5] Factoring the resultant polynomial ...")
    t_start = time.time()
    factors = factor_list(real_pol.as_expr(), v, t0, t1, t2, t3, t4, t5, t6, t7)
    print(f"    Done in {time.time()-t_start:.1f}s")
    print(f"    Number of factors: {len(factors[1])}")
    for i, (fac, mult) in enumerate(factors[1]):
        deg = sp.degree(fac, v)
        print(f"    Factor {i}: degree_v={deg}, multiplicity={mult}")
    save("factors", factors)

    # ── Stage 6: Extract degree-16 factor and export ─────────
    print("\n[6] Extracting degree-16 factor ...")
    poly16 = None
    for fac, mult in factors[1]:
        if sp.degree(fac, v) == 16:
            poly16 = fac
            break

    if poly16 is None:
        print("  WARNING: No degree-16 factor found! Degrees present:")
        for fac, mult in factors[1]:
            print(f"    degree={sp.degree(fac, v)}")
    else:
        print("  Found degree-16 factor!")
        # Extract coefficients as functions of t0..t7
        coeffs = sp.Poly(poly16, v).all_coeffs()  # [c16, c15, ..., c0]
        print(f"  Coefficients (c16..c0) simplified:")
        coeff_strings = []
        for i, c in enumerate(coeffs):
            cs = str(sp.simplify(c))
            coeff_strings.append(cs)
            print(f"    c[{16-i}] = {cs[:120]}{'...' if len(cs)>120 else ''}")

        # Save human-readable output
        out_path = os.path.join(os.path.dirname(__file__), "poly16_coeffs.txt")
        with open(out_path, "w") as f:
            f.write("# Degree-16 IK polynomial coefficients for GoFa5\n")
            f.write("# p(v) = sum_{k=0}^{16} c[k] * v^k\n")
            f.write("# where v = tan(theta4/2) and t0..t7 are Study parameters of TCP\n\n")
            for i, c in enumerate(reversed(coeffs)):
                f.write(f"c[{i}] = {c}\n\n")
        print(f"\n  Coefficients saved to {out_path}")

        save("poly16", poly16)
        save("poly16_coeffs", coeffs)

    print("\nDone!")

if __name__ == "__main__":
    main()