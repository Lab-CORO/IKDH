"""
_format.py  -  Shared radian formatting for RoboDK <-> YAML conversion scripts.

Not a public API; imported directly by sibling scripts run as
`python3 robodk/<script>.py`.
"""

import math

# (value, canonical string) pairs, restricted to the "pi", "pi/N" forms that
# include/robots.h's parsePiExpr() can read back  -  no "N*pi/D" multiples.
# Angles are folded to their (-pi, pi] representative (e.g. 3*pi/2 -> -pi/2)
# rather than spelled out with a numerator, since they're equivalent rotations.
_PI_FRACS = [
    (math.pi,     "pi"),    (-math.pi,     "-pi"),
    (math.pi / 2, "pi/2"),  (-math.pi / 2, "-pi/2"),
    (math.pi / 3, "pi/3"),  (-math.pi / 3, "-pi/3"),
    (math.pi / 4, "pi/4"),  (-math.pi / 4, "-pi/4"),
    (0.0, "0"),
]


def fmt_angle(rad, tol=1e-6):
    """Format a radian value as a compact pi-fraction when possible, else %.8g."""
    folded = math.remainder(rad, 2 * math.pi)  # fold into (-pi, pi]
    for val, s in _PI_FRACS:
        if abs(folded - val) < tol:
            return s
    return f"{rad:.8g}"
