# IKDH — Inverse Kinematics via Denavit-Hartenberg

A C++ solver for the closed-form inverse kinematics of general 6-DOF serial robots,
based on the algebraic method of Husty and Pfurner.

## Algorithm

The solver uses the **Study-quadric / dual-quaternion formulation** introduced by
Husty (2007) and systematically extended by Capco, Quam, Pfurner, Schröcker and
Sinn (2021). The end-effector constraint is encoded as a system of eight quadratic
equations on the Study parameters; the system is reduced to a univariate degree-16
polynomial whose real roots yield all joint-angle solutions.

### References

- M. Husty, M. Pfurner, H.-P. Schröcker, K. Brunnthaler,
  *Algebraic methods in mechanism analysis and synthesis*,
  Robotica **25** (2007), 661–675.
- J. Capco, M. Quam, J. Pfurner, H.-P. Schröcker, J. Sinn,
  *Robots, computer algebra and eight connected components*,
  Proc. ISSAC 2021.

The HuPf solver source (``src/hupf/``) is redistributed with permission from the
original authors and integrated here without modification except for bug fixes
documented in the commit history.

### Extensions to the original HuPf method

The original HuPf solver is an exact algebraic method, but two structural
limitations prevent it from returning all solutions for certain robot
configurations. This implementation adds a numerical post-processing layer that
recovers the missing roots.

**Problem 1 — Rank-deficient kinematic surface.**
When the target pose has a special orientation (e.g. pure rotation about a single
world axis), the 7 × 7 matrix that defines the kinematic surface becomes
rank-deficient. The resulting univariate polynomial degenerates and some real roots
collapse or disappear. This is not a bug in the implementation; it is an intrinsic
algebraic singularity of the Study-quadric encoding for that pose.

**Problem 2 — Roots at infinity.**
The HuPf right-chain parametrisation uses the substitution u₄ = 1/tan(θ₄/2).
A joint angle θ₄ = 0° maps to u₄ = ∞, which is not a root of any finite
polynomial. Up to four of the eight solutions of a 6R robot can involve θ₄ = 0°
and are therefore inaccessible to the algebraic solver.

**Fix — Perturbation + Newton-Raphson refinement.**
The solver wraps the HuPf core with the following strategy:

1. *Primary solve.* The HuPf polynomial is solved for the target pose and all
   returned roots are Newton-refined (damped least-squares, Tikhonov λ = 10⁻⁶)
   to correct the loss of precision near degenerate poses.

2. *Perturbation sweep.* Small rotations ε ∈ {1°, 5°, 10°} are applied to the
   target in twelve independent directions — six *right-multiplied* on the full
   4 × 4 matrix (local tool rotations) and six applied to the 3 × 3 rotation
   block only, with the position unchanged (world-frame orientation change
   equivalent to an RPW Euler-angle perturbation). The perturbed poses are
   non-degenerate and their HuPf roots are well-defined, including roots near
   θ₄ ≈ 0°.

3. *Shoulder and wrist flips.* For each seed returned by the perturbed solve,
   three additional starting points are formed by shifting J1 by +180°
   (shoulder flip), J4 by +180° (wrist flip), and both simultaneously. For a
   6R serial robot the eight solutions split into groups related by exactly
   these ±180° symmetries, so each flip exposes a distinct solution basin.

4. *Newton refinement with angle normalisation.* Every candidate is refined back
   onto the original target by Newton-Raphson iteration. Joint angles are wrapped
   to (−180°, 180°] after each step to prevent winding (accumulation of
   multiples of 360° that otherwise require hundreds of extra iterations).

5. *Deduplication.* Converged solutions are accepted only if their circular
   distance (modulo 360°) from every already-found solution exceeds a threshold.

The combined procedure reliably recovers all eight solutions, including those
at algebraic singularities of the Study-quadric formulation, with Frobenius
FK errors below 10⁻⁹ for all solutions.

### Limitations

The algorithm requires the robot to have at least one non-zero link length among
joints 4–5 (parameters *a*₄ or *a*₅ in 1-indexed DH notation). Robots with a
**strict spherical wrist** (*a*₄ = *a*₅ = 0, *d*₅ = 0) are not supported; the
Pieper decomposition is more appropriate for that class.

## Build

Requires CMake ≥ 3.14 and a C++11-capable compiler.

```bash
cmake -B build
cmake --build build
```

The build produces:
- `build/demo` — demonstration executable
- `build/libikdh.a` — static library for linking into other projects

## Usage

### C++ API

```cpp
#include <ikdh.h>

// 1. Define the robot
IKDH::DHTable dh;
// a, d in meters; alpha, theta in radians
dh.a     = { 0.0,  0.444,  0.110,  0.0,   0.080,  0.0  };
dh.d     = { 0.265, 0.0,   0.0,    0.470,  0.0,   0.101};
dh.alpha = {-M_PI/2, 0.0, -M_PI/2, M_PI/2,-M_PI/2, 0.0 };
dh.theta = { 0.0, -M_PI/2, 0.0,    0.0,    0.0,    M_PI};
for (int i = 0; i < 6; ++i) dh.revolute[i] = true;

// 2. Specify joint limits (degrees)
IKDH::JointLimits limits({
    {-180.0,  180.0},
    {-180.0,  180.0},
    {-225.0,   85.0},
    {-180.0,  180.0},
    {-180.0,  180.0},
    {-180.0,  180.0},
});

// 3. Solve
IKDH::Solver    solver(dh, limits);
IKDH::Transform ee   = IKDH::forwardKin(dh, {10.0, -30.0, 60.0, -20.0, 45.0, 15.0});
auto            sols = solver.solve(ee);

for (auto& q : sols) {
    for (double v : q) printf("  %.3f", v);
    printf("\n");
}
```

Solutions are returned in degrees, guaranteed to lie within the specified joint limits.
Angles are mapped into the limits by ±360° equivalence; configurations that cannot be
mapped are discarded.

### DH convention

```
T_i = Rz(θ_i + θ_offset_i) · Tz(d_i) · Tx(a_i) · Rx(α_i)
```

### RoboDK verification

`tools/robodk_verify.py` connects to a running RoboDK instance, sends each IK
solution to the robot model and reports the Cartesian position error.

**Prerequisites**

```bash
pip install robodk
cmake --build build   # build/demo must exist
```

**Steps**

1. Open RoboDK with your robot loaded in the scene.
2. Edit the configuration block at the top of the script:

```python
ROBOT_LABEL  = 'ABB GoFa CRB 15000-5'  # must match the name in demo.cpp
ROBOT_RDK    = 'ABB CRB 15000'          # exact name in the RoboDK scene
JOINT_LIMITS = [(-180, 180), (-180, 180), (-225, 85), ...]  # degrees, per joint
DEMO_BIN     = 'build/demo'
INTERVAL     = 5.0                       # seconds between configurations
```

3. Run from the repository root:

```bash
python3 tools/robodk_verify.py
```

The script executes `build/demo`, extracts the solutions for the selected robot,
normalises angles into the joint limits, and streams each configuration to RoboDK
with a pause of `INTERVAL` seconds. The internal FK error and the Cartesian
position deviation measured by RoboDK (Δ) are printed for each solution.

## Repository structure

```
include/
  ikdh.h              public API header
src/
  ikdh.cpp            solver wrapper with Newton post-processing
  hupf/               HuPf algebraic IK solver (Husty–Pfurner)
examples/
  robots.h            DH tables and joint limits for common robots
  demo.cpp            round-trip IK demonstration
  diag.cpp            diagnostic tool (solution-basin analysis, GoFa5)
  poly.cpp            characteristic-polynomial inspector
tools/
  robodk_dh.py        extract Standard DH parameters from a RoboDK robot
  robodk_verify.py    stream IK solutions to RoboDK and verify FK consistency
```
