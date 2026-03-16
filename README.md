# IKDH — Inverse Kinematics via Denavit-Hartenberg

IK solver for general 6R serial robots based on the algebraic method of Husty and Pfurner.
Returns **all solutions** within joint limits for any reachable pose.

---

## Quick start (Python)

```bash
cmake -B build -S . && cmake --build build
```

```python
import sys; sys.path.insert(0, 'build')
import ikdh

robot  = ikdh.load_robot("robots/gofa5.yaml")
solver = ikdh.Solver(robot.dh, robot.limits)

ee   = ikdh.pose_from_xyzrpw(500.0, 0.0, 500.0, 0.0, 90.0, 0.0)  # x y z mm, rx ry rz deg
sols = solver.solve(ee)   # list of (6,) numpy arrays, in degrees

for q in sols:
    print(q)
```

Pre-defined robots: `ur5e`, `gofa5`, `gofa12`, `crb15000_10`, `fanuc_crx_5ia`, `fanuc_crx_10ia`, `doosan_a0509`.

---

## Path planning — warm-start IK

For trajectory tracking, use `solve_from_seed()` instead of `solve()`.
It refines the previous solution directly via Newton-Raphson, bypassing the algebraic solver (~100× faster).

```python
poses = interpolate(pose_A, pose_B, n=100)

# Seed all branches from the first pose
branches = [[q] for q in solver.solve(poses[0])]

# Track each branch with warm-start IK
for pose in poses[1:]:
    for branch in branches:
        if branch[-1] is None:
            branch.append(None)
            continue
        result = solver.solve_from_seed(pose, branch[-1])
        branch.append(result[0] if result else None)
```

See [`examples/python/path_planning.ipynb`](examples/python/path_planning.ipynb) for a full example with branch tracking and visualization.

---

A C++ API is also available — see [`include/ikdh.h`](include/ikdh.h) and the examples in [`examples/cpp/`](examples/cpp/).

---

## Adding a robot

**From RoboDK** (robot must be open in RoboDK):

```bash
pip install robodk
python3 tools/robodk_dh.py "Robot Name"
# Saves robots/<name>.yaml automatically
```

**Manually** — create `robots/<name>.yaml`:

```yaml
name: My Robot
dh:
  a:     [0.0, 0.444, 0.110, 0.0, 0.080, 0.0]   # metres
  d:     [0.265, 0.0, 0.0, 0.470, 0.0, 0.101]    # metres
  alpha: [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]        # radians
  theta: [0, -pi/2, 0, 0, 0, pi]                  # radians (joint offsets)
limits:
  J1: [-180.0, 180.0]   # degrees
  J2: [-180.0, 180.0]
  J3: [-225.0,  85.0]
  J4: [-180.0, 180.0]
  J5: [-180.0, 180.0]
  J6: [-180.0, 180.0]
```

DH convention: `T_i = Rz(θ_i + θ_offset_i) · Tz(d_i) · Tx(a_i) · Rx(α_i)`

---

## Repository structure

```
robots/                   robot YAML files (DH + joint limits)
include/
  ikdh.h                  public C++ API
  robots.h                YAML loader (no external deps)
src/
  ikdh.cpp                solver (HuPf core + Newton post-processing)
  ikdh_bindings.cpp       Python bindings (pybind11)
  hupf/                   HuPf algebraic IK core
examples/
  cpp/                    minimal C++ examples
  python/
    path_planning.ipynb   branch tracking along a Cartesian trajectory
    singularity.py        singularity analysis
tools/
  robodk_dh.py            export DH from RoboDK → robots/*.yaml
  robodk_verify.py        verify solutions against RoboDK FK
```

> Run all executables from the repository root (e.g. `./build/demo`).

---

## Algorithm

The HuPf method encodes the IK constraint as quadratic equations on Study parameters,
reduced to a degree-16 univariate polynomial. A post-processing layer recovers solutions
missed by the algebraic core:

- **Perturbation sweep** (1°/5°/10°, 12 directions) — breaks algebraic singularities
- **Joint flips** (J1/J4/J6 ± 180°) — explores all solution basins
- **Newton-Raphson refinement** (damped LS, analytical Jacobian) — converges each seed to the exact pose
- **Halton sampling** (100 quasi-random seeds, last resort) — handles robots with no real algebraic roots

### References

- Husty, Pfurner, Schröcker, Brunnthaler — *Algebraic methods in mechanism analysis*, Robotica **25** (2007)
- Capco, Quam, Pfurner, Schröcker, Sinn — *Robots, computer algebra and eight connected components*, ISSAC 2021
