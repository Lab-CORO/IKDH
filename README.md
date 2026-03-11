# IKDH — Inverse Kinematics via Denavit-Hartenberg

IK solver for general 6R serial robots, based on the algebraic
method of Husty and Pfurner. Returns all the solutions for every geometry, within the joint limits. DH parameters and joint limits can be exported from RoboDK via the included Python script.

## Using IKDH as a library

Requires CMake ≥ 3.14 and a C++11 compiler.

```cmake
include(FetchContent)
FetchContent_Declare(ikdh
    GIT_REPOSITORY https://github.com/Lab-CORO/IKDH.git
    GIT_TAG master
)
FetchContent_MakeAvailable(ikdh)

target_link_libraries(your_target PRIVATE ikdh)
```

## Python bindings

Build the Python module with:

```bash
cmake -B build -DIKDH_BUILD_PYTHON=ON
cmake --build build --target ikdh_py
```

This produces `build/ikdh.cpython-<version>.so`, usable directly:

```python
import sys; sys.path.insert(0, 'build')
import ikdh

robot  = ikdh.load_robot("robots/gofa5.yaml")
solver = ikdh.Solver(robot.dh, robot.limits)
ee     = ikdh.pose_from_xyzrpw(500.0, 0.0, 500.0, 0.0, 90.0, 0.0)
sols   = solver.solve(ee)   # list of numpy arrays of shape (6,), in degrees
```

## Adding a robot

**Option 1 — from RoboDK**, with your robot open in RoboDK:

```bash
pip install robodk
python3 tools/robodk_dh.py "Robot Name"
# Connects to the running RoboDK session, extracts the Standard DH parameters,
# and saves robots/<name>.yaml automatically.
```

**Option 2 — manually**, create `robots/<name>.yaml`:

```yaml
name: My Robot
dh:
  a:     [0.0, 0.444, 0.110, 0.0, 0.080, 0.0]   # metres
  d:     [0.265, 0.0, 0.0, 0.470, 0.0, 0.101]    # metres
  alpha: [-pi/2, 0, -pi/2, pi/2, -pi/2, 0]        # radians (pi expressions ok)
  theta: [0, -pi/2, 0, 0, 0, pi]                  # radians (offsets)
limits:
  J1: [-180.0, 180.0]   # degrees
  J2: [-180.0, 180.0]
  J3: [-225.0,  85.0]
  J4: [-180.0, 180.0]
  J5: [-180.0, 180.0]
  J6: [-180.0, 180.0]
```

DH convention: `T_i = Rz(θ_i + θ_offset_i) · Tz(d_i) · Tx(a_i) · Rx(α_i)`

Pre-defined robots are in `robots/`: `ur5e`, `gofa5`, `gofa12`, `crb15000_10`,
`fanuc_crx_5ia`, `fanuc_crx_10ia`, `doosan_a0509`.

## C++ API

```cpp
#include <ikdh.h>
#include <robots.h>   // minimal YAML loader, no external deps

auto robot = Robots::loadRobot("robots/gofa5.yaml");
IKDH::Solver solver(robot.dh, robot.limits);

// Option A — from a Cartesian pose (RoboDK convention: Rz*Ry*Rx)
IKDH::Transform ee = IKDH::poseFromXYZRPW(500.0, 0.0, 500.0,  // x y z (mm)
                                            0.0,  90.0, 0.0);  // rx ry rz (deg)

// Option B — from a known joint configuration (degrees)
IKDH::Transform ee = IKDH::forwardKin(robot.dh, {0, 45, 0, 0, 45, 0});

auto sols = solver.solve(ee);

for (auto& q : sols) {
    for (double v : q) printf("  %.3f", v);
    printf("\n");
}
```

Solutions are returned in degrees, within the specified joint limits.

## Verifying against RoboDK

`tools/robodk_verify.py` connects to a running RoboDK session, runs `build/demo`
to get the IK solutions, then sends each configuration to the robot model and
compares the TCP position reported by RoboDK's own FK. This gives an independent
check that the DH parameters and joint angles are correct.

```bash
python3 tools/robodk_verify.py robots/gofa12.yaml
```

The robot name in the YAML must match the robot name in the RoboDK scene exactly.

## Repository structure

```
robots/                   robot YAML files (DH + joint limits)
include/
  ikdh.h                  public C++ API
  robots.h                minimal YAML loader (no external deps)
src/
  ikdh.cpp                solver (HuPf core + Newton post-processing)
  ikdh_bindings.cpp       Python bindings (pybind11)
  hupf/                   HuPf algebraic IK (Husty-Pfurner)
examples/
  cpp/
    gofa5_DH.cpp          minimal example: DH defined inline
    gofa5_yaml.cpp        minimal example: DH loaded from YAML
    demo.cpp              all robots round-trip demo
  python/
    gofa5_DH.py           same as gofa5_DH.cpp in Python
    gofa5_yaml.py         same as gofa5_yaml.cpp in Python
    demo.py               same as demo.cpp in Python
tools/
  robodk_dh.py            export robot DH from RoboDK → robots/*.yaml
  robodk_verify.py        stream solutions to RoboDK, report FK error
```

> All executables must be run from the repository root so that paths like
> `"robots/gofa5.yaml"` resolve correctly (e.g. `./build/demo`).

## Algorithm notes

The HuPf method encodes the IK constraint as a system of quadratic equations on
Study parameters, reduced to a univariate degree-16 polynomial. This implementation
adds a post-processing layer to recover solutions missed by the algebraic core:

- **Perturbation sweep** (1°/5°/10° in 12 directions) bypasses algebraic singularities
- **Joint flips** (J1/J4/J6 ± 180°) explore all solution basins
- **Newton-Raphson refinement** (damped LS, analytical Jacobian) converges each seed back to the exact pose
- **Halton joint-space sampling** (last resort, 100 quasi-random seeds) handles robots where the HuPf polynomial returns no real roots

### References

- Husty, Pfurner, Schröcker, Brunnthaler — *Algebraic methods in mechanism analysis*, Robotica **25** (2007)
- Capco, Quam, Pfurner, Schröcker, Sinn — *Robots, computer algebra and eight connected components*, ISSAC 2021
