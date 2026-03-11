#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <ikdh.h>
#include <robots.h>

namespace py = pybind11;
using namespace IKDH;

// ── Helpers ────────────────────────────────────────────────────────────────────

static Transform numpy_to_transform(py::array_t<double, py::array::c_style | py::array::forcecast> arr)
{
    if (arr.size() != 16)
        throw std::runtime_error("Expected a 4x4 matrix (16 elements)");
    Transform t;
    const double* data = arr.data();
    for (int i = 0; i < 16; ++i) t[i] = data[i];
    return t;
}

static py::array_t<double> transform_to_numpy(const Transform& t)
{
    auto arr = py::array_t<double>({4, 4});
    double* data = arr.mutable_data();
    for (int i = 0; i < 16; ++i) data[i] = t[i];
    return arr;
}

// ── Module ─────────────────────────────────────────────────────────────────────

PYBIND11_MODULE(ikdh, m)
{
    m.doc() = "IKDH — Inverse kinematics for 6R serial robots (HuPf algebraic method)";

    // ── DHTable ──────────────────────────────────────────────────────────────
    py::class_<DHTable>(m, "DHTable")
        .def(py::init<>())
        .def(py::init([](std::vector<double> a, std::vector<double> d,
                         std::vector<double> alpha, std::vector<double> theta) {
            if (a.size() != 6 || d.size() != 6 || alpha.size() != 6 || theta.size() != 6)
                throw std::runtime_error("All DH arrays must have exactly 6 elements");
            DHTable dh{};
            for (int i = 0; i < 6; ++i) {
                dh.a[i]       = a[i];
                dh.d[i]       = d[i];
                dh.alpha[i]   = alpha[i];
                dh.theta[i]   = theta[i];
                dh.revolute[i] = true;
            }
            return dh;
        }), py::arg("a"), py::arg("d"), py::arg("alpha"), py::arg("theta"),
            "DH parameters. Units: a/d in metres, alpha/theta in radians.")
        .def_property("a",
            [](const DHTable& dh) { return std::vector<double>(dh.a, dh.a + 6); },
            [](DHTable& dh, std::vector<double> v) { for (int i=0;i<6;++i) dh.a[i]=v[i]; })
        .def_property("d",
            [](const DHTable& dh) { return std::vector<double>(dh.d, dh.d + 6); },
            [](DHTable& dh, std::vector<double> v) { for (int i=0;i<6;++i) dh.d[i]=v[i]; })
        .def_property("alpha",
            [](const DHTable& dh) { return std::vector<double>(dh.alpha, dh.alpha + 6); },
            [](DHTable& dh, std::vector<double> v) { for (int i=0;i<6;++i) dh.alpha[i]=v[i]; })
        .def_property("theta",
            [](const DHTable& dh) { return std::vector<double>(dh.theta, dh.theta + 6); },
            [](DHTable& dh, std::vector<double> v) { for (int i=0;i<6;++i) dh.theta[i]=v[i]; })
        .def("__repr__", [](const DHTable& dh) {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "DHTable(a=[%.4f,...], d=[%.4f,...], alpha=[%.4f,...], theta=[%.4f,...])",
                dh.a[0], dh.d[0], dh.alpha[0], dh.theta[0]);
            return std::string(buf);
        });

    // ── JointLimits ──────────────────────────────────────────────────────────
    py::class_<JointLimits>(m, "JointLimits")
        .def(py::init<>(), "Default: ±180° for all joints")
        .def(py::init([](std::vector<std::pair<double,double>> pairs) {
            if (pairs.size() != 6)
                throw std::runtime_error("Expected exactly 6 (lo, hi) pairs");
            JointLimits jl;
            for (int i = 0; i < 6; ++i) { jl.lo[i] = pairs[i].first; jl.hi[i] = pairs[i].second; }
            return jl;
        }), py::arg("limits"),
            "List of 6 (lo, hi) tuples in degrees.")
        .def_property("lo",
            [](const JointLimits& jl) { return std::vector<double>(jl.lo, jl.lo + 6); },
            [](JointLimits& jl, std::vector<double> v) { for (int i=0;i<6;++i) jl.lo[i]=v[i]; })
        .def_property("hi",
            [](const JointLimits& jl) { return std::vector<double>(jl.hi, jl.hi + 6); },
            [](JointLimits& jl, std::vector<double> v) { for (int i=0;i<6;++i) jl.hi[i]=v[i]; });

    // ── Solver ───────────────────────────────────────────────────────────────
    py::class_<Solver>(m, "Solver")
        .def(py::init<const DHTable&, const JointLimits&>(),
             py::arg("dh"), py::arg("limits") = JointLimits{},
             "Construct once per robot, then call solve() for each pose.")
        .def("solve",
            [](const Solver& solver, py::array_t<double, py::array::c_style | py::array::forcecast> ee) {
                auto sols = solver.solve(numpy_to_transform(ee));
                // Return list of numpy arrays shape (6,)
                py::list result;
                for (const auto& q : sols) {
                    auto arr = py::array_t<double>(6);
                    double* data = arr.mutable_data();
                    for (int i = 0; i < 6; ++i) data[i] = q[i];
                    result.append(arr);
                }
                return result;
            },
            py::arg("ee"),
            "Solve IK for a 4x4 end-effector transform (numpy array).\n"
            "Returns a list of numpy arrays of shape (6,), one per solution, in degrees.");

    // ── Free functions ────────────────────────────────────────────────────────
    m.def("forward_kin",
        [](const DHTable& dh, py::array_t<double, py::array::c_style | py::array::forcecast> q) {
            if (q.size() != 6)
                throw std::runtime_error("Expected 6 joint values");
            JointConfig jc;
            const double* data = q.data();
            for (int i = 0; i < 6; ++i) jc[i] = data[i];
            return transform_to_numpy(forwardKin(dh, jc));
        },
        py::arg("dh"), py::arg("q"),
        "Forward kinematics. q in degrees. Returns 4x4 numpy array.");

    m.def("pose_from_xyzrpw",
        [](double x_mm, double y_mm, double z_mm,
           double rx_deg, double ry_deg, double rz_deg) {
            return transform_to_numpy(poseFromXYZRPW(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg));
        },
        py::arg("x_mm"), py::arg("y_mm"), py::arg("z_mm"),
        py::arg("rx_deg"), py::arg("ry_deg"), py::arg("rz_deg"),
        "Build a 4x4 transform from a RoboDK-style pose (mm, degrees, Rz*Ry*Rx).");

    m.def("fk_error",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> A,
           py::array_t<double, py::array::c_style | py::array::forcecast> B) {
            return fkError(numpy_to_transform(A), numpy_to_transform(B));
        },
        py::arg("A"), py::arg("B"),
        "Sum of squared element-wise differences between two 4x4 transforms.");

    // ── Robot loader ──────────────────────────────────────────────────────────
    py::class_<Robots::Robot>(m, "Robot")
        .def_readonly("name",   &Robots::Robot::name)
        .def_readonly("dh",     &Robots::Robot::dh)
        .def_readonly("limits", &Robots::Robot::limits);

    m.def("load_robot", &Robots::loadRobot, py::arg("yaml_path"),
          "Load a robot from a YAML file (robots/*.yaml).\n"
          "Returns a Robot object with .name, .dh, and .limits attributes.");
}
