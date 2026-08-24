#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <ikdh.h>
#include <robot.h>

namespace py = pybind11;
using namespace IKDH;

// Helpers

static Transform numpyToTransform(py::array_t<double, py::array::c_style | py::array::forcecast> arr)
{
    if (arr.size() != 16)
        throw std::runtime_error("Expected a 4x4 matrix (16 elements)");
    Transform t;
    const double* data = arr.data();
    for (int i = 0; i < 16; ++i) t[i] = data[i];
    return t;
}

static py::array_t<double> transformToNumpy(const Transform& t)
{
    auto arr = py::array_t<double>({4, 4});
    double* data = arr.mutable_data();
    for (int i = 0; i < 16; ++i) data[i] = t[i];
    return arr;
}

// Convert IK solutions to a Python list of (6,) numpy arrays.
static py::list jointConfigsToList(const std::vector<JointConfig>& sols)
{
    py::list result;
    for (const auto& q : sols) {
        auto arr = py::array_t<double>(6);
        double* data = arr.mutable_data();
        for (int i = 0; i < 6; ++i) data[i] = q[i];
        result.append(arr);
    }
    return result;
}

static std::array<double, 6> vecToArray6(const std::vector<double>& v)
{
    if (v.size() != 6) throw std::runtime_error("Expected exactly 6 elements");
    std::array<double, 6> a;
    for (int i = 0; i < 6; ++i) a[i] = v[i];
    return a;
}

static JointConfig numpyToJointConfig(py::array_t<double, py::array::c_style | py::array::forcecast> arr)
{
    if (arr.size() != 6) throw std::runtime_error("Expected 6 joint values");
    JointConfig q;
    const double* data = arr.data();
    for (int i = 0; i < 6; ++i) q[i] = data[i];
    return q;
}

static py::array_t<double> matrix6ToNumpy(const Matrix6& M)
{
    auto arr = py::array_t<double>({6, 6});
    double* data = arr.mutable_data();
    for (int r = 0; r < 6; ++r)
        for (int c = 0; c < 6; ++c) data[r*6 + c] = M[r][c];
    return arr;
}

// Module

PYBIND11_MODULE(_ikdh, m)
{
    m.doc() = "IKDH  -  Lightweight inverse kinematics for 6R serial robots (Jacobian + Halton multi-start Newton refinement)";

    // DHTable
    py::class_<DHTable>(m, "DHTable")
        .def(py::init<>())
        .def(py::init([](std::vector<double> a, std::vector<double> d,
                         std::vector<double> alpha, std::vector<double> theta) {
            if (a.size() != 6 || d.size() != 6 || alpha.size() != 6 || theta.size() != 6)
                throw std::runtime_error("All DH arrays must have exactly 6 elements");
            DHTable dh{};
            for (int i = 0; i < 6; ++i) {
                dh.a[i]        = a[i];
                dh.d[i]        = d[i];
                dh.alpha[i]    = alpha[i];
                dh.theta[i]    = theta[i];
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

    // JointLimits
    py::class_<JointLimits>(m, "JointLimits")
        .def(py::init<>(), "Default: ±180° for all joints.")
        .def(py::init([](std::vector<std::pair<double,double>> pairs) {
            if (pairs.size() != 6)
                throw std::runtime_error("Expected exactly 6 (lo, hi) pairs");
            JointLimits jl;
            for (int i = 0; i < 6; ++i) { jl.lo[i] = pairs[i].first; jl.hi[i] = pairs[i].second; }
            return jl;
        }), py::arg("limits"), "List of 6 (lo, hi) tuples in degrees.")
        .def_property("lo",
            [](const JointLimits& jl) { return std::vector<double>(jl.lo, jl.lo + 6); },
            [](JointLimits& jl, std::vector<double> v) { for (int i=0;i<6;++i) jl.lo[i]=v[i]; })
        .def_property("hi",
            [](const JointLimits& jl) { return std::vector<double>(jl.hi, jl.hi + 6); },
            [](JointLimits& jl, std::vector<double> v) { for (int i=0;i<6;++i) jl.hi[i]=v[i]; });

    // Solver
    py::class_<Solver>(m, "Solver")
        .def(py::init<const DHTable&, const JointLimits&>(),
             py::arg("dh"), py::arg("limits") = JointLimits{},
             "Construct once per robot.")
        .def("solve",
            [](const Solver& solver,
               py::array_t<double, py::array::c_style | py::array::forcecast> ee,
               bool expand_wraps,
               int n_seeds) {
                auto sols = solver.solve(numpyToTransform(ee), expand_wraps, n_seeds);
                return jointConfigsToList(sols);
            },
            py::arg("ee"), py::arg("expand_wraps") = false, py::arg("n_seeds") = 64,
            "Return all IK solutions for a 4x4 end-effector transform (numpy array).\n"
            "Each solution is a (6,) array in degrees. If expand_wraps=True, also\n"
            "includes ±360° equivalents within joint limits. n_seeds sets how many\n"
            "Halton quasi-random seeds cover the joint-space cube before Newton\n"
            "refinement; higher values trade speed for a better chance of finding\n"
            "distant solution branches.")
        .def("solve_from_seed",
            [](const Solver& solver,
               py::array_t<double, py::array::c_style | py::array::forcecast> ee,
               py::array_t<double, py::array::c_style | py::array::forcecast> seed,
               int max_iter) {
                if (seed.size() != 6)
                    throw std::runtime_error("seed must have exactly 6 elements");
                JointConfig jc;
                const double* s = seed.data();
                for (int i = 0; i < 6; ++i) jc[i] = s[i];
                auto sols = solver.solveFromSeed(numpyToTransform(ee), jc, max_iter);
                return jointConfigsToList(sols);
            },
            py::arg("ee"), py::arg("seed"), py::arg("max_iter") = 100,
            "Warm-start IK via damped Newton-Raphson from seed (degrees, shape (6,)).\n"
            "Returns a list with one (6,) solution if Newton converges within limits,\n"
            "otherwise an empty list. ~100x faster than solve()  -  use for path planning.");

    // Free functions
    m.def("forward_kin",
        [](const DHTable& dh, py::array_t<double, py::array::c_style | py::array::forcecast> q) {
            if (q.size() != 6)
                throw std::runtime_error("Expected 6 joint values");
            JointConfig jc;
            const double* data = q.data();
            for (int i = 0; i < 6; ++i) jc[i] = data[i];
            return transformToNumpy(forwardKin(dh, jc));
        },
        py::arg("dh"), py::arg("q"),
        "Forward kinematics. q in degrees. Returns a (4, 4) numpy array.");

    m.def("pose_from_xyzrpw",
        [](double x_mm, double y_mm, double z_mm,
           double rx_deg, double ry_deg, double rz_deg) {
            return transformToNumpy(poseFromXYZRPW(x_mm, y_mm, z_mm, rx_deg, ry_deg, rz_deg));
        },
        py::arg("x_mm"), py::arg("y_mm"), py::arg("z_mm"),
        py::arg("rx_deg"), py::arg("ry_deg"), py::arg("rz_deg"),
        "Build a (4, 4) transform from a RoboDK pose (mm, degrees, Rz·Ry·Rx).");

    m.def("fk_error",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> A,
           py::array_t<double, py::array::c_style | py::array::forcecast> B) {
            return fkError(numpyToTransform(A), numpyToTransform(B));
        },
        py::arg("A"), py::arg("B"),
        "Σ(A_ij − B_ij)² over all 16 elements of two 4×4 transforms.");

    // Robot  -  central object built around a DH table: FK, IK, Jacobian /
    // manipulability, and moveL / moveJ path planning.
    py::class_<Robot>(m, "Robot")
        .def(py::init<>())
        .def("define_dh",
            [](Robot& r, std::vector<double> a, std::vector<double> d,
               std::vector<double> alpha, std::vector<double> theta, const JointLimits& limits) {
                r.defineDH(vecToArray6(a), vecToArray6(d), vecToArray6(alpha), vecToArray6(theta), limits);
            },
            py::arg("a"), py::arg("d"), py::arg("alpha"), py::arg("theta"),
            py::arg("limits") = JointLimits{},
            "Set the DH parameters (a/d in metres, alpha/theta in radians) and joint limits (degrees).")
        .def_static("load_yaml", &Robot::loadYAML, py::arg("yaml_path"),
            "Load DH parameters and joint limits from a robot YAML file.")
        .def("forward_kinematics",
            [](const Robot& r, py::array_t<double, py::array::c_style | py::array::forcecast> q) {
                return transformToNumpy(r.forwardKinematics(numpyToJointConfig(q)));
            },
            py::arg("q"), "Forward kinematics. q in degrees. Returns a (4, 4) numpy array.")
        .def("inverse_kinematics",
            [](const Robot& r, py::array_t<double, py::array::c_style | py::array::forcecast> pose,
               int n_seeds, bool expand_wraps) {
                return jointConfigsToList(r.inverseKinematics(numpyToTransform(pose), n_seeds, expand_wraps));
            },
            py::arg("pose"), py::arg("n_seeds") = 64, py::arg("expand_wraps") = false,
            "All IK solutions for a 4x4 end-effector pose (numpy array). n_seeds sets how "
            "many Halton seeds cover the joint-space cube before Newton refinement.")
        .def("jacobian",
            [](const Robot& r, py::array_t<double, py::array::c_style | py::array::forcecast> q) {
                return matrix6ToNumpy(r.jacobian(numpyToJointConfig(q)));
            },
            py::arg("q"),
            "Geometric Jacobian at q (degrees). Returns a (6, 6) numpy array, rows "
            "vx,vy,vz,wx,wy,wz and columns joints.")
        .def("jacobian_determinant",
            [](const Robot& r, py::array_t<double, py::array::c_style | py::array::forcecast> q) {
                return r.jacobianDeterminant(numpyToJointConfig(q));
            },
            py::arg("q"), "Determinant of the geometric Jacobian at q (degrees).")
        .def("manipulability",
            [](const Robot& r, py::array_t<double, py::array::c_style | py::array::forcecast> q) {
                return r.manipulability(numpyToJointConfig(q));
            },
            py::arg("q"), "Yoshikawa manipulability sqrt(det(J * J^T)) at q (degrees).")
        .def("move_l",
            [](const Robot& r,
               py::array_t<double, py::array::c_style | py::array::forcecast> pose_a,
               py::array_t<double, py::array::c_style | py::array::forcecast> pose_b,
               py::array_t<double, py::array::c_style | py::array::forcecast> q0,
               int n_points) {
                return jointConfigsToList(r.moveL(numpyToTransform(pose_a), numpyToTransform(pose_b),
                                                   numpyToJointConfig(q0), n_points));
            },
            py::arg("pose_a"), py::arg("pose_b"), py::arg("q0"), py::arg("n_points"),
            "Cartesian path from pose_a to pose_b (4x4 numpy arrays) via inverse-Jacobian "
            "steps, n_points waypoints including both ends, starting from joint config q0. "
            "Stops early (returning what was found so far) if a step's error grows too "
            "large or leaves the joint limits.")
        .def("move_j",
            [](const Robot& r,
               py::array_t<double, py::array::c_style | py::array::forcecast> q_a,
               py::array_t<double, py::array::c_style | py::array::forcecast> q_b,
               int n_points) {
                return jointConfigsToList(r.moveJ(numpyToJointConfig(q_a), numpyToJointConfig(q_b), n_points));
            },
            py::arg("q_a"), py::arg("q_b"), py::arg("n_points"),
            "Joint-space path from q_a to q_b (degrees) via linear interpolation, n_points "
            "waypoints including both ends.")
        .def_readonly("name",   &Robot::name)
        .def_readonly("dh",     &Robot::dh)
        .def_readonly("limits", &Robot::limits);

    m.def("load_robot", &Robot::loadYAML, py::arg("yaml_path"),
          "Load a robot from a YAML file. Returns a Robot with .name, .dh, .limits.");
}
