#include <emscripten/bind.h>
#include <ikdh.h>
#include <memory>

using namespace emscripten;
using namespace IKDH;

static Solver* g_solver = nullptr;

void setRobot(val a, val d, val alpha, val theta, val lo, val hi)
{
    DHTable dh;
    JointLimits limits;
    for (int i = 0; i < 6; ++i) {
        dh.a[i]     = a[i].as<double>();
        dh.d[i]     = d[i].as<double>();
        dh.alpha[i] = alpha[i].as<double>();
        dh.theta[i] = theta[i].as<double>();
        dh.revolute[i] = true;
        limits.lo[i] = lo[i].as<double>();
        limits.hi[i] = hi[i].as<double>();
    }
    delete g_solver;
    g_solver = new Solver(dh, limits);
}

val solve(double x, double y, double z, double rx, double ry, double rz)
{
    val result = val::array();
    if (!g_solver) return result;
    Transform ee = poseFromXYZRPW(x, y, z, rx, ry, rz);
    auto solutions = g_solver->solve(ee);
    for (size_t i = 0; i < solutions.size(); ++i) {
        val sol = val::array();
        for (int j = 0; j < 6; ++j)
            sol.set(j, solutions[i][j]);
        result.set(i, sol);
    }
    return result;
}

EMSCRIPTEN_BINDINGS(ikdh)
{
    function("setRobot", &setRobot);
    function("solve", &solve);
}
