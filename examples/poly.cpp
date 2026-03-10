#include <ikdh.h>
#include <robots.h>

#include <cmath>
#include <cstdio>

static void printPoly(const char* label, const std::vector<double>& c)
{
    printf("%s  (degree %d)\n", label, (int)c.size() - 1);
    for (int k = 0; k < (int)c.size(); ++k)
        printf("  %+.6e", c[k]);
    printf("\n");
}

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

int main()
{
    auto robot = Robots::loadRobot(ROBOTS_DIR "/gofa5.yaml");
    IKDH::Solver solver(robot.dh, robot.limits);
    auto& dh = robot.dh;

    IKDH::Transform T1 = IKDH::forwardKin(dh, {0,  0,  0,  0,  0,  0});
    IKDH::Transform T2 = IKDH::forwardKin(dh, {10,-30, 60,-20, 45, 15});
    IKDH::Transform T3 = IKDH::forwardKin(dh, {45,-60, 30, 90,-45, 90});

    solver.solve(T1); printPoly("home  {0,0,0,0,0,0}",          solver.lastPolynomial());
    solver.solve(T2); printPoly("pose2 {10,-30,60,-20,45,15}",   solver.lastPolynomial());
    solver.solve(T3); printPoly("pose3 {45,-60,30,90,-45,90}",   solver.lastPolynomial());

    return 0;
}
