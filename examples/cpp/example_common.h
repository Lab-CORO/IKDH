#pragma once

#include <ikdh.h>

#include <cstdio>

// Solve IK for `ee` and print all solutions with their FK error.
inline void printSolutions(const IKDH::Solver& solver, const IKDH::DHTable& dh,
                           const IKDH::Transform& ee, bool expand_wraps = false)
{
    auto sols = solver.solve(ee, expand_wraps);

    printf("%zu solution(s) found%s\n", sols.size(),
           expand_wraps ? " (with wrap expansion)" : "");
    for (size_t i = 0; i < sols.size(); ++i) {
        IKDH::Transform check = IKDH::forwardKin(dh, sols[i]);
        double err = IKDH::fkError(ee, check);
        printf("  [%2zu]", i);
        for (double v : sols[i]) printf("  %7.3f", v);
        printf("   FK err = %.1e\n", err);
    }
    printf("\n");
}
