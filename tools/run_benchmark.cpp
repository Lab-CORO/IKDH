#include <ikdh.h>
#include <robots.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <string>
#include <vector>

#ifndef ROBOTS_DIR
#  define ROBOTS_DIR "robots"
#endif

// ── Directory listing ─────────────────────────────────────────────────────────
#if defined(_WIN32)
#  include <windows.h>
static std::vector<std::string> listYaml(const std::string& dir)
{
    std::vector<std::string> out;
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA((dir + "\\*.yaml").c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) return out;
    do { out.push_back(dir + "\\" + fd.cFileName); } while (FindNextFileA(h, &fd));
    FindClose(h);
    std::sort(out.begin(), out.end());
    return out;
}
#else
#  include <dirent.h>
static std::vector<std::string> listYaml(const std::string& dir)
{
    std::vector<std::string> out;
    DIR* d = opendir(dir.c_str());
    if (!d) { fprintf(stderr, "Cannot open directory: %s\n", dir.c_str()); return out; }
    struct dirent* e;
    while ((e = readdir(d))) {
        const char* name = e->d_name;
        size_t len = strlen(name);
        if (len > 5 && strcmp(name + len - 5, ".yaml") == 0)
            out.push_back(dir + "/" + name);
    }
    closedir(d);
    std::sort(out.begin(), out.end());
    return out;
}
#endif

// ── Diverse joint configurations used as FK seeds ────────────────────────────
// Angles are scaled to [-1, 1] and mapped to each robot's actual joint limits,
// so the poses are always reachable regardless of the robot's range of motion.
static const double SEEDS[][6] = {
    {  0.00,  0.50, -0.50,  0.00,  0.50,  0.00 },  // mid-range
    {  0.33,  0.25, -0.67,  0.11, -0.25,  0.06 },  // generic pose 1
    { -0.33,  0.60, -0.33,  0.50, -0.50,  0.33 },  // generic pose 2
    {  0.50, -0.25,  0.25, -0.33,  0.40, -0.17 },  // generic pose 3
    { -0.17,  0.75, -0.50, -0.25,  0.60,  0.50 },  // generic pose 4
    {  0.00,  0.25,  0.00,  0.00,  0.25,  0.00 },  // near home
    { -0.50,  0.50, -0.75,  0.50, -0.33, -0.50 },  // generic pose 5
    {  0.67,  0.33, -0.40, -0.60,  0.20,  0.25 },  // generic pose 6
};
static const int N_SEEDS = (int)(sizeof(SEEDS) / sizeof(SEEDS[0]));

static IKDH::JointConfig scaledSeed(const Robots::Robot& robot, int s)
{
    IKDH::JointConfig q;
    for (int j = 0; j < 6; ++j) {
        double lo = robot.limits.lo[j], hi = robot.limits.hi[j];
        q[j] = lo + (SEEDS[s][j] * 0.5 + 0.5) * (hi - lo);
    }
    return q;
}

int main(int argc, char* argv[])
{
    const int REPS = 50;   // solves per pose

    // Accept an optional robots directory as the first argument.
    // Falls back to ROBOTS_DIR ("robots"), then "../robots" if the first doesn't exist.
    std::string dir = (argc > 1) ? argv[1] : ROBOTS_DIR;
    if (argc == 1 && listYaml(dir).empty())
        dir = "../robots";

    auto files = listYaml(dir);
    if (files.empty()) {
        fprintf(stderr, "No YAML files found in %s\n", dir.c_str());
        return 1;
    }

    printf("IKDH benchmark — %d poses × %d solves each\n\n", N_SEEDS, REPS);
    printf("  %-32s  %9s  %9s  %9s  %7s\n",
           "Robot", "avg ms", "min ms", "max ms", "avg sol");
    printf("  %s\n", std::string(66, '-').c_str());

    double grand_avg = 0.0;
    int    n_robots  = 0;

    for (const auto& path : files) {
        Robots::Robot robot;
        try { robot = Robots::loadRobot(path); }
        catch (...) { fprintf(stderr, "  [skip] %s\n", path.c_str()); continue; }

        IKDH::Solver solver(robot.dh, robot.limits);

        // Build target poses from FK so they are always reachable.
        std::vector<IKDH::Transform> poses;
        for (int s = 0; s < N_SEEDS; ++s)
            poses.push_back(IKDH::forwardKin(robot.dh, scaledSeed(robot, s)));

        // Warm-up (one pass, not timed).
        for (auto& ee : poses) solver.solve(ee);

        std::vector<double> times;
        times.reserve(N_SEEDS);
        double total_sols = 0.0;

        for (auto& ee : poses) {
            auto t0 = std::chrono::high_resolution_clock::now();
            int s = 0;
            for (int i = 0; i < REPS; ++i) s += (int)solver.solve(ee).size();
            auto t1 = std::chrono::high_resolution_clock::now();

            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count() / REPS;
            times.push_back(ms);
            total_sols += (double)s / REPS;
        }

        double avg = 0.0;
        for (double t : times) avg += t;
        avg /= (double)times.size();

        double mn = *std::min_element(times.begin(), times.end());
        double mx = *std::max_element(times.begin(), times.end());
        double avg_sol = total_sols / N_SEEDS;

        grand_avg += avg;
        ++n_robots;

        printf("  %-32s  %7.2f ms  %7.2f ms  %7.2f ms  %5.1f\n",
               robot.name.c_str(), avg, mn, mx, avg_sol);
    }

    printf("  %s\n", std::string(66, '-').c_str());
    printf("  %-32s  %7.2f ms\n", "Grand average", grand_avg / n_robots);

    return 0;
}
