#pragma once

#include <ikdh.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Minimal YAML robot loader — no external dependencies.
// Expected format: robots/<name>.yaml  (see robots/ directory for examples)

namespace Robots {

namespace detail {

// Parse a pi-expression like "pi/2", "-pi/2", "pi", "0", "1.5", etc.
inline double parsePiExpr(const std::string& s)
{
    std::string t = s;
    // strip whitespace
    t.erase(0, t.find_first_not_of(" \t\r\n"));
    auto e = t.find_last_not_of(" \t\r\n");
    if (e != std::string::npos) t = t.substr(0, e + 1);

    bool neg = false;
    if (!t.empty() && t[0] == '-') { neg = true; t = t.substr(1); }

    double val;
    if (t == "pi") {
        val = M_PI;
    } else if (t.substr(0, 3) == "pi/") {
        val = M_PI / std::stod(t.substr(3));
    } else if (t.substr(0, 3) == "pi*") {
        val = M_PI * std::stod(t.substr(3));
    } else {
        val = std::stod(t);
    }
    return neg ? -val : val;
}

// Parse "[v1, v2, v3, v4, v5, v6]" into a 6-element vector.
inline std::vector<double> parseArray(const std::string& line)
{
    std::vector<double> out;
    auto lb = line.find('[');
    auto rb = line.find(']');
    if (lb == std::string::npos || rb == std::string::npos)
        throw std::runtime_error("Expected '[...]' in: " + line);

    std::string inner = line.substr(lb + 1, rb - lb - 1);
    std::istringstream ss(inner);
    std::string token;
    while (std::getline(ss, token, ','))
        out.push_back(parsePiExpr(token));
    return out;
}

// Strip inline YAML comment (#...) and trailing whitespace.
inline std::string stripComment(const std::string& s)
{
    auto pos = s.find('#');
    std::string t = (pos != std::string::npos) ? s.substr(0, pos) : s;
    auto e = t.find_last_not_of(" \t\r\n");
    return (e != std::string::npos) ? t.substr(0, e + 1) : "";
}

} // namespace detail

struct Robot {
    std::string      name;
    IKDH::DHTable    dh;
    IKDH::JointLimits limits;
};

inline Robot loadRobot(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("Cannot open robot file: " + path);

    Robot robot;
    std::vector<double> a, d, alpha, theta;
    std::string line;
    std::string section;  // "dh" or "limits"
    int limCount = 0;

    while (std::getline(f, line)) {
        std::string s = detail::stripComment(line);
        if (s.empty()) continue;

        // Detect section headers
        if (s.find("dh:") != std::string::npos && s.find(':') != std::string::npos) {
            section = "dh"; continue;
        }
        if (s.find("limits:") != std::string::npos) {
            section = "limits"; continue;
        }

        // name field
        if (s.find("name:") != std::string::npos) {
            auto pos = s.find(':');
            std::string v = s.substr(pos + 1);
            v.erase(0, v.find_first_not_of(" \t"));
            robot.name = v;
            continue;
        }

        if (section == "dh" && s.find('[') != std::string::npos) {
            // Match key: strip leading whitespace and check prefix
            std::string trimmed = s;
            trimmed.erase(0, trimmed.find_first_not_of(" \t"));
            if      (trimmed.substr(0, 6) == "alpha:") alpha = detail::parseArray(s);
            else if (trimmed.substr(0, 6) == "theta:") theta = detail::parseArray(s);
            else if (trimmed.substr(0, 2) == "a:")     a     = detail::parseArray(s);
            else if (trimmed.substr(0, 2) == "d:")     d     = detail::parseArray(s);
        }

        if (section == "limits" && s.find('[') != std::string::npos && limCount < 6) {
            auto row = detail::parseArray(s);
            if (row.size() >= 2) {
                robot.limits.lo[limCount] = row[0];
                robot.limits.hi[limCount] = row[1];
                ++limCount;
            }
        }
    }

    if (a.size() != 6 || d.size() != 6 || alpha.size() != 6 || theta.size() != 6)
        throw std::runtime_error("DH table incomplete in: " + path);
    if (limCount != 6)
        throw std::runtime_error("Joint limits incomplete in: " + path);

    for (int i = 0; i < 6; ++i) {
        robot.dh.a[i]       = a[i];
        robot.dh.d[i]       = d[i];
        robot.dh.alpha[i]   = alpha[i];
        robot.dh.theta[i]   = theta[i];
        robot.dh.revolute[i] = true;
    }

    return robot;
}

} // namespace Robots
