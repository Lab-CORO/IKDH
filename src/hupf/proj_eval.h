#pragma once

namespace LibHUPF
{

// Evaluate a function that is LINEAR in one variable al, given projective
// coordinates [p:q] (al = p/q), returning f(p/q)*q  (cleared denominator).
// Degree-1 multilinear decomposition: evaluates f at al=0 and al=1.
template<typename F>
inline double proj_eval_1(F f, double p, double q)
{
    double f0 = f(0.0);
    double f1 = f(1.0);
    return f0 * q + (f1 - f0) * p;
}

// Evaluate a function that is MULTILINEAR in three variables (degree ≤ 1 in each)
// given projective coordinates [p1:q1], [p2:q2], [p3:q3].
// Returns f(p1/q1, p2/q2, p3/q3) * (q1*q2*q3)  (cleared denominators).
// Uses 8-point multilinear decomposition on the vertices of {0,1}^3.
template<typename F>
inline double proj_eval_3(F f,
    double p1, double q1, double p2, double q2, double p3, double q3)
{
    double f000 = f(0,0,0);
    double f100 = f(1,0,0);
    double f010 = f(0,1,0);
    double f001 = f(0,0,1);
    double f110 = f(1,1,0);
    double f101 = f(1,0,1);
    double f011 = f(0,1,1);
    double f111 = f(1,1,1);

    double A000 = f000;
    double A100 = f100 - f000;
    double A010 = f010 - f000;
    double A001 = f001 - f000;
    double A110 = f110 - f100 - f010 + f000;
    double A101 = f101 - f100 - f001 + f000;
    double A011 = f011 - f010 - f001 + f000;
    double A111 = f111 - f110 - f101 - f011 + f100 + f010 + f001 - f000;

    double m000 = q1*q2*q3;
    double m100 = p1*q2*q3;
    double m010 = q1*p2*q3;
    double m001 = q1*q2*p3;
    double m110 = p1*p2*q3;
    double m101 = p1*q2*p3;
    double m011 = q1*p2*p3;
    double m111 = p1*p2*p3;

    return A000*m000 + A100*m100 + A010*m010 + A001*m001
         + A110*m110 + A101*m101 + A011*m011 + A111*m111;
}

// Evaluate a function of degree (2,2,2) in (al4, al5, al6) given projective
// coordinates [p4:q4], [p5:q5], [p6:q6].
// Returns f(p4/q4, p5/q5, p6/q6) * (q4^2 * q5^2 * q6^2)  (cleared denominators).
// Uses nested degree-2 Lagrange interpolation at {-1,0,1} in each dimension.
// Handles wrist_sol.h expressions which are degree 2 in al6 due to pow(linear,2) terms.
template<typename F>
inline double proj_eval_222(F f,
    double p4, double q4, double p5, double q5, double p6, double q6)
{
    // Step 1: evaluate at al4 ∈ {-1,0,1}, al5 ∈ {-1,0,1}, al6 ∈ {-1,0,1}
    double g[3][3][3];
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j)
            for (int k = 0; k <= 2; ++k)
                g[i][j][k] = f(i - 1.0, j - 1.0, k - 1.0);

    // Step 2: project over al6 (degree 2)  -  cleared by q6^2
    double h[3][3];
    for (int i = 0; i <= 2; ++i)
        for (int j = 0; j <= 2; ++j) {
            double C0 = g[i][j][1];
            double C1 = (g[i][j][2] - g[i][j][0]) * 0.5;
            double C2 = (g[i][j][2] + g[i][j][0]) * 0.5 - g[i][j][1];
            h[i][j] = C0*(q6*q6) + C1*(p6*q6) + C2*(p6*p6);
        }

    // Step 3: project over al5 (degree 2)  -  cleared by q5^2
    double k[3];
    for (int i = 0; i <= 2; ++i) {
        double C0 = h[i][1];
        double C1 = (h[i][2] - h[i][0]) * 0.5;
        double C2 = (h[i][2] + h[i][0]) * 0.5 - h[i][1];
        k[i] = C0*(q5*q5) + C1*(p5*q5) + C2*(p5*p5);
    }

    // Step 4: project over al4 (degree 2)  -  cleared by q4^2
    double C0 = k[1];
    double C1 = (k[2] - k[0]) * 0.5;
    double C2 = (k[2] + k[0]) * 0.5 - k[1];
    return C0*(q4*q4) + C1*(p4*q4) + C2*(p4*p4);
}

} // namespace LibHUPF
