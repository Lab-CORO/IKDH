#pragma once
#include <hupf/Input.h>
#include <hupf/proj_eval.h>

namespace LibHUPF
{
  std::vector<double> new_wrist_sol1_v6(std::vector<double> &rValue, Input& a)
  {
    std::vector<double> v6;
    v6.resize(2);
    Polynomial v6_Quadratic(2);
    Matrix e = a.studyParams();

    double e0 = e.get(0, 0);
    double e1 = e.get(1, 0);
    double e2 = e.get(2, 0);
    double e3 = e.get(3, 0);

    double me0 = rValue[0];
    double me1 = rValue[1];
    double me2 = rValue[2];
    double me3 = rValue[3];

    double p4=a.alp[3], q4=a.alq[3];
    double p5=a.alp[4], q5=a.alq[4];
    double p6=a.alp[5], q6=a.alq[5];

    v6_Quadratic.set(0, proj_eval_222([&](double al4, double al5, double al6){
      return (-4 * me1 * e0 + 4 * me0 * e1 + 4 * me3 * e2 - 4 * me2 * e3 - 4 * al6 * me0 * e0 - 4 * al6 * me1 * e1 - 4 * al6 * me2 * e2 - 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5
                     + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
                     + al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
                     - pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
                     + al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2) - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2)
                     - pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                     - al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                     + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                     - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5 + 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0);
    }, p4,q4,p5,q5,p6,q6));

    v6_Quadratic.set(1, proj_eval_222([&](double al4, double al5, double al6){
      return (8 * me1 * e0 - 8 * me0 * e1 - 8 * me3 * e2 + 8 * me2 * e3 + 8 * al6 * me0 * e0 + 8 * al6 * me1 * e1 + 8 * al6 * me2 * e2 + 8 * al6 * me3 * e3) * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * al5 + (8 * me1 * e0 - 8 * me0 * e1 - 8 * me3 * e2 + 8 * me2 * e3 + 8 * al6 * me0 * e0 + 8 * al6 * me1 * e1 + 8 * al6 * me2 * e2 + 8 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) + 4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + 4 * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0);
    }, p4,q4,p5,q5,p6,q6));

    v6_Quadratic.set(2, proj_eval_222([&](double al4, double al5, double al6){
      return -pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
                     + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
                     + al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 *al5 + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5
                     - pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
                     - pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2) - 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
                     + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                     - al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                     + al4 *al4  * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2) - 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0);
    }, p4,q4,p5,q5,p6,q6));

    v6[0] = (-v6_Quadratic.get(1) + pow(v6_Quadratic.get(1) * v6_Quadratic.get(1) - 4 * v6_Quadratic.get(2) * v6_Quadratic.get(0), 0.5)) / (2 * v6_Quadratic.get(2));
    v6[1] = (-v6_Quadratic.get(1) - pow(v6_Quadratic.get(1) * v6_Quadratic.get(1) - 4 * v6_Quadratic.get(2) * v6_Quadratic.get(0), 0.5)) / (2 * v6_Quadratic.get(2));

    return v6;

  };

  double new_wrist_sol1_v5(std::vector<double> &rValue, double v6, Input& a)
  {
    double v5;
    Matrix e = a.studyParams();

    double e0 = e.get(0, 0);
    double e1 = e.get(1, 0);
    double e2 = e.get(2, 0);
    double e3 = e.get(3, 0);

    double me0 = rValue[0];
    double me1 = rValue[1];
    double me2 = rValue[2];
    double me3 = rValue[3];

    double p4=a.alp[3], q4=a.alq[3];
    double p5=a.alp[4], q5=a.alq[4];
    double p6=a.alp[5], q6=a.alq[5];

    // v6 (argument) is the computed joint variable  -  captured by [&] inside lambdas.
    // Both numer and denom are scaled by the same factor q4^2*q5^2*q6^2, which cancels.
    double proj_numer = proj_eval_222([&](double al4, double al5, double al6){
      return -((-4 * me1 * e0 + 4 * me0 * e1 + 4 * me3 * e2 - 4 * me2 * e3 - 4 * al6 * me0 * e0 - 4 * al6 * me1 * e1 - 4 * al6 * me2 * e2 - 4 * al6 * me3 * e3) * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * al5 - 2 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - 2 * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
             + (-al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
                - al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
                - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
                + al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5 + 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
                + pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
                + pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)) * v6);
    }, p4,q4,p5,q5,p6,q6);

    double proj_denom = proj_eval_222([&](double al4, double al5, double al6){
      return (pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
          - al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
          - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
          - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
          + pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
          + 2 * al4 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
          - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
          + al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
          + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
          + 2 * al4 * pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5
          + 2 * al4 * al5 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
          + 2 * al4 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2));
    }, p4,q4,p5,q5,p6,q6);

    v5 = proj_numer / proj_denom;

    return v5;
  };

  double new_wrist_sol1_v4(std::vector<double> &rValue, double v6, Input& a)
  {
    double v4;
    Matrix e = a.studyParams();

    double e0 = e.get(0, 0);
    double e1 = e.get(1, 0);
    double e2 = e.get(2, 0);
    double e3 = e.get(3, 0);

    double me0 = rValue[0];
    double me1 = rValue[1];
    double me2 = rValue[2];
    double me3 = rValue[3];

    double p4=a.alp[3], q4=a.alq[3];
    double p5=a.alp[4], q5=a.alq[4];
    double p6=a.alp[5], q6=a.alq[5];

    double proj_numer = proj_eval_222([&](double al4, double al5, double al6){
      return -((-4 * me1 * e0 + 4 * me0 * e1 + 4 * me3 * e2 - 4 * me2 * e3 - 4 * al6 * me0 * e0 - 4 * al6 * me1 * e1 - 4 * al6 * me2 * e2 - 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - 2 * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) - 2 * al4 * al5 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al5 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - 2 * al4 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * al5 - 2 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
             + (-al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
                - al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
                - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
                - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
                + al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5 + 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
                + pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
                + pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)) * v6);
    }, p4,q4,p5,q5,p6,q6);

    double proj_denom = proj_eval_222([&](double al4, double al5, double al6){
      return ((al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
           + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
           - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
           - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
           + 2 * al4 * al5 * al5 * (2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * al5 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + 2 * al4 * (2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
           + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
           + al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
           - pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
           - pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)));
    }, p4,q4,p5,q5,p6,q6);

    v4 = proj_numer / proj_denom;

    return v4;
  };

}
