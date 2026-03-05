#pragma once
#include <hupf/Input.h>

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

    double al4 = a.al[3];
    double al5 = a.al[4];
    double al6 = a.al[5];

    v6_Quadratic.set(0, (-4 * me1 * e0 + 4 * me0 * e1 + 4 * me3 * e2 - 4 * me2 * e3 - 4 * al6 * me0 * e0 - 4 * al6 * me1 * e1 - 4 * al6 * me2 * e2 - 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5
                     + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
                     + al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
                     - pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
                     + al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2) - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2)
                     - pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                     - al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                     + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                     - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5 + 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0));

    v6_Quadratic.set(1,(8 * me1 * e0 - 8 * me0 * e1 - 8 * me3 * e2 + 8 * me2 * e3 + 8 * al6 * me0 * e0 + 8 * al6 * me1 * e1 + 8 * al6 * me2 * e2 + 8 * al6 * me3 * e3) * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * al5 + (8 * me1 * e0 - 8 * me0 * e1 - 8 * me3 * e2 + 8 * me2 * e3 + 8 * al6 * me0 * e0 + 8 * al6 * me1 * e1 + 8 * al6 * me2 * e2 + 8 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) + 4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + 4 * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0));

    v6_Quadratic.set(2,-pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
                     + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
                     + al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 *al5 + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5
                     - pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
                     - pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2) - 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
                     + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
                     - al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
                     + al4 *al4  * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2) - 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0));

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

    double al4 = a.al[3];
    double al5 = a.al[4];
    double al6 = a.al[5];

    v5 = -((-4 * me1 * e0 + 4 * me0 * e1 + 4 * me3 * e2 - 4 * me2 * e3 - 4 * al6 * me0 * e0 - 4 * al6 * me1 * e1 - 4 * al6 * me2 * e2 - 4 * al6 * me3 * e3) * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * al5 - 2 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - 2 * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
           + (-al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
              - al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
              - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
              - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
              + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
              + al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
              - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5 + 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
              + pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
              + pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)) * v6)
         /
         (pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
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

    double al4 = a.al[3];
    double al5 = a.al[4];
    double al6 = a.al[5];

    v4 = -((-4 * me1 * e0 + 4 * me0 * e1 + 4 * me3 * e2 - 4 * me2 * e3 - 4 * al6 * me0 * e0 - 4 * al6 * me1 * e1 - 4 * al6 * me2 * e2 - 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - 2 * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) - 2 * al4 * al5 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al5 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - 2 * al4 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * al5 - 2 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
           + (-al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
              - al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
              - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * al4 * al4 * al5 * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * al4 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
              - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
              - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
              + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
              + al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
              - (4 * me1 * e0 - 4 * me0 * e1 - 4 * me3 * e2 + 4 * me2 * e3 + 4 * al6 * me0 * e0 + 4 * al6 * me1 * e1 + 4 * al6 * me2 * e2 + 4 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) * al5 + 2 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
              + pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
              + pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)) * v6)
         /
         ((al4 * al4 * pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2) * al5 * al5
           + al4 * al4 * al5 * al5 * pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)
           - pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al4 * al4
           - al4 * al4 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
           + 2 * al4 * al5 * al5 * (2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * al5 * al5 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0) + 2 * al4 * (2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3) * (2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2) + 2 * al4 * (2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3) * (2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0)
           + pow((2 * me0 * e0 + 2 * me1 * e1 + 2 * me2 * e2 + 2 * me3 * e3 + 2 * al6 * me0 * e1 - 2 * al6 * me1 * e0 - 2 * al6 * me2 * e3 + 2 * al6 * me3 * e2), 2) * al5 * al5
           + al5 * al5 * pow((2 * me1 * e2 - 2 * me0 * e3 + 2 * me3 * e0 - 2 * me2 * e1 + 2 * al6 * me2 * e0 + 2 * al6 * me3 * e1 - 2 * al6 * me0 * e2 - 2 * al6 * me1 * e3), 2)
           - pow((2 * me1 * e0 - 2 * me0 * e1 - 2 * me3 * e2 + 2 * me2 * e3 + 2 * al6 * me0 * e0 + 2 * al6 * me1 * e1 + 2 * al6 * me2 * e2 + 2 * al6 * me3 * e3), 2)
           - pow((2 * me2 * e0 + 2 * me3 * e1 - 2 * me0 * e2 - 2 * me1 * e3 + 2 * al6 * me0 * e3 - 2 * al6 * me1 * e2 + 2 * al6 * me2 * e1 - 2 * al6 * me3 * e0), 2)));

    return v4;
  };

}