#pragma once
#include <hupf/Input.h>

namespace LibHUPF
{
  std::vector<Polynomial> h1_v4q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((-t0 * al6 * a6 - 2 * t4 + t2 * al4 * d5 - 2 * t4 * al4 * al6 + t3 * d4 + t3 * d6 + t0 * al4 * a6 - t2 * al4 * d4 - t0 * al4 * a4 + t0 * a5 * al5 + t0 * al6 * a4 + t2 * al6 * d5 + t2 * al6 * d4 + t2 * al6 * d6 + t2 * d6 * al4 - t1 * a5 * al5 * al4 - t1 * al6 * al4 * a4 + t1 * al6 * a5 * al5 + t1 * al6 * a6 * al4 - t3 * al6 * al4 * d5 + t3 * al6 * al4 * d4 - t3 * al6 * al4 * d6 - t1 * a4 + t1 * a6 + t3 * d5 + t0 * al6 * a5 * al5 * al4 + 2 * t5 * al4 - 2 * t5 * al6),(t1*al6*d4-t3*a5*al5-2*t6*al4+2*t7+t0*d5-t3*al6*a5*al5*al4-t2*al6*a6*al4+t2*a4-t0*al6*al4*d5+t1*al4*d5+t1*al6*d5+t1*al6*d6+t3*al4*a4-t3*al6*a4+t3*al6*a6-t3*a6*al4+2*t7*al6*al4+t2*a5*al5*al4-t2*al6*a5*al5+t2*al6*al4*a4-t1*al4*d4+t1*al4*d6+t0*al6*al4*d4-t0*al6*al4*d6+2*t6*al6+t0*d6+t0*d4-t2*a6)));
    result.push_back(Polynomial((t1 * a5 * al5 + t1 * al4 * a6 - 2 * t5 - 2 * t4 * al4 - t2 * d6 - t1 * al6 * a6 - t0 * al6 * a5 * al5 + t3 * al6 * d5 - 2 * t5 * al4 * al6 - t3 * al4 * d4 + t3 * al6 * d4 + t3 * al4 * d5 + t0 * a5 * al5 * al4 + t0 * al6 * al4 * a4 - t0 * al6 * a6 * al4 + t2 * al6 * al4 * d5 - t2 * al6 * al4 * d4 - t1 * al4 * a4 + t2 * al6 * al4 * d6 + t1 * al6 * a4 + t3 * al6 * d6 + t3 * d6 * al4 + t0 * a4 - t0 * a6 - t2 * d5 - t2 * d4 + 2 * t4 * al6 + t1 * al6 * a5 * al5 * al4),(-t3*al6*a5*al5-2*t6-t2*al6*a6-t3*a6+2*t7*al6+t2*a5*al5+t3*a4+t1*d5+t2*al6*a5*al5*al4-t0*al4*d5-t0*al6*d5-t0*al6*d6-t2*al4*a4+t2*al6*a4+t2*a6*al4-2*t6*al6*al4-2*t7*al4-t1*al6*al4*d5+t3*a5*al5*al4+t3*al6*al4*a4-t3*al6*a6*al4+t1*d6+t1*d4+t1*al6*al4*d4-t1*al6*al4*d6+t0*al4*d4-t0*al6*d4-t0*al4*d6)));
    result.push_back(Polynomial((-t3 * al6 * a5 * al5 - 2 * t6 - t2 * al6 * a6 - t3 * a6 + 2 * t7 * al6 + t2 * a5 * al5 + t3 * a4 + t1 * d5 + t2 * al6 * a5 * al5 * al4 - t0 * al4 * d5 - t0 * al6 * d5 - t0 * al6 * d6 - t2 * al4 * a4 + t2 * al6 * a4 + t2 * a6 * al4 - 2 * t6 * al6 * al4 - 2 * t7 * al4 - t1 * al6 * al4 * d5 + t3 * a5 * al5 * al4 + t3 * al6 * al4 * a4 - t3 * al6 * a6 * al4 + t1 * d6 + t1 * d4 + t1 * al6 * al4 * d4 - t1 * al6 * al4 * d6 + t0 * al4 * d4 - t0 * al6 * d4 - t0 * al4 * d6),(-t1*a5*al5-t1*al4*a6+2*t5+2*t4*al4+t2*d6+t1*al6*a6+t0*al6*a5*al5-t3*al6*d5+2*t5*al4*al6+t3*al4*d4-t3*al6*d4-t3*al4*d5-t0*a5*al5*al4-t0*al6*al4*a4+t0*al6*a6*al4-t2*al6*al4*d5+t2*al6*al4*d4+t1*al4*a4-t2*al6*al4*d6-t1*al6*a4-t3*al6*d6-t3*d6*al4-t0*a4+t0*a6+t2*d5+t2*d4-2*t4*al6-t1*al6*a5*al5*al4)));
    result.push_back(Polynomial((-t1 * al6 * d4 + t3 * a5 * al5 + 2 * t6 * al4 - 2 * t7 - t0 * d5 + t3 * al6 * a5 * al5 * al4 + t2 * al6 * a6 * al4 - t2 * a4 + t0 * al6 * al4 * d5 - t1 * al4 * d5 - t1 * al6 * d5 - t1 * al6 * d6 - t3 * al4 * a4 + t3 * al6 * a4 - t3 * al6 * a6 + t3 * a6 * al4 - 2 * t7 * al6 * al4 - t2 * a5 * al5 * al4 + t2 * al6 * a5 * al5 - t2 * al6 * al4 * a4 + t1 * al4 * d4 - t1 * al4 * d6 - t0 * al6 * al4 * d4 + t0 * al6 * al4 * d6 - 2 * t6 * al6 - t0 * d6 - t0 * d4 + t2 * a6),(-t0*al6*a6-2*t4+t2*al4*d5-2*t4*al4*al6+t3*d4+t3*d6+t0*al4*a6-t2*al4*d4-t0*al4*a4+t0*a5*al5+t0*al6*a4+t2*al6*d5+t2*al6*d4+t2*al6*d6+t2*d6*al4-t1*a5*al5*al4-t1*al6*al4*a4+t1*al6*a5*al5+t1*al6*a6*al4-t3*al6*al4*d5+t3*al6*al4*d4-t3*al6*al4*d6-t1*a4+t1*a6+t3*d5+t0*al6*a5*al5*al4+2*t5*al4-2*t5*al6)));
    result.push_back(Polynomial((-2 * t0 - 2 * t0 * al6 * al4 + 2 * t1 * al4 - 2 * t1 * al6),(-2*t2*al4+2*t2*al6+2*t3+2*t3*al6*al4)));
    result.push_back(Polynomial((-2 * t1 - 2 * t1 * al6 * al4 - 2 * t0 * al4 + 2 * t0 * al6),(-2*t3*al4+2*t3*al6-2*t2-2*t2*al6*al4)));
    result.push_back(Polynomial((-2 * t3 * al4 + 2 * t3 * al6 - 2 * t2 - 2 * t2 * al6 * al4),(2*t0*al4-2*t0*al6+2*t1+2*t1*al6*al4)));
    result.push_back(Polynomial((-2 * t3 - 2 * t3 * al6 * al4 + 2 * t2 * al4 - 2 * t2 * al6),(-2*t0-2*t0*al6*al4+2*t1*al4-2*t1*al6)));

    return result;
  };

  std::vector<Polynomial> h2_v4q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((2 * t4 * al5 * al6 - t1 * a5 - 2 * t5 * al5 - t3 * al5 * al4 * d5 + t1 * al6 * al5 * a4 - t1 * al6 * a6 * al5 - t1 * al5 * al4 * a4 + t3 * al5 * al4 * d4 - t2 * al6 * al5 * al4 * d5 + t2 * al6 * al5 * al4 * d4 - t0 * al6 * al5 * al4 * a6 + t2 * al5 * d5 + t2 * al6 * al5 * al4 * d6 - t3 * al6 * al5 * d5 + t1 * a6 * al5 * al4 + t0 * al6 * al5 * al4 * a4 - 2 * t5 * al6 * al5 * al4 - t3 * al6 * al5 * d4 + t3 * d6 * al5 * al4 + t0 * al6 * a5-t1*al6*a5*al4+t3*al6*al5*d6-2*t4*al5*al4-t0*al5*a6+t2*al5*d4-t2*al5*d6-t0*a5*al4+t0*al5*a4),(2*t6*al5+t3*a6*al5+2*t7*al5*al4-t3*al5*a4-2*t7*al6*al5-t1*al5*d6+t3*a5*al4-t0*al6*al5*d5+t1*al6*d6*al5*al4-t3*al6*a5+t2*al5*al4*a4-t0*al6*al5*d4+t0*al5*al4*d6+t2*al5*al6*a6-t2*al5*al4*a6+2*t6*al5*al4*al6+t0*al5*al4*d4+t2*a5+t0*al6*al5*d6-t0*al5*al4*d5+t3*al6*a6*al5*al4-t1*al6*al5*al4*d5+t1*al5*d5+t1*al5*d4+t2*al6*a5*al4+t1*al6*al5*al4*d4-t3*al6*al5*al4*a4-t2*al6*al5*a4)));
    result.push_back(Polynomial((2*t4*al5+t3*al5*d4+t1*al6*a5-t1*a6*al5+t0*a5+2*t5*al6*al5-t2*al6*al5*d6-t2*al5*d6*al4-t3*d6*al5+t2*al6*al5*d5-t3*al6*al5*al4*d5+t0*al5*al4*a4-t0*al6*al5*a4-t1*a5*al4+t3*al6*al5*al4*d4+t0*al6*a5*al4+t0*al6*a6*al5+t1*al5*a4+t2*al6*al5*d4+2*t4*al5*al6*al4+t3*al6*al5*al4*d6-t0*al5*a6*al4-2*t5*al5*al4+t3*al5*d5+t1*al6*al5*al4*a4+t2*al5*al4*d5-t2*al5*al4*d4-t1*al6*al5*al4*a6),(t3*a5+t1*al6*d6*al5+t2*al6*a5+2*t7*al5-2*t6*al5*al4+t3*al6*a5*al4+t3*al5*al4*a4+2*t6*al6*al5-t2*al5*al6*a6*al4-t0*al5*d4+t0*al5*d6+t2*al6*al5*al4*a4-t1*al5*al4*d5+t0*al6*al5*al4*d5-t3*al6*al5*a4-t0*al6*al5*d6*al4-t1*al6*al5*d5-t2*a6*al5-t0*al6*al5*al4*d4+t3*al5*al6*a6-t3*al5*al4*a6+t1*al5*al4*d6+2*t7*al5*al4*al6-t1*al5*al6*d4+t1*al5*al4*d4-t0*al5*d5-t2*a5*al4+t2*al5*a4)));
    result.push_back(Polynomial((t3 * a5 + t1 * al6 * d6 * al5 + t2 * al6 * a5 + 2 * t7 * al5 - 2 * t6 * al5 * al4 + t3 * al6 * a5 * al4 + t3 * al5 * al4 * a4 + 2 * t6 * al6 * al5 - t2 * al5 * al6 * a6 * al4 - t0 * al5 * d4 + t0 * al5 * d6 + t2 * al6 * al5 * al4 * a4 - t1 * al5 * al4 * d5 + t0 * al6 * al5 * al4 * d5 - t3 * al6 * al5 * a4 - t0 * al6 * al5 * d6 * al4 - t1 * al6 * al5 * d5-t2*a6*al5-t0*al6*al5*al4*d4+t3*al5*al6*a6-t3*al5*al4*a6+t1*al5*al4*d6+2*t7*al5*al4*al6-t1*al5*al6*d4+t1*al5*al4*d4-t0*al5*d5-t2*a5*al4+t2*al5*a4),(-2*t4*al5-t3*al5*d4-t1*al6*a5+t1*a6*al5-t0*a5-2*t5*al6*al5+t2*al6*al5*d6+t2*al5*d6*al4+t3*d6*al5-t2*al6*al5*d5+t3*al6*al5*al4*d5-t0*al5*al4*a4+t0*al6*al5*a4+t1*a5*al4-t3*al6*al5*al4*d4-t0*al6*a5*al4-t0*al6*a6*al5-t1*al5*a4-t2*al6*al5*d4-2*t4*al5*al6*al4-t3*al6*al5*al4*d6+t0*al5*a6*al4+2*t5*al5*al4-t3*al5*d5-t1*al6*al5*al4*a4-t2*al5*al4*d5+t2*al5*al4*d4+t1*al6*al5*al4*a6)));
    result.push_back(Polynomial((-2 * t6 * al5 - t3 * a6 * al5 - 2 * t7 * al5 * al4 + t3 * al5 * a4 + 2 * t7 * al6 * al5 + t1 * al5 * d6 - t3 * a5 * al4 + t0 * al6 * al5 * d5 - t1 * al6 * d6 * al5 * al4 + t3 * al6 * a5 - t2 * al5 * al4 * a4 + t0 * al6 * al5 * d4 - t0 * al5 * al4 * d6 - t2 * al5 * al6 * a6 + t2 * al5 * al4 * a6 - 2 * t6 * al5 * al4 * al6 - t0 * al5 * al4 * d4 - t2 * a5 - t0 * al6 * al5 * d6 + t0 * al5 * al4 * d5 - t3 * al6 * a6 * al5 * al4 + t1 * al6 * al5 * al4 * d5 - t1 * al5 * d5 - t1 * al5 * d4-t2*al6*a5*al4-t1*al6*al5*al4*d4+t3*al6*al5*al4*a4+t2*al6*al5*a4),(2*t4*al5*al6-t1*a5-2*t5*al5-t3*al5*al4*d5+t1*al6*al5*a4-t1*al6*a6*al5-t1*al5*al4*a4+t3*al5*al4*d4-t2*al6*al5*al4*d5+t2*al6*al5*al4*d4-t0*al6*al5*al4*a6+t2*al5*d5+t2*al6*al5*al4*d6-t3*al6*al5*d5+t1*a6*al5*al4+t0*al6*al5*al4*a4-2*t5*al6*al5*al4-t3*al6*al5*d4+t3*d6*al5*al4+t0*al6*a5-t1*al6*a5*al4+t3*al6*al5*d6-2*t4*al5*al4-t0*al5*a6+t2*al5*d4-t2*al5*d6-t0*a5*al4+t0*al5*a4)));
    result.push_back(Polynomial((-2 * t0 * al5 * al4 + 2 * t0 * al6 * al5 - 2 * t1 * al5 - 2 * t1 * al6 * al5 * al4),(2*t2*al5+2*t2*al6*al5*al4+2*t3*al5*al4-2*t3*al6*al5)));
    result.push_back(Polynomial((-2*t1*al5*al4+2*t1*al6*al5+2*t0*al5+2*t0*al6*al5*al4),(2*t3*al5+2*t3*al6*al5*al4-2*t2*al5*al4+2*t2*al6*al5)));
    result.push_back(Polynomial((2 * t3 * al5 + 2 * t3 * al6 * al5 * al4 - 2 * t2 * al5 * al4 + 2 * t2 * al6 * al5),(-2*t0*al5-2*t0*al6*al5*al4+2*t1*al5*al4-2*t1*al6*al5)));
    result.push_back(Polynomial((-2 * t3 * al5 * al4 + 2 * t3 * al6 * al5 - 2 * t2 * al5 - 2 * t2 * al6 * al5 * al4),(-2*t0*al5*al4+2*t0*al6*al5-2*t1*al5-2*t1*al6*al5*al4)));

    return result;
  };

  std::vector<Polynomial> h3_v4q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((-2 * t6 * al5 - t3 * a6 * al5 + 2 * t7 * al5 * al4 - t3 * al5 * a4 + 2 * t7 * al6 * al5 + t1 * al5 * d6 + t3 * a5 * al4 + t0 * al6 * al5 * d5 + t1 * al6 * d6 * al5 * al4 + t3 * al6 * a5 - t2 * al5 * al4 * a4 + t0 * al6 * al5 * d4 + t0 * al5 * al4 * d6 - t2 * al5 * al6 * a6 - t2 * al5 * al4 * a6 + 2 * t6 * al5 * al4 * al6 + t0 * al5 * al4 * d4 - t2 * a5 - t0 * al6 * al5 * d6 - t0 * al5 * al4 * d5 + t3 * al6 * a6 * al5 * al4 - t1 * al6 * al5 * al4 * d5 - t1 * al5 * d5 - t1 * al5 * d4 + t2 * al6 * a5 * al4 + t1 * al6 * al5 * al4 * d4 + t3 * al6 * al5 * al4 * a4-t2*al6*al5*a4),(2*t4*al5*al6-t1*a5-2*t5*al5+t3*al5*al4*d5-t1*al6*al5*a4-t1*al6*a6*al5-t1*al5*al4*a4-t3*al5*al4*d4+t2*al6*al5*al4*d5-t2*al6*al5*al4*d4+t0*al6*al5*al4*a6+t2*al5*d5-t2*al6*al5*al4*d6-t3*al6*al5*d5-t1*a6*al5*al4+t0*al6*al5*al4*a4+2*t5*al6*al5*al4-t3*al6*al5*d4-t3*d6*al5*al4+t0*al6*a5+t1*al6*a5*al4+t3*al6*al5*d6+2*t4*al5*al4-t0*al5*a6+t2*al5*d4-t2*al5*d6+t0*a5*al4-t0*al5*a4)));
    result.push_back(Polynomial((-t3 * a5 - t1 * al6 * d6 * al5 - t2 * al6 * a5 - 2 * t7 * al5 - 2 * t6 * al5 * al4 + t3 * al6 * a5 * al4 - t3 * al5 * al4 * a4 - 2 * t6 * al6 * al5 - t2 * al5 * al6 * a6 * al4 + t0 * al5 * d4 - t0 * al5 * d6 - t2 * al6 * al5 * al4 * a4 - t1 * al5 * al4 * d5 + t0 * al6 * al5 * al4 * d5 - t3 * al6 * al5 * a4 - t0 * al6 * al5 * d6 * al4 + t1 * al6 * al5 * d5 + t2 * a6 * al5 - t0 * al6 * al5 * al4 * d4 - t3 * al5 * al6 * a6 - t3 * al5 * al4 * a6 + t1 * al5 * al4 * d6 + 2 * t7 * al5 * al4 * al6 + t1 * al5 * al6 * d4 + t1 * al5 * al4 * d4 + t0 * al5 * d5 - t2 * a5 * al4 + t2 * al5 * a4),(2*t4*al5+t3*al5*d4+t1*al6*a5-t1*a6*al5+t0*a5+2*t5*al6*al5-t2*al6*al5*d6+t2*al5*d6*al4-t3*d6*al5+t2*al6*al5*d5+t3*al6*al5*al4*d5+t0*al5*al4*a4+t0*al6*al5*a4+t1*a5*al4-t3*al6*al5*al4*d4-t0*al6*a5*al4+t0*al6*a6*al5-t1*al5*a4+t2*al6*al5*d4-2*t4*al5*al6*al4-t3*al6*al5*al4*d6+t0*al5*a6*al4+2*t5*al5*al4+t3*al5*d5+t1*al6*al5*al4*a4-t2*al5*al4*d5+t2*al5*al4*d4+t1*al6*al5*al4*a6)));
    result.push_back(Polynomial((t3 * al5 * d4 + 2 * t4 * al5 + t2 * al5 * al4 * d4 + t1 * a5 * al4 + t1 * al6 * al5 * al4 * a4 - t1 * al5 * a6 + t3 * al6 * al5 * al4 * d5 + t1 * al6 * al5 * al4 * a6 - t1 * al5 * a4 + t2 * d6 * al5 * al4 + t1 * al6 * a5 - t2 * al5 * al4 * d5 + t2 * al6 * al5 * d5 - t3 * al6 * al5 * al4 * d6 + 2 * t5 * al6 * al5 + t0 * a5 - t3 * d6 * al5 + t0 * al5 * al4 * a4 + t0 * al6 * al5 * a4 - t3 * al6 * al5 * al4 * d4 - t0 * al6 * a5 * al4 + t0 * al6 * a6 * al5 + t2 * al6 * al5 * d4-t2*al6*al5*d6-2*t4*al6*al5*al4+t0*al5*a6*al4+2*t5*al5*al4+t3*al5*d5),(t3*a5+t1*al6*d6*al5+t2*al6*a5+2*t7*al5+2*t6*al5*al4-t3*al6*a5*al4+t3*al5*al4*a4+2*t6*al6*al5+t2*al5*al6*a6*al4-t0*al5*d4+t0*al5*d6+t2*al6*al5*al4*a4+t1*al5*al4*d5-t0*al6*al5*al4*d5+t3*al6*al5*a4+t0*al6*al5*d6*al4-t1*al6*al5*d5-t2*a6*al5+t0*al6*al5*al4*d4+t3*al5*al6*a6+t3*al5*al4*a6-t1*al5*al4*d6-2*t7*al5*al4*al6-t1*al5*al6*d4-t1*al5*al4*d4-t0*al5*d5+t2*a5*al4-t2*al5*a4)));
    result.push_back(Polynomial((-t0 * al6 * a5 - t0 * al6 * al5 * al4 * a4+t1*a5+t2*d6*al5+t1*al5*al4*a4-t2*al5*d5-2*t5*al6*al5*al4+2*t5*al5+t2*al6*al5*al4*d6+t2*al6*al5*al4*d4+t3*al6*al5*d4+t1*al5*a6*al4+t0*al5*a4-t3*al5*al4*d5-t2*al6*al5*al4*d5-t1*al6*a5*al4+t3*al5*al4*d4+t1*al6*al5*a4+t1*al6*a6*al5-t2*al5*d4+t3*d6*al5*al4+t3*al6*al5*d5-2*t4*al5*al4-t0*a5*al4+t0*al5*a6-t3*al6*al5*d6-2*t4*al6*al5-t0*al6*al5*al4*a6),(2*t7*al6*al5+2*t7*al5*al4-t2*a5-t2*al6*a6*al5+t1*al6*al5*al4*d4+t3*al6*a5-t0*al6*d6*al5-2*t6*al5-t1*al6*al5*al4*d5+t3*a5*al4-t1*al5*d5+t1*al6*d6*al5*al4+t0*al5*al4*d4+t0*al6*al5*d4-t2*al5*al4*a6+t0*al5*al4*d6+2*t6*al5*al4*al6-t2*al6*al5*a4+t2*al6*a5*al4+t3*al5*al6*a6*al4-t2*al5*al4*a4-t3*a6*al5+t3*al6*al5*al4*a4-t3*al5*a4-t0*al5*al4*d5-t1*al5*d4+t1*al5*d6+t0*al6*al5*d5)));
    result.push_back(Polynomial((-2 * t2 * al5 + 2 * t2 * al6 * al5 * al4 + 2 * t3 * al5 * al4 + 2 * t3 * al6 * al5),(2*t0*al5*al4+2*t0*al6*al5-2*t1*al5+2*t1*al6*al5*al4)));
    result.push_back(Polynomial((-2 * t3 * al5 + 2 * t3 * al6 * al5 * al4 - 2 * t2 * al5 * al4 - 2 * t2 * al6 * al5),(2*t1*al5*al4+2*t1*al6*al5+2*t0*al5-2*t0*al6*al5*al4)));
    result.push_back(Polynomial((2 * t1 * al5 * al4 + 2 * t1 * al6 * al5 + 2 * t0 * al5 - 2 * t0 * al6 * al5 * al4),(2*t2*al5*al4+2*t2*al6*al5+2*t3*al5-2*t3*al6*al5*al4)));
    result.push_back(Polynomial((2 * t1 * al5 - 2 * t1 * al6 * al5 * al4 - 2 * t0 * al5 * al4 - 2 * t0 * al6 * al5),(-2*t2*al5+2*t2*al6*al5*al4+2*t3*al5*al4+2*t3*al6*al5)));

    return result;
  };

  std::vector<Polynomial> h4_v4q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((2 * t6 * al6 + 2 * t7 - t2 * a6 + t1 * al6 * d6 + t0 * d4 - t1 * al4 * d5 + t3 * al6 * a6-t2*a4+t3*al6*a4+t1*al4*d4+2*t6*al4+t0*d6+t3*al4*a6+t0*d5+t0*al6*al4*d6-t1*d6*al4+t1*al6*d5-t3*a5*al5+t0*al6*al4*d5-2*t7*al4*al6+t1*al6*d4-t0*al6*al4*d4+t2*al6*a6*al4+t2*al6*al4*a4-t2*al6*a5*al5+t3*al4*a4+t3*al6*a5*al5*al4-t2*a5*al5*al4),(t2*al4*d5-2*t4*al6*al4+t0*al4*a4-t1*a6+2*t4+t0*al6*a4-t0*a5*al5-t3*d4-t3*d5-t3*d6+t3*al6*al4*d4-t3*al6*al4*d6-t2*al6*d5+2*t5*al6+t1*al6*a6*al4+t2*al4*d6-t2*al6*d4-t2*al4*d4+t1*al6*al4*a4-t1*a5*al5*al4-t1*a4+2*t5*al4+t0*al6*a5*al5*al4+t0*al6*a6-t3*al6*al4*d5+t0*a6*al4-t1*al6*a5*al5-t2*al6*d6)));
    result.push_back(Polynomial((-t2 * al6 * a6 + 2 * t7 * al6 - t0 * al6 * d6 - 2 * t6 - t2 * al6 * a4 - t2 * al4 * a6 - t0 * al6 * d4 - t2 * al4 * a4 + t1 * al6 * al4 * d5 - t3 * a6 + t2 * a5 * al5 + t1 * d5 - t1 * al6 * al4 * d4 + t0 * d6 * al4+t1*d4-t0*al4*d4+2*t6*al4*al6-t0*al6*d5-t3*a4+t1*d6-t2*al6*a5*al5*al4+t1*al6*al4*d6-t3*a5*al5*al4+2*t7*al4-t3*al6*a5*al5+t3*al6*al4*a4+t0*al4*d5+t3*al6*a6*al4),(2*t5+t2*d6-2*t4*al6+t1*al6*a5*al5*al4-t3*al6*d4-t3*al4*d4+t3*al4*d6+t2*d5+t0*a6+t0*al6*a5*al5+t0*a5*al5*al4-t3*al6*d5-t3*al6*d6+t1*al4*a4+t1*al6*a4+t2*d4-2*t5*al6*al4+t3*al4*d5-t0*al6*al4*a4-t0*al6*a6*al4-2*t4*al4+t0*a4+t1*al6*a6+t2*al6*al4*d5+t1*a6*al4+t2*al6*al4*d6-t2*al6*al4*d4-t1*a5*al5)));
    result.push_back(Polynomial((2 * t5 + t2 * d6 - 2 * t4 * al6 + t1 * al6 * a5 * al5 * al4 - t3 * al6 * d4 - t3 * al4 * d4 + t3 * al4 * d6 + t2 * d5 + t0 * a6 + t0 * al6 * a5 * al5 + t0 * a5 * al5 * al4 - t3 * al6 * d5 - t3 * al6 * d6+t1*al4*a4+t1*al6*a4+t2*d4-2*t5*al6*al4+t3*al4*d5-t0*al6*al4*a4-t0*al6*a6*al4-2*t4*al4+t0*a4+t1*al6*a6+t2*al6*al4*d5+t1*a6*al4+t2*al6*al4*d6-t2*al6*al4*d4-t1*a5*al5),(t2*al6*a6-2*t7*al6+t0*al6*d6+2*t6+t2*al6*a4+t2*al4*a6+t0*al6*d4+t2*al4*a4-t1*al6*al4*d5+t3*a6-t2*a5*al5-t1*d5+t1*al6*al4*d4-t0*d6*al4-t1*d4+t0*al4*d4-2*t6*al4*al6+t0*al6*d5+t3*a4-t1*d6+t2*al6*a5*al5*al4-t1*al6*al4*d6+t3*a5*al5*al4-2*t7*al4+t3*al6*a5*al5-t3*al6*al4*a4-t0*al4*d5-t3*al6*a6*al4)));
    result.push_back(Polynomial((-t2 * al4 * d5 + 2 * t4 * al6 * al4 - t0 * al4 * a4 + t1 * a6 - 2 * t4 - t0 * al6 * a4 + t0 * a5 * al5 + t3 * d4 + t3 * d5 + t3 * d6 - t3 * al6 * al4 * d4 + t3 * al6 * al4 * d6 + t2 * al6 * d5 - 2 * t5 * al6 - t1 * al6 * a6 * al4 - t2 * al4 * d6 + t2 * al6 * d4 + t2 * al4 * d4 - t1 * al6 * al4 * a4 + t1 * a5 * al5 * al4 + t1 * a4 - 2 * t5 * al4-t0*al6*a5*al5*al4-t0*al6*a6+t3*al6*al4*d5-t0*a6*al4+t1*al6*a5*al5+t2*al6*d6),(2*t6*al6+2*t7-t2*a6+t1*al6*d6+t0*d4-t1*al4*d5+t3*al6*a6-t2*a4+t3*al6*a4+t1*al4*d4+2*t6*al4+t0*d6+t3*al4*a6+t0*d5+t0*al6*al4*d6-t1*d6*al4+t1*al6*d5-t3*a5*al5+t0*al6*al4*d5-2*t7*al4*al6+t1*al6*d4-t0*al6*al4*d4+t2*al6*a6*al4+t2*al6*al4*a4-t2*al6*a5*al5+t3*al4*a4+t3*al6*a5*al5*al4-t2*a5*al5*al4)));
    result.push_back(Polynomial((2 * t2 * al4 + 2 * t2 * al6 + 2 * t3 - 2 * t3 * al6 * al4),(2*t0-2*t0*al6*al4+2*t1*al4+2*t1*al6)));
    result.push_back(Polynomial((2 * t3 * al4 + 2 * t3 * al6 - 2 * t2 + 2 * t2 * al6 * al4),(2*t1-2*t1*al6*al4-2*t0*al4-2*t0*al6)));
    result.push_back(Polynomial((2 * t1 - 2 * t1 * al6 * al4 - 2 * t0 * al4 - 2 * t0 * al6),(2*t2-2*t2*al6*al4-2*t3*al4-2*t3*al6)));
    result.push_back(Polynomial((-2 * t1 * al4 - 2 * t1 * al6 - 2 * t0 + 2 * t0 * al6 * al4),(2*t2*al4+2*t2*al6+2*t3-2*t3*al6*al4)));

    return result;
  };

  std::vector<Polynomial> h1_v5q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((t1*d4*al5+t1*al5*d5-t1*d6*al5-t3*al6*a5-t3*al6*a4+t3*a6*al5-t1*d6*al4+t3*a6*al4+t2*a4+2*t6*al4+2*t6*al5+t1*d4*al4-t1*al4*d5-t0*al6*d4*al5-2*t7*al6*al4-2*t7*al6*al5+t2*a5-t0*al6*al5*d5+t0*al6*al4*d5-t0*al6*d4*al4+t0*al6*d6*al4+t0*al6*d6*al5-t2*a4*al5*al4-t2*a5*al5*al4+t2*al6*a6*al4+t2*al6*a6*al5+t3*al6*a4*al5*al4+t3*al6*a5*al5*al4),(-t0*al6*a4+t2*d4*al5-t0*a6*al5+t2*d5*al5-t1*a5+t1*a4-2*t5*al5+2*t5*al4+t0*al6*a5+t0*a6*al4+t2*d5*al4-t2*d4*al4-t2*d6*al5+t2*d6*al4+2*t4*al6*al5-2*t4*al6*al4-t0*al6*al5*al4*a4+t1*al5*al4*a4+t0*al6*al5*a5*al4-t1*al5*a5*al4-t1*al6*a6*al5+t1*al6*a6*al4-t3*al6*d5*al5-t3*al6*d5*al4-t3*al6*d4*al5+t3*al6*d4*al4+t3*al6*d6*al5-t3*al6*d6*al4)));
    result.push_back(Polynomial((t0*al4*d5+2*t6*al6*al4-t0*al5*d5+t2*al6*a5-t0*d4*al5+t0*d6*al4-t0*d4*al4-t2*a6*al4-t2*a6*al5+t2*al6*a4-t1*al6*d4*al5-t1*al6*al5*d5+t1*al6*al4*d5-t1*al6*d4*al4+t3*a4+t1*al6*d6*al4+t1*al6*d6*al5+t3*a5-t3*a4*al5*al4-t3*a5*al5*al4+t3*al6*a6*al5+t3*al6*a6*al4+t0*d6*al5+2*t6*al6*al5-t2*al6*a4*al5*al4-t2*al6*a5*al5*al4+2*t7*al4+2*t7*al5),(t3*d5*al5+2*t5*al6*al5+t2*al6*d5*al5+2*t4*al5+t0*a5-2*t4*al4-t1*al6*a4+t1*al6*a5-t1*a6*al5+t1*a6*al4+t3*d5*al4+t3*d4*al5-t3*d4*al4-t3*d6*al5+t3*d6*al4-2*t5*al6*al4-t0*a4-t0*al5*al4*a4-t1*al6*al5*al4*a4+t0*al5*a5*al4+t1*al6*al5*a5*al4+t0*al6*a6*al5-t0*al6*a6*al4+t2*al6*d5*al4+t2*al6*d4*al5-t2*al6*d4*al4-t2*al6*d6*al5+t2*al6*d6*al4)));
    result.push_back(Polynomial((-t3*d5*al5-2*t5*al6*al5-t2*al6*d5*al5-2*t4*al5-t0*a5-2*t4*al4-t1*al6*a4-t1*al6*a5+t1*a6*al5+t1*a6*al4+t3*d5*al4-t3*d4*al5-t3*d4*al4+t3*d6*al5+t3*d6*al4-2*t5*al6*al4-t0*a4+t0*al5*al4*a4+t1*al6*al5*al4*a4+t0*al5*a5*al4+t1*al6*al5*a5*al4-t0*al6*a6*al5-t0*al6*a6*al4+t2*al6*d5*al4-t2*al6*d4*al5-t2*al6*d4*al4+t2*al6*d6*al5+t2*al6*d6*al4),(-t0*al4*d5-2*t6*al6*al4-t0*al5*d5+t2*al6*a5-t0*d4*al5-t0*d6*al4+t0*d4*al4+t2*a6*al4-t2*a6*al5-t2*al6*a4-t1*al6*d4*al5-t1*al6*al5*d5-t1*al6*al4*d5+t1*al6*d4*al4-t3*a4-t1*al6*d6*al4+t1*al6*d6*al5+t3*a5-t3*a4*al5*al4+t3*a5*al5*al4+t3*al6*a6*al5-t3*al6*a6*al4+t0*d6*al5+2*t6*al6*al5-t2*al6*a4*al5*al4+t2*al6*a5*al5*al4-2*t7*al4+2*t7*al5)));
    result.push_back(Polynomial((t0*al6*a4+t2*d4*al5-t0*a6*al5+t2*d5*al5-t1*a5-t1*a4-2*t5*al5-2*t5*al4+t0*al6*a5-t0*a6*al4-t2*d5*al4+t2*d4*al4-t2*d6*al5-t2*d6*al4+2*t4*al6*al5+2*t4*al6*al4-t0*al6*al5*al4*a4+t1*al5*al4*a4-t0*al6*al5*a5*al4+t1*al5*a5*al4-t1*al6*a6*al5-t1*al6*a6*al4-t3*al6*d5*al5+t3*al6*d5*al4-t3*al6*d4*al5-t3*al6*d4*al4+t3*al6*d6*al5+t3*al6*d6*al4),(-t1*d4*al5-t1*al5*d5+t1*d6*al5+t3*al6*a5-t3*al6*a4-t3*a6*al5-t1*d6*al4+t3*a6*al4+t2*a4+2*t6*al4-2*t6*al5+t1*d4*al4-t1*al4*d5+t0*al6*d4*al5-2*t7*al6*al4+2*t7*al6*al5-t2*a5+t0*al6*al5*d5+t0*al6*al4*d5-t0*al6*d4*al4+t0*al6*d6*al4-t0*al6*d6*al5+t2*a4*al5*al4-t2*a5*al5*al4+t2*al6*a6*al4-t2*al6*a6*al5-t3*al6*a4*al5*al4+t3*al6*a5*al5*al4)));
    result.push_back(Polynomial((2*t2*al4+2*t2*al5-2*t3*al6*al4-2*t3*al6*al5),(2*t0*al6*al5-2*t0*al6*al4-2*t1*al5+2*t1*al4)));
    result.push_back(Polynomial((2*t3*al4+2*t3*al5+2*t2*al6*al4+2*t2*al6*al5),(2*t1*al6*al5-2*t1*al6*al4+2*t0*al5-2*t0*al4)));
    result.push_back(Polynomial((-2*t0*al4-2*t0*al5-2*t1*al6*al4-2*t1*al6*al5),(2*t2*al6*al5-2*t2*al6*al4+2*t3*al5-2*t3*al4)));
    result.push_back(Polynomial((-2*t1*al4-2*t1*al5+2*t0*al6*al4+2*t0*al6*al5),(2*t3*al6*al5-2*t3*al6*al4-2*t2*al5+2*t2*al4)));

    return result;
  };

  std::vector<Polynomial> h2_v5q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((-t3*al5*a4-t3*al5*a5+2*t6*al6+2*t7+t1*al6*d4-t3*al4*a5-2*t7*al5*al4+t1*al6*d6-t3*al4*a4+t1*al6*d5+t3*al6*a6+t0*d4+t0*d6-t2*a6-t0*d4*al5*al4+t0*al5*al4*d5-t0*d6*al5*al4-t2*al6*al5*a5-t2*al6*al4*a4-t2*al6*al5*a4-t2*al6*al4*a5+t2*a6*al5*al4-2*t6*al6*al5*al4+t0*d5-t1*al6*d4*al5*al4+t1*al6*al5*al4*d5-t1*al6*d6*al5*al4-t3*al6*a6*al5*al4),(-t0*al6*a6+t1*a6*al5*al4+t2*al6*d4+t3*d5-2*t5*al6+t3*d6-2*t4+t1*a6+t2*al6*d4*al5*al4+t0*al4*a4+t2*al6*d5+t2*al6*d6-2*t4*al5*al4+t3*d4-t0*al6*a6*al5*al4+t1*al6*al4*a4-t2*al6*al5*d5*al4+t2*al6*d6*al5*al4-t3*al5*d5*al4+t3*d4*al5*al4+t3*d6*al5*al4-2*t5*al6*al5*al4-t0*al5*a4+t0*al5*a5-t0*al4*a5-t1*al6*al5*a4-t1*al6*al4*a5+t1*al6*al5*a5)));
    result.push_back(Polynomial((t2*al4*a5-2*t6-t1*d4*al5*al4+t2*al5*a4-t2*al6*a6-t0*al6*d5+2*t6*al5*al4-t0*al6*d4+2*t7*al6+t1*d4+t1*d6+t2*al4*a4-t0*al6*d6+t2*al5*a5+t1*al5*al4*d5-t1*d6*al5*al4-t3*al6*al5*a5-t3*al6*al4*a5-t3*al6*al4*a4-t3*al6*al5*a4-t3*a6+t3*a6*al5*al4-2*t7*al6*al5*al4+t1*d5+t0*al6*d4*al5*al4-t0*al6*al5*al4*d5+t0*al6*d6*al5*al4+t2*al6*a6*al5*al4),(-2*t5-t2*d6-t2*d4-t0*a6-t2*d5-t1*al6*a6-t1*al6*a6*al5*al4+t1*al4*a4+t3*al6*d5+t3*al6*d4+2*t4*al6+t3*al6*d6-2*t5*al5*al4-t0*al6*al4*a4-t0*a6*al5*al4+t3*al6*d4*al5*al4-t3*al6*al5*d5*al4+t3*al6*d6*al5*al4+t2*al5*d5*al4-t2*d4*al5*al4-t2*d6*al5*al4+2*t4*al6*al5*al4+t0*al6*al5*a4+t0*al6*al4*a5-t0*al6*al5*a5-t1*al5*a4-t1*al4*a5+t1*al5*a5)));
    result.push_back(Polynomial((2*t5+t2*d6+t2*d4+t0*a6+t2*d5+t1*al6*a6-t1*al6*a6*al5*al4-t1*al4*a4-t3*al6*d5-t3*al6*d4-2*t4*al6-t3*al6*d6-2*t5*al5*al4+t0*al6*al4*a4-t0*a6*al5*al4+t3*al6*d4*al5*al4-t3*al6*al5*d5*al4+t3*al6*d6*al5*al4+t2*al5*d5*al4-t2*d4*al5*al4-t2*d6*al5*al4+2*t4*al6*al5*al4+t0*al6*al5*a4+t0*al6*al4*a5+t0*al6*al5*a5-t1*al5*a4-t1*al4*a5-t1*al5*a5),(-t2*al4*a5-2*t6+t1*d4*al5*al4-t2*al5*a4-t2*al6*a6-t0*al6*d5-2*t6*al5*al4-t0*al6*d4+2*t7*al6+t1*d4+t1*d6+t2*al4*a4-t0*al6*d6+t2*al5*a5-t1*al5*al4*d5+t1*d6*al5*al4-t3*al6*al5*a5+t3*al6*al4*a5-t3*al6*al4*a4+t3*al6*al5*a4-t3*a6-t3*a6*al5*al4+2*t7*al6*al5*al4+t1*d5-t0*al6*d4*al5*al4+t0*al6*al5*al4*d5-t0*al6*d6*al5*al4-t2*al6*a6*al5*al4)));
    result.push_back(Polynomial((-t0*al6*a6-t1*a6*al5*al4+t2*al6*d4+t3*d5-2*t5*al6+t3*d6-2*t4+t1*a6-t2*al6*d4*al5*al4+t0*al4*a4+t2*al6*d5+t2*al6*d6+2*t4*al5*al4+t3*d4+t0*al6*a6*al5*al4+t1*al6*al4*a4+t2*al6*al5*d5*al4-t2*al6*d6*al5*al4+t3*al5*d5*al4-t3*d4*al5*al4-t3*d6*al5*al4+2*t5*al6*al5*al4+t0*al5*a4+t0*al5*a5+t0*al4*a5+t1*al6*al5*a4+t1*al6*al4*a5+t1*al6*al5*a5),(-t3*al5*a4+t3*al5*a5-2*t6*al6-2*t7-t1*al6*d4-t3*al4*a5-2*t7*al5*al4-t1*al6*d6+t3*al4*a4-t1*al6*d5-t3*al6*a6-t0*d4-t0*d6+t2*a6-t0*d4*al5*al4+t0*al5*al4*d5-t0*d6*al5*al4+t2*al6*al5*a5+t2*al6*al4*a4-t2*al6*al5*a4-t2*al6*al4*a5+t2*a6*al5*al4-2*t6*al6*al5*al4-t0*d5-t1*al6*d4*al5*al4+t1*al6*al5*al4*d5-t1*al6*d6*al5*al4-t3*al6*a6*al5*al4)));
    result.push_back(Polynomial((-2*t2*al6*al5*al4+2*t2*al6-2*t3*al5*al4+2*t3),(-2*t0*al5*al4-2*t0-2*t1*al6*al5*al4-2*t1*al6)));
    result.push_back(Polynomial((-2*t3*al6*al5*al4+2*t3*al6+2*t2*al5*al4-2*t2),(-2*t1*al5*al4-2*t1+2*t0*al6*al5*al4+2*t0*al6)));
    result.push_back(Polynomial((2*t0*al6*al5*al4-2*t0*al6-2*t1*al5*al4+2*t1),(-2*t2*al5*al4-2*t2+2*t3*al6*al5*al4+2*t3*al6)));
    result.push_back(Polynomial((2*t1*al6*al5*al4-2*t1*al6+2*t0*al5*al4-2*t0),(-2*t3*al5*al4-2*t3-2*t2*al6*al5*al4-2*t2*al6)));

    return result;
  };

  std::vector<Polynomial> h3_v5q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((t1*al6*a6*al5+2*t5*al5+t0*a6*al5+t3*al6*d4*al5+t3*al6*al5*d5-t3*al6*d6*al5-t2*al5*d5-t2*d4*al5+t2*d6*al5-2*t4*al6*al5+t0*a6*al4+t2*d6*al4-2*t4*al6*al4+t2*al4*d5-t0*al6*a4-t2*d4*al4-t0*al6*a5-t3*al6*al4*d5+t3*al6*d4*al4+t1*a4-t3*al6*d6*al4-t1*al5*a4*al4+t1*a5-t1*al5*a5*al4+t1*al6*a6*al4+2*t5*al4+t0*al6*al5*a4*al4+t0*al6*al5*a5*al4),(-2*t7*al6*al5+t1*d6*al4-t1*d6*al5-t1*d4*al4+t1*d4*al5+t1*d5*al4+t1*al5*d5-t3*a6*al4-t3*al6*a5+t3*a6*al5+2*t7*al6*al4-t2*a4+t2*a5+2*t6*al5+t0*al6*d6*al5-2*t6*al4-t0*al6*d6*al4+t3*al6*al5*al4*a4+t3*al6*a4-t3*al6*al5*a5*al4-t0*al6*al5*d5-t2*al5*al4*a4+t2*al5*a5*al4-t2*al6*a6*al4+t2*al6*a6*al5-t0*al6*d4*al5-t0*al6*d5*al4+t0*al6*d4*al4)));
    result.push_back(Polynomial((-2*t4*al5+t1*a6*al5-t3*al5*d5-t3*d4*al4+t0*a5*al5*al4-t1*al6*a4+t2*al6*d6*al4+t1*al6*a5*al5*al4-2*t5*al6*al4+t3*al4*d5-t2*al6*d4*al4+t1*a6*al4-t0*al6*a6*al5-t1*al6*a5-t0*a4+t2*al6*d6*al5+t2*al6*al4*d5-t2*al6*al5*d5-t3*d4*al5-t0*a5+t3*d6*al4+t0*a4*al5*al4-t2*al6*d4*al5-t0*al6*a6*al4-2*t4*al4+t1*al6*a4*al5*al4+t3*d6*al5-2*t5*al6*al5),(-t0*d6*al4+t2*al6*al5*a5*al4+t2*a6*al4+t0*d4*al4-t0*d5*al4-2*t6*al6*al4-t1*al6*d5*al5-t2*al6*al5*al4*a4-t1*al6*d5*al4-t3*a4-t2*a6*al5-t0*d5*al5+2*t6*al6*al5+2*t7*al5+t3*a5+t3*al6*a6*al5-t1*al6*d6*al4-2*t7*al4-t2*al6*a4+t0*d6*al5+t3*al5*a5*al4-t1*al6*d4*al5-t3*al6*a6*al4-t0*d4*al5+t2*al6*a5+t1*al6*d4*al4-t3*al5*al4*a4+t1*al6*d6*al5)));
    result.push_back(Polynomial((-t0*d6*al4+t2*al6*al5*a5*al4+t2*a6*al4+t0*d4*al4-t0*d5*al4-2*t6*al6*al4+t1*al6*d5*al5+t2*al6*al5*al4*a4-t1*al6*d5*al4-t3*a4+t2*a6*al5+t0*d5*al5-2*t6*al6*al5-2*t7*al5-t3*a5-t3*al6*a6*al5-t1*al6*d6*al4-2*t7*al4-t2*al6*a4-t0*d6*al5+t3*al5*a5*al4+t1*al6*d4*al5-t3*al6*a6*al4+t0*d4*al5-t2*al6*a5+t1*al6*d4*al4+t3*al5*al4*a4-t1*al6*d6*al5),(-2*t4*al5+t1*a6*al5-t3*al5*d5+t3*d4*al4-t0*a5*al5*al4+t1*al6*a4-t2*al6*d6*al4-t1*al6*a5*al5*al4+2*t5*al6*al4-t3*al4*d5+t2*al6*d4*al4-t1*a6*al4-t0*al6*a6*al5-t1*al6*a5+t0*a4+t2*al6*d6*al5-t2*al6*al4*d5-t2*al6*al5*d5-t3*d4*al5-t0*a5-t3*d6*al4+t0*a4*al5*al4-t2*al6*d4*al5+t0*al6*a6*al4+2*t4*al4+t1*al6*a4*al5*al4+t3*d6*al5-2*t5*al6*al5)));
    result.push_back(Polynomial((-2*t7*al6*al5-t1*d6*al4-t1*d6*al5+t1*d4*al4+t1*d4*al5-t1*d5*al4+t1*al5*d5+t3*a6*al4-t3*al6*a5+t3*a6*al5-2*t7*al6*al4+t2*a4+t2*a5+2*t6*al5+t0*al6*d6*al5+2*t6*al4+t0*al6*d6*al4+t3*al6*al5*al4*a4-t3*al6*a4+t3*al6*al5*a5*al4-t0*al6*al5*d5-t2*al5*al4*a4-t2*al5*a5*al4+t2*al6*a6*al4+t2*al6*a6*al5-t0*al6*d4*al5+t0*al6*d5*al4-t0*al6*d4*al4),(-t1*al6*a6*al5-2*t5*al5-t0*a6*al5-t3*al6*d4*al5-t3*al6*al5*d5+t3*al6*d6*al5+t2*al5*d5+t2*d4*al5-t2*d6*al5+2*t4*al6*al5+t0*a6*al4+t2*d6*al4-2*t4*al6*al4+t2*al4*d5-t0*al6*a4-t2*d4*al4+t0*al6*a5-t3*al6*al4*d5+t3*al6*d4*al4+t1*a4-t3*al6*d6*al4+t1*al5*a4*al4-t1*a5-t1*al5*a5*al4+t1*al6*a6*al4+2*t5*al4-t0*al6*al5*a4*al4+t0*al6*al5*a5*al4)));
    result.push_back(Polynomial((-2*t0*al6*al4-2*t0*al6*al5+2*t1*al4+2*t1*al5),(2*t2*al5-2*t2*al4-2*t3*al6*al5+2*t3*al6*al4)));
    result.push_back(Polynomial((-2*t0*al4-2*t0*al5-2*t1*al6*al4-2*t1*al6*al5),(2*t2*al6*al5-2*t2*al6*al4+2*t3*al5-2*t3*al4)));
    result.push_back(Polynomial((-2*t2*al6*al4-2*t2*al6*al5-2*t3*al4-2*t3*al5),(-2*t0*al5+2*t0*al4-2*t1*al6*al5+2*t1*al6*al4)));
    result.push_back(Polynomial((2*t2*al4+2*t2*al5-2*t3*al6*al4-2*t3*al6*al5),(2*t0*al6*al5-2*t0*al6*al4-2*t1*al5+2*t1*al4)));

    return result;
  };

  std::vector<Polynomial> h4_v5q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((-t2*al6*d6-t3*d6-t0*a5*al5-t1*a6+2*t5*al6-t0*a4*al5-t3*d4-t2*al6*d4-t1*al6*a5*al5+2*t4-t3*d5+t3*d6*al5*al4+t2*al6*d4*al5*al4+t3*d4*al5*al4-t3*al5*al4*d5-t0*al6*a6*al5*al4-t1*al6*a4*al5-t1*al6*al4*a5-t0*al4*a4+t2*al6*d6*al5*al4-t1*al6*al4*a4-t0*al4*a5+t1*a6*al5*al4-2*t5*al6*al5*al4-2*t4*al5*al4-t2*al6*al5*al4*d5-t2*al6*d5+t0*al6*a6),(2*t6*al6+t0*d5-t2*a6-t2*al6*al5*a5+t1*al6*d5+t3*al6*a6+t1*al6*d6+t0*d4+2*t7-t3*al5*a5+t1*al6*d4+t3*al6*a6*al5*al4-t2*a6*al5*al4+t0*d6+2*t6*al6*al5*al4+t2*al6*al4*a5-t1*al6*d5*al5*al4-t0*d5*al5*al4+t3*al4*a5+t0*d4*al5*al4-t2*al6*al4*a4-t3*al4*a4+t1*al6*d6*al5*al4+t3*a4*al5+t2*al6*a4*al5+t1*al6*d4*al5*al4+t0*d6*al5*al4+2*t7*al5*al4)));
    result.push_back(Polynomial((2*t5-2*t4*al6+t3*al6*d4*al5*al4-t2*d4*al5*al4-t1*al6*a6*al5*al4+2*t4*al6*al5*al4+t0*a6+t2*d5+t2*d4+t0*al6*al5*a5-t3*al6*d4+t0*al6*al4*a5-t0*a6*al5*al4+t2*al5*al4*d5-t2*d6*al5*al4+t0*al6*al4*a4+t3*al6*d6*al5*al4-t3*al6*al5*al4*d5-t1*al4*a4-t1*al5*a4+t0*al6*al5*a4-t3*al6*d5-t1*al4*a5-2*t5*al5*al4+t2*d6-t3*al6*d6-t1*al5*a5+t1*al6*a6),(-t0*al6*d6-2*t6-t3*al6*a5*al5-t0*al6*d4+t1*d5+2*t7*al6-t2*al6*a6-2*t6*al5*al4+2*t7*al6*al5*al4+t2*al4*a4+t1*d6*al5*al4-t0*al6*d5-t2*al6*a6*al5*al4+t1*d4*al5*al4+t1*d4+t1*d6-t3*a6*al5*al4-t3*al6*al4*a4-t2*al4*a5+t2*a5*al5+t0*al6*al5*d5*al4-t0*al6*d4*al5*al4-t1*al5*d5*al4-t0*al6*d6*al5*al4+t3*al6*al4*a5+t3*al6*a4*al5-t2*a4*al5-t3*a6)));
    result.push_back(Polynomial((t0*al6*d6+2*t6+t3*al6*a5*al5+t0*al6*d4-t1*d5-2*t7*al6+t2*al6*a6-2*t6*al5*al4+2*t7*al6*al5*al4-t2*al4*a4+t1*d6*al5*al4+t0*al6*d5-t2*al6*a6*al5*al4+t1*d4*al5*al4-t1*d4-t1*d6-t3*a6*al5*al4+t3*al6*al4*a4-t2*al4*a5-t2*a5*al5+t0*al6*al5*d5*al4-t0*al6*d4*al5*al4-t1*al5*d5*al4-t0*al6*d6*al5*al4+t3*al6*al4*a5+t3*al6*a4*al5-t2*a4*al5+t3*a6),(2*t5-2*t4*al6-t3*al6*d4*al5*al4+t2*d4*al5*al4+t1*al6*a6*al5*al4-2*t4*al6*al5*al4+t0*a6+t2*d5+t2*d4+t0*al6*al5*a5-t3*al6*d4-t0*al6*al4*a5+t0*a6*al5*al4-t2*al5*al4*d5+t2*d6*al5*al4+t0*al6*al4*a4-t3*al6*d6*al5*al4+t3*al6*al5*al4*d5-t1*al4*a4+t1*al5*a4-t0*al6*al5*a4-t3*al6*d5+t1*al4*a5+2*t5*al5*al4+t2*d6-t3*al6*d6-t1*al5*a5+t1*al6*a6)));
    result.push_back(Polynomial((2*t6*al6+t0*d5-t2*a6-t2*al6*al5*a5+t1*al6*d5+t3*al6*a6+t1*al6*d6+t0*d4+2*t7-t3*al5*a5+t1*al6*d4-t3*al6*a6*al5*al4+t2*a6*al5*al4+t0*d6-2*t6*al6*al5*al4-t2*al6*al4*a5+t1*al6*d5*al5*al4+t0*d5*al5*al4-t3*al4*a5-t0*d4*al5*al4-t2*al6*al4*a4-t3*al4*a4-t1*al6*d6*al5*al4-t3*a4*al5-t2*al6*a4*al5-t1*al6*d4*al5*al4-t0*d6*al5*al4-2*t7*al5*al4),(t2*al6*d6+t3*d6+t0*a5*al5+t1*a6-2*t5*al6-t0*a4*al5+t3*d4+t2*al6*d4+t1*al6*a5*al5-2*t4+t3*d5+t3*d6*al5*al4+t2*al6*d4*al5*al4+t3*d4*al5*al4-t3*al5*al4*d5-t0*al6*a6*al5*al4-t1*al6*a4*al5-t1*al6*al4*a5+t0*al4*a4+t2*al6*d6*al5*al4+t1*al6*al4*a4-t0*al4*a5+t1*a6*al5*al4-2*t5*al6*al5*al4-2*t4*al5*al4-t2*al6*al5*al4*d5+t2*al6*d5-t0*al6*a6)));
    result.push_back(Polynomial((-2*t0*al5*al4+2*t0-2*t1*al6*al5*al4+2*t1*al6),(2*t2*al6*al5*al4+2*t2*al6+2*t3*al5*al4+2*t3)));
    result.push_back(Polynomial((2*t0*al6*al5*al4-2*t0*al6-2*t1*al5*al4+2*t1),(-2*t2*al5*al4-2*t2+2*t3*al6*al5*al4+2*t3*al6)));
    result.push_back(Polynomial((-2*t2*al5*al4+2*t2+2*t3*al6*al5*al4-2*t3*al6),(-2*t0*al6*al5*al4-2*t0*al6+2*t1*al5*al4+2*t1)));
    result.push_back(Polynomial((-2*t2*al6*al5*al4+2*t2*al6-2*t3*al5*al4+2*t3),(-2*t0*al5*al4-2*t0-2*t1*al6*al5*al4-2*t1*al6)));

    return result;
  };

  std::vector<Polynomial> h1_v6q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((t3*d5+t0*al6*a4*al4*al5+t3*d6+t0*a4*al4-2*t4*al5*al6-2*t5*al6+t0*al5*a6-t2*d5*al5+2*t5*al5-t2*d4*al5+t2*al6*d4+t2*al6*d5+t2*al5*d6-2*t4+t2*al6*d6-t0*al5*a5+t0*al6*a5-t0*al6*a6-t1*al6*al5*a5-t1*a4*al4*al5+t3*d4+t1*al6*a4*al4+t1*al6*al5*a6+t3*al6*d5*al5+t3*al6*d4*al5-t3*al5*al6*d6-t1*a5+t1*a6),(t1*d5*al5-t0*al6*d5*al5+t3*al6*a5+2*t6*al5-t2*a5+t0*d4+t0*d5+t1*d4*al5+t1*al6*d4+t1*al6*d5-t3*a4*al4+t3*al5*a5+2*t7-t0*al6*d4*al5-t2*a4*al4*al5+t2*al6*al5*a5-t2*al6*a4*al4+t3*al6*a4*al4*al5+t1*al6*d6-t1*al5*d6-2*t7*al5*al6+t2*al6*al5*a6+t0*al6*al5*d6+t3*al6*a6+t3*al5*a6+t0*d6+2*t6*al6-t2*a6)));
    result.push_back(Polynomial((t1 * a4 * al4 + t3 * al6 * d5 + t3 * al6 * d4 - 2 * t5 - t0 * al6 * a4 * al4 - t2 * d5 - t2 * d4 - 2 * t4 * al5 - t3 * d5 * al5 - 2 * t5 * al5 * al6 + t1 * al5 * a6 - t1 * al6 * a6 - t1 * al5 * a5 + t1 * al6 * a5 - t3 * d4 * al5 + t0 * a4 * al4 * al5 + t0 * a5 - t0 * a6 + t3 * al6 * d6 + t3 * al5 * d6 + t0 * al6 * al5 * a5 - t0 * al6 * al5 * a6 - t2 * al6 * d5 * al5 - t2 * al6 * d4 * al5 + t2 * al6 * al5 * d6 + t1 * al6 * a4 * al4 * al5 - t2 * d6 + 2 * t4 * al6),(-2*t6-t0*al6*d5+t1*d4+t1*d5+2*t7*al5-t0*al6*d4+t2*a4*al4-t3*al6*a4*al4-t3*a5-t1*al6*d5*al5-t0*d4*al5-t0*d5*al5-t2*al5*a5-t2*al6*a5-t1*al6*d4*al5-t3*a4*al4*al5+t3*al5*al6*a5-t2*al6*a4*al4*al5-t2*al6*a6+2*t6*al5*al6-t2*al5*a6+t3*al5*al6*a6+t1*al6*al5*d6-t0*al6*d6+t0*al5*d6-t3*a6+t1*d6+2*t7*al6)));
    result.push_back(Polynomial((-2*t6-t0*al6*d5+t1*d4+t1*d5-2*t7*al5-t0*al6*d4+t2*a4*al4-t3*al6*a4*al4+t3*a5+t1*al6*d5*al5+t0*d4*al5+t0*d5*al5-t2*al5*a5+t2*al6*a5+t1*al6*d4*al5+t3*a4*al4*al5+t3*al5*al6*a5+t2*al6*a4*al4*al5-t2*al6*a6-2*t6*al5*al6+t2*al5*a6-t3*al5*al6*a6-t1*al6*al5*d6-t0*al6*d6-t0*al5*d6-t3*a6+t1*d6+2*t7*al6),(-t1*a4*al4-t3*al6*d5-t3*al6*d4+2*t5+t0*al6*a4*al4+t2*d5+t2*d4-2*t4*al5-t3*d5*al5-2*t5*al5*al6+t1*al5*a6+t1*al6*a6+t1*al5*a5+t1*al6*a5-t3*d4*al5+t0*a4*al4*al5+t0*a5+t0*a6-t3*al6*d6+t3*al5*d6-t0*al6*al5*a5-t0*al6*al5*a6-t2*al6*d5*al5-t2*al6*d4*al5+t2*al6*al5*d6+t1*al6*a4*al4*al5+t2*d6-2*t4*al6)));
    result.push_back(Polynomial((t1*d5*al5-t0*al6*d5*al5+t3*al6*a5+2*t6*al5-t2*a5-t0*d4-t0*d5+t1*d4*al5-t1*al6*d4-t1*al6*d5+t3*a4*al4-t3*al5*a5-2*t7-t0*al6*d4*al5-t2*a4*al4*al5-t2*al6*al5*a5+t2*al6*a4*al4+t3*al6*a4*al4*al5-t1*al6*d6-t1*al5*d6-2*t7*al5*al6+t2*al6*al5*a6+t0*al6*al5*d6-t3*al6*a6+t3*al5*a6-t0*d6-2*t6*al6+t2*a6),(t3*d5-t0*al6*a4*al4*al5+t3*d6+t0*a4*al4+2*t4*al5*al6-2*t5*al6-t0*al5*a6+t2*d5*al5-2*t5*al5+t2*d4*al5+t2*al6*d4+t2*al6*d5-t2*al5*d6-2*t4+t2*al6*d6-t0*al5*a5-t0*al6*a5-t0*al6*a6-t1*al6*al5*a5+t1*a4*al4*al5+t3*d4+t1*al6*a4*al4-t1*al6*al5*a6-t3*al6*d5*al5-t3*al6*d4*al5+t3*al5*al6*d6+t1*a5+t1*a6)));
    result.push_back(Polynomial((-2 * t0 - 2 * t0 * al6 * al5 + 2 * t1 * al5 - 2 * t1 * al6),(2*t2*al5+2*t2*al6+2*t3-2*t3*al5*al6)));
    result.push_back(Polynomial((-2 * t1 - 2 * t1 * al6 * al5 - 2 * t0 * al5 + 2 * t0 * al6),(2*t3*al5+2*t3*al6-2*t2+2*t2*al6*al5)));
    result.push_back(Polynomial((-2 * t2 - 2 * t2 * al6 * al5 - 2 * t3 * al5 + 2 * t3 * al6),(-2*t0*al5-2*t0*al6+2*t1-2*t1*al6*al5)));
    result.push_back(Polynomial((-2 * t3 - 2 * t3 * al5 * al6 + 2 * t2 * al5 - 2 * t2 * al6),(-2*t1*al5-2*t1*al6-2*t0+2*t0*al6*al5)));

    return result;
  };

  std::vector<Polynomial> h2_v6q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((-2*t5*al4+t3*al6*al4*d5+t0*al4*a5-t1*a4-t0*al4*a6+t3*al4*al5*d6-2*t5*al4*al5*al6+t2*al4*d4+2*t4*al4*al6-t2*al4*d5-t0*al6*al4*al5*a6+t1*al6*al4*a5-t1*al6*al4*a6+t1*al4*al5*a6-t0*a4*al5-t1*al4*al5*a5-t1*al6*a4*al5+t3*al6*al4*d6-t3*al6*al4*d4-t3*al4*d5*al5+t3*al4*d4*al5+t0*al6*al4*al5*a5+t2*al4*al5*al6*d6-t2*al4*d6-2*t4*al4*al5-t2*al6*al4*d5*al5+t2*al6*al4*d4*al5+t0*al6*a4),(t3*al6*a4-t2*a4-t0*al4*d5*al5+t0*al6*al4*d4-t0*al6*al4*d5-t2*al4*al5*a5+2*t7*al4*al5-t2*al6*al4*a5+t0*al4*d4*al5+t3*a4*al5-t3*al4*a5-t1*al4*d4+t1*al4*d5-t0*al4*al6*d6-t2*al4*al6*a6+t0*al4*al5*d6-t2*al4*al5*a6+2*t6*al4*al5*al6+t1*al4*d6+2*t7*al4*al6-t3*al4*a6-t1*al6*al4*d5*al5-2*t6*al4+t2*al6*a4*al5+t3*al6*al4*al5*a5+t1*al6*al4*d4*al5+t3*al4*al5*al6*a6+t1*al4*al5*al6*d6)));
    result.push_back(Polynomial((t1*al6*a4-2*t5*al4*al5-t3*al4*d5+t3*al4*d4+t0*a4-t1*al4*a6+2*t5*al4*al6+2*t4*al4+t0*al6*a4*al5-t0*al6*al4*a5+t0*al6*al4*a6-t0*al4*al5*a6-t2*al4*d4*al5+t2*al4*d5*al5+t3*al6*al4*al5*d6-t2*al6*al4*d6-t2*al4*al5*d6+2*t4*al4*al5*al6-t3*al4*d6-t1*al6*al4*al5*a6-t2*al6*al4*d5+t2*al6*al4*d4+t3*al6*al4*d4*al5+t1*al6*al4*al5*a5+t0*al4*al5*a5-t1*a4*al5+t1*al4*a5-t3*al6*al4*d5*al5),(-2*t7*al4-t2*al6*a4+t0*al4*d4-t0*al4*d5-t3*a4-t0*al6*al4*d4*al5-t3*al4*al5*a5-2*t6*al4*al5-t1*al6*al4*d5-t2*al4*al5*al6*a5-t3*al6*al4*a5+t3*al6*a4*al5-t0*al4*d6-t2*a4*al5+t2*al4*a5+t1*al4*d4*al5-t1*al4*d5*al5+t1*al6*al4*d4+t0*al6*al4*d5*al5+t1*al4*al5*d6-t1*al4*al6*d6-t3*al4*al5*a6-t3*al4*al6*a6+2*t7*al4*al5*al6-2*t6*al4*al6+t2*al4*a6-t0*al4*al5*al6*d6-t2*al4*al5*al6*a6)));
    result.push_back(Polynomial((2*t7*al4+t2*al6*a4-t0*al4*d4+t0*al4*d5+t3*a4-t0*al6*al4*d4*al5+t3*al4*al5*a5-2*t6*al4*al5+t1*al6*al4*d5+t2*al4*al5*al6*a5-t3*al6*al4*a5+t3*al6*a4*al5+t0*al4*d6-t2*a4*al5+t2*al4*a5+t1*al4*d4*al5-t1*al4*d5*al5-t1*al6*al4*d4+t0*al6*al4*d5*al5+t1*al4*al5*d6+t1*al4*al6*d6-t3*al4*al5*a6+t3*al4*al6*a6+2*t7*al4*al5*al6+2*t6*al4*al6-t2*al4*a6-t0*al4*al5*al6*d6-t2*al4*al5*al6*a6),(t1*al6*a4+2*t5*al4*al5-t3*al4*d5+t3*al4*d4+t0*a4-t1*al4*a6+2*t5*al4*al6+2*t4*al4-t0*al6*a4*al5+t0*al6*al4*a5+t0*al6*al4*a6+t0*al4*al5*a6+t2*al4*d4*al5-t2*al4*d5*al5-t3*al6*al4*al5*d6-t2*al6*al4*d6+t2*al4*al5*d6-2*t4*al4*al5*al6-t3*al4*d6+t1*al6*al4*al5*a6-t2*al6*al4*d5+t2*al6*al4*d4-t3*al6*al4*d4*al5+t1*al6*al4*al5*a5+t0*al4*al5*a5+t1*a4*al5-t1*al4*a5+t3*al6*al4*d5*al5)));
    result.push_back(Polynomial((t3 * al6 * a4 - t2 * a4 + t0 * al4 * d5 * al5 + t0 * al6 * al4 * d4 - t0 * al6 * al4 * d5 - t2 * al4 * al5 * a5 - 2 * t7 * al4 * al5 + t2 * al6 * al4 * a5 - t0 * al4 * d4 * al5 - t3 * a4 * al5 + t3 * al4 * a5 - t1 * al4 * d4 + t1 * al4 * d5 - t0 * al4 * al6 * d6 - t2 * al4 * al6 * a6 - t0 * al4 * al5 * d6 + t2 * al4 * al5 * a6 - 2 * t6 * al4 * al5 * al6 + t1 * al4 * d6 + 2 * t7 * al4 * al6 - t3 * al4 * a6 + t1 * al6 * al4 * d5 * al5 - 2 * t6 * al4 - t2 * al6 * a4 * al5 + t3 * al6 * al4 * al5 * a5 - t1 * al6 * al4 * d4 * al5 - t3 * al4 * al5 * al6 * a6 - t1 * al4 * al5 * al6 * d6),(t1*al6*al4*a5+t1*al4*al5*a5+t3*al4*al5*d6-t1*al6*a4*al5-2*t4*al4*al6+t0*al4*a6-t0*al6*a4+t0*al4*a5-t0*a4*al5+t2*al4*d6-t2*al4*d4+t2*al4*d5-2*t4*al4*al5+t1*al4*al5*a6+t1*al4*al6*a6+t1*a4+2*t5*al4-t2*al6*al4*d5*al5+t2*al6*al4*d4*al5-t0*al6*al4*al5*a5-t3*al4*d5*al5+t3*al4*d4*al5-2*t5*al4*al5*al6-t3*al4*al6*d6-t0*al4*al5*al6*a6+t2*al4*al5*al6*d6+t3*al6*al4*d4-t3*al6*al4*d5)));
    result.push_back(Polynomial((-2 * t0 * al4 * al5 + 2 * t0 * al4 * al6 - 2 * t1 * al4 - 2 * t1 * al4 * al5 * al6),(-2*t2*al4+2*t2*al4*al5*al6+2*t3*al4*al5+2*t3*al4*al6)));
    result.push_back(Polynomial((-2 * t1 * al4 * al5 + 2 * t1 * al4 * al6 + 2 * t0 * al4 + 2 * t0 * al4 * al5 * al6),(-2*t3*al4+2*t3*al4*al5*al6-2*t2*al4*al5-2*t2*al4*al6)));
    result.push_back(Polynomial((-2*t2*al4*al5+2*t2*al4*al6+2*t3*al4+2*t3*al4*al5*al6),(2*t0*al4-2*t0*al4*al5*al6+2*t1*al4*al5+2*t1*al4*al6)));
    result.push_back(Polynomial((-2*t3*al4*al5+2*t3*al4*al6-2*t2*al4-2*t2*al4*al5*al6),(2*t1*al4-2*t1*al4*al5*al6-2*t0*al4*al5-2*t0*al4*al6)));

    return result;
  };

  std::vector<Polynomial> h3_v6q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((t3 * al6 * a4 - t2 * a4 + t0 * al4 * d5 * al5 + t0 * al6 * al4 * d4 - t0 * al6 * al4 * d5 - t2 * al4 * al5 * a5 - 2 * t7 * al4 * al5 + t2 * al6 * al4 * a5 - t0 * al4 * d4 * al5 - t3 * a4 * al5 + t3 * al4 * a5 - t1 * al4 * d4 + t1 * al4 * d5 - t0 * al4 * al6 * d6 - t2 * al4 * al6 * a6 - t0 * al4 * al5 * d6 + t2 * al4 * al5 * a6 - 2 * t6 * al4 * al5 * al6 + t1 * al4 * d6 + 2 * t7 * al4 * al6 - t3 * al4 * a6 + t1 * al6 * al4 * d5 * al5 - 2 * t6 * al4 - t2 * al6 * a4 * al5 + t3 * al6 * al4 * al5 * a5 - t1 * al6 * al4 * d4 * al5 - t3 * al4 * al5 * al6 * a6 - t1 * al4 * al5 * al6 * d6),(t1*al6*al4*a5+t1*al4*al5*a5+t3*al4*al5*d6-t1*al6*a4*al5-2*t4*al4*al6+t0*al4*a6-t0*al6*a4+t0*al4*a5-t0*a4*al5+t2*al4*d6-t2*al4*d4+t2*al4*d5-2*t4*al4*al5+t1*al4*al5*a6+t1*al4*al6*a6+t1*a4+2*t5*al4-t2*al6*al4*d5*al5+t2*al6*al4*d4*al5-t0*al6*al4*al5*a5-t3*al4*d5*al5+t3*al4*d4*al5-2*t5*al4*al5*al6-t3*al4*al6*d6-t0*al4*al5*al6*a6+t2*al4*al5*al6*d6+t3*al6*al4*d4-t3*al6*al4*d5)));
    result.push_back(Polynomial((-2*t7*al4-t2*al6*a4+t0*al4*d4-t0*al4*d5-t3*a4+t0*al6*al4*d4*al5-t3*al4*al5*a5+2*t6*al4*al5-t1*al6*al4*d5-t2*al4*al5*al6*a5+t3*al6*al4*a5-t3*al6*a4*al5-t0*al4*d6+t2*a4*al5-t2*al4*a5-t1*al4*d4*al5+t1*al4*d5*al5+t1*al6*al4*d4-t0*al6*al4*d5*al5-t1*al4*al5*d6-t1*al4*al6*d6+t3*al4*al5*a6-t3*al4*al6*a6-2*t7*al4*al5*al6-2*t6*al4*al6+t2*al4*a6+t0*al4*al5*al6*d6+t2*al4*al5*al6*a6),(-t0*al4*al5*a5+t2*al4*d5*al5-2*t4*al4+t3*al4*d5+t0*al6*a4*al5-t2*al4*al6*d4-t1*al4*al5*al6*a6+t3*al4*al5*al6*d6-t2*al4*al5*d6-t0*al4*al6*a6-t0*al4*al5*a6+t2*al4*al6*d6+2*t4*al4*al5*al6-2*t5*al4*al5-t1*al6*a4-t3*al6*al4*d5*al5-t1*a4*al5+t1*al4*a5-t0*al4*al6*a5+t2*al4*al6*d5-t3*al4*d4-t1*al4*al5*al6*a5-t2*al4*d4*al5+t3*al6*al4*d4*al5-t0*a4+t1*al4*a6+t3*al4*d6-2*t5*al4*al6)));
    result.push_back(Polynomial((t0*al4*al5*a5+t2*al4*d5*al5+2*t4*al4-t3*al4*d5+t0*al6*a4*al5+t2*al4*al6*d4-t1*al4*al5*al6*a6+t3*al4*al5*al6*d6-t2*al4*al5*d6+t0*al4*al6*a6-t0*al4*al5*a6-t2*al4*al6*d6+2*t4*al4*al5*al6-2*t5*al4*al5+t1*al6*a4-t3*al6*al4*d5*al5-t1*a4*al5+t1*al4*a5-t0*al4*al6*a5-t2*al4*al6*d5+t3*al4*d4+t1*al4*al5*al6*a5-t2*al4*d4*al5+t3*al6*al4*d4*al5+t0*a4-t1*al4*a6-t3*al4*d6+2*t5*al4*al6),(-2*t7*al4-t2*al6*a4+t0*al4*d4-t0*al4*d5-t3*a4-t0*al6*al4*d4*al5-t3*al4*al5*a5-2*t6*al4*al5-t1*al6*al4*d5-t2*al4*al5*al6*a5-t3*al6*al4*a5+t3*al6*a4*al5-t0*al4*d6-t2*a4*al5+t2*al4*a5+t1*al4*d4*al5-t1*al4*d5*al5+t1*al6*al4*d4+t0*al6*al4*d5*al5+t1*al4*al5*d6-t1*al4*al6*d6-t3*al4*al5*a6-t3*al4*al6*a6+2*t7*al4*al5*al6-2*t6*al4*al6+t2*al4*a6-t0*al4*al5*al6*d6-t2*al4*al5*al6*a6)));
    result.push_back(Polynomial((-t1*al6*al4*a5+t1*al4*al5*a5-t3*al4*al5*d6+t1*al6*a4*al5-2*t4*al4*al6+t0*al4*a6-t0*al6*a4-t0*al4*a5+t0*a4*al5+t2*al4*d6-t2*al4*d4+t2*al4*d5+2*t4*al4*al5-t1*al4*al5*a6+t1*al4*al6*a6+t1*a4+2*t5*al4+t2*al6*al4*d5*al5-t2*al6*al4*d4*al5-t0*al6*al4*al5*a5+t3*al4*d5*al5-t3*al4*d4*al5+2*t5*al4*al5*al6-t3*al4*al6*d6+t0*al4*al5*al6*a6-t2*al4*al5*al6*d6+t3*al6*al4*d4-t3*al6*al4*d5),(-t3*al6*a4+t2*a4+t0*al4*d5*al5-t0*al6*al4*d4+t0*al6*al4*d5+t2*al4*al5*a5-2*t7*al4*al5+t2*al6*al4*a5-t0*al4*d4*al5-t3*a4*al5+t3*al4*a5+t1*al4*d4-t1*al4*d5+t0*al4*al6*d6+t2*al4*al6*a6-t0*al4*al5*d6+t2*al4*al5*a6-2*t6*al4*al5*al6-t1*al4*d6-2*t7*al4*al6+t3*al4*a6+t1*al6*al4*d5*al5+2*t6*al4-t2*al6*a4*al5-t3*al6*al4*al5*a5-t1*al6*al4*d4*al5-t3*al4*al5*al6*a6-t1*al4*al5*al6*d6)));
    result.push_back(Polynomial((-2 * t3 * al4 * al5 + 2 * t3 * al4 * al6 - 2 * t2 * al4 - 2 * t2 * al4 * al5 * al6),(2*t1*al4-2*t1*al4*al5*al6-2*t0*al4*al5-2*t0*al4*al6)));
    result.push_back(Polynomial((-2 * t3 * al4 - 2 * t3 * al4 * al5 * al6 + 2 * t2 * al4 * al5 - 2 * t2 * al4 * al6),(-2*t1*al4*al5-2*t1*al4*al6-2*t0*al4+2*t0*al4*al5*al6)));
    result.push_back(Polynomial((-2 * t1 * al4 * al5 + 2 * t1 * al4 * al6 + 2 * t0 * al4 + 2 * t0 * al4 * al5 * al6),(-2*t3*al4+2*t3*al4*al5*al6-2*t2*al4*al5-2*t2*al4*al6)));
    result.push_back(Polynomial((2 * t1 * al4 + 2 * t1 * al4 * al5 * al6 + 2 * t0 * al4 * al5 - 2 * t0 * al4 * al6),(-2*t3*al4*al5-2*t3*al4*al6+2*t2*al4-2*t2*al4*al5*al6)));

    return result;
  };

  std::vector<Polynomial> h4_v6q(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double t4 = t.get(4,0);
    double t5 = t.get(5,0);
    double t6 = t.get(6,0);
    double t7 = t.get(7,0);
    double al4 = a.al[3], al5 = a.al[4], al6 = a.al[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial((2 * t7 - t2 * al6 * a4 * al4 + t0 * d5 - t2 * a6 + t2 * al6 * al5 * a5 + 2 * t7 * al5 * al6 - t1 * d5 * al5 - t1 * d4 * al5 + t1 * al6 * d4 + t1 * al6 * d5 - t3 * al6 * a5 + t1 * al6 * d6 + t1 * al5 * d6 - t0 * al5 * al6 * d6 + t3 * al5 * a5 - t3 * a4 * al4 + t0 * d4 + t0 * al6 * d5 * al5 - t3 * al5 * a6 + t2 * a5 + t0 * al6 * d4 * al5 + t3 * al6 * a6 - t2 * al6 * al5 * a6 - t3 * al6 * a4 * al4 * al5 + t2 * a4 * al4 * al5+t0*d6-2*t6*al5+2*t6*al6),(2*t4-t2*al6*d4-t2*al6*d5-t2*d4*al5-t1*a4*al4*al5-t0*a4*al4+t0*al5*a5+t0*al6*a5-t3*d5+2*t5*al5-t1*a5+t0*al6*a4*al4*al5+t3*al6*d4*al5+t3*al6*d5*al5-t3*d4-t3*al5*al6*d6+t1*al5*al6*a6-2*t4*al5*al6+t0*al6*a6+t0*al5*a6-t2*al6*d6+t2*al5*d6-t2*d5*al5-t1*a6-t3*d6+2*t5*al6+t1*al6*al5*a5-t1*al6*a4*al4)));
    result.push_back(Polynomial((-2*t7*al5+2*t7*al6-2*t6+t2*a4*al4+t1*d4+t1*d5-t0*al6*d5-t3*al6*a4*al4+t1*al6*d5*al5-t3*a6-t0*al6*d4+t3*a5+t2*al6*a5-t0*al6*d6-t0*al5*d6-2*t6*al5*al6+t3*a4*al4*al5+t3*al6*al5*a5+t0*d5*al5+t0*d4*al5-t2*al6*a6+t1*al6*d4*al5+t2*al5*a6-t3*al6*al5*a6-t1*al6*al5*d6-t2*al5*a5+t1*d6+t2*al6*a4*al4*al5),(-2*t5*al5*al6-t1*a4*al4+t3*al5*d6-t3*al6*d6+t1*al6*a6+t1*al5*a6-t3*d5*al5-t3*d4*al5-t3*al6*d5-t3*al6*d4+t1*al5*a5+t1*al6*a5-2*t4*al6-t2*al6*d5*al5+t2*d6-t2*al6*d4*al5+t0*a6+t0*al6*a4*al4-t0*al5*al6*a5+t0*a5-2*t4*al5-t0*al5*al6*a6+t0*a4*al4*al5+t2*d4+t2*al5*al6*d6+t2*d5+2*t5+t1*al6*a4*al4*al5)));
    result.push_back(Polynomial((2*t5*al5*al6-t1*a4*al4-t3*al5*d6-t3*al6*d6+t1*al6*a6-t1*al5*a6+t3*d5*al5+t3*d4*al5-t3*al6*d5-t3*al6*d4+t1*al5*a5-t1*al6*a5-2*t4*al6+t2*al6*d5*al5+t2*d6+t2*al6*d4*al5+t0*a6+t0*al6*a4*al4-t0*al5*al6*a5-t0*a5+2*t4*al5+t0*al5*al6*a6-t0*a4*al4*al5+t2*d4-t2*al5*al6*d6+t2*d5+2*t5-t1*al6*a4*al4*al5),(-2*t7*al5-2*t7*al6+2*t6-t2*a4*al4-t1*d4-t1*d5+t0*al6*d5+t3*al6*a4*al4+t1*al6*d5*al5+t3*a6+t0*al6*d4+t3*a5+t2*al6*a5+t0*al6*d6-t0*al5*d6-2*t6*al5*al6+t3*a4*al4*al5-t3*al6*al5*a5+t0*d5*al5+t0*d4*al5+t2*al6*a6+t1*al6*d4*al5+t2*al5*a6-t3*al6*al5*a6-t1*al6*al5*d6+t2*al5*a5-t1*d6+t2*al6*a4*al4*al5)));
    result.push_back(Polynomial((-t0*al6*a6+t3*d4+t0*a4*al4+t3*d5+2*t5*al5+t2*al6*d4+t2*al6*d5-2*t4+t3*al6*d4*al5-t2*d5*al5-2*t4*al5*al6+t2*al6*d6+t0*al6*a5+t0*al6*a4*al4*al5-t1*al6*al5*a5+t3*d6-t1*a5-t0*al5*a5+t0*al5*a6+t3*al6*d5*al5+t1*a6-t2*d4*al5-t3*al6*al5*d6-t1*a4*al4*al5-2*t5*al6+t2*al5*d6+t1*al6*al5*a6+t1*al6*a4*al4),(2*t7-t2*al6*a4*al4+t0*d5-t2*a6+t2*al6*al5*a5-2*t7*al5*al6+t1*d5*al5+t1*d4*al5+t1*al6*d4+t1*al6*d5+t3*al6*a5+t1*al6*d6-t1*al5*d6+t0*al5*al6*d6+t3*al5*a5-t3*a4*al4+t0*d4-t0*al6*d5*al5+t3*al5*a6-t2*a5-t0*al6*d4*al5+t3*al6*a6+t2*al6*al5*a6+t3*al6*a4*al4*al5-t2*a4*al4*al5+t0*d6+2*t6*al5+2*t6*al6)));
    result.push_back(Polynomial((-2 * t2 * al5 + 2 * t2 * al6 + 2 * t3 + 2 * t3 * al6 * al5),(2*t0-2*t0*al6*al5+2*t1*al5+2*t1*al6)));
    result.push_back(Polynomial((-2 * t2 - 2 * t2 * al5 * al6 - 2 * t3 * al5 + 2 * t3 * al6),(-2*t0*al5-2*t0*al6+2*t1-2*t1*al6*al5)));
    result.push_back(Polynomial((2 * t0 * al5 - 2 * t0 * al6 + 2 * t1 + 2 * t1 * al6 * al5),(2*t2-2*t2*al5*al6-2*t3*al5-2*t3*al6)));
    result.push_back(Polynomial((-2 * t0 - 2 * t0 * al6 * al5 + 2 * t1 * al5 - 2 * t1 * al6),(2*t2*al5+2*t2*al6+2*t3-2*t3*al6*al5)));

    return result;
  };
}