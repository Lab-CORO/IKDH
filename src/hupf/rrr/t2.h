#pragma once
#include <hupf/Input.h>

namespace LibHUPF
{
  std::vector<Polynomial> h1_twq(Input& a)
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
    double al6 = a.al[5];
    double d4 = a.d[3], d6 = a.d[5];
    double a6 = a.a[5];

    result.push_back(Polynomial((t0 * al6 * a6 - t1 * a6 - t2 * (al6 * d4 + al6 * d6) - t3 * (d4 + d6) + 2 * t4 + 2 * t5 * al6)));
    result.push_back(Polynomial((t1 * al6 * a6 + t0 * a6 - t3 * (al6 * d4 + al6 * d6) + t2 * (d4 + d6) + 2 * t5 - 2 * t4 * al6)));
    result.push_back(Polynomial((t2 * al6 * a6 + t3 * a6 + t0 * (al6 * d4 + al6 * d6) - t1 * (d4 + d6) + 2 * t6 - 2 * t7 * al6)));
    result.push_back(Polynomial((t3 * al6 * a6 - t2 * a6 + t1 * (al6 * d4 + al6 * d6) + t0 * (d4 + d6) + 2 * t7 + 2 * t6 * al6)));
    result.push_back(Polynomial((2 * t0 + 2 * t1 * al6)));
    result.push_back(Polynomial((2 * t1 - 2 * t0 * al6)));
    result.push_back(Polynomial((2 * t2 - 2 * t3 * al6)));
    result.push_back(Polynomial((2 * t3 + 2 * t2 * al6)));

    return result;
  };

  std::vector<Polynomial> h2_twq(Input& a)
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
    double al6 = a.al[5];
    double d4 = a.d[3], d6 = a.d[5];
    double a6 = a.a[5];

    result.push_back(Polynomial((-t0 * a6 - t1 * al6 * a6 - t2 * (-d4 + d6) - t3 * (al6 * d4 - al6 * d6) + 2 * t4 * al6 - 2 * t5)));
    result.push_back(Polynomial((-t1 * a6 + t0 * al6 * a6 - t3 * (-d4 + d6) + t2 * (al6 * d4 - al6 * d6) + 2 * t5 * al6 + 2 * t4)));
    result.push_back(Polynomial((-t2 * a6 + t3 * al6 * a6 + t0 * (-d4 + d6) - t1 * (al6 * d4 - al6 * d6) + 2 * t6 * al6 + 2 * t7)));
    result.push_back(Polynomial((-t3 * a6 - t2 * al6 * a6 + t1 * (-d4 + d6) + t0 * (al6 * d4 - al6 * d6) + 2 * t7 * al6 - 2 * t6)));
    result.push_back(Polynomial((2 * t0 * al6 - 2 * t1)));
    result.push_back(Polynomial((2 * t0 + 2 * t1 * al6)));
    result.push_back(Polynomial((2 * t3 + 2 * t2 * al6)));
    result.push_back(Polynomial((2 * t3 * al6 - 2 * t2)));

    return result;
  };

  std::vector<Polynomial> h3_twq(Input& a)
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
    double al6 = a.al[5];
    double d4 = a.d[3], d6 = a.d[5];
    double a6 = a.a[5];

    result.push_back(Polynomial((t0 * (al6 * d4 - al6 * d6) - t1 * (d4 - d6) - t2 * al6 * a6 - t3 * a6 - 2 * t6 + 2 * t7 * al6)));
    result.push_back(Polynomial((t1 * (al6 * d4 - al6 * d6) + t0 * (d4 - d6) - t3 * al6 * a6 + t2 * a6 - 2 * t7 - 2 * t6 * al6)));
    result.push_back(Polynomial((t2 * (al6 * d4 - al6 * d6) + t3 * (d4 - d6) + t0 * al6 * a6 - t1 * a6 + 2 * t4 + 2 * t5 * al6)));
    result.push_back(Polynomial((t3 * (al6 * d4 - al6 * d6) - t2 * (d4 - d6) + t1 * al6 * a6 + t0 * a6 + 2 * t5 - 2 * t4 * al6)));
    result.push_back(Polynomial((2 * t3 * al6 - 2 * t2)));
    result.push_back(Polynomial((-2 * t3 - 2 * t2 * al6)));
    result.push_back(Polynomial((2 * t0 + 2 * t1 * al6)));
    result.push_back(Polynomial((2 * t1 - 2 * t0 * al6)));

    return result;
  };

  std::vector<Polynomial> h4_twq(Input& a)
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
    double al6 = a.al[5];
    double d4 = a.d[3], d6 = a.d[5];
    double a6 = a.a[5];

    result.push_back(Polynomial((t0 * (-d4 - d6) - t1 * (al6 * d4 + al6 * d6) + t2 * a6 - t3 * al6 * a6 - 2 * t6 * al6 - 2 * t7)));
    result.push_back(Polynomial((t1 * (-d4 - d6) + t0 * (al6 * d4 + al6 * d6) + t3 * a6 + t2 * al6 * a6 - 2 * t7 * al6 + 2 * t6)));
    result.push_back(Polynomial((t2 * (-d4 - d6) + t3 * (al6 * d4 + al6 * d6) - t0 * a6 - t1 * al6 * a6 + 2 * t4 * al6 - 2 * t5)));
    result.push_back(Polynomial((t3 * (-d4 - d6) - t2 * (al6 * d4 + al6 * d6) - t1 * a6 + t0 * al6 * a6 + 2 * t5 * al6 + 2 * t4)));
    result.push_back(Polynomial((-2 * t3 - 2 * t2 * al6)));
    result.push_back(Polynomial((2 * t2 - 2 * t3 * al6)));
    result.push_back(Polynomial((2 * t0 * al6 - 2 * t1)));
    result.push_back(Polynomial((2 * t0 + 2 * t1 * al6)));

    return result;
  };

  std::vector<Polynomial> h1_tpq(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double al6 = a.al[5];

    result.push_back(Polynomial((2 * t0 * al6 - 2 * t1)));
    result.push_back(Polynomial((2 * t1 * al6 + 2 * t0)));
    result.push_back(Polynomial((2 * t2 * al6 + 2 * t3)));
    result.push_back(Polynomial((2 * t3 * al6 - 2 * t2)));

    return result;
  };

  std::vector<Polynomial> h2_tpq(Input& a)
  {
    std::vector<Polynomial> result;
    Matrix t = a.studyParams();
    double t0 = t.get(0,0);
    double t1 = t.get(1,0);
    double t2 = t.get(2,0);
    double t3 = t.get(3,0);
    double al6 = a.al[5];

    result.push_back(Polynomial((2 * t3 * al6 - 2 * t2)));
    result.push_back(Polynomial((-2 * t3 - 2 * t2 * al6)));
    result.push_back(Polynomial((2 * t1 * al6 + 2 * t0)));
    result.push_back(Polynomial((2 * t1 - 2 * t0 * al6)));

    return result;
  };

  std::vector<Polynomial> h3_tpq(Input& a)
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
    double al6 = a.al[5];
    double d4 = a.d[3], d6 = a.d[5];
    double a6 = a.a[5];

    result.push_back(Polynomial((t0 * al6 * a6 - t1 * a6 - t2 * (al6 * d4 + al6 * d6) - t3 * (d4 + d6) + 2 * t4 + 2 * t5 * al6)));
    result.push_back(Polynomial((t1 * al6 * a6 + t0 * a6 - t3 * (al6 * d4 + al6 * d6) + t2 * (d4 + d6) + 2 * t5 - 2 * t4 * al6)));
    result.push_back(Polynomial((t2 * al6 * a6 + t3 * a6 + t0 * (al6 * d4 + al6 * d6) - t1 * (d4 + d6) + 2 * t6 - 2 * t7 * al6)));
    result.push_back(Polynomial((t3 * al6 * a6 - t2 * a6 + t1 * (al6 * d4 + al6 * d6) + t0 * (d4 + d6) + 2 * t7 + 2 * t6 * al6)));
    result.push_back(Polynomial((2 * t1 * al6 + 2 * t0)));
    result.push_back(Polynomial((2 * t1 - 2 * t0 * al6)));
    result.push_back(Polynomial((2 * t2 - 2 * t3 * al6)));
    result.push_back(Polynomial((2 * t2 * al6 + 2 * t3)));

    return result;
  };

  std::vector<Polynomial> h4_tpq(Input& a)
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
    double al6 = a.al[5];
    double d4 = a.d[3], d6 = a.d[5];
    double a6 = a.a[5];

    result.push_back(Polynomial((t0 * (-d4 - d6) - t1 * (al6 * d4 + al6 * d6) + t2 * a6 - t3 * al6 * a6 - 2 * t6 * al6 - 2 * t7)));
    result.push_back(Polynomial((t1 * (-d4 - d6) + t0 * (al6 * d4 + al6 * d6) + t3 * a6 + t2 * al6 * a6 - 2 * t7 * al6 + 2 * t6)));
    result.push_back(Polynomial((t2 * (-d4 - d6) + t3 * (al6 * d4 + al6 * d6) - t0 * a6 - t1 * al6 * a6 + 2 * t4 * al6 - 2 * t5)));
    result.push_back(Polynomial((t3 * (-d4 - d6) - t2 * (al6 * d4 + al6 * d6) - t1 * a6 + t0 * al6 * a6 + 2 * t5 * al6 + 2 * t4)));
    result.push_back(Polynomial((-2 * t3 - 2 * t2 * al6)));
    result.push_back(Polynomial((2 * t2 - 2 * t3 * al6)));
    result.push_back(Polynomial((2 * t0 * al6 - 2 * t1)));
    result.push_back(Polynomial((2 * t1 * al6 + 2 * t0)));

    return result;
  };  
}