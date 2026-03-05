#pragma once
#include <hupf/Input.h>

namespace LibHUPF
{
  std::vector<Polynomial> h1_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(((8 - 8 * al1 * al3) * a2 * al2 - 8 * a1 * al1 - 8 * al1 * a3 - 8 * a1 * al3 - 8 * al3 * a3),(-8*al1*al3*d2-8*d2-8*al1*al3*d3-8*d3)));
    result.push_back(Polynomial(((8 * al3 + 8 * al1) * a2 * al2 - 8 * a1 * al1 * al3 - 8 * al1 * al3 * a3 + 8 * a1 + 8 * a3),(8*al1*d2-8*al3*d2+8*al1*d3-8*al3*d3)));
    result.push_back(Polynomial((8 * al3 * d3 - 8 * al1 * d2 + 8 * al3 * d2 - 8 * al1 * d3),((8*al3+8*al1)*a2*al2-8*a1*al1*al3-8*al1*al3*a3+8*a1+8*a3)));
    result.push_back(Polynomial((8 * d2 + 8 * al1 * al3 * d2 + 8 * al1 * al3 * d3 + 8 * d3),((8-8*al1*al3)*a2*al2-8*a1*al1-8*al1*a3-8*a1*al3-8*al3*a3)));
    result.push_back(Polynomial((-16 + 16 * al3 * al1)));
    result.push_back(Polynomial((-16 * al3 - 16 * al1)));
    result.push_back(Polynomial(0.0, (-16 * al1 - 16 * al3)));
    result.push_back(Polynomial(0.0, (-16 + 16 * al3 * al1)));

    return result;
  };

  std::vector<Polynomial> h2_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(((-8 * al1 - 8 * al3) * a2 + (8 * a1 + 8 * a3 - 8 * a1 * al1 * al3 - 8 * al1 * al3 * a3) * al2),(-8*al3*d2+8*al1*d2+8*al3*d3-8*al1*d3)*al2));
    result.push_back(Polynomial(((8 - 8 * al1 * al3) * a2 + (8 * a1 * al3 + 8 * al3 * a3 + 8 * a1 * al1 + 8 * al1 * a3) * al2),(8*d2+8*al1*al3*d2-8*d3-8*al1*al3*d3)*al2));
    result.push_back(Polynomial(((8 * al1 * al3 * d3 - 8 * d2 - 8 * al1 * al3 * d2 + 8 * d3) * al2),((8-8*al1*al3)*a2+(8*a1*al3+8*al3*a3+8*a1*al1+8*al1*a3)*al2)));
    result.push_back(Polynomial(((8 * al3 * d2 - 8 * al3 * d3 - 8 * al1 * d2 + 8 * al1 * d3) * al2),((-8*al1-8*al3)*a2+(8*a1+8*a3-8*a1*al1*al3-8*al1*al3*a3)*al2)));
    result.push_back(Polynomial((2 * (-8 * al1 - 8 * al3) * al2)));
    result.push_back(Polynomial((2 * (8 - 8 * al3 * al1) * al2)));
    result.push_back(Polynomial(0.0, (2 * (8 - 8 * al3 * al1) * al2)));
    result.push_back(Polynomial(0.0, (2 * (-8 * al3 - 8 * al1) * al2)));

    return result;
  };

  std::vector<Polynomial> h3_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(((-8 * al1 * d2 - 8 * al3 * d2 + 8 * al1 * d3 + 8 * al3 * d3) * al2),((-8*al1+8*al3)*a2+(8*a1-8*al1*al3*a3-8*a3+8*a1*al1*al3)*al2)));
    result.push_back(Polynomial(((-8 * al1 * al3 * d2 + 8 * d2 + 8 * al1 * al3 * d3 - 8 * d3) * al2),((-8-8*al1*al3)*a2+(8*a1*al3+8*al1*a3-8*al3*a3-8*a1*al1)*al2)));
    result.push_back(Polynomial(((8 + 8 * al1 * al3) * a2 + (8 * a1 * al1 - 8 * a1 * al3 - 8 * al1 * a3 + 8 * al3 * a3) * al2),(-8*al1*al3*d2+8*d2+8*al1*al3*d3-8*d3)*al2));
    result.push_back(Polynomial(((-8 * al3 + 8 * al1) * a2 + (-8 * a1 * al1 * al3 - 8 * a1 + 8 * al1 * al3 * a3 + 8 * a3) * al2),(-8*al1*d2-8*al3*d2+8*al1*d3+8*al3*d3)*al2));
    result.push_back(Polynomial(0.0, (2 * (-8 * al1 + 8 * al3) * al2) ));
    result.push_back(Polynomial(0.0, (2 * (-8 - 8 * al3 * al1) * al2) ));
    result.push_back(Polynomial((2 * (8 + 8 * al3 * al1) * al2) ));
    result.push_back(Polynomial((2 * (-8 * al3 + 8 * al1) * al2) ));

    return result;
  };

  std::vector<Polynomial> h4_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-8 * d3 - 8 * d2 + 8 * al1 * al3 * d2 + 8 * al1 * al3 * d3),(8*a1*al1+(-8-8*al1*al3)*a2*al2+8*al3*a3-8*al1*a3-8*a1*al3)));
    result.push_back(Polynomial((-8 * al1 * d2 - 8 * al3 * d2 - 8 * al1 * d3 - 8 * al3 * d3),(8*a1*al1*al3-8*a3-8*al1*al3*a3+8*a1+(-8*al3+8*al1)*a2*al2)));
    result.push_back(Polynomial(((-8 * al1 + 8 * al3) * a2 * al2 - 8 * a1 - 8 * a1 * al1 * al3 + 8 * a3 + 8 * al1 * al3 * a3),(-8*al1*d2-8*al3*d2-8*al1*d3-8*al3*d3)));
    result.push_back(Polynomial(((8 + 8 * al1 * al3) * a2 * al2 + 8 * a1 * al3 - 8 * a1 * al1 - 8 * al3 * a3 + 8 * al1 * a3),(-8*d3-8*d2+8*al1*al3*d2+8*al1*al3*d3)));
    result.push_back(Polynomial(0.0, (16 + 16 * al3 * al1) ));
    result.push_back(Polynomial(0.0, (16 * al3 - 16 * al1) ));
    result.push_back(Polynomial((16 * al1 - 16 * al3) ));
    result.push_back(Polynomial((-16 - 16 * al3 * al1) ));

    return result;
  };

  std::vector<Polynomial> h1_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-2 * al3 * al1 * d2 + 2 * al3 * al2 * d2 - 2 * al3 * d3 * al2 - 2 * al3 * d3 * al1),(2 * al3 * a2 + 2 * al3 * al1 * al2 * a2 - 2 * al3 * al1 * a1 * al2 - 2 * al3 * a1  + 2 * a3 * al1 - 2 * al2 * a3 )));
    result.push_back(Polynomial((+2 * al1 * d2 - 2 * al2 * d2 + 2 * d3 * al2 + 2 * d3 * al1),(-2 * a2 - 2 * al1 * al2 * a2 + 2 * al1 * a1 * al2 + 2 * a1  + 2 * al3 * a3 * al1 - 2 * a3 * al3 * al2 )));
    result.push_back(Polynomial((-2 * a1 - 2 * a2 + 2 * a2 * al1 * al2 + 2 * a1 * al1 * al2 - 2 * a3 * al3 * al2 - 2 * a3 * al3 * al1),(2 * d2 * al1 + 2 * d2 * al2  + 2 * d3 * al1 - 2 * d3 * al2 )));
    result.push_back(Polynomial((+2 * al3 * a1 + 2 * al3 * a2 - 2 * al3 * al1 * al2 * a2 - 2 * al3 * al1 * al2 * a1 - 2 * al2 * a3 - 2 * al1 * a3),(-2 * al1 * d2 * al3 - 2 * d2 * al3 * al2  - 2 * d3 * al1 * al3 + 2 * al3 * d3 * al2 )));
    result.push_back(Polynomial(0.0, (-4 * al3 * (al1 - al2)) ));
    result.push_back(Polynomial(0.0, (4 * al1 - 4 * al2) ));
    result.push_back(Polynomial((-4 * al2 - 4 * al1) ));
    result.push_back(Polynomial((4 * al3 * (al1 + al2)) ));

    return result;
  };

  std::vector<Polynomial> h2_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-2 * d2 - 2 * al1 * al2 * d2 + 2 * d3 * al1 * al2 - 2 * d3),(2 * al2 * a2 + 2 * al1 * a1 - 2 * al1 * a2 - 2 * al2 * a1  - 2 * al3 * a3 * al1 * al2 - 2 * al3 * a3 )));
    result.push_back(Polynomial((-2 * al3 * d2 - 2 * al3 * al1 * al2 * d2 + 2 * al3 * d3 * al1 * al2 - 2 * al3 * d3),(2 * a2 * al3 * al2 + 2 * al3 * al1 * a1 - 2 * al3 * al1 * a2 - 2 * a1 * al3 * al2  + 2 * a3 * al1 * al2 + 2 * a3 )));
    result.push_back(Polynomial((+2 * a2 * al3 * al1 + 2 * a1 * al3 * al1 + 2 * a1 * al3 * al2 + 2 * a2 * al3 * al2 - 2 * a3 * al1 * al2 + 2 * a3),(2 * d2 * al3 - 2 * al1 * d2 * al3 * al2  + 2 * al3 * d3 * al1 * al2 + 2 * d3 * al3 )));
    result.push_back(Polynomial((+2 * al1 * a2 + 2 * al1 * a1 + 2 * al2 * a1 + 2 * al2 * a2 + 2 * al3 * al1 * al2 * a3 - 2 * al3 * a3),(2 * d2 - 2 * al1 * d2 * al2  + 2 * d3 * al1 * al2 + 2 * d3 )));
    result.push_back(Polynomial(0.0, (-4 * al1 * al2 - 4) ));
    result.push_back(Polynomial(0.0, (-4 * al3 * (1 + al1 * al2)) ));
    result.push_back(Polynomial((4 * al3 * (-1 + al1 * al2)) ));
    result.push_back(Polynomial((4 * al1 * al2 - 4) ));

    return result;
  };

  std::vector<Polynomial> h3_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((+2 * al3 * a1 + 2 * al3 * a2 - 2 * al3 * al1 * al2 * a2 - 2 * al3 * al1 * al2 * a1 - 2 * al2 * a3 - 2 * al1 * a3),(-2 * al1 * d2 * al3 - 2 * d2 * al3 * al2  - 2 * d3 * al1 * al3 + 2 * al3 * d3 * al2 )));
    result.push_back(Polynomial((-2 * a1 - 2 * a2 + 2 * a2 * al1 * al2 + 2 * a1 * al1 * al2 - 2 * a3 * al3 * al2 - 2 * a3 * al3 * al1),(2 * d2 * al1 + 2 * d2 * al2  + 2 * d3 * al1 - 2 * d3 * al2 )));
    result.push_back(Polynomial((+2 * al2 * d2 - 2 * al1 * d2 - 2 * d3 * al2 - 2 * d3 * al1),(2 * a2 - 2 * a1 + 2 * al1 * al2 * a2 - 2 * al1 * a1 * al2  + 2 * a3 * al3 * al2 - 2 * al3 * a3 * al1)));
    result.push_back(Polynomial((- 2 * al3 * al2 * d2 + 2 * al3 * al1 * d2 + 2 * al3 * d3 * al2 + 2 * al3 * d3 * al1),(-2 * al3 * a2 + 2 * al3 * a1 - 2 * al3 * al1 * al2 * a2 + 2 * al3 * al1 * a1 * al2  + 2 * al2 * a3 - 2 * a3 * al1)));
    result.push_back(Polynomial((4 * al3 * (al1 + al2)) ));
    result.push_back(Polynomial((-4 * al2 - 4 * al1) ));
    result.push_back(Polynomial(0.0, (4 * al2 - 4 * al1) ));
    result.push_back(Polynomial(0.0, (4 * al3 * (al1 - al2)) ));

    return result;
  };

  std::vector<Polynomial> h4_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-2 * al2 * a1 - 2 * al1 * a1 - 2 * al2 * a2 - 2 * al1 * a2 + 2 * al3 * a3 - 2 * al3 * al1 * al2 * a3),(-2 * d2 + 2 * al1 * d2 * al2  - 2 * d3 * al1 * al2 - 2 * d3)));
    result.push_back(Polynomial((-2 * a1 * al3 * al2 - 2 * a1 * al3 * al1 - 2 * a2 * al3 * al2 - 2 * a2 * al3 * al1 - 2 * a3 + 2 * a3 * al1 * al2),(-2 * d2 * al3 + 2 * al1 * d2 * al3 * al2  - 2 * al3 * d3 * al1 * al2 - 2 * d3 * al3)));
    result.push_back(Polynomial((-2 * al3 * d2 - 2 * al3 * al1 * al2 * d2 + 2 * al3 * d3 * al1 * al2 - 2 * al3 * d3),(2 * a2 * al3 * al2 + 2 * al3 * al1 * a1 - 2 * al3 * al1 * a2 - 2 * a1 * al3 * al2  + 2 * a3 * al1 * al2 + 2 * a3 )));
    result.push_back(Polynomial((-2 * d2 - 2 * al1 * al2 * d2 + 2 * d3 * al1 * al2 - 2 * d3),(2 * al2 * a2 + 2 * al1 * a1 - 2 * al1 * a2 - 2 * al2 * a1  - 2 * al3 * a3 * al1 * al2 - 2 * al3 * a3 )));
    result.push_back(Polynomial((4 - 4 * al1 * al2) ));
    result.push_back(Polynomial(-(4 * al3 * (-1 + al1 * al2)) ));
    result.push_back(Polynomial(0.0, -(4 * al3 * (1 + al1 * al2)) ));
    result.push_back(Polynomial(0.0, (-4 * al1 * al2 - 4) ));

    return result;
  };

  std::vector<Polynomial> h1_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-al3 * a2 - al3 * a3 - al2 * a3 + a1 * al1 - al2 * a2 - a1 * al1 * al2 * al3),(-d2-al3*al2*d2+al3*al2*d3-d3)));
    result.push_back(Polynomial((a1 * al1 * al3 - al3 * al2 * a2 + a3 + a2 + a1 * al1 * al2 - al3 * al2 * a3),(-al3*d3-al2*d3-al3*d2+al2*d2)));
    result.push_back(Polynomial((-al2 * d3 + al3 * d3 + al2 * d2 + al3 * d2),(a3-a2-al3*al2*a2+al3*al2*a3-a1*al1*al2+a1*al1*al3)));
    result.push_back(Polynomial((d3 + al3 * al2 * d3 - al3 * al2 * d2 + d2),(al3*a2+a1*al1*al2*al3+al2*a3-al2*a2+a1*al1-al3*a3)));
    result.push_back(Polynomial((-2 + 2 * al2 * al3)));
    result.push_back(Polynomial((-2 * al2 - 2 * al3)));
    result.push_back(Polynomial(0.0, (2 * al2 - 2 * al3)));
    result.push_back(Polynomial(0.0, (-2 * al2 * al3 - 2)));

    return result;
  };

  std::vector<Polynomial> h2_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-al1 * al3 * al2 * a3 - a1 * al2 - al1 * al3 * al2 * a2 - a1 * al3 + al1 * a2 + al1 * a3),(-al1*al3*d2-al1*al2*d3+al1*al2*d2-al1*al3*d3)));
    result.push_back(Polynomial((a1 + al1 * al2 * a3 + al1 * al3 * a2 - a1 * al2 * al3 + al1 * al3 * a3 + al1 * al2 * a2),(al1*al3*al2*d2+al1*d2-al1*al3*al2*d3+al1*d3)));
    result.push_back(Polynomial((al1 * d2 - al1 * al3 * al2 * d2 + al1 * d3 + al1 * al3 * al2 * d3),(-a1*al2*al3-al1*al3*a3-a1+al1*al3*a2+al1*al2*a3-al1*al2*a2)));
    result.push_back(Polynomial((-al1 * al2 * d2 - al1 * al3 * d3 - al1 * al3 * d2 + al1 * al2 * d3),(-al1*a3+a1*al3-a1*al2+al1*a2+al1*al3*al2*a2-al1*al3*al2*a3)));
    result.push_back(Polynomial((-2 * al1 * al3 - 2 * al1 * al2) ));
    result.push_back(Polynomial((2 * al1 - 2 * al1 * al2 * al3) ));
    result.push_back(Polynomial(0.0, (-2 * al1 - 2 * al1 * al2 * al3) ));
    result.push_back(Polynomial(0.0, (-2 * al1 * al2 + 2 * al1 * al3) ));

    return result;
  };

  std::vector<Polynomial> h3_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-al1 * al2 * d3 + al1 * al3 * d3 + al1 * al2 * d2 + al1 * al3 * d2),(al1*a3+a1*al2-a1*al3-al1*al3*al2*a2-al1*a2+al1*al3*al2*a3)));
    result.push_back(Polynomial((-al1 * d2 + al1 * al3 * al2 * d2 - al1 * d3 - al1 * al3 * al2 * d3),(a1*al2*al3+al1*al2*a2+al1*al3*a3+a1-al1*al2*a3-al1*al3*a2)));
    result.push_back(Polynomial((a1 + al1 * al2 * a3 + al1 * al3 * a2 - a1 * al2 * al3 + al1 * al3 * a3 + al1 * al2 * a2),(al1*al3*al2*d2+al1*d2-al1*al3*al2*d3+al1*d3)));
    result.push_back(Polynomial((-al1 * al3 * al2 * a3 - a1 * al2 - al1 * al3 * al2 * a2 - a1 * al3 + al1 * a2 + al1 * a3),(-al1*al3*d2-al1*al2*d3+al1*al2*d2-al1*al3*d3)));
    result.push_back(Polynomial(0.0, (2 * al1 * al2 - 2 * al1 * al3) ));
    result.push_back(Polynomial(0.0, (2 * al1 * al2 * al3 + 2 * al1) ));
    result.push_back(Polynomial((2 * al1 - 2 * al1 * al2 * al3) ));
    result.push_back(Polynomial((-2 * al1 * al3 - 2 * al1 * al2) ));

    return result;
  };

  std::vector<Polynomial> h4_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double al1 = a.al[0], al2 = a.al[1], al3 = a.al[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial((-d2 + al3 * al2 * d2 - al3 * al2 * d3 - d3),(al2*a2-al3*a2-a1*al1+al3*a3-al2*a3-a1*al1*al2*al3)));
    result.push_back(Polynomial((-al3 * d2 - al2 * d2 - al3 * d3 + al2 * d3),(al3*al2*a2-a3-a1*al1*al3-al3*al2*a3+a2+a1*al1*al2)));
    result.push_back(Polynomial((a1 * al1 * al3 - al3 * al2 * a2 + a3 + a2 + a1 * al1 * al2 - al3 * al2 * a3),(-al3*d3-al2*d3-al3*d2+al2*d2)));
    result.push_back(Polynomial((-al3 * a2 - al3 * a3 - al2 * a3 + a1 * al1 - al2 * a2 - a1 * al1 * al2 * al3),(-d2-al3*al2*d2+al3*al2*d3-d3)));
    result.push_back(Polynomial(0.0, (2 * al2 * al3 + 2) ));
    result.push_back(Polynomial(0.0, (2 * al3 - 2 * al2) ));
    result.push_back(Polynomial((-2 * al2 - 2 * al3) ));
    result.push_back(Polynomial((2 * al3 * al2 - 2) ));

    return result;
  };

}
