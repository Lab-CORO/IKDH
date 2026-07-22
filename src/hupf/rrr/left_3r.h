#pragma once
#include <hupf/proj_eval.h>
#include <hupf/Input.h>

namespace LibHUPF
{
  std::vector<Polynomial> h1_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 - 8 * al1 * al3) * a2 * al2 - 8 * a1 * al1 - 8 * al1 * a3 - 8 * a1 * al3 - 8 * al3 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-8*al1*al3*d2-8*d2-8*al1*al3*d3-8*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 * al3 + 8 * al1) * a2 * al2 - 8 * a1 * al1 * al3 - 8 * al1 * al3 * a3 + 8 * a1 + 8 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (8*al1*d2-8*al3*d2+8*al1*d3-8*al3*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (8 * al3 * d3 - 8 * al1 * d2 + 8 * al3 * d2 - 8 * al1 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return ((8*al3+8*al1)*a2*al2-8*a1*al1*al3-8*al1*al3*a3+8*a1+8*a3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (8 * d2 + 8 * al1 * al3 * d2 + 8 * al1 * al3 * d3 + 8 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return ((8-8*al1*al3)*a2*al2-8*a1*al1-8*al1*a3-8*a1*al3-8*al3*a3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-16 + 16 * al3 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-16 * al3 - 16 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-16 * al1 - 16 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-16 + 16 * al3 * al1);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h2_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((-8 * al1 - 8 * al3) * a2 + (8 * a1 + 8 * a3 - 8 * a1 * al1 * al3 - 8 * al1 * al3 * a3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-8*al3*d2+8*al1*d2+8*al3*d3-8*al1*d3)*al2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 - 8 * al1 * al3) * a2 + (8 * a1 * al3 + 8 * al3 * a3 + 8 * a1 * al1 + 8 * al1 * a3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (8*d2+8*al1*al3*d2-8*d3-8*al1*al3*d3)*al2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 * al1 * al3 * d3 - 8 * d2 - 8 * al1 * al3 * d2 + 8 * d3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return ((8-8*al1*al3)*a2+(8*a1*al3+8*al3*a3+8*a1*al1+8*al1*a3)*al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 * al3 * d2 - 8 * al3 * d3 - 8 * al1 * d2 + 8 * al1 * d3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return ((-8*al1-8*al3)*a2+(8*a1+8*a3-8*a1*al1*al3-8*al1*al3*a3)*al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * (-8 * al1 - 8 * al3) * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * (8 - 8 * al3 * al1) * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * (8 - 8 * al3 * al1) * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * (-8 * al3 - 8 * al1) * al2);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h3_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((-8 * al1 * d2 - 8 * al3 * d2 + 8 * al1 * d3 + 8 * al3 * d3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return ((-8*al1+8*al3)*a2+(8*a1-8*al1*al3*a3-8*a3+8*a1*al1*al3)*al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((-8 * al1 * al3 * d2 + 8 * d2 + 8 * al1 * al3 * d3 - 8 * d3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return ((-8-8*al1*al3)*a2+(8*a1*al3+8*al1*a3-8*al3*a3-8*a1*al1)*al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 + 8 * al1 * al3) * a2 + (8 * a1 * al1 - 8 * a1 * al3 - 8 * al1 * a3 + 8 * al3 * a3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-8*al1*al3*d2+8*d2+8*al1*al3*d3-8*d3)*al2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((-8 * al3 + 8 * al1) * a2 + (-8 * a1 * al1 * al3 - 8 * a1 + 8 * al1 * al3 * a3 + 8 * a3) * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-8*al1*d2-8*al3*d2+8*al1*d3+8*al3*d3)*al2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * (-8 * al1 + 8 * al3) * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * (-8 - 8 * al3 * al1) * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * (8 + 8 * al3 * al1) * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * (-8 * al3 + 8 * al1) * al2);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h4_tc_v1(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-8 * d3 - 8 * d2 + 8 * al1 * al3 * d2 + 8 * al1 * al3 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (8*a1*al1+(-8-8*al1*al3)*a2*al2+8*al3*a3-8*al1*a3-8*a1*al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-8 * al1 * d2 - 8 * al3 * d2 - 8 * al1 * d3 - 8 * al3 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (8*a1*al1*al3-8*a3-8*al1*al3*a3+8*a1+(-8*al3+8*al1)*a2*al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((-8 * al1 + 8 * al3) * a2 * al2 - 8 * a1 - 8 * a1 * al1 * al3 + 8 * a3 + 8 * al1 * al3 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-8*al1*d2-8*al3*d2-8*al1*d3-8*al3*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return ((8 + 8 * al1 * al3) * a2 * al2 + 8 * a1 * al3 - 8 * a1 * al1 - 8 * al3 * a3 + 8 * al1 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-8*d3-8*d2+8*al1*al3*d2+8*al1*al3*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (16 + 16 * al3 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (16 * al3 - 16 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (16 * al1 - 16 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-16 - 16 * al3 * al1);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h1_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al3 * al1 * d2 + 2 * al3 * al2 * d2 - 2 * al3 * d3 * al2 - 2 * al3 * d3 * al1);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al3 * a2 + 2 * al3 * al1 * al2 * a2 - 2 * al3 * al1 * a1 * al2 - 2 * al3 * a1  + 2 * a3 * al1 - 2 * al2 * a3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (+2 * al1 * d2 - 2 * al2 * d2 + 2 * d3 * al2 + 2 * d3 * al1);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * a2 - 2 * al1 * al2 * a2 + 2 * al1 * a1 * al2 + 2 * a1  + 2 * al3 * a3 * al1 - 2 * a3 * al3 * al2 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * a1 - 2 * a2 + 2 * a2 * al1 * al2 + 2 * a1 * al1 * al2 - 2 * a3 * al3 * al2 - 2 * a3 * al3 * al1);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * d2 * al1 + 2 * d2 * al2  + 2 * d3 * al1 - 2 * d3 * al2 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (+2 * al3 * a1 + 2 * al3 * a2 - 2 * al3 * al1 * al2 * a2 - 2 * al3 * al1 * al2 * a1 - 2 * al2 * a3 - 2 * al1 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * al1 * d2 * al3 - 2 * d2 * al3 * al2  - 2 * d3 * al1 * al3 + 2 * al3 * d3 * al2 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-4 * al3 * (al1 - al2));},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (4 * al1 - 4 * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-4 * al2 - 4 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (4 * al3 * (al1 + al2));},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h2_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * d2 - 2 * al1 * al2 * d2 + 2 * d3 * al1 * al2 - 2 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al2 * a2 + 2 * al1 * a1 - 2 * al1 * a2 - 2 * al2 * a1  - 2 * al3 * a3 * al1 * al2 - 2 * al3 * a3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al3 * d2 - 2 * al3 * al1 * al2 * d2 + 2 * al3 * d3 * al1 * al2 - 2 * al3 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * a2 * al3 * al2 + 2 * al3 * al1 * a1 - 2 * al3 * al1 * a2 - 2 * a1 * al3 * al2  + 2 * a3 * al1 * al2 + 2 * a3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (+2 * a2 * al3 * al1 + 2 * a1 * al3 * al1 + 2 * a1 * al3 * al2 + 2 * a2 * al3 * al2 - 2 * a3 * al1 * al2 + 2 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * d2 * al3 - 2 * al1 * d2 * al3 * al2  + 2 * al3 * d3 * al1 * al2 + 2 * d3 * al3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (+2 * al1 * a2 + 2 * al1 * a1 + 2 * al2 * a1 + 2 * al2 * a2 + 2 * al3 * al1 * al2 * a3 - 2 * al3 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * d2 - 2 * al1 * d2 * al2  + 2 * d3 * al1 * al2 + 2 * d3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-4 * al1 * al2 - 4);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-4 * al3 * (1 + al1 * al2));},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (4 * al3 * (-1 + al1 * al2));},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (4 * al1 * al2 - 4);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h3_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (+2 * al3 * a1 + 2 * al3 * a2 - 2 * al3 * al1 * al2 * a2 - 2 * al3 * al1 * al2 * a1 - 2 * al2 * a3 - 2 * al1 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * al1 * d2 * al3 - 2 * d2 * al3 * al2  - 2 * d3 * al1 * al3 + 2 * al3 * d3 * al2 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * a1 - 2 * a2 + 2 * a2 * al1 * al2 + 2 * a1 * al1 * al2 - 2 * a3 * al3 * al2 - 2 * a3 * al3 * al1);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * d2 * al1 + 2 * d2 * al2  + 2 * d3 * al1 - 2 * d3 * al2 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (+2 * al2 * d2 - 2 * al1 * d2 - 2 * d3 * al2 - 2 * d3 * al1);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * a2 - 2 * a1 + 2 * al1 * al2 * a2 - 2 * al1 * a1 * al2  + 2 * a3 * al3 * al2 - 2 * al3 * a3 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (- 2 * al3 * al2 * d2 + 2 * al3 * al1 * d2 + 2 * al3 * d3 * al2 + 2 * al3 * d3 * al1);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * al3 * a2 + 2 * al3 * a1 - 2 * al3 * al1 * al2 * a2 + 2 * al3 * al1 * a1 * al2  + 2 * al2 * a3 - 2 * a3 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (4 * al3 * (al1 + al2));},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-4 * al2 - 4 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (4 * al2 - 4 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (4 * al3 * (al1 - al2));},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h4_tc_v2(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al2 * a1 - 2 * al1 * a1 - 2 * al2 * a2 - 2 * al1 * a2 + 2 * al3 * a3 - 2 * al3 * al1 * al2 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * d2 + 2 * al1 * d2 * al2  - 2 * d3 * al1 * al2 - 2 * d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * a1 * al3 * al2 - 2 * a1 * al3 * al1 - 2 * a2 * al3 * al2 - 2 * a2 * al3 * al1 - 2 * a3 + 2 * a3 * al1 * al2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * d2 * al3 + 2 * al1 * d2 * al3 * al2  - 2 * al3 * d3 * al1 * al2 - 2 * d3 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al3 * d2 - 2 * al3 * al1 * al2 * d2 + 2 * al3 * d3 * al1 * al2 - 2 * al3 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * a2 * al3 * al2 + 2 * al3 * al1 * a1 - 2 * al3 * al1 * a2 - 2 * a1 * al3 * al2  + 2 * a3 * al1 * al2 + 2 * a3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * d2 - 2 * al1 * al2 * d2 + 2 * d3 * al1 * al2 - 2 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al2 * a2 + 2 * al1 * a1 - 2 * al1 * a2 - 2 * al2 * a1  - 2 * al3 * a3 * al1 * al2 - 2 * al3 * a3 );},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (4 - 4 * al1 * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return -(4 * al3 * (-1 + al1 * al2));},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return -(4 * al3 * (1 + al1 * al2));},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-4 * al1 * al2 - 4);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h1_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al3 * a2 - al3 * a3 - al2 * a3 + a1 * al1 - al2 * a2 - a1 * al1 * al2 * al3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-d2-al3*al2*d2+al3*al2*d3-d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (a1 * al1 * al3 - al3 * al2 * a2 + a3 + a2 + a1 * al1 * al2 - al3 * al2 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-al3*d3-al2*d3-al3*d2+al2*d2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al2 * d3 + al3 * d3 + al2 * d2 + al3 * d2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (a3-a2-al3*al2*a2+al3*al2*a3-a1*al1*al2+a1*al1*al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (d3 + al3 * al2 * d3 - al3 * al2 * d2 + d2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (al3*a2+a1*al1*al2*al3+al2*a3-al2*a2+a1*al1-al3*a3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 + 2 * al2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al2 - 2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al2 - 2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * al2 * al3 - 2);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h2_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al1 * al3 * al2 * a3 - a1 * al2 - al1 * al3 * al2 * a2 - a1 * al3 + al1 * a2 + al1 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-al1*al3*d2-al1*al2*d3+al1*al2*d2-al1*al3*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (a1 + al1 * al2 * a3 + al1 * al3 * a2 - a1 * al2 * al3 + al1 * al3 * a3 + al1 * al2 * a2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (al1*al3*al2*d2+al1*d2-al1*al3*al2*d3+al1*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (al1 * d2 - al1 * al3 * al2 * d2 + al1 * d3 + al1 * al3 * al2 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-a1*al2*al3-al1*al3*a3-a1+al1*al3*a2+al1*al2*a3-al1*al2*a2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al1 * al2 * d2 - al1 * al3 * d3 - al1 * al3 * d2 + al1 * al2 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-al1*a3+a1*al3-a1*al2+al1*a2+al1*al3*al2*a2-al1*al3*al2*a3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al1 * al3 - 2 * al1 * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * al1 - 2 * al1 * al2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * al1 - 2 * al1 * al2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-2 * al1 * al2 + 2 * al1 * al3);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h3_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al1 * al2 * d3 + al1 * al3 * d3 + al1 * al2 * d2 + al1 * al3 * d2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (al1*a3+a1*al2-a1*al3-al1*al3*al2*a2-al1*a2+al1*al3*al2*a3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al1 * d2 + al1 * al3 * al2 * d2 - al1 * d3 - al1 * al3 * al2 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (a1*al2*al3+al1*al2*a2+al1*al3*a3+a1-al1*al2*a3-al1*al3*a2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (a1 + al1 * al2 * a3 + al1 * al3 * a2 - a1 * al2 * al3 + al1 * al3 * a3 + al1 * al2 * a2);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (al1*al3*al2*d2+al1*d2-al1*al3*al2*d3+al1*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al1 * al3 * al2 * a3 - a1 * al2 - al1 * al3 * al2 * a2 - a1 * al3 + al1 * a2 + al1 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-al1*al3*d2-al1*al2*d3+al1*al2*d2-al1*al3*d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al1 * al2 - 2 * al1 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al1 * al2 * al3 + 2 * al1);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * al1 - 2 * al1 * al2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al1 * al3 - 2 * al1 * al2);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

  std::vector<Polynomial> h4_tc_v3(Input& a)
  {
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], d3 = a.d[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-d2 + al3 * al2 * d2 - al3 * al2 * d3 - d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (al2*a2-al3*a2-a1*al1+al3*a3-al2*a3-a1*al1*al2*al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al3 * d2 - al2 * d2 - al3 * d3 + al2 * d3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (al3*al2*a2-a3-a1*al1*al3-al3*al2*a3+a2+a1*al1*al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (a1 * al1 * al3 - al3 * al2 * a2 + a3 + a2 + a1 * al1 * al2 - al3 * al2 * a3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-al3*d3-al2*d3-al3*d2+al2*d2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-al3 * a2 - al3 * a3 - al2 * a3 + a1 * al1 - al2 * a2 - a1 * al1 * al2 * al3);},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (-d2-al3*al2*d2+al3*al2*d3-d3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al2 * al3 + 2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return 0.0;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double al1,double al2,double al3){return (2 * al3 - 2 * al2);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (-2 * al2 - 2 * al3);},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double al1,double al2,double al3){return (2 * al3 * al2 - 2);},p1,q1,p2,q2,p3,q3)));

    return result;
  };

}
