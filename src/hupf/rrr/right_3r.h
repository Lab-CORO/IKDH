#pragma once
#include <hupf/proj_eval.h>
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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t0 * al6 * a6 - 2 * t4 + t2 * al4 * d5 - 2 * t4 * al4 * al6 + t3 * d4 + t3 * d6 + t0 * al4 * a6 - t2 * al4 * d4 - t0 * al4 * a4 + t0 * a5 * al5 + t0 * al6 * a4 + t2 * al6 * d5 + t2 * al6 * d4 + t2 * al6 * d6 + t2 * d6 * al4 - t1 * a5 * al5 * al4 - t1 * al6 * al4 * a4 + t1 * al6 * a5 * al5 + t1 * al6 * a6 * al4 - t3 * al6 * al4 * d5 + t3 * al6 * al4 * d4 - t3 * al6 * al4 * d6 - t1 * a4 + t1 * a6 + t3 * d5 + t0 * al6 * a5 * al5 * al4 + 2 * t5 * al4 - 2 * t5 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t1*al6*d4-t3*a5*al5-2*t6*al4+2*t7+t0*d5-t3*al6*a5*al5*al4-t2*al6*a6*al4+t2*a4-t0*al6*al4*d5+t1*al4*d5+t1*al6*d5+t1*al6*d6+t3*al4*a4-t3*al6*a4+t3*al6*a6-t3*a6*al4+2*t7*al6*al4+t2*a5*al5*al4-t2*al6*a5*al5+t2*al6*al4*a4-t1*al4*d4+t1*al4*d6+t0*al6*al4*d4-t0*al6*al4*d6+2*t6*al6+t0*d6+t0*d4-t2*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t1 * a5 * al5 + t1 * al4 * a6 - 2 * t5 - 2 * t4 * al4 - t2 * d6 - t1 * al6 * a6 - t0 * al6 * a5 * al5 + t3 * al6 * d5 - 2 * t5 * al4 * al6 - t3 * al4 * d4 + t3 * al6 * d4 + t3 * al4 * d5 + t0 * a5 * al5 * al4 + t0 * al6 * al4 * a4 - t0 * al6 * a6 * al4 + t2 * al6 * al4 * d5 - t2 * al6 * al4 * d4 - t1 * al4 * a4 + t2 * al6 * al4 * d6 + t1 * al6 * a4 + t3 * al6 * d6 + t3 * d6 * al4 + t0 * a4 - t0 * a6 - t2 * d5 - t2 * d4 + 2 * t4 * al6 + t1 * al6 * a5 * al5 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t3*al6*a5*al5-2*t6-t2*al6*a6-t3*a6+2*t7*al6+t2*a5*al5+t3*a4+t1*d5+t2*al6*a5*al5*al4-t0*al4*d5-t0*al6*d5-t0*al6*d6-t2*al4*a4+t2*al6*a4+t2*a6*al4-2*t6*al6*al4-2*t7*al4-t1*al6*al4*d5+t3*a5*al5*al4+t3*al6*al4*a4-t3*al6*a6*al4+t1*d6+t1*d4+t1*al6*al4*d4-t1*al6*al4*d6+t0*al4*d4-t0*al6*d4-t0*al4*d6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t3 * al6 * a5 * al5 - 2 * t6 - t2 * al6 * a6 - t3 * a6 + 2 * t7 * al6 + t2 * a5 * al5 + t3 * a4 + t1 * d5 + t2 * al6 * a5 * al5 * al4 - t0 * al4 * d5 - t0 * al6 * d5 - t0 * al6 * d6 - t2 * al4 * a4 + t2 * al6 * a4 + t2 * a6 * al4 - 2 * t6 * al6 * al4 - 2 * t7 * al4 - t1 * al6 * al4 * d5 + t3 * a5 * al5 * al4 + t3 * al6 * al4 * a4 - t3 * al6 * a6 * al4 + t1 * d6 + t1 * d4 + t1 * al6 * al4 * d4 - t1 * al6 * al4 * d6 + t0 * al4 * d4 - t0 * al6 * d4 - t0 * al4 * d6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t1*a5*al5-t1*al4*a6+2*t5+2*t4*al4+t2*d6+t1*al6*a6+t0*al6*a5*al5-t3*al6*d5+2*t5*al4*al6+t3*al4*d4-t3*al6*d4-t3*al4*d5-t0*a5*al5*al4-t0*al6*al4*a4+t0*al6*a6*al4-t2*al6*al4*d5+t2*al6*al4*d4+t1*al4*a4-t2*al6*al4*d6-t1*al6*a4-t3*al6*d6-t3*d6*al4-t0*a4+t0*a6+t2*d5+t2*d4-2*t4*al6-t1*al6*a5*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t1 * al6 * d4 + t3 * a5 * al5 + 2 * t6 * al4 - 2 * t7 - t0 * d5 + t3 * al6 * a5 * al5 * al4 + t2 * al6 * a6 * al4 - t2 * a4 + t0 * al6 * al4 * d5 - t1 * al4 * d5 - t1 * al6 * d5 - t1 * al6 * d6 - t3 * al4 * a4 + t3 * al6 * a4 - t3 * al6 * a6 + t3 * a6 * al4 - 2 * t7 * al6 * al4 - t2 * a5 * al5 * al4 + t2 * al6 * a5 * al5 - t2 * al6 * al4 * a4 + t1 * al4 * d4 - t1 * al4 * d6 - t0 * al6 * al4 * d4 + t0 * al6 * al4 * d6 - 2 * t6 * al6 - t0 * d6 - t0 * d4 + t2 * a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*al6*a6-2*t4+t2*al4*d5-2*t4*al4*al6+t3*d4+t3*d6+t0*al4*a6-t2*al4*d4-t0*al4*a4+t0*a5*al5+t0*al6*a4+t2*al6*d5+t2*al6*d4+t2*al6*d6+t2*d6*al4-t1*a5*al5*al4-t1*al6*al4*a4+t1*al6*a5*al5+t1*al6*a6*al4-t3*al6*al4*d5+t3*al6*al4*d4-t3*al6*al4*d6-t1*a4+t1*a6+t3*d5+t0*al6*a5*al5*al4+2*t5*al4-2*t5*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t0 - 2 * t0 * al6 * al4 + 2 * t1 * al4 - 2 * t1 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al4+2*t2*al6+2*t3+2*t3*al6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t1 - 2 * t1 * al6 * al4 - 2 * t0 * al4 + 2 * t0 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al4+2*t3*al6-2*t2-2*t2*al6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 * al4 + 2 * t3 * al6 - 2 * t2 - 2 * t2 * al6 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al4-2*t0*al6+2*t1+2*t1*al6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 - 2 * t3 * al6 * al4 + 2 * t2 * al4 - 2 * t2 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0-2*t0*al6*al4+2*t1*al4-2*t1*al6);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t4 * al5 * al6 - t1 * a5 - 2 * t5 * al5 - t3 * al5 * al4 * d5 + t1 * al6 * al5 * a4 - t1 * al6 * a6 * al5 - t1 * al5 * al4 * a4 + t3 * al5 * al4 * d4 - t2 * al6 * al5 * al4 * d5 + t2 * al6 * al5 * al4 * d4 - t0 * al6 * al5 * al4 * a6 + t2 * al5 * d5 + t2 * al6 * al5 * al4 * d6 - t3 * al6 * al5 * d5 + t1 * a6 * al5 * al4 + t0 * al6 * al5 * al4 * a4 - 2 * t5 * al6 * al5 * al4 - t3 * al6 * al5 * d4 + t3 * d6 * al5 * al4 + t0 * al6 * a5-t1*al6*a5*al4+t3*al6*al5*d6-2*t4*al5*al4-t0*al5*a6+t2*al5*d4-t2*al5*d6-t0*a5*al4+t0*al5*a4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t6*al5+t3*a6*al5+2*t7*al5*al4-t3*al5*a4-2*t7*al6*al5-t1*al5*d6+t3*a5*al4-t0*al6*al5*d5+t1*al6*d6*al5*al4-t3*al6*a5+t2*al5*al4*a4-t0*al6*al5*d4+t0*al5*al4*d6+t2*al5*al6*a6-t2*al5*al4*a6+2*t6*al5*al4*al6+t0*al5*al4*d4+t2*a5+t0*al6*al5*d6-t0*al5*al4*d5+t3*al6*a6*al5*al4-t1*al6*al5*al4*d5+t1*al5*d5+t1*al5*d4+t2*al6*a5*al4+t1*al6*al5*al4*d4-t3*al6*al5*al4*a4-t2*al6*al5*a4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t4*al5+t3*al5*d4+t1*al6*a5-t1*a6*al5+t0*a5+2*t5*al6*al5-t2*al6*al5*d6-t2*al5*d6*al4-t3*d6*al5+t2*al6*al5*d5-t3*al6*al5*al4*d5+t0*al5*al4*a4-t0*al6*al5*a4-t1*a5*al4+t3*al6*al5*al4*d4+t0*al6*a5*al4+t0*al6*a6*al5+t1*al5*a4+t2*al6*al5*d4+2*t4*al5*al6*al4+t3*al6*al5*al4*d6-t0*al5*a6*al4-2*t5*al5*al4+t3*al5*d5+t1*al6*al5*al4*a4+t2*al5*al4*d5-t2*al5*al4*d4-t1*al6*al5*al4*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t3*a5+t1*al6*d6*al5+t2*al6*a5+2*t7*al5-2*t6*al5*al4+t3*al6*a5*al4+t3*al5*al4*a4+2*t6*al6*al5-t2*al5*al6*a6*al4-t0*al5*d4+t0*al5*d6+t2*al6*al5*al4*a4-t1*al5*al4*d5+t0*al6*al5*al4*d5-t3*al6*al5*a4-t0*al6*al5*d6*al4-t1*al6*al5*d5-t2*a6*al5-t0*al6*al5*al4*d4+t3*al5*al6*a6-t3*al5*al4*a6+t1*al5*al4*d6+2*t7*al5*al4*al6-t1*al5*al6*d4+t1*al5*al4*d4-t0*al5*d5-t2*a5*al4+t2*al5*a4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t3 * a5 + t1 * al6 * d6 * al5 + t2 * al6 * a5 + 2 * t7 * al5 - 2 * t6 * al5 * al4 + t3 * al6 * a5 * al4 + t3 * al5 * al4 * a4 + 2 * t6 * al6 * al5 - t2 * al5 * al6 * a6 * al4 - t0 * al5 * d4 + t0 * al5 * d6 + t2 * al6 * al5 * al4 * a4 - t1 * al5 * al4 * d5 + t0 * al6 * al5 * al4 * d5 - t3 * al6 * al5 * a4 - t0 * al6 * al5 * d6 * al4 - t1 * al6 * al5 * d5-t2*a6*al5-t0*al6*al5*al4*d4+t3*al5*al6*a6-t3*al5*al4*a6+t1*al5*al4*d6+2*t7*al5*al4*al6-t1*al5*al6*d4+t1*al5*al4*d4-t0*al5*d5-t2*a5*al4+t2*al5*a4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t4*al5-t3*al5*d4-t1*al6*a5+t1*a6*al5-t0*a5-2*t5*al6*al5+t2*al6*al5*d6+t2*al5*d6*al4+t3*d6*al5-t2*al6*al5*d5+t3*al6*al5*al4*d5-t0*al5*al4*a4+t0*al6*al5*a4+t1*a5*al4-t3*al6*al5*al4*d4-t0*al6*a5*al4-t0*al6*a6*al5-t1*al5*a4-t2*al6*al5*d4-2*t4*al5*al6*al4-t3*al6*al5*al4*d6+t0*al5*a6*al4+2*t5*al5*al4-t3*al5*d5-t1*al6*al5*al4*a4-t2*al5*al4*d5+t2*al5*al4*d4+t1*al6*al5*al4*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t6 * al5 - t3 * a6 * al5 - 2 * t7 * al5 * al4 + t3 * al5 * a4 + 2 * t7 * al6 * al5 + t1 * al5 * d6 - t3 * a5 * al4 + t0 * al6 * al5 * d5 - t1 * al6 * d6 * al5 * al4 + t3 * al6 * a5 - t2 * al5 * al4 * a4 + t0 * al6 * al5 * d4 - t0 * al5 * al4 * d6 - t2 * al5 * al6 * a6 + t2 * al5 * al4 * a6 - 2 * t6 * al5 * al4 * al6 - t0 * al5 * al4 * d4 - t2 * a5 - t0 * al6 * al5 * d6 + t0 * al5 * al4 * d5 - t3 * al6 * a6 * al5 * al4 + t1 * al6 * al5 * al4 * d5 - t1 * al5 * d5 - t1 * al5 * d4-t2*al6*a5*al4-t1*al6*al5*al4*d4+t3*al6*al5*al4*a4+t2*al6*al5*a4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t4*al5*al6-t1*a5-2*t5*al5-t3*al5*al4*d5+t1*al6*al5*a4-t1*al6*a6*al5-t1*al5*al4*a4+t3*al5*al4*d4-t2*al6*al5*al4*d5+t2*al6*al5*al4*d4-t0*al6*al5*al4*a6+t2*al5*d5+t2*al6*al5*al4*d6-t3*al6*al5*d5+t1*a6*al5*al4+t0*al6*al5*al4*a4-2*t5*al6*al5*al4-t3*al6*al5*d4+t3*d6*al5*al4+t0*al6*a5-t1*al6*a5*al4+t3*al6*al5*d6-2*t4*al5*al4-t0*al5*a6+t2*al5*d4-t2*al5*d6-t0*a5*al4+t0*al5*a4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t0 * al5 * al4 + 2 * t0 * al6 * al5 - 2 * t1 * al5 - 2 * t1 * al6 * al5 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al5+2*t2*al6*al5*al4+2*t3*al5*al4-2*t3*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t1*al5*al4+2*t1*al6*al5+2*t0*al5+2*t0*al6*al5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t3*al5+2*t3*al6*al5*al4-2*t2*al5*al4+2*t2*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t3 * al5 + 2 * t3 * al6 * al5 * al4 - 2 * t2 * al5 * al4 + 2 * t2 * al6 * al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5-2*t0*al6*al5*al4+2*t1*al5*al4-2*t1*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 * al5 * al4 + 2 * t3 * al6 * al5 - 2 * t2 * al5 - 2 * t2 * al6 * al5 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5*al4+2*t0*al6*al5-2*t1*al5-2*t1*al6*al5*al4);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t6 * al5 - t3 * a6 * al5 + 2 * t7 * al5 * al4 - t3 * al5 * a4 + 2 * t7 * al6 * al5 + t1 * al5 * d6 + t3 * a5 * al4 + t0 * al6 * al5 * d5 + t1 * al6 * d6 * al5 * al4 + t3 * al6 * a5 - t2 * al5 * al4 * a4 + t0 * al6 * al5 * d4 + t0 * al5 * al4 * d6 - t2 * al5 * al6 * a6 - t2 * al5 * al4 * a6 + 2 * t6 * al5 * al4 * al6 + t0 * al5 * al4 * d4 - t2 * a5 - t0 * al6 * al5 * d6 - t0 * al5 * al4 * d5 + t3 * al6 * a6 * al5 * al4 - t1 * al6 * al5 * al4 * d5 - t1 * al5 * d5 - t1 * al5 * d4 + t2 * al6 * a5 * al4 + t1 * al6 * al5 * al4 * d4 + t3 * al6 * al5 * al4 * a4-t2*al6*al5*a4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t4*al5*al6-t1*a5-2*t5*al5+t3*al5*al4*d5-t1*al6*al5*a4-t1*al6*a6*al5-t1*al5*al4*a4-t3*al5*al4*d4+t2*al6*al5*al4*d5-t2*al6*al5*al4*d4+t0*al6*al5*al4*a6+t2*al5*d5-t2*al6*al5*al4*d6-t3*al6*al5*d5-t1*a6*al5*al4+t0*al6*al5*al4*a4+2*t5*al6*al5*al4-t3*al6*al5*d4-t3*d6*al5*al4+t0*al6*a5+t1*al6*a5*al4+t3*al6*al5*d6+2*t4*al5*al4-t0*al5*a6+t2*al5*d4-t2*al5*d6+t0*a5*al4-t0*al5*a4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t3 * a5 - t1 * al6 * d6 * al5 - t2 * al6 * a5 - 2 * t7 * al5 - 2 * t6 * al5 * al4 + t3 * al6 * a5 * al4 - t3 * al5 * al4 * a4 - 2 * t6 * al6 * al5 - t2 * al5 * al6 * a6 * al4 + t0 * al5 * d4 - t0 * al5 * d6 - t2 * al6 * al5 * al4 * a4 - t1 * al5 * al4 * d5 + t0 * al6 * al5 * al4 * d5 - t3 * al6 * al5 * a4 - t0 * al6 * al5 * d6 * al4 + t1 * al6 * al5 * d5 + t2 * a6 * al5 - t0 * al6 * al5 * al4 * d4 - t3 * al5 * al6 * a6 - t3 * al5 * al4 * a6 + t1 * al5 * al4 * d6 + 2 * t7 * al5 * al4 * al6 + t1 * al5 * al6 * d4 + t1 * al5 * al4 * d4 + t0 * al5 * d5 - t2 * a5 * al4 + t2 * al5 * a4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t4*al5+t3*al5*d4+t1*al6*a5-t1*a6*al5+t0*a5+2*t5*al6*al5-t2*al6*al5*d6+t2*al5*d6*al4-t3*d6*al5+t2*al6*al5*d5+t3*al6*al5*al4*d5+t0*al5*al4*a4+t0*al6*al5*a4+t1*a5*al4-t3*al6*al5*al4*d4-t0*al6*a5*al4+t0*al6*a6*al5-t1*al5*a4+t2*al6*al5*d4-2*t4*al5*al6*al4-t3*al6*al5*al4*d6+t0*al5*a6*al4+2*t5*al5*al4+t3*al5*d5+t1*al6*al5*al4*a4-t2*al5*al4*d5+t2*al5*al4*d4+t1*al6*al5*al4*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t3 * al5 * d4 + 2 * t4 * al5 + t2 * al5 * al4 * d4 + t1 * a5 * al4 + t1 * al6 * al5 * al4 * a4 - t1 * al5 * a6 + t3 * al6 * al5 * al4 * d5 + t1 * al6 * al5 * al4 * a6 - t1 * al5 * a4 + t2 * d6 * al5 * al4 + t1 * al6 * a5 - t2 * al5 * al4 * d5 + t2 * al6 * al5 * d5 - t3 * al6 * al5 * al4 * d6 + 2 * t5 * al6 * al5 + t0 * a5 - t3 * d6 * al5 + t0 * al5 * al4 * a4 + t0 * al6 * al5 * a4 - t3 * al6 * al5 * al4 * d4 - t0 * al6 * a5 * al4 + t0 * al6 * a6 * al5 + t2 * al6 * al5 * d4-t2*al6*al5*d6-2*t4*al6*al5*al4+t0*al5*a6*al4+2*t5*al5*al4+t3*al5*d5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t3*a5+t1*al6*d6*al5+t2*al6*a5+2*t7*al5+2*t6*al5*al4-t3*al6*a5*al4+t3*al5*al4*a4+2*t6*al6*al5+t2*al5*al6*a6*al4-t0*al5*d4+t0*al5*d6+t2*al6*al5*al4*a4+t1*al5*al4*d5-t0*al6*al5*al4*d5+t3*al6*al5*a4+t0*al6*al5*d6*al4-t1*al6*al5*d5-t2*a6*al5+t0*al6*al5*al4*d4+t3*al5*al6*a6+t3*al5*al4*a6-t1*al5*al4*d6-2*t7*al5*al4*al6-t1*al5*al6*d4-t1*al5*al4*d4-t0*al5*d5+t2*a5*al4-t2*al5*a4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t0 * al6 * a5 - t0 * al6 * al5 * al4 * a4+t1*a5+t2*d6*al5+t1*al5*al4*a4-t2*al5*d5-2*t5*al6*al5*al4+2*t5*al5+t2*al6*al5*al4*d6+t2*al6*al5*al4*d4+t3*al6*al5*d4+t1*al5*a6*al4+t0*al5*a4-t3*al5*al4*d5-t2*al6*al5*al4*d5-t1*al6*a5*al4+t3*al5*al4*d4+t1*al6*al5*a4+t1*al6*a6*al5-t2*al5*d4+t3*d6*al5*al4+t3*al6*al5*d5-2*t4*al5*al4-t0*a5*al4+t0*al5*a6-t3*al6*al5*d6-2*t4*al6*al5-t0*al6*al5*al4*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t7*al6*al5+2*t7*al5*al4-t2*a5-t2*al6*a6*al5+t1*al6*al5*al4*d4+t3*al6*a5-t0*al6*d6*al5-2*t6*al5-t1*al6*al5*al4*d5+t3*a5*al4-t1*al5*d5+t1*al6*d6*al5*al4+t0*al5*al4*d4+t0*al6*al5*d4-t2*al5*al4*a6+t0*al5*al4*d6+2*t6*al5*al4*al6-t2*al6*al5*a4+t2*al6*a5*al4+t3*al5*al6*a6*al4-t2*al5*al4*a4-t3*a6*al5+t3*al6*al5*al4*a4-t3*al5*a4-t0*al5*al4*d5-t1*al5*d4+t1*al5*d6+t0*al6*al5*d5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t2 * al5 + 2 * t2 * al6 * al5 * al4 + 2 * t3 * al5 * al4 + 2 * t3 * al6 * al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al5*al4+2*t0*al6*al5-2*t1*al5+2*t1*al6*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 * al5 + 2 * t3 * al6 * al5 * al4 - 2 * t2 * al5 * al4 - 2 * t2 * al6 * al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t1*al5*al4+2*t1*al6*al5+2*t0*al5-2*t0*al6*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t1 * al5 * al4 + 2 * t1 * al6 * al5 + 2 * t0 * al5 - 2 * t0 * al6 * al5 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al5*al4+2*t2*al6*al5+2*t3*al5-2*t3*al6*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t1 * al5 - 2 * t1 * al6 * al5 * al4 - 2 * t0 * al5 * al4 - 2 * t0 * al6 * al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al5+2*t2*al6*al5*al4+2*t3*al5*al4+2*t3*al6*al5);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t6 * al6 + 2 * t7 - t2 * a6 + t1 * al6 * d6 + t0 * d4 - t1 * al4 * d5 + t3 * al6 * a6-t2*a4+t3*al6*a4+t1*al4*d4+2*t6*al4+t0*d6+t3*al4*a6+t0*d5+t0*al6*al4*d6-t1*d6*al4+t1*al6*d5-t3*a5*al5+t0*al6*al4*d5-2*t7*al4*al6+t1*al6*d4-t0*al6*al4*d4+t2*al6*a6*al4+t2*al6*al4*a4-t2*al6*a5*al5+t3*al4*a4+t3*al6*a5*al5*al4-t2*a5*al5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t2*al4*d5-2*t4*al6*al4+t0*al4*a4-t1*a6+2*t4+t0*al6*a4-t0*a5*al5-t3*d4-t3*d5-t3*d6+t3*al6*al4*d4-t3*al6*al4*d6-t2*al6*d5+2*t5*al6+t1*al6*a6*al4+t2*al4*d6-t2*al6*d4-t2*al4*d4+t1*al6*al4*a4-t1*a5*al5*al4-t1*a4+2*t5*al4+t0*al6*a5*al5*al4+t0*al6*a6-t3*al6*al4*d5+t0*a6*al4-t1*al6*a5*al5-t2*al6*d6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t2 * al6 * a6 + 2 * t7 * al6 - t0 * al6 * d6 - 2 * t6 - t2 * al6 * a4 - t2 * al4 * a6 - t0 * al6 * d4 - t2 * al4 * a4 + t1 * al6 * al4 * d5 - t3 * a6 + t2 * a5 * al5 + t1 * d5 - t1 * al6 * al4 * d4 + t0 * d6 * al4+t1*d4-t0*al4*d4+2*t6*al4*al6-t0*al6*d5-t3*a4+t1*d6-t2*al6*a5*al5*al4+t1*al6*al4*d6-t3*a5*al5*al4+2*t7*al4-t3*al6*a5*al5+t3*al6*al4*a4+t0*al4*d5+t3*al6*a6*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t5+t2*d6-2*t4*al6+t1*al6*a5*al5*al4-t3*al6*d4-t3*al4*d4+t3*al4*d6+t2*d5+t0*a6+t0*al6*a5*al5+t0*a5*al5*al4-t3*al6*d5-t3*al6*d6+t1*al4*a4+t1*al6*a4+t2*d4-2*t5*al6*al4+t3*al4*d5-t0*al6*al4*a4-t0*al6*a6*al4-2*t4*al4+t0*a4+t1*al6*a6+t2*al6*al4*d5+t1*a6*al4+t2*al6*al4*d6-t2*al6*al4*d4-t1*a5*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t5 + t2 * d6 - 2 * t4 * al6 + t1 * al6 * a5 * al5 * al4 - t3 * al6 * d4 - t3 * al4 * d4 + t3 * al4 * d6 + t2 * d5 + t0 * a6 + t0 * al6 * a5 * al5 + t0 * a5 * al5 * al4 - t3 * al6 * d5 - t3 * al6 * d6+t1*al4*a4+t1*al6*a4+t2*d4-2*t5*al6*al4+t3*al4*d5-t0*al6*al4*a4-t0*al6*a6*al4-2*t4*al4+t0*a4+t1*al6*a6+t2*al6*al4*d5+t1*a6*al4+t2*al6*al4*d6-t2*al6*al4*d4-t1*a5*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t2*al6*a6-2*t7*al6+t0*al6*d6+2*t6+t2*al6*a4+t2*al4*a6+t0*al6*d4+t2*al4*a4-t1*al6*al4*d5+t3*a6-t2*a5*al5-t1*d5+t1*al6*al4*d4-t0*d6*al4-t1*d4+t0*al4*d4-2*t6*al4*al6+t0*al6*d5+t3*a4-t1*d6+t2*al6*a5*al5*al4-t1*al6*al4*d6+t3*a5*al5*al4-2*t7*al4+t3*al6*a5*al5-t3*al6*al4*a4-t0*al4*d5-t3*al6*a6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t2 * al4 * d5 + 2 * t4 * al6 * al4 - t0 * al4 * a4 + t1 * a6 - 2 * t4 - t0 * al6 * a4 + t0 * a5 * al5 + t3 * d4 + t3 * d5 + t3 * d6 - t3 * al6 * al4 * d4 + t3 * al6 * al4 * d6 + t2 * al6 * d5 - 2 * t5 * al6 - t1 * al6 * a6 * al4 - t2 * al4 * d6 + t2 * al6 * d4 + t2 * al4 * d4 - t1 * al6 * al4 * a4 + t1 * a5 * al5 * al4 + t1 * a4 - 2 * t5 * al4-t0*al6*a5*al5*al4-t0*al6*a6+t3*al6*al4*d5-t0*a6*al4+t1*al6*a5*al5+t2*al6*d6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t6*al6+2*t7-t2*a6+t1*al6*d6+t0*d4-t1*al4*d5+t3*al6*a6-t2*a4+t3*al6*a4+t1*al4*d4+2*t6*al4+t0*d6+t3*al4*a6+t0*d5+t0*al6*al4*d6-t1*d6*al4+t1*al6*d5-t3*a5*al5+t0*al6*al4*d5-2*t7*al4*al6+t1*al6*d4-t0*al6*al4*d4+t2*al6*a6*al4+t2*al6*al4*a4-t2*al6*a5*al5+t3*al4*a4+t3*al6*a5*al5*al4-t2*a5*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t2 * al4 + 2 * t2 * al6 + 2 * t3 - 2 * t3 * al6 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0-2*t0*al6*al4+2*t1*al4+2*t1*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t3 * al4 + 2 * t3 * al6 - 2 * t2 + 2 * t2 * al6 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t1-2*t1*al6*al4-2*t0*al4-2*t0*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t1 - 2 * t1 * al6 * al4 - 2 * t0 * al4 - 2 * t0 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2-2*t2*al6*al4-2*t3*al4-2*t3*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t1 * al4 - 2 * t1 * al6 - 2 * t0 + 2 * t0 * al6 * al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al4+2*t2*al6+2*t3-2*t3*al6*al4);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t1*d4*al5+t1*al5*d5-t1*d6*al5-t3*al6*a5-t3*al6*a4+t3*a6*al5-t1*d6*al4+t3*a6*al4+t2*a4+2*t6*al4+2*t6*al5+t1*d4*al4-t1*al4*d5-t0*al6*d4*al5-2*t7*al6*al4-2*t7*al6*al5+t2*a5-t0*al6*al5*d5+t0*al6*al4*d5-t0*al6*d4*al4+t0*al6*d6*al4+t0*al6*d6*al5-t2*a4*al5*al4-t2*a5*al5*al4+t2*al6*a6*al4+t2*al6*a6*al5+t3*al6*a4*al5*al4+t3*al6*a5*al5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*al6*a4+t2*d4*al5-t0*a6*al5+t2*d5*al5-t1*a5+t1*a4-2*t5*al5+2*t5*al4+t0*al6*a5+t0*a6*al4+t2*d5*al4-t2*d4*al4-t2*d6*al5+t2*d6*al4+2*t4*al6*al5-2*t4*al6*al4-t0*al6*al5*al4*a4+t1*al5*al4*a4+t0*al6*al5*a5*al4-t1*al5*a5*al4-t1*al6*a6*al5+t1*al6*a6*al4-t3*al6*d5*al5-t3*al6*d5*al4-t3*al6*d4*al5+t3*al6*d4*al4+t3*al6*d6*al5-t3*al6*d6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t0*al4*d5+2*t6*al6*al4-t0*al5*d5+t2*al6*a5-t0*d4*al5+t0*d6*al4-t0*d4*al4-t2*a6*al4-t2*a6*al5+t2*al6*a4-t1*al6*d4*al5-t1*al6*al5*d5+t1*al6*al4*d5-t1*al6*d4*al4+t3*a4+t1*al6*d6*al4+t1*al6*d6*al5+t3*a5-t3*a4*al5*al4-t3*a5*al5*al4+t3*al6*a6*al5+t3*al6*a6*al4+t0*d6*al5+2*t6*al6*al5-t2*al6*a4*al5*al4-t2*al6*a5*al5*al4+2*t7*al4+2*t7*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t3*d5*al5+2*t5*al6*al5+t2*al6*d5*al5+2*t4*al5+t0*a5-2*t4*al4-t1*al6*a4+t1*al6*a5-t1*a6*al5+t1*a6*al4+t3*d5*al4+t3*d4*al5-t3*d4*al4-t3*d6*al5+t3*d6*al4-2*t5*al6*al4-t0*a4-t0*al5*al4*a4-t1*al6*al5*al4*a4+t0*al5*a5*al4+t1*al6*al5*a5*al4+t0*al6*a6*al5-t0*al6*a6*al4+t2*al6*d5*al4+t2*al6*d4*al5-t2*al6*d4*al4-t2*al6*d6*al5+t2*al6*d6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t3*d5*al5-2*t5*al6*al5-t2*al6*d5*al5-2*t4*al5-t0*a5-2*t4*al4-t1*al6*a4-t1*al6*a5+t1*a6*al5+t1*a6*al4+t3*d5*al4-t3*d4*al5-t3*d4*al4+t3*d6*al5+t3*d6*al4-2*t5*al6*al4-t0*a4+t0*al5*al4*a4+t1*al6*al5*al4*a4+t0*al5*a5*al4+t1*al6*al5*a5*al4-t0*al6*a6*al5-t0*al6*a6*al4+t2*al6*d5*al4-t2*al6*d4*al5-t2*al6*d4*al4+t2*al6*d6*al5+t2*al6*d6*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*al4*d5-2*t6*al6*al4-t0*al5*d5+t2*al6*a5-t0*d4*al5-t0*d6*al4+t0*d4*al4+t2*a6*al4-t2*a6*al5-t2*al6*a4-t1*al6*d4*al5-t1*al6*al5*d5-t1*al6*al4*d5+t1*al6*d4*al4-t3*a4-t1*al6*d6*al4+t1*al6*d6*al5+t3*a5-t3*a4*al5*al4+t3*a5*al5*al4+t3*al6*a6*al5-t3*al6*a6*al4+t0*d6*al5+2*t6*al6*al5-t2*al6*a4*al5*al4+t2*al6*a5*al5*al4-2*t7*al4+2*t7*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t0*al6*a4+t2*d4*al5-t0*a6*al5+t2*d5*al5-t1*a5-t1*a4-2*t5*al5-2*t5*al4+t0*al6*a5-t0*a6*al4-t2*d5*al4+t2*d4*al4-t2*d6*al5-t2*d6*al4+2*t4*al6*al5+2*t4*al6*al4-t0*al6*al5*al4*a4+t1*al5*al4*a4-t0*al6*al5*a5*al4+t1*al5*a5*al4-t1*al6*a6*al5-t1*al6*a6*al4-t3*al6*d5*al5+t3*al6*d5*al4-t3*al6*d4*al5-t3*al6*d4*al4+t3*al6*d6*al5+t3*al6*d6*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t1*d4*al5-t1*al5*d5+t1*d6*al5+t3*al6*a5-t3*al6*a4-t3*a6*al5-t1*d6*al4+t3*a6*al4+t2*a4+2*t6*al4-2*t6*al5+t1*d4*al4-t1*al4*d5+t0*al6*d4*al5-2*t7*al6*al4+2*t7*al6*al5-t2*a5+t0*al6*al5*d5+t0*al6*al4*d5-t0*al6*d4*al4+t0*al6*d6*al4-t0*al6*d6*al5+t2*a4*al5*al4-t2*a5*al5*al4+t2*al6*a6*al4-t2*al6*a6*al5-t3*al6*a4*al5*al4+t3*al6*a5*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al4+2*t2*al5-2*t3*al6*al4-2*t3*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al6*al5-2*t0*al6*al4-2*t1*al5+2*t1*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t3*al4+2*t3*al5+2*t2*al6*al4+2*t2*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t1*al6*al5-2*t1*al6*al4+2*t0*al5-2*t0*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al4-2*t0*al5-2*t1*al6*al4-2*t1*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al6*al5-2*t2*al6*al4+2*t3*al5-2*t3*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t1*al4-2*t1*al5+2*t0*al6*al4+2*t0*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t3*al6*al5-2*t3*al6*al4-2*t2*al5+2*t2*al4);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t3*al5*a4-t3*al5*a5+2*t6*al6+2*t7+t1*al6*d4-t3*al4*a5-2*t7*al5*al4+t1*al6*d6-t3*al4*a4+t1*al6*d5+t3*al6*a6+t0*d4+t0*d6-t2*a6-t0*d4*al5*al4+t0*al5*al4*d5-t0*d6*al5*al4-t2*al6*al5*a5-t2*al6*al4*a4-t2*al6*al5*a4-t2*al6*al4*a5+t2*a6*al5*al4-2*t6*al6*al5*al4+t0*d5-t1*al6*d4*al5*al4+t1*al6*al5*al4*d5-t1*al6*d6*al5*al4-t3*al6*a6*al5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*al6*a6+t1*a6*al5*al4+t2*al6*d4+t3*d5-2*t5*al6+t3*d6-2*t4+t1*a6+t2*al6*d4*al5*al4+t0*al4*a4+t2*al6*d5+t2*al6*d6-2*t4*al5*al4+t3*d4-t0*al6*a6*al5*al4+t1*al6*al4*a4-t2*al6*al5*d5*al4+t2*al6*d6*al5*al4-t3*al5*d5*al4+t3*d4*al5*al4+t3*d6*al5*al4-2*t5*al6*al5*al4-t0*al5*a4+t0*al5*a5-t0*al4*a5-t1*al6*al5*a4-t1*al6*al4*a5+t1*al6*al5*a5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t2*al4*a5-2*t6-t1*d4*al5*al4+t2*al5*a4-t2*al6*a6-t0*al6*d5+2*t6*al5*al4-t0*al6*d4+2*t7*al6+t1*d4+t1*d6+t2*al4*a4-t0*al6*d6+t2*al5*a5+t1*al5*al4*d5-t1*d6*al5*al4-t3*al6*al5*a5-t3*al6*al4*a5-t3*al6*al4*a4-t3*al6*al5*a4-t3*a6+t3*a6*al5*al4-2*t7*al6*al5*al4+t1*d5+t0*al6*d4*al5*al4-t0*al6*al5*al4*d5+t0*al6*d6*al5*al4+t2*al6*a6*al5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t5-t2*d6-t2*d4-t0*a6-t2*d5-t1*al6*a6-t1*al6*a6*al5*al4+t1*al4*a4+t3*al6*d5+t3*al6*d4+2*t4*al6+t3*al6*d6-2*t5*al5*al4-t0*al6*al4*a4-t0*a6*al5*al4+t3*al6*d4*al5*al4-t3*al6*al5*d5*al4+t3*al6*d6*al5*al4+t2*al5*d5*al4-t2*d4*al5*al4-t2*d6*al5*al4+2*t4*al6*al5*al4+t0*al6*al5*a4+t0*al6*al4*a5-t0*al6*al5*a5-t1*al5*a4-t1*al4*a5+t1*al5*a5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t5+t2*d6+t2*d4+t0*a6+t2*d5+t1*al6*a6-t1*al6*a6*al5*al4-t1*al4*a4-t3*al6*d5-t3*al6*d4-2*t4*al6-t3*al6*d6-2*t5*al5*al4+t0*al6*al4*a4-t0*a6*al5*al4+t3*al6*d4*al5*al4-t3*al6*al5*d5*al4+t3*al6*d6*al5*al4+t2*al5*d5*al4-t2*d4*al5*al4-t2*d6*al5*al4+2*t4*al6*al5*al4+t0*al6*al5*a4+t0*al6*al4*a5+t0*al6*al5*a5-t1*al5*a4-t1*al4*a5-t1*al5*a5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t2*al4*a5-2*t6+t1*d4*al5*al4-t2*al5*a4-t2*al6*a6-t0*al6*d5-2*t6*al5*al4-t0*al6*d4+2*t7*al6+t1*d4+t1*d6+t2*al4*a4-t0*al6*d6+t2*al5*a5-t1*al5*al4*d5+t1*d6*al5*al4-t3*al6*al5*a5+t3*al6*al4*a5-t3*al6*al4*a4+t3*al6*al5*a4-t3*a6-t3*a6*al5*al4+2*t7*al6*al5*al4+t1*d5-t0*al6*d4*al5*al4+t0*al6*al5*al4*d5-t0*al6*d6*al5*al4-t2*al6*a6*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t0*al6*a6-t1*a6*al5*al4+t2*al6*d4+t3*d5-2*t5*al6+t3*d6-2*t4+t1*a6-t2*al6*d4*al5*al4+t0*al4*a4+t2*al6*d5+t2*al6*d6+2*t4*al5*al4+t3*d4+t0*al6*a6*al5*al4+t1*al6*al4*a4+t2*al6*al5*d5*al4-t2*al6*d6*al5*al4+t3*al5*d5*al4-t3*d4*al5*al4-t3*d6*al5*al4+2*t5*al6*al5*al4+t0*al5*a4+t0*al5*a5+t0*al4*a5+t1*al6*al5*a4+t1*al6*al4*a5+t1*al6*al5*a5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t3*al5*a4+t3*al5*a5-2*t6*al6-2*t7-t1*al6*d4-t3*al4*a5-2*t7*al5*al4-t1*al6*d6+t3*al4*a4-t1*al6*d5-t3*al6*a6-t0*d4-t0*d6+t2*a6-t0*d4*al5*al4+t0*al5*al4*d5-t0*d6*al5*al4+t2*al6*al5*a5+t2*al6*al4*a4-t2*al6*al5*a4-t2*al6*al4*a5+t2*a6*al5*al4-2*t6*al6*al5*al4-t0*d5-t1*al6*d4*al5*al4+t1*al6*al5*al4*d5-t1*al6*d6*al5*al4-t3*al6*a6*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al6*al5*al4+2*t2*al6-2*t3*al5*al4+2*t3);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5*al4-2*t0-2*t1*al6*al5*al4-2*t1*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al6*al5*al4+2*t3*al6+2*t2*al5*al4-2*t2);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t1*al5*al4-2*t1+2*t0*al6*al5*al4+2*t0*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al6*al5*al4-2*t0*al6-2*t1*al5*al4+2*t1);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al5*al4-2*t2+2*t3*al6*al5*al4+2*t3*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t1*al6*al5*al4-2*t1*al6+2*t0*al5*al4-2*t0);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al5*al4-2*t3-2*t2*al6*al5*al4-2*t2*al6);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t1*al6*a6*al5+2*t5*al5+t0*a6*al5+t3*al6*d4*al5+t3*al6*al5*d5-t3*al6*d6*al5-t2*al5*d5-t2*d4*al5+t2*d6*al5-2*t4*al6*al5+t0*a6*al4+t2*d6*al4-2*t4*al6*al4+t2*al4*d5-t0*al6*a4-t2*d4*al4-t0*al6*a5-t3*al6*al4*d5+t3*al6*d4*al4+t1*a4-t3*al6*d6*al4-t1*al5*a4*al4+t1*a5-t1*al5*a5*al4+t1*al6*a6*al4+2*t5*al4+t0*al6*al5*a4*al4+t0*al6*al5*a5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al6*al5+t1*d6*al4-t1*d6*al5-t1*d4*al4+t1*d4*al5+t1*d5*al4+t1*al5*d5-t3*a6*al4-t3*al6*a5+t3*a6*al5+2*t7*al6*al4-t2*a4+t2*a5+2*t6*al5+t0*al6*d6*al5-2*t6*al4-t0*al6*d6*al4+t3*al6*al5*al4*a4+t3*al6*a4-t3*al6*al5*a5*al4-t0*al6*al5*d5-t2*al5*al4*a4+t2*al5*a5*al4-t2*al6*a6*al4+t2*al6*a6*al5-t0*al6*d4*al5-t0*al6*d5*al4+t0*al6*d4*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t4*al5+t1*a6*al5-t3*al5*d5-t3*d4*al4+t0*a5*al5*al4-t1*al6*a4+t2*al6*d6*al4+t1*al6*a5*al5*al4-2*t5*al6*al4+t3*al4*d5-t2*al6*d4*al4+t1*a6*al4-t0*al6*a6*al5-t1*al6*a5-t0*a4+t2*al6*d6*al5+t2*al6*al4*d5-t2*al6*al5*d5-t3*d4*al5-t0*a5+t3*d6*al4+t0*a4*al5*al4-t2*al6*d4*al5-t0*al6*a6*al4-2*t4*al4+t1*al6*a4*al5*al4+t3*d6*al5-2*t5*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*d6*al4+t2*al6*al5*a5*al4+t2*a6*al4+t0*d4*al4-t0*d5*al4-2*t6*al6*al4-t1*al6*d5*al5-t2*al6*al5*al4*a4-t1*al6*d5*al4-t3*a4-t2*a6*al5-t0*d5*al5+2*t6*al6*al5+2*t7*al5+t3*a5+t3*al6*a6*al5-t1*al6*d6*al4-2*t7*al4-t2*al6*a4+t0*d6*al5+t3*al5*a5*al4-t1*al6*d4*al5-t3*al6*a6*al4-t0*d4*al5+t2*al6*a5+t1*al6*d4*al4-t3*al5*al4*a4+t1*al6*d6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t0*d6*al4+t2*al6*al5*a5*al4+t2*a6*al4+t0*d4*al4-t0*d5*al4-2*t6*al6*al4+t1*al6*d5*al5+t2*al6*al5*al4*a4-t1*al6*d5*al4-t3*a4+t2*a6*al5+t0*d5*al5-2*t6*al6*al5-2*t7*al5-t3*a5-t3*al6*a6*al5-t1*al6*d6*al4-2*t7*al4-t2*al6*a4-t0*d6*al5+t3*al5*a5*al4+t1*al6*d4*al5-t3*al6*a6*al4+t0*d4*al5-t2*al6*a5+t1*al6*d4*al4+t3*al5*al4*a4-t1*al6*d6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t4*al5+t1*a6*al5-t3*al5*d5+t3*d4*al4-t0*a5*al5*al4+t1*al6*a4-t2*al6*d6*al4-t1*al6*a5*al5*al4+2*t5*al6*al4-t3*al4*d5+t2*al6*d4*al4-t1*a6*al4-t0*al6*a6*al5-t1*al6*a5+t0*a4+t2*al6*d6*al5-t2*al6*al4*d5-t2*al6*al5*d5-t3*d4*al5-t0*a5-t3*d6*al4+t0*a4*al5*al4-t2*al6*d4*al5+t0*al6*a6*al4+2*t4*al4+t1*al6*a4*al5*al4+t3*d6*al5-2*t5*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al6*al5-t1*d6*al4-t1*d6*al5+t1*d4*al4+t1*d4*al5-t1*d5*al4+t1*al5*d5+t3*a6*al4-t3*al6*a5+t3*a6*al5-2*t7*al6*al4+t2*a4+t2*a5+2*t6*al5+t0*al6*d6*al5+2*t6*al4+t0*al6*d6*al4+t3*al6*al5*al4*a4-t3*al6*a4+t3*al6*al5*a5*al4-t0*al6*al5*d5-t2*al5*al4*a4-t2*al5*a5*al4+t2*al6*a6*al4+t2*al6*a6*al5-t0*al6*d4*al5+t0*al6*d5*al4-t0*al6*d4*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t1*al6*a6*al5-2*t5*al5-t0*a6*al5-t3*al6*d4*al5-t3*al6*al5*d5+t3*al6*d6*al5+t2*al5*d5+t2*d4*al5-t2*d6*al5+2*t4*al6*al5+t0*a6*al4+t2*d6*al4-2*t4*al6*al4+t2*al4*d5-t0*al6*a4-t2*d4*al4+t0*al6*a5-t3*al6*al4*d5+t3*al6*d4*al4+t1*a4-t3*al6*d6*al4+t1*al5*a4*al4-t1*a5-t1*al5*a5*al4+t1*al6*a6*al4+2*t5*al4-t0*al6*al5*a4*al4+t0*al6*al5*a5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al6*al4-2*t0*al6*al5+2*t1*al4+2*t1*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al5-2*t2*al4-2*t3*al6*al5+2*t3*al6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al4-2*t0*al5-2*t1*al6*al4-2*t1*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al6*al5-2*t2*al6*al4+2*t3*al5-2*t3*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al6*al4-2*t2*al6*al5-2*t3*al4-2*t3*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5+2*t0*al4-2*t1*al6*al5+2*t1*al6*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al4+2*t2*al5-2*t3*al6*al4-2*t3*al6*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al6*al5-2*t0*al6*al4-2*t1*al5+2*t1*al4);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t2*al6*d6-t3*d6-t0*a5*al5-t1*a6+2*t5*al6-t0*a4*al5-t3*d4-t2*al6*d4-t1*al6*a5*al5+2*t4-t3*d5+t3*d6*al5*al4+t2*al6*d4*al5*al4+t3*d4*al5*al4-t3*al5*al4*d5-t0*al6*a6*al5*al4-t1*al6*a4*al5-t1*al6*al4*a5-t0*al4*a4+t2*al6*d6*al5*al4-t1*al6*al4*a4-t0*al4*a5+t1*a6*al5*al4-2*t5*al6*al5*al4-2*t4*al5*al4-t2*al6*al5*al4*d5-t2*al6*d5+t0*al6*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t6*al6+t0*d5-t2*a6-t2*al6*al5*a5+t1*al6*d5+t3*al6*a6+t1*al6*d6+t0*d4+2*t7-t3*al5*a5+t1*al6*d4+t3*al6*a6*al5*al4-t2*a6*al5*al4+t0*d6+2*t6*al6*al5*al4+t2*al6*al4*a5-t1*al6*d5*al5*al4-t0*d5*al5*al4+t3*al4*a5+t0*d4*al5*al4-t2*al6*al4*a4-t3*al4*a4+t1*al6*d6*al5*al4+t3*a4*al5+t2*al6*a4*al5+t1*al6*d4*al5*al4+t0*d6*al5*al4+2*t7*al5*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t5-2*t4*al6+t3*al6*d4*al5*al4-t2*d4*al5*al4-t1*al6*a6*al5*al4+2*t4*al6*al5*al4+t0*a6+t2*d5+t2*d4+t0*al6*al5*a5-t3*al6*d4+t0*al6*al4*a5-t0*a6*al5*al4+t2*al5*al4*d5-t2*d6*al5*al4+t0*al6*al4*a4+t3*al6*d6*al5*al4-t3*al6*al5*al4*d5-t1*al4*a4-t1*al5*a4+t0*al6*al5*a4-t3*al6*d5-t1*al4*a5-2*t5*al5*al4+t2*d6-t3*al6*d6-t1*al5*a5+t1*al6*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*al6*d6-2*t6-t3*al6*a5*al5-t0*al6*d4+t1*d5+2*t7*al6-t2*al6*a6-2*t6*al5*al4+2*t7*al6*al5*al4+t2*al4*a4+t1*d6*al5*al4-t0*al6*d5-t2*al6*a6*al5*al4+t1*d4*al5*al4+t1*d4+t1*d6-t3*a6*al5*al4-t3*al6*al4*a4-t2*al4*a5+t2*a5*al5+t0*al6*al5*d5*al4-t0*al6*d4*al5*al4-t1*al5*d5*al4-t0*al6*d6*al5*al4+t3*al6*al4*a5+t3*al6*a4*al5-t2*a4*al5-t3*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t0*al6*d6+2*t6+t3*al6*a5*al5+t0*al6*d4-t1*d5-2*t7*al6+t2*al6*a6-2*t6*al5*al4+2*t7*al6*al5*al4-t2*al4*a4+t1*d6*al5*al4+t0*al6*d5-t2*al6*a6*al5*al4+t1*d4*al5*al4-t1*d4-t1*d6-t3*a6*al5*al4+t3*al6*al4*a4-t2*al4*a5-t2*a5*al5+t0*al6*al5*d5*al4-t0*al6*d4*al5*al4-t1*al5*d5*al4-t0*al6*d6*al5*al4+t3*al6*al4*a5+t3*al6*a4*al5-t2*a4*al5+t3*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t5-2*t4*al6-t3*al6*d4*al5*al4+t2*d4*al5*al4+t1*al6*a6*al5*al4-2*t4*al6*al5*al4+t0*a6+t2*d5+t2*d4+t0*al6*al5*a5-t3*al6*d4-t0*al6*al4*a5+t0*a6*al5*al4-t2*al5*al4*d5+t2*d6*al5*al4+t0*al6*al4*a4-t3*al6*d6*al5*al4+t3*al6*al5*al4*d5-t1*al4*a4+t1*al5*a4-t0*al6*al5*a4-t3*al6*d5+t1*al4*a5+2*t5*al5*al4+t2*d6-t3*al6*d6-t1*al5*a5+t1*al6*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t6*al6+t0*d5-t2*a6-t2*al6*al5*a5+t1*al6*d5+t3*al6*a6+t1*al6*d6+t0*d4+2*t7-t3*al5*a5+t1*al6*d4-t3*al6*a6*al5*al4+t2*a6*al5*al4+t0*d6-2*t6*al6*al5*al4-t2*al6*al4*a5+t1*al6*d5*al5*al4+t0*d5*al5*al4-t3*al4*a5-t0*d4*al5*al4-t2*al6*al4*a4-t3*al4*a4-t1*al6*d6*al5*al4-t3*a4*al5-t2*al6*a4*al5-t1*al6*d4*al5*al4-t0*d6*al5*al4-2*t7*al5*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t2*al6*d6+t3*d6+t0*a5*al5+t1*a6-2*t5*al6-t0*a4*al5+t3*d4+t2*al6*d4+t1*al6*a5*al5-2*t4+t3*d5+t3*d6*al5*al4+t2*al6*d4*al5*al4+t3*d4*al5*al4-t3*al5*al4*d5-t0*al6*a6*al5*al4-t1*al6*a4*al5-t1*al6*al4*a5+t0*al4*a4+t2*al6*d6*al5*al4+t1*al6*al4*a4-t0*al4*a5+t1*a6*al5*al4-2*t5*al6*al5*al4-2*t4*al5*al4-t2*al6*al5*al4*d5+t2*al6*d5-t0*al6*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5*al4+2*t0-2*t1*al6*al5*al4+2*t1*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al6*al5*al4+2*t2*al6+2*t3*al5*al4+2*t3);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al6*al5*al4-2*t0*al6-2*t1*al5*al4+2*t1);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al5*al4-2*t2+2*t3*al6*al5*al4+2*t3*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al5*al4+2*t2+2*t3*al6*al5*al4-2*t3*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al6*al5*al4-2*t0*al6+2*t1*al5*al4+2*t1);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al6*al5*al4+2*t2*al6-2*t3*al5*al4+2*t3);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5*al4-2*t0-2*t1*al6*al5*al4-2*t1*al6);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t3*d5+t0*al6*a4*al4*al5+t3*d6+t0*a4*al4-2*t4*al5*al6-2*t5*al6+t0*al5*a6-t2*d5*al5+2*t5*al5-t2*d4*al5+t2*al6*d4+t2*al6*d5+t2*al5*d6-2*t4+t2*al6*d6-t0*al5*a5+t0*al6*a5-t0*al6*a6-t1*al6*al5*a5-t1*a4*al4*al5+t3*d4+t1*al6*a4*al4+t1*al6*al5*a6+t3*al6*d5*al5+t3*al6*d4*al5-t3*al5*al6*d6-t1*a5+t1*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t1*d5*al5-t0*al6*d5*al5+t3*al6*a5+2*t6*al5-t2*a5+t0*d4+t0*d5+t1*d4*al5+t1*al6*d4+t1*al6*d5-t3*a4*al4+t3*al5*a5+2*t7-t0*al6*d4*al5-t2*a4*al4*al5+t2*al6*al5*a5-t2*al6*a4*al4+t3*al6*a4*al4*al5+t1*al6*d6-t1*al5*d6-2*t7*al5*al6+t2*al6*al5*a6+t0*al6*al5*d6+t3*al6*a6+t3*al5*a6+t0*d6+2*t6*al6-t2*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t1 * a4 * al4 + t3 * al6 * d5 + t3 * al6 * d4 - 2 * t5 - t0 * al6 * a4 * al4 - t2 * d5 - t2 * d4 - 2 * t4 * al5 - t3 * d5 * al5 - 2 * t5 * al5 * al6 + t1 * al5 * a6 - t1 * al6 * a6 - t1 * al5 * a5 + t1 * al6 * a5 - t3 * d4 * al5 + t0 * a4 * al4 * al5 + t0 * a5 - t0 * a6 + t3 * al6 * d6 + t3 * al5 * d6 + t0 * al6 * al5 * a5 - t0 * al6 * al5 * a6 - t2 * al6 * d5 * al5 - t2 * al6 * d4 * al5 + t2 * al6 * al5 * d6 + t1 * al6 * a4 * al4 * al5 - t2 * d6 + 2 * t4 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t6-t0*al6*d5+t1*d4+t1*d5+2*t7*al5-t0*al6*d4+t2*a4*al4-t3*al6*a4*al4-t3*a5-t1*al6*d5*al5-t0*d4*al5-t0*d5*al5-t2*al5*a5-t2*al6*a5-t1*al6*d4*al5-t3*a4*al4*al5+t3*al5*al6*a5-t2*al6*a4*al4*al5-t2*al6*a6+2*t6*al5*al6-t2*al5*a6+t3*al5*al6*a6+t1*al6*al5*d6-t0*al6*d6+t0*al5*d6-t3*a6+t1*d6+2*t7*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t6-t0*al6*d5+t1*d4+t1*d5-2*t7*al5-t0*al6*d4+t2*a4*al4-t3*al6*a4*al4+t3*a5+t1*al6*d5*al5+t0*d4*al5+t0*d5*al5-t2*al5*a5+t2*al6*a5+t1*al6*d4*al5+t3*a4*al4*al5+t3*al5*al6*a5+t2*al6*a4*al4*al5-t2*al6*a6-2*t6*al5*al6+t2*al5*a6-t3*al5*al6*a6-t1*al6*al5*d6-t0*al6*d6-t0*al5*d6-t3*a6+t1*d6+2*t7*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t1*a4*al4-t3*al6*d5-t3*al6*d4+2*t5+t0*al6*a4*al4+t2*d5+t2*d4-2*t4*al5-t3*d5*al5-2*t5*al5*al6+t1*al5*a6+t1*al6*a6+t1*al5*a5+t1*al6*a5-t3*d4*al5+t0*a4*al4*al5+t0*a5+t0*a6-t3*al6*d6+t3*al5*d6-t0*al6*al5*a5-t0*al6*al5*a6-t2*al6*d5*al5-t2*al6*d4*al5+t2*al6*al5*d6+t1*al6*a4*al4*al5+t2*d6-2*t4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t1*d5*al5-t0*al6*d5*al5+t3*al6*a5+2*t6*al5-t2*a5-t0*d4-t0*d5+t1*d4*al5-t1*al6*d4-t1*al6*d5+t3*a4*al4-t3*al5*a5-2*t7-t0*al6*d4*al5-t2*a4*al4*al5-t2*al6*al5*a5+t2*al6*a4*al4+t3*al6*a4*al4*al5-t1*al6*d6-t1*al5*d6-2*t7*al5*al6+t2*al6*al5*a6+t0*al6*al5*d6-t3*al6*a6+t3*al5*a6-t0*d6-2*t6*al6+t2*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t3*d5-t0*al6*a4*al4*al5+t3*d6+t0*a4*al4+2*t4*al5*al6-2*t5*al6-t0*al5*a6+t2*d5*al5-2*t5*al5+t2*d4*al5+t2*al6*d4+t2*al6*d5-t2*al5*d6-2*t4+t2*al6*d6-t0*al5*a5-t0*al6*a5-t0*al6*a6-t1*al6*al5*a5+t1*a4*al4*al5+t3*d4+t1*al6*a4*al4-t1*al6*al5*a6-t3*al6*d5*al5-t3*al6*d4*al5+t3*al5*al6*d6+t1*a5+t1*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t0 - 2 * t0 * al6 * al5 + 2 * t1 * al5 - 2 * t1 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al5+2*t2*al6+2*t3-2*t3*al5*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t1 - 2 * t1 * al6 * al5 - 2 * t0 * al5 + 2 * t0 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t3*al5+2*t3*al6-2*t2+2*t2*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t2 - 2 * t2 * al6 * al5 - 2 * t3 * al5 + 2 * t3 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5-2*t0*al6+2*t1-2*t1*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 - 2 * t3 * al5 * al6 + 2 * t2 * al5 - 2 * t2 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t1*al5-2*t1*al6-2*t0+2*t0*al6*al5);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t5*al4+t3*al6*al4*d5+t0*al4*a5-t1*a4-t0*al4*a6+t3*al4*al5*d6-2*t5*al4*al5*al6+t2*al4*d4+2*t4*al4*al6-t2*al4*d5-t0*al6*al4*al5*a6+t1*al6*al4*a5-t1*al6*al4*a6+t1*al4*al5*a6-t0*a4*al5-t1*al4*al5*a5-t1*al6*a4*al5+t3*al6*al4*d6-t3*al6*al4*d4-t3*al4*d5*al5+t3*al4*d4*al5+t0*al6*al4*al5*a5+t2*al4*al5*al6*d6-t2*al4*d6-2*t4*al4*al5-t2*al6*al4*d5*al5+t2*al6*al4*d4*al5+t0*al6*a4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t3*al6*a4-t2*a4-t0*al4*d5*al5+t0*al6*al4*d4-t0*al6*al4*d5-t2*al4*al5*a5+2*t7*al4*al5-t2*al6*al4*a5+t0*al4*d4*al5+t3*a4*al5-t3*al4*a5-t1*al4*d4+t1*al4*d5-t0*al4*al6*d6-t2*al4*al6*a6+t0*al4*al5*d6-t2*al4*al5*a6+2*t6*al4*al5*al6+t1*al4*d6+2*t7*al4*al6-t3*al4*a6-t1*al6*al4*d5*al5-2*t6*al4+t2*al6*a4*al5+t3*al6*al4*al5*a5+t1*al6*al4*d4*al5+t3*al4*al5*al6*a6+t1*al4*al5*al6*d6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t1*al6*a4-2*t5*al4*al5-t3*al4*d5+t3*al4*d4+t0*a4-t1*al4*a6+2*t5*al4*al6+2*t4*al4+t0*al6*a4*al5-t0*al6*al4*a5+t0*al6*al4*a6-t0*al4*al5*a6-t2*al4*d4*al5+t2*al4*d5*al5+t3*al6*al4*al5*d6-t2*al6*al4*d6-t2*al4*al5*d6+2*t4*al4*al5*al6-t3*al4*d6-t1*al6*al4*al5*a6-t2*al6*al4*d5+t2*al6*al4*d4+t3*al6*al4*d4*al5+t1*al6*al4*al5*a5+t0*al4*al5*a5-t1*a4*al5+t1*al4*a5-t3*al6*al4*d5*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al4-t2*al6*a4+t0*al4*d4-t0*al4*d5-t3*a4-t0*al6*al4*d4*al5-t3*al4*al5*a5-2*t6*al4*al5-t1*al6*al4*d5-t2*al4*al5*al6*a5-t3*al6*al4*a5+t3*al6*a4*al5-t0*al4*d6-t2*a4*al5+t2*al4*a5+t1*al4*d4*al5-t1*al4*d5*al5+t1*al6*al4*d4+t0*al6*al4*d5*al5+t1*al4*al5*d6-t1*al4*al6*d6-t3*al4*al5*a6-t3*al4*al6*a6+2*t7*al4*al5*al6-2*t6*al4*al6+t2*al4*a6-t0*al4*al5*al6*d6-t2*al4*al5*al6*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t7*al4+t2*al6*a4-t0*al4*d4+t0*al4*d5+t3*a4-t0*al6*al4*d4*al5+t3*al4*al5*a5-2*t6*al4*al5+t1*al6*al4*d5+t2*al4*al5*al6*a5-t3*al6*al4*a5+t3*al6*a4*al5+t0*al4*d6-t2*a4*al5+t2*al4*a5+t1*al4*d4*al5-t1*al4*d5*al5-t1*al6*al4*d4+t0*al6*al4*d5*al5+t1*al4*al5*d6+t1*al4*al6*d6-t3*al4*al5*a6+t3*al4*al6*a6+2*t7*al4*al5*al6+2*t6*al4*al6-t2*al4*a6-t0*al4*al5*al6*d6-t2*al4*al5*al6*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t1*al6*a4+2*t5*al4*al5-t3*al4*d5+t3*al4*d4+t0*a4-t1*al4*a6+2*t5*al4*al6+2*t4*al4-t0*al6*a4*al5+t0*al6*al4*a5+t0*al6*al4*a6+t0*al4*al5*a6+t2*al4*d4*al5-t2*al4*d5*al5-t3*al6*al4*al5*d6-t2*al6*al4*d6+t2*al4*al5*d6-2*t4*al4*al5*al6-t3*al4*d6+t1*al6*al4*al5*a6-t2*al6*al4*d5+t2*al6*al4*d4-t3*al6*al4*d4*al5+t1*al6*al4*al5*a5+t0*al4*al5*a5+t1*a4*al5-t1*al4*a5+t3*al6*al4*d5*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t3 * al6 * a4 - t2 * a4 + t0 * al4 * d5 * al5 + t0 * al6 * al4 * d4 - t0 * al6 * al4 * d5 - t2 * al4 * al5 * a5 - 2 * t7 * al4 * al5 + t2 * al6 * al4 * a5 - t0 * al4 * d4 * al5 - t3 * a4 * al5 + t3 * al4 * a5 - t1 * al4 * d4 + t1 * al4 * d5 - t0 * al4 * al6 * d6 - t2 * al4 * al6 * a6 - t0 * al4 * al5 * d6 + t2 * al4 * al5 * a6 - 2 * t6 * al4 * al5 * al6 + t1 * al4 * d6 + 2 * t7 * al4 * al6 - t3 * al4 * a6 + t1 * al6 * al4 * d5 * al5 - 2 * t6 * al4 - t2 * al6 * a4 * al5 + t3 * al6 * al4 * al5 * a5 - t1 * al6 * al4 * d4 * al5 - t3 * al4 * al5 * al6 * a6 - t1 * al4 * al5 * al6 * d6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t1*al6*al4*a5+t1*al4*al5*a5+t3*al4*al5*d6-t1*al6*a4*al5-2*t4*al4*al6+t0*al4*a6-t0*al6*a4+t0*al4*a5-t0*a4*al5+t2*al4*d6-t2*al4*d4+t2*al4*d5-2*t4*al4*al5+t1*al4*al5*a6+t1*al4*al6*a6+t1*a4+2*t5*al4-t2*al6*al4*d5*al5+t2*al6*al4*d4*al5-t0*al6*al4*al5*a5-t3*al4*d5*al5+t3*al4*d4*al5-2*t5*al4*al5*al6-t3*al4*al6*d6-t0*al4*al5*al6*a6+t2*al4*al5*al6*d6+t3*al6*al4*d4-t3*al6*al4*d5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t0 * al4 * al5 + 2 * t0 * al4 * al6 - 2 * t1 * al4 - 2 * t1 * al4 * al5 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al4+2*t2*al4*al5*al6+2*t3*al4*al5+2*t3*al4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t1 * al4 * al5 + 2 * t1 * al4 * al6 + 2 * t0 * al4 + 2 * t0 * al4 * al5 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al4+2*t3*al4*al5*al6-2*t2*al4*al5-2*t2*al4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t2*al4*al5+2*t2*al4*al6+2*t3*al4+2*t3*al4*al5*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0*al4-2*t0*al4*al5*al6+2*t1*al4*al5+2*t1*al4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al4*al5+2*t3*al4*al6-2*t2*al4-2*t2*al4*al5*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t1*al4-2*t1*al4*al5*al6-2*t0*al4*al5-2*t0*al4*al6);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t3 * al6 * a4 - t2 * a4 + t0 * al4 * d5 * al5 + t0 * al6 * al4 * d4 - t0 * al6 * al4 * d5 - t2 * al4 * al5 * a5 - 2 * t7 * al4 * al5 + t2 * al6 * al4 * a5 - t0 * al4 * d4 * al5 - t3 * a4 * al5 + t3 * al4 * a5 - t1 * al4 * d4 + t1 * al4 * d5 - t0 * al4 * al6 * d6 - t2 * al4 * al6 * a6 - t0 * al4 * al5 * d6 + t2 * al4 * al5 * a6 - 2 * t6 * al4 * al5 * al6 + t1 * al4 * d6 + 2 * t7 * al4 * al6 - t3 * al4 * a6 + t1 * al6 * al4 * d5 * al5 - 2 * t6 * al4 - t2 * al6 * a4 * al5 + t3 * al6 * al4 * al5 * a5 - t1 * al6 * al4 * d4 * al5 - t3 * al4 * al5 * al6 * a6 - t1 * al4 * al5 * al6 * d6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (t1*al6*al4*a5+t1*al4*al5*a5+t3*al4*al5*d6-t1*al6*a4*al5-2*t4*al4*al6+t0*al4*a6-t0*al6*a4+t0*al4*a5-t0*a4*al5+t2*al4*d6-t2*al4*d4+t2*al4*d5-2*t4*al4*al5+t1*al4*al5*a6+t1*al4*al6*a6+t1*a4+2*t5*al4-t2*al6*al4*d5*al5+t2*al6*al4*d4*al5-t0*al6*al4*al5*a5-t3*al4*d5*al5+t3*al4*d4*al5-2*t5*al4*al5*al6-t3*al4*al6*d6-t0*al4*al5*al6*a6+t2*al4*al5*al6*d6+t3*al6*al4*d4-t3*al6*al4*d5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al4-t2*al6*a4+t0*al4*d4-t0*al4*d5-t3*a4+t0*al6*al4*d4*al5-t3*al4*al5*a5+2*t6*al4*al5-t1*al6*al4*d5-t2*al4*al5*al6*a5+t3*al6*al4*a5-t3*al6*a4*al5-t0*al4*d6+t2*a4*al5-t2*al4*a5-t1*al4*d4*al5+t1*al4*d5*al5+t1*al6*al4*d4-t0*al6*al4*d5*al5-t1*al4*al5*d6-t1*al4*al6*d6+t3*al4*al5*a6-t3*al4*al6*a6-2*t7*al4*al5*al6-2*t6*al4*al6+t2*al4*a6+t0*al4*al5*al6*d6+t2*al4*al5*al6*a6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t0*al4*al5*a5+t2*al4*d5*al5-2*t4*al4+t3*al4*d5+t0*al6*a4*al5-t2*al4*al6*d4-t1*al4*al5*al6*a6+t3*al4*al5*al6*d6-t2*al4*al5*d6-t0*al4*al6*a6-t0*al4*al5*a6+t2*al4*al6*d6+2*t4*al4*al5*al6-2*t5*al4*al5-t1*al6*a4-t3*al6*al4*d5*al5-t1*a4*al5+t1*al4*a5-t0*al4*al6*a5+t2*al4*al6*d5-t3*al4*d4-t1*al4*al5*al6*a5-t2*al4*d4*al5+t3*al6*al4*d4*al5-t0*a4+t1*al4*a6+t3*al4*d6-2*t5*al4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (t0*al4*al5*a5+t2*al4*d5*al5+2*t4*al4-t3*al4*d5+t0*al6*a4*al5+t2*al4*al6*d4-t1*al4*al5*al6*a6+t3*al4*al5*al6*d6-t2*al4*al5*d6+t0*al4*al6*a6-t0*al4*al5*a6-t2*al4*al6*d6+2*t4*al4*al5*al6-2*t5*al4*al5+t1*al6*a4-t3*al6*al4*d5*al5-t1*a4*al5+t1*al4*a5-t0*al4*al6*a5-t2*al4*al6*d5+t3*al4*d4+t1*al4*al5*al6*a5-t2*al4*d4*al5+t3*al6*al4*d4*al5+t0*a4-t1*al4*a6-t3*al4*d6+2*t5*al4*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al4-t2*al6*a4+t0*al4*d4-t0*al4*d5-t3*a4-t0*al6*al4*d4*al5-t3*al4*al5*a5-2*t6*al4*al5-t1*al6*al4*d5-t2*al4*al5*al6*a5-t3*al6*al4*a5+t3*al6*a4*al5-t0*al4*d6-t2*a4*al5+t2*al4*a5+t1*al4*d4*al5-t1*al4*d5*al5+t1*al6*al4*d4+t0*al6*al4*d5*al5+t1*al4*al5*d6-t1*al4*al6*d6-t3*al4*al5*a6-t3*al4*al6*a6+2*t7*al4*al5*al6-2*t6*al4*al6+t2*al4*a6-t0*al4*al5*al6*d6-t2*al4*al5*al6*a6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t1*al6*al4*a5+t1*al4*al5*a5-t3*al4*al5*d6+t1*al6*a4*al5-2*t4*al4*al6+t0*al4*a6-t0*al6*a4-t0*al4*a5+t0*a4*al5+t2*al4*d6-t2*al4*d4+t2*al4*d5+2*t4*al4*al5-t1*al4*al5*a6+t1*al4*al6*a6+t1*a4+2*t5*al4+t2*al6*al4*d5*al5-t2*al6*al4*d4*al5-t0*al6*al4*al5*a5+t3*al4*d5*al5-t3*al4*d4*al5+2*t5*al4*al5*al6-t3*al4*al6*d6+t0*al4*al5*al6*a6-t2*al4*al5*al6*d6+t3*al6*al4*d4-t3*al6*al4*d5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-t3*al6*a4+t2*a4+t0*al4*d5*al5-t0*al6*al4*d4+t0*al6*al4*d5+t2*al4*al5*a5-2*t7*al4*al5+t2*al6*al4*a5-t0*al4*d4*al5-t3*a4*al5+t3*al4*a5+t1*al4*d4-t1*al4*d5+t0*al4*al6*d6+t2*al4*al6*a6-t0*al4*al5*d6+t2*al4*al5*a6-2*t6*al4*al5*al6-t1*al4*d6-2*t7*al4*al6+t3*al4*a6+t1*al6*al4*d5*al5+2*t6*al4-t2*al6*a4*al5-t3*al6*al4*al5*a5-t1*al6*al4*d4*al5-t3*al4*al5*al6*a6-t1*al4*al5*al6*d6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 * al4 * al5 + 2 * t3 * al4 * al6 - 2 * t2 * al4 - 2 * t2 * al4 * al5 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t1*al4-2*t1*al4*al5*al6-2*t0*al4*al5-2*t0*al4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t3 * al4 - 2 * t3 * al4 * al5 * al6 + 2 * t2 * al4 * al5 - 2 * t2 * al4 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t1*al4*al5-2*t1*al4*al6-2*t0*al4+2*t0*al4*al5*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t1 * al4 * al5 + 2 * t1 * al4 * al6 + 2 * t0 * al4 + 2 * t0 * al4 * al5 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al4+2*t3*al4*al5*al6-2*t2*al4*al5-2*t2*al4*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t1 * al4 + 2 * t1 * al4 * al5 * al6 + 2 * t0 * al4 * al5 - 2 * t0 * al4 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t3*al4*al5-2*t3*al4*al6+2*t2*al4-2*t2*al4*al5*al6);},p4,q4,p5,q5,p6,q6)));

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
    double p4=a.alp[3],q4=a.alq[3], p5=a.alp[4],q5=a.alq[4], p6=a.alp[5],q6=a.alq[5];
    double d4 = a.d[3], d5 = a.d[4], d6 = a.d[5];
    double a4 = a.a[3], a5 = a.a[4], a6 = a.a[5];

    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t7 - t2 * al6 * a4 * al4 + t0 * d5 - t2 * a6 + t2 * al6 * al5 * a5 + 2 * t7 * al5 * al6 - t1 * d5 * al5 - t1 * d4 * al5 + t1 * al6 * d4 + t1 * al6 * d5 - t3 * al6 * a5 + t1 * al6 * d6 + t1 * al5 * d6 - t0 * al5 * al6 * d6 + t3 * al5 * a5 - t3 * a4 * al4 + t0 * d4 + t0 * al6 * d5 * al5 - t3 * al5 * a6 + t2 * a5 + t0 * al6 * d4 * al5 + t3 * al6 * a6 - t2 * al6 * al5 * a6 - t3 * al6 * a4 * al4 * al5 + t2 * a4 * al4 * al5+t0*d6-2*t6*al5+2*t6*al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t4-t2*al6*d4-t2*al6*d5-t2*d4*al5-t1*a4*al4*al5-t0*a4*al4+t0*al5*a5+t0*al6*a5-t3*d5+2*t5*al5-t1*a5+t0*al6*a4*al4*al5+t3*al6*d4*al5+t3*al6*d5*al5-t3*d4-t3*al5*al6*d6+t1*al5*al6*a6-2*t4*al5*al6+t0*al6*a6+t0*al5*a6-t2*al6*d6+t2*al5*d6-t2*d5*al5-t1*a6-t3*d6+2*t5*al6+t1*al6*al5*a5-t1*al6*a4*al4);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al5+2*t7*al6-2*t6+t2*a4*al4+t1*d4+t1*d5-t0*al6*d5-t3*al6*a4*al4+t1*al6*d5*al5-t3*a6-t0*al6*d4+t3*a5+t2*al6*a5-t0*al6*d6-t0*al5*d6-2*t6*al5*al6+t3*a4*al4*al5+t3*al6*al5*a5+t0*d5*al5+t0*d4*al5-t2*al6*a6+t1*al6*d4*al5+t2*al5*a6-t3*al6*al5*a6-t1*al6*al5*d6-t2*al5*a5+t1*d6+t2*al6*a4*al4*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t5*al5*al6-t1*a4*al4+t3*al5*d6-t3*al6*d6+t1*al6*a6+t1*al5*a6-t3*d5*al5-t3*d4*al5-t3*al6*d5-t3*al6*d4+t1*al5*a5+t1*al6*a5-2*t4*al6-t2*al6*d5*al5+t2*d6-t2*al6*d4*al5+t0*a6+t0*al6*a4*al4-t0*al5*al6*a5+t0*a5-2*t4*al5-t0*al5*al6*a6+t0*a4*al4*al5+t2*d4+t2*al5*al6*d6+t2*d5+2*t5+t1*al6*a4*al4*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2*t5*al5*al6-t1*a4*al4-t3*al5*d6-t3*al6*d6+t1*al6*a6-t1*al5*a6+t3*d5*al5+t3*d4*al5-t3*al6*d5-t3*al6*d4+t1*al5*a5-t1*al6*a5-2*t4*al6+t2*al6*d5*al5+t2*d6+t2*al6*d4*al5+t0*a6+t0*al6*a4*al4-t0*al5*al6*a5-t0*a5+2*t4*al5+t0*al5*al6*a6-t0*a4*al4*al5+t2*d4-t2*al5*al6*d6+t2*d5+2*t5-t1*al6*a4*al4*al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t7*al5-2*t7*al6+2*t6-t2*a4*al4-t1*d4-t1*d5+t0*al6*d5+t3*al6*a4*al4+t1*al6*d5*al5+t3*a6+t0*al6*d4+t3*a5+t2*al6*a5+t0*al6*d6-t0*al5*d6-2*t6*al5*al6+t3*a4*al4*al5-t3*al6*al5*a5+t0*d5*al5+t0*d4*al5+t2*al6*a6+t1*al6*d4*al5+t2*al5*a6-t3*al6*al5*a6-t1*al6*al5*d6+t2*al5*a5-t1*d6+t2*al6*a4*al4*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-t0*al6*a6+t3*d4+t0*a4*al4+t3*d5+2*t5*al5+t2*al6*d4+t2*al6*d5-2*t4+t3*al6*d4*al5-t2*d5*al5-2*t4*al5*al6+t2*al6*d6+t0*al6*a5+t0*al6*a4*al4*al5-t1*al6*al5*a5+t3*d6-t1*a5-t0*al5*a5+t0*al5*a6+t3*al6*d5*al5+t1*a6-t2*d4*al5-t3*al6*al5*d6-t1*a4*al4*al5-2*t5*al6+t2*al5*d6+t1*al6*al5*a6+t1*al6*a4*al4);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t7-t2*al6*a4*al4+t0*d5-t2*a6+t2*al6*al5*a5-2*t7*al5*al6+t1*d5*al5+t1*d4*al5+t1*al6*d4+t1*al6*d5+t3*al6*a5+t1*al6*d6-t1*al5*d6+t0*al5*al6*d6+t3*al5*a5-t3*a4*al4+t0*d4-t0*al6*d5*al5+t3*al5*a6-t2*a5-t0*al6*d4*al5+t3*al6*a6+t2*al6*al5*a6+t3*al6*a4*al4*al5-t2*a4*al4*al5+t0*d6+2*t6*al5+2*t6*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t2 * al5 + 2 * t2 * al6 + 2 * t3 + 2 * t3 * al6 * al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t0-2*t0*al6*al5+2*t1*al5+2*t1*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t2 - 2 * t2 * al5 * al6 - 2 * t3 * al5 + 2 * t3 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (-2*t0*al5-2*t0*al6+2*t1-2*t1*al6*al5);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (2 * t0 * al5 - 2 * t0 * al6 + 2 * t1 + 2 * t1 * al6 * al5);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2-2*t2*al5*al6-2*t3*al5-2*t3*al6);},p4,q4,p5,q5,p6,q6)));
    result.push_back(Polynomial(proj_eval_3([&](double al4,double al5,double al6){return (-2 * t0 - 2 * t0 * al6 * al5 + 2 * t1 * al5 - 2 * t1 * al6);},p4,q4,p5,q5,p6,q6), proj_eval_3([&](double al4,double al5,double al6){return (2*t2*al5+2*t2*al6+2*t3-2*t3*al6*al5);},p4,q4,p5,q5,p6,q6)));

    return result;
  };
}