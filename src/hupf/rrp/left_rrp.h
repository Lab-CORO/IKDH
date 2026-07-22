#pragma once
#include <hupf/proj_eval.h>
#include <hupf/Input.h>

namespace LibHUPF
{
namespace RRP
{

// jcapco Hyperplanes Td3: => At the moment we change sign for translation because apparently somewhere in this C code pfurner's convention is involved
std::vector<Polynomial> h1_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l1*l2*l3+a1*l1-l2*l3*v3*d2-l2*a2-l2*a3-l3*a2-l3*a3-v3*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l3*v3-v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l1*l2+a1*l1*l3-l2*l3*a2-l2*l3*a3+l2*v3*d2-l3*v3*d2+a2+a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l2*v3-l3*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l1*l2*v3+a1*l1*v3*l3-l2*v3*l3*a2+l2*v3*l3*a3+l2*d2-v3*a2+v3*a3+l3*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l2+l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l1*l2*l3*v3+a1*l1*v3-l2*l3*d2-l2*v3*a2+l2*v3*a3+l3*v3*a2-l3*v3*a3+d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l3+1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3-2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2-2*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*v3-2*v3*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3*v3-2*v3;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

std::vector<Polynomial> h2_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2+a1*l3+l2*l3*a2*l1+l2*l3*l1*a3-l2*l1*v3*d2+l3*l1*v3*d2-a2*l1-l1*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l1*v3+l3*l1*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l3-a1-l2*l3*v3*l1*d2-l2*l1*a2-l2*l1*a3-l3*l1*a2-l3*l1*a3-v3*l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l3*v3*l1-v3*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l3*v3+a1*v3+l2*l3*l1*d2+l2*v3*l1*a2-l2*v3*l1*a3-l3*v3*l1*a2+l3*v3*l1*a3-l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l2*l3*l1-l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*v3-a1*v3*l3-l2*v3*l3*a2*l1+l2*v3*l3*l1*a3+l2*l1*d2-v3*a2*l1+v3*l1*a3+l3*l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l2*l1+l3*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l1+2*l1*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*l1-2*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*v3*l1+2*v3*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*v3*l1-2*v3*l1*l3;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

std::vector<Polynomial> h3_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*v3+a1*v3*l3+l2*v3*l3*a2*l1-l2*v3*l3*l1*a3-l2*l1*d2+v3*a2*l1-v3*l1*a3-l3*l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l1-l3*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*l3*v3-a1*v3-l2*l3*l1*d2-l2*v3*l1*a2+l2*v3*l1*a3+l3*v3*l1*a2-l3*v3*l1*a3+l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l3*l1+l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l3-a1-l2*l3*v3*l1*d2-l2*l1*a2-l2*l1*a3-l3*l1*a2-l3*l1*a3-v3*l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l3*v3*l1-v3*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2+a1*l3+l2*l3*a2*l1+l2*l3*l1*a3-l2*l1*v3*d2+l3*l1*v3*d2-a2*l1-l1*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l1*v3+l3*l1*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*v3*l1+2*v3*l1*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3*v3*l1-2*v3*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*l1-2*l1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l1+2*l1*l3;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

std::vector<Polynomial> h4_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double d2 = a.d[1], v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l1*l2*l3*v3-a1*l1*v3+l2*l3*d2+l2*v3*a2-l2*v3*a3-l3*v3*a2+l3*v3*a3-d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l2*l3-1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l1*l2*v3-a1*l1*v3*l3+l2*v3*l3*a2-l2*v3*l3*a3-l2*d2+v3*a2-v3*a3-l3*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2-l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l1*l2+a1*l1*l3-l2*l3*a2-l2*l3*a3+l2*v3*d2-l3*v3*d2+a2+a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l2*v3-l3*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l1*l2*l3+a1*l1-l2*l3*v3*d2-l2*a2-l2*a3-l3*a2-l3*a3-v3*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l2*l3*v3-v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*v3+2*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*v3+2*v3*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2-2*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3-2;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

//jcapco todo above: investigate where sign change becomes relevant in the C code!

// jcapco Hyperplanes Tv1: 

std::vector<Polynomial> h1_tc_v1 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double v3 = a.v[2];
    // double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*l3+l1+l2+l3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*l3*v3-l1*v3-l2*v3+l3*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l1*l2+l1*l3+l2*l3-1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l2*v3-l1*v3*l3-l2*v3*l3-v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*v3+l1*v3*l3+l2*v3*l3+v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l2+l1*l3+l2*l3-1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l1*l2*l3*v3+l1*v3+l2*v3-l3*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*l3+l1+l2+l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h2_tc_v1 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*l3*v3-l1*v3+l2*v3-l3*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l2*l3-l1+l2+l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l1*l2*v3-l1*v3*l3+l2*v3*l3+v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l2-l1*l3+l2*l3-1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l1*l2+l1*l3-l2*l3+1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l2*v3-l1*v3*l3+l2*v3*l3+v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*l3+l1-l2-l3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l2*l3*v3-l1*v3+l2*v3-l3*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h3_tc_v1 (Input& a)
{
  std::vector<Polynomial> result;
  double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
  //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
  double d2 = a.d[1], v3 = a.v[2];
  double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

  if (abs(a1)<1E-3) //rref condition
  {
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l3+a2+l2*l3*a3+l2*v3*d2+3*l3*v3*d2+a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l3*v3-a2*v3-l2*l3*v3*a3+l2*d2-3*l3*d2+v3*a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a2*l2+a2*l3+l2*v3*l3*d2-l2*a3-3*v3*d2+l3*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a2*l2*v3-a2*v3*l3+l2*v3*a3+l2*l3*d2+v3*l3*a3+3*d2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a2*l2*v3+a2*v3*l3-l2*v3*a3-l2*l3*d2-v3*l3*a3-3*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a2*l2+a2*l3+l2*v3*l3*d2-l2*a3-3*v3*d2+l3*a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a2*l2*l3*v3+a2*v3+l2*l3*v3*a3-l2*d2+3*l3*d2-v3*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l3+a2+l2*l3*a3+l2*v3*d2+3*l3*v3*d2+a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2+2*l3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*v3+2*v3*l3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3-2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*v3-2*v3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3*v3+2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3-2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*v3-2*v3*l3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l2+2*l3;},p1,q1,p2,q2,p3,q3)));
  }
  else
  {
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l3+a1+l2*l2*l3*l1*a2-l2*l2*a2-l2*l3*v3*l1*d2+l2*l3*a3+l2*v3*d2+l2*l1*a3+l3*v3*d2-l3*l1*a2-l3*l1*a3+v3*l1*d2+a2+a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l3*v3-a1*v3+l2*l2*l3*v3*l1*a2+l2*l2*v3*a2-l2*l3*v3*a3+l2*l3*l1*d2+l2*v3*l1*a3+l2*d2-l3*v3*l1*a2+l3*v3*l1*a3-l3*d2-v3*a2+v3*a3+l1*d2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l2+a1*l3-l2*l2*l3*a2-l2*l2*a2*l1+l2*v3*l3*d2+l2*v3*l1*d2+l2*l3*l1*a3-l2*a3+v3*l3*l1*d2-v3*d2+l3*a2+l3*a3+a2*l1+l1*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*v3-a1*v3*l3+l2*l2*v3*l3*a2-l2*l2*v3*a2*l1+l2*v3*l3*l1*a3+l2*v3*a3+l2*l3*d2-l2*l1*d2-v3*l3*a2+v3*l3*a3+v3*a2*l1-v3*l1*a3+l3*l1*d2+d2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*v3+a1*v3*l3-l2*l2*v3*l3*a2+l2*l2*v3*a2*l1-l2*v3*l3*l1*a3-l2*v3*a3-l2*l3*d2+l2*l1*d2+v3*l3*a2-v3*l3*a3-v3*a2*l1+v3*l1*a3-l3*l1*d2-d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a1*l2+a1*l3-l2*l2*l3*a2-l2*l2*a2*l1+l2*v3*l3*d2+l2*v3*l1*d2+l2*l3*l1*a3-l2*a3+v3*l3*l1*d2-v3*d2+l3*a2+l3*a3+a2*l1+l1*a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*l3*v3+a1*v3-l2*l2*l3*v3*l1*a2-l2*l2*v3*a2+l2*l3*v3*a3-l2*l3*l1*d2-l2*v3*l1*a3-l2*d2+l3*v3*l1*a2-l3*v3*l1*a3+l3*d2+v3*a2-v3*a3-l1*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l3+a1+l2*l2*l3*l1*a2-l2*l2*a2-l2*l3*v3*l1*d2+l2*l3*a3+l2*v3*d2+l2*l1*a3+l3*v3*d2-l3*l1*a2-l3*l1*a3+v3*l1*d2+a2+a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*l1-2*l2+2*l3+2*l1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*v3*l1+2*l2*v3+2*l3*v3-2*v3*l1;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3-2*l2*l1+2*l3*l1-2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*v3-2*l2*v3*l1-2*l3*v3*l1-2*v3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3*v3+2*l2*v3*l1+2*l3*v3*l1+2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3-2*l2*l1+2*l3*l1-2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l3*l1*v3-2*l2*v3-2*l3*v3+2*l1*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l3*l1-2*l2+2*l3+2*l1;},p1,q1,p2,q2,p3,q3)));
  }
  return result;
 };

std::vector<Polynomial> h4_tc_v1 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

  if (abs(a1)<1E-3) //rref condition
  {
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l2*l3*v3-a2*l3*v3-l2*l2*l3*v3*a3+l2*l2*d2-3*l2*l3*d2+l2*v3*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a2*l2*l2*l3+a2*l3-l2*l2*l3*a3-l2*l2*v3*d2-3*l2*l3*v3*d2-l2*a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a2*l2*l2*v3+a2*v3+l2*l2*v3*a3+l2*l2*l3*d2+l2*v3*l3*a3+3*l2*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l2-a2-l2*l2*v3*l3*d2+l2*l2*a3+3*l2*v3*d2-l2*l3*a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a2*l2*l2+a2+l2*l2*v3*l3*d2-l2*l2*a3-3*l2*v3*d2+l2*l3*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a2*l2*l2*v3+a2*v3+l2*l2*v3*a3+l2*l2*l3*d2+l2*v3*l3*a3+3*l2*d2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l2*l3-a2*l3+l2*l2*l3*a3+l2*l2*v3*d2+3*l2*l3*v3*d2+l2*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a2*l2*l2*l3*v3-a2*l3*v3-l2*l2*l3*v3*a3+l2*l2*d2-3*l2*l3*d2+l2*v3*a3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*v3+2*l2*v3*l3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2-2*l2*l3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3*v3-2*l2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3+2*l2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l2*l3-2*l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3*v3-2*l2*v3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l2+2*l2*l3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*v3+2*l2*v3*l3;},p1,q1,p2,q2,p3,q3)));
  }
  else
  {
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*l2*l3*v3+a1*l2*v3+l2*l2*l3*v3*a2-l2*l2*l3*v3*a3-l2*l2*l3*l1*d2+l2*l2*v3*l1*a2-l2*l2*v3*l1*a3+l2*l2*d2-l2*l3*v3*l1*a3-l2*l3*d2+l2*v3*a3-l2*l1*d2-l3*v3*a2-v3*l1*a2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l2*l3+a1*l2-l2*l2*l3*v3*l1*d2-l2*l2*l3*a2-l2*l2*l3*a3-l2*l2*v3*d2+l2*l2*l1*a2+l2*l2*l1*a3-l2*l3*v3*d2-l2*l3*l1*a3+l2*v3*l1*d2-l2*a3+l3*a2-l1*a2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l2*v3+a1*l2*v3*l3+l2*l2*v3*l3*l1*a2-l2*l2*v3*l3*l1*a3-l2*l2*v3*a2+l2*l2*v3*a3+l2*l2*l3*d2+l2*l2*d2*l1+l2*v3*l3*a3+l2*v3*l1*a3-l2*l3*d2*l1+l2*d2-v3*l3*l1*a2+v3*a2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*l2+a1*l2*l3-l2*l2*v3*l3*d2+l2*l2*v3*d2*l1+l2*l2*l3*l1*a2+l2*l2*l3*l1*a3+l2*l2*a2+l2*l2*a3+l2*v3*l3*d2*l1+l2*v3*d2-l2*l3*a3+l2*l1*a3-l3*l1*a2-a2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l2-a1*l2*l3+l2*l2*v3*l3*d2-l2*l2*v3*l1*d2-l2*l2*l3*l1*a2-l2*l2*l3*l1*a3-l2*l2*a2-l2*l2*a3-l2*v3*l3*l1*d2-l2*v3*d2+l2*l3*a3-l2*l1*a3+l3*l1*a2+a2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return a1*l2*l2*v3+a1*l2*v3*l3+l2*l2*v3*l3*l1*a2-l2*l2*v3*l3*l1*a3-l2*l2*v3*a2+l2*l2*v3*a3+l2*l2*l3*d2+l2*l2*l1*d2+l2*v3*l3*a3+l2*v3*l1*a3-l2*l3*l1*d2+l2*d2-v3*l3*l1*a2+v3*a2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*l2*l3-a1*l2+l2*l2*l3*v3*l1*d2+l2*l2*l3*a2+l2*l2*l3*a3+l2*l2*v3*d2-l2*l2*a2*l1-l2*l2*a3*l1+l2*l3*v3*d2+l2*l3*a3*l1-l2*v3*l1*d2+l2*a3-l3*a2+a2*l1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -a1*l2*l2*l3*v3+a1*l2*v3+l2*l2*l3*v3*a2-l2*l2*l3*v3*a3-l2*l2*l3*l1*d2+l2*l2*v3*a2*l1-l2*l2*v3*a3*l1+l2*l2*d2-l2*l3*v3*a3*l1-l2*l3*d2+l2*v3*a3-l2*l1*d2-l3*v3*a2-v3*a2*l1;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l2*l3*l1*v3+2*l2*l2*v3+2*l2*l3*v3+2*l2*l1*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3*l1+2*l2*l2-2*l2*l3+2*l2*l1;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3*v3+2*l2*l2*v3*l1+2*l2*l3*v3*l1-2*l2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3-2*l2*l2*l1+2*l2*l3*l1+2*l2;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l2*l3+2*l2*l2*l1-2*l2*l3*l1-2*l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l2*l2*l3*v3+2*l2*l2*v3*l1+2*l2*l3*v3*l1-2*l2*v3;},p1,q1,p2,q2,p3,q3)));
    result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l2*l3*l1-2*l2*l2+2*l2*l3-2*l2*l1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l2*l2*l3*v3*l1+2*l2*l2*v3+2*l2*l3*v3+2*l2*v3*l1;},p1,q1,p2,q2,p3,q3)));
  }
  return result;
};

// jcapco Hyperplanes Tv2: 

std::vector<Polynomial> h1_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*l2*v3-l1*v3+l3*v3-l2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l3*l2-l1+l3+l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*v3+l1*v3*l2-l3*v3*l2-v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l3-l1*l2+l3*l2-1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l1*l3+l1*l2+l3*l2-1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*v3+l1*v3*l2+l3*v3*l2+v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*l2+l1+l3+l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*v3*l2-l1*v3-l3*v3+v3*l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h2_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    double v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l1*l3*l2-l1-l3-l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l3*v3*l2+l1*v3+l3*v3-v3*l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l3-l1*l2-l3*l2+1;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l3*v3-l1*v3*l2-l3*v3*l2-v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*v3+l1*v3*l2-l3*v3*l2-v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -l1*l3-l1*l2+l3*l2-1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l1*l3*l2*v3-l1*v3+l3*v3-l2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return l1*l3*l2-l1+l3+l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h3_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l3*v3-4*a2*l1*l3*v3*l2-l1*l1*l1*l3*v3*l2*a3-2*l1*l1*l1*l3*d2-l1*l1*l1*v3*a3-l1*l1*l3*v3*a3-2*l1*l1*l3*l2*d2+l1*l1*v3*l2*a3+4*l1*l1*d2+l1*l3*v3*l2*a3+2*l1*l3*d2+l1*v3*a3+4*l1*l2*d2+l3*v3*a3+2*l3*l2*d2-v3*l2*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l2-4*a2*l1-2*l1*l1*l1*l3*v3*d2-l1*l1*l1*l3*l2*a3+l1*l1*l1*a3-2*l1*l1*l3*v3*l2*d2+l1*l1*l3*a3+l1*l1*l2*a3-2*l1*l3*v3*d2+l1*l3*l2*a3-l1*a3-2*l3*v3*l2*d2-l3*a3-l2*a3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -4*a2*l1*l1*v3+4*a2*l1*l2*v3-l1*l1*l1*l3*v3*a3+l1*l1*l1*l2*v3*a3+2*l1*l1*l1*d2+l1*l1*l3*l2*v3*a3+4*l1*l1*l3*d2+2*l1*l1*l2*d2+l1*l1*v3*a3+4*l1*l3*l2*d2+l1*l3*v3*a3-l1*l2*v3*a3-2*l1*d2-l3*l2*v3*a3-2*l2*d2-v3*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l3*l2-4*a2*l1*l3+l1*l1*l1*l3*a3+l1*l1*l1*l2*a3+2*l1*l1*l1*v3*d2+l1*l1*l3*l2*a3+2*l1*l1*l2*v3*d2-l1*l1*a3-l1*l3*a3-l1*l2*a3+2*l1*v3*d2-l3*l2*a3+2*l2*v3*d2+a3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -4*a2*l1*l1+4*a2*l1*l2-l1*l1*l1*l3*a3-2*l1*l1*l1*v3*d2-l1*l1*l1*l2*a3+4*l1*l1*l3*v3*d2+l1*l1*l3*l2*a3-2*l1*l1*v3*l2*d2-l1*l1*a3+4*l1*l3*v3*l2*d2+l1*l3*a3+2*l1*v3*d2+l1*l2*a3-l3*l2*a3+2*v3*l2*d2+a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l3*v3*l2-4*a2*l1*l3*v3-l1*l1*l1*l3*v3*a3+l1*l1*l1*v3*l2*a3+2*l1*l1*l1*d2-l1*l1*l3*v3*l2*a3-l1*l1*v3*a3+2*l1*l1*l2*d2+l1*l3*v3*a3-l1*v3*l2*a3+2*l1*d2+l3*v3*l2*a3+v3*a3+2*l2*d2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l3-4*a2*l1*l3*l2+2*l1*l1*l1*l3*v3*d2+l1*l1*l1*l3*l2*a3-l1*l1*l1*a3+2*l1*l1*l3*v3*l2*d2+l1*l1*l3*a3+4*l1*l1*v3*d2+l1*l1*l2*a3-2*l1*l3*v3*d2-l1*l3*l2*a3+4*l1*v3*l2*d2+l1*a3-2*l3*v3*l2*d2-l3*a3-l2*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*v3*l2-4*a2*l1*v3-l1*l1*l1*l3*v3*l2*a3-2*l1*l1*l1*l3*d2-l1*l1*l1*v3*a3+l1*l1*l3*v3*a3-2*l1*l1*l3*l2*d2-l1*l1*v3*l2*a3+l1*l3*v3*l2*a3-2*l1*l3*d2+l1*v3*a3-l3*v3*a3-2*l3*l2*d2+v3*l2*a3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*v3+2*l1*l1*l1*v3*l2+2*l1*l1*l3*v3*l2+2*l1*l1*v3+2*l1*l3*v3-2*l1*v3*l2-2*l3*v3*l2-2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l1*l1*l1*l3+2*l1*l1*l1*l2+2*l1*l1*l3*l2-2*l1*l1-2*l1*l3-2*l1*l2-2*l3*l2+2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l1*l1*l1*l3*v3*l2+2*l1*l1*l1*v3+2*l1*l1*l3*v3-2*l1*l1*v3*l2-2*l1*l3*v3*l2-2*l1*v3-2*l3*v3+2*v3*l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l1*l1*l1*l3*l2-2*l1*l1*l1-2*l1*l1*l3-2*l1*l1*l2-2*l1*l3*l2+2*l1+2*l3+2*l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*l2+2*l1*l1*l1-2*l1*l1*l3-2*l1*l1*l2+2*l1*l3*l2-2*l1+2*l3+2*l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l1*l1*l1*l3*v3*l2+2*l1*l1*l1*v3-2*l1*l1*l3*v3+2*l1*l1*v3*l2-2*l1*l3*v3*l2-2*l1*v3+2*l3*v3-2*v3*l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3-2*l1*l1*l1*l2+2*l1*l1*l3*l2-2*l1*l1+2*l1*l3+2*l1*l2-2*l3*l2+2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*v3+2*l1*l1*l1*v3*l2-2*l1*l1*l3*v3*l2-2*l1*l1*v3+2*l1*l3*v3-2*l1*v3*l2+2*l3*v3*l2+2*v3;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

std::vector<Polynomial> h4_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l2+4*a2*l1+l1*l1*l1*l3*l2*a3-2*l1*l1*l1*v3*l2*d2-l1*l1*l1*a3+4*l1*l1*l3*v3*l2*d2+l1*l1*l3*a3+2*l1*l1*v3*d2+l1*l1*l2*a3-4*l1*l3*v3*d2-l1*l3*l2*a3+2*l1*v3*l2*d2+l1*a3-l3*a3-2*v3*d2-l2*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -4*a2*l1*l1*l3*v3-4*a2*l1*l3*v3*l2-l1*l1*l1*l3*v3*l2*a3-l1*l1*l1*v3*a3-2*l1*l1*l1*l2*d2+l1*l1*l3*v3*a3-l1*l1*v3*l2*a3+2*l1*l1*d2+l1*l3*v3*l2*a3+l1*v3*a3-2*l1*l2*d2-l3*v3*a3+v3*l2*a3+2*d2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l3*l2+4*a2*l1*l3-2*l1*l1*l1*l3*l2*v3*d2-l1*l1*l1*l3*a3-l1*l1*l1*l2*a3+l1*l1*l3*l2*a3+2*l1*l1*l3*v3*d2-4*l1*l1*l2*v3*d2-l1*l1*a3+2*l1*l3*l2*v3*d2+l1*l3*a3+l1*l2*a3+4*l1*v3*d2-l3*l2*a3-2*l3*v3*d2+a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*v3+4*a2*l1*l2*v3-2*l1*l1*l1*l3*l2*d2-l1*l1*l1*l3*v3*a3+l1*l1*l1*l2*v3*a3-l1*l1*l3*l2*v3*a3+2*l1*l1*l3*d2-l1*l1*v3*a3-2*l1*l3*l2*d2+l1*l3*v3*a3-l1*l2*v3*a3+l3*l2*v3*a3+2*l3*d2+v3*a3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*l3*v3*l2+4*a2*l1*l3*v3+l1*l1*l1*l3*v3*a3+2*l1*l1*l1*l3*l2*d2-l1*l1*l1*v3*l2*a3-l1*l1*l3*v3*l2*a3-2*l1*l1*l3*d2-l1*l1*v3*a3-4*l1*l1*l2*d2-l1*l3*v3*a3-2*l1*l3*l2*d2+l1*v3*l2*a3+4*l1*d2+l3*v3*l2*a3+2*l3*d2+v3*a3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1+4*a2*l1*l2-2*l1*l1*l1*l3*v3*l2*d2-l1*l1*l1*l3*a3-l1*l1*l1*l2*a3+2*l1*l1*l3*v3*d2-l1*l1*l3*l2*a3+l1*l1*a3-2*l1*l3*v3*l2*d2+l1*l3*a3+l1*l2*a3+2*l3*v3*d2+l3*l2*a3-a3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 4*a2*l1*l1*v3*l2+4*a2*l1*v3+l1*l1*l1*l3*v3*l2*a3+l1*l1*l1*v3*a3+2*l1*l1*l1*l2*d2+l1*l1*l3*v3*a3+4*l1*l1*l3*l2*d2-l1*l1*v3*l2*a3-2*l1*l1*d2-l1*l3*v3*l2*a3-4*l1*l3*d2-l1*v3*a3-2*l1*l2*d2-l3*v3*a3+v3*l2*a3+2*d2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -4*a2*l1*l1*l3-4*a2*l1*l3*l2+l1*l1*l1*l3*l2*a3-2*l1*l1*l1*v3*l2*d2-l1*l1*l1*a3-l1*l1*l3*a3+2*l1*l1*v3*d2-l1*l1*l2*a3-l1*l3*l2*a3-2*l1*v3*l2*d2+l1*a3+l3*a3+2*v3*d2+l2*a3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3-2*l1*l1*l1*l2+2*l1*l1*l3*l2-2*l1*l1+2*l1*l3+2*l1*l2-2*l3*l2+2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*v3+2*l1*l1*l1*v3*l2-2*l1*l1*l3*v3*l2-2*l1*l1*v3+2*l1*l3*v3-2*l1*v3*l2+2*l3*v3*l2+2*v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*l2+2*l1*l1*l1-2*l1*l1*l3-2*l1*l1*l2+2*l1*l3*l2-2*l1+2*l3+2*l2;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return 2*l1*l1*l1*l3*v3*l2+2*l1*l1*l1*v3-2*l1*l1*l3*v3+2*l1*l1*v3*l2-2*l1*l3*v3*l2-2*l1*v3+2*l3*v3-2*v3*l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*l2*v3-2*l1*l1*l1*v3-2*l1*l1*l3*v3+2*l1*l1*l2*v3+2*l1*l3*l2*v3+2*l1*v3+2*l3*v3-2*l2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3*l2+2*l1*l1*l1+2*l1*l1*l3+2*l1*l1*l2+2*l1*l3*l2-2*l1-2*l3-2*l2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return 2*l1*l1*l1*l3*v3-2*l1*l1*l1*v3*l2-2*l1*l1*l3*v3*l2-2*l1*l1*v3-2*l1*l3*v3+2*l1*v3*l2+2*l3*v3*l2+2*v3;},p1,q1,p2,q2,p3,q3), proj_eval_3([&](double l1,double l2,double l3){return -2*l1*l1*l1*l3-2*l1*l1*l1*l2-2*l1*l1*l3*l2+2*l1*l1+2*l1*l3+2*l1*l2+2*l3*l2-2;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

// jcapco Hyperplanes RRP special case (planar): 

std::vector<Polynomial> h1_tsp (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2-l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*l3+1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*l3*v3-v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*v3+v3*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h2_tsp (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l2*v3-v3*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l2*l3*v3+v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*l3+1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2-l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h3_tsp (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (-a2*l2-a2*l3-l2*a3-l3*a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (-a2*l2*l3+a2-l2*l3*a3+a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (-a2*l2*l3*v3-a2*v3+l2*l3*v3*a3+v3*a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (-a2*l2*v3+a2*v3*l3+l2*v3*a3-v3*l3*a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*l3+1;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l2+l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*v3+v3*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l2*l3*v3+v3;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

std::vector<Polynomial> h4_tsp (Input& a)
{
    std::vector<Polynomial> result;
    double p1=a.alp[0],q1=a.alq[0], p2=a.alp[1],q2=a.alq[1], p3=a.alp[2],q3=a.alq[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (a2*l2*v3-a2*v3*l3-l2*v3*a3+v3*l3*a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (a2*l2*l3*v3+a2*v3-l2*l3*v3*a3-v3*a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (-a2*l2*l3+a2-l2*l3*a3+a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return (-a2*l2-a2*l3-l2*a3-l3*a3)/2;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*l3*v3-v3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l2*v3-v3*l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return l2+l3;},p1,q1,p2,q2,p3,q3)));
result.push_back(Polynomial(proj_eval_3([&](double l1,double l2,double l3){return -l2*l3+1;},p1,q1,p2,q2,p3,q3)));

 return result;
 };

}
}
