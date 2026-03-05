#pragma once
#include <hupf/Input.h>

namespace LibHUPF
{
namespace RRP
{

// jcapco Hyperplanes Td3: => At the moment we change sign for translation because apparently somewhere in this C code pfurner's convention is involved
std::vector<Polynomial> h1_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-a1*l1*l2*l3+a1*l1-l2*l3*v3*d2-l2*a2-l2*a3-l3*a2-l3*a3-v3*d2,l2*l3*v3-v3));
result.push_back(Polynomial(a1*l1*l2+a1*l1*l3-l2*l3*a2-l2*l3*a3+l2*v3*d2-l3*v3*d2+a2+a3,-l2*v3-l3*v3));
result.push_back(Polynomial(-a1*l1*l2*v3+a1*l1*v3*l3-l2*v3*l3*a2+l2*v3*l3*a3+l2*d2-v3*a2+v3*a3+l3*d2,-l2+l3));
result.push_back(Polynomial(a1*l1*l2*l3*v3+a1*l1*v3-l2*l3*d2-l2*v3*a2+l2*v3*a3+l3*v3*a2-l3*v3*a3+d2,l2*l3+1));
result.push_back(Polynomial(2*l2*l3-2));
result.push_back(Polynomial(-2*l2-2*l3));
result.push_back(Polynomial(2*l2*v3-2*v3*l3));
result.push_back(Polynomial(-2*l2*l3*v3-2*v3));

 return result;
 };

std::vector<Polynomial> h2_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(a1*l2+a1*l3+l2*l3*a2*l1+l2*l3*l1*a3-l2*l1*v3*d2+l3*l1*v3*d2-a2*l1-l1*a3,l2*l1*v3+l3*l1*v3));
result.push_back(Polynomial(a1*l2*l3-a1-l2*l3*v3*l1*d2-l2*l1*a2-l2*l1*a3-l3*l1*a2-l3*l1*a3-v3*l1*d2,l2*l3*v3*l1-v3*l1));
result.push_back(Polynomial(a1*l2*l3*v3+a1*v3+l2*l3*l1*d2+l2*v3*l1*a2-l2*v3*l1*a3-l3*v3*l1*a2+l3*v3*l1*a3-l1*d2,-l2*l3*l1-l1));
result.push_back(Polynomial(a1*l2*v3-a1*v3*l3-l2*v3*l3*a2*l1+l2*v3*l3*l1*a3+l2*l1*d2-v3*a2*l1+v3*l1*a3+l3*l1*d2,-l2*l1+l3*l1));
result.push_back(Polynomial(2*l2*l1+2*l1*l3));
result.push_back(Polynomial(2*l2*l3*l1-2*l1));
result.push_back(Polynomial(2*l2*l3*v3*l1+2*v3*l1));
result.push_back(Polynomial(2*l2*v3*l1-2*v3*l1*l3));

 return result;
 };

std::vector<Polynomial> h3_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-a1*l2*v3+a1*v3*l3+l2*v3*l3*a2*l1-l2*v3*l3*l1*a3-l2*l1*d2+v3*a2*l1-v3*l1*a3-l3*l1*d2,l2*l1-l3*l1));
result.push_back(Polynomial(-a1*l2*l3*v3-a1*v3-l2*l3*l1*d2-l2*v3*l1*a2+l2*v3*l1*a3+l3*v3*l1*a2-l3*v3*l1*a3+l1*d2,l2*l3*l1+l1));
result.push_back(Polynomial(a1*l2*l3-a1-l2*l3*v3*l1*d2-l2*l1*a2-l2*l1*a3-l3*l1*a2-l3*l1*a3-v3*l1*d2,l2*l3*v3*l1-v3*l1));
result.push_back(Polynomial(a1*l2+a1*l3+l2*l3*a2*l1+l2*l3*l1*a3-l2*l1*v3*d2+l3*l1*v3*d2-a2*l1-l1*a3,l2*l1*v3+l3*l1*v3));
result.push_back(Polynomial(-2*l2*v3*l1+2*v3*l1*l3));
result.push_back(Polynomial(-2*l2*l3*v3*l1-2*v3*l1));
result.push_back(Polynomial(2*l2*l3*l1-2*l1));
result.push_back(Polynomial(2*l2*l1+2*l1*l3));

 return result;
 };

std::vector<Polynomial> h4_tc_d3 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double d2 = a.d[1], v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-a1*l1*l2*l3*v3-a1*l1*v3+l2*l3*d2+l2*v3*a2-l2*v3*a3-l3*v3*a2+l3*v3*a3-d2,-l2*l3-1));
result.push_back(Polynomial(a1*l1*l2*v3-a1*l1*v3*l3+l2*v3*l3*a2-l2*v3*l3*a3-l2*d2+v3*a2-v3*a3-l3*d2,l2-l3));
result.push_back(Polynomial(a1*l1*l2+a1*l1*l3-l2*l3*a2-l2*l3*a3+l2*v3*d2-l3*v3*d2+a2+a3,-l2*v3-l3*v3));
result.push_back(Polynomial(-a1*l1*l2*l3+a1*l1-l2*l3*v3*d2-l2*a2-l2*a3-l3*a2-l3*a3-v3*d2,l2*l3*v3-v3));
result.push_back(Polynomial(2*l2*l3*v3+2*v3));
result.push_back(Polynomial(-2*l2*v3+2*v3*l3));
result.push_back(Polynomial(-2*l2-2*l3));
result.push_back(Polynomial(2*l2*l3-2));

 return result;
 };

//jcapco todo above: investigate where sign change becomes relevant in the C code!

// jcapco Hyperplanes Tv1: 

std::vector<Polynomial> h1_tc_v1 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double v3 = a.v[2];
    // double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-l1*l2*l3+l1+l2+l3,-l1*l2*l3*v3-l1*v3-l2*v3+l3*v3));
result.push_back(Polynomial(l1*l2+l1*l3+l2*l3-1,l1*l2*v3-l1*v3*l3-l2*v3*l3-v3));
result.push_back(Polynomial(-l1*l2*v3+l1*v3*l3+l2*v3*l3+v3,l1*l2+l1*l3+l2*l3-1));
result.push_back(Polynomial(l1*l2*l3*v3+l1*v3+l2*v3-l3*v3,-l1*l2*l3+l1+l2+l3));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h2_tc_v1 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-l1*l2*l3*v3-l1*v3+l2*v3-l3*v3,l1*l2*l3-l1+l2+l3));
result.push_back(Polynomial(l1*l2*v3-l1*v3*l3+l2*v3*l3+v3,-l1*l2-l1*l3+l2*l3-1));
result.push_back(Polynomial(l1*l2+l1*l3-l2*l3+1,l1*l2*v3-l1*v3*l3+l2*v3*l3+v3));
result.push_back(Polynomial(-l1*l2*l3+l1-l2-l3,-l1*l2*l3*v3-l1*v3+l2*v3-l3*v3));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h3_tc_v1 (Input& a)
{
  std::vector<Polynomial> result;
  double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
  //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
  double d2 = a.d[1], v3 = a.v[2];
  double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

  if (abs(a1)<1E-3) //rref condition
  {
    result.push_back(Polynomial(a2*l2*l3+a2+l2*l3*a3+l2*v3*d2+3*l3*v3*d2+a3,a2*l2*l3*v3-a2*v3-l2*l3*v3*a3+l2*d2-3*l3*d2+v3*a3));
    result.push_back(Polynomial(-a2*l2+a2*l3+l2*v3*l3*d2-l2*a3-3*v3*d2+l3*a3,-a2*l2*v3-a2*v3*l3+l2*v3*a3+l2*l3*d2+v3*l3*a3+3*d2));
    result.push_back(Polynomial(a2*l2*v3+a2*v3*l3-l2*v3*a3-l2*l3*d2-v3*l3*a3-3*d2,-a2*l2+a2*l3+l2*v3*l3*d2-l2*a3-3*v3*d2+l3*a3));
    result.push_back(Polynomial(-a2*l2*l3*v3+a2*v3+l2*l3*v3*a3-l2*d2+3*l3*d2-v3*a3,a2*l2*l3+a2+l2*l3*a3+l2*v3*d2+3*l3*v3*d2+a3));
    result.push_back(Polynomial(-2*l2+2*l3,2*l2*v3+2*v3*l3));
    result.push_back(Polynomial(-2*l2*l3-2,2*l2*l3*v3-2*v3));
    result.push_back(Polynomial(-2*l2*l3*v3+2*v3,-2*l2*l3-2));
    result.push_back(Polynomial(-2*l2*v3-2*v3*l3,-2*l2+2*l3));
  }
  else
  {
    result.push_back(Polynomial(a1*l2*l3+a1+l2*l2*l3*l1*a2-l2*l2*a2-l2*l3*v3*l1*d2+l2*l3*a3+l2*v3*d2+l2*l1*a3+l3*v3*d2-l3*l1*a2-l3*l1*a3+v3*l1*d2+a2+a3,a1*l2*l3*v3-a1*v3+l2*l2*l3*v3*l1*a2+l2*l2*v3*a2-l2*l3*v3*a3+l2*l3*l1*d2+l2*v3*l1*a3+l2*d2-l3*v3*l1*a2+l3*v3*l1*a3-l3*d2-v3*a2+v3*a3+l1*d2));
    result.push_back(Polynomial(-a1*l2+a1*l3-l2*l2*l3*a2-l2*l2*a2*l1+l2*v3*l3*d2+l2*v3*l1*d2+l2*l3*l1*a3-l2*a3+v3*l3*l1*d2-v3*d2+l3*a2+l3*a3+a2*l1+l1*a3,-a1*l2*v3-a1*v3*l3+l2*l2*v3*l3*a2-l2*l2*v3*a2*l1+l2*v3*l3*l1*a3+l2*v3*a3+l2*l3*d2-l2*l1*d2-v3*l3*a2+v3*l3*a3+v3*a2*l1-v3*l1*a3+l3*l1*d2+d2));
    result.push_back(Polynomial(a1*l2*v3+a1*v3*l3-l2*l2*v3*l3*a2+l2*l2*v3*a2*l1-l2*v3*l3*l1*a3-l2*v3*a3-l2*l3*d2+l2*l1*d2+v3*l3*a2-v3*l3*a3-v3*a2*l1+v3*l1*a3-l3*l1*d2-d2,-a1*l2+a1*l3-l2*l2*l3*a2-l2*l2*a2*l1+l2*v3*l3*d2+l2*v3*l1*d2+l2*l3*l1*a3-l2*a3+v3*l3*l1*d2-v3*d2+l3*a2+l3*a3+a2*l1+l1*a3));
    result.push_back(Polynomial(-a1*l2*l3*v3+a1*v3-l2*l2*l3*v3*l1*a2-l2*l2*v3*a2+l2*l3*v3*a3-l2*l3*l1*d2-l2*v3*l1*a3-l2*d2+l3*v3*l1*a2-l3*v3*l1*a3+l3*d2+v3*a2-v3*a3-l1*d2,a1*l2*l3+a1+l2*l2*l3*l1*a2-l2*l2*a2-l2*l3*v3*l1*d2+l2*l3*a3+l2*v3*d2+l2*l1*a3+l3*v3*d2-l3*l1*a2-l3*l1*a3+v3*l1*d2+a2+a3));
    result.push_back(Polynomial(2*l2*l3*l1-2*l2+2*l3+2*l1,2*l2*l3*v3*l1+2*l2*v3+2*l3*v3-2*v3*l1));
    result.push_back(Polynomial(-2*l2*l3-2*l2*l1+2*l3*l1-2,2*l2*l3*v3-2*l2*v3*l1-2*l3*v3*l1-2*v3));
    result.push_back(Polynomial(-2*l2*l3*v3+2*l2*v3*l1+2*l3*v3*l1+2*v3,-2*l2*l3-2*l2*l1+2*l3*l1-2));
    result.push_back(Polynomial(-2*l2*l3*l1*v3-2*l2*v3-2*l3*v3+2*l1*v3,2*l2*l3*l1-2*l2+2*l3+2*l1));
  }
  return result;
 };

std::vector<Polynomial> h4_tc_v1 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];

  if (abs(a1)<1E-3) //rref condition
  {
    result.push_back(Polynomial(a2*l2*l2*l3*v3-a2*l3*v3-l2*l2*l3*v3*a3+l2*l2*d2-3*l2*l3*d2+l2*v3*a3,-a2*l2*l2*l3+a2*l3-l2*l2*l3*a3-l2*l2*v3*d2-3*l2*l3*v3*d2-l2*a3));
    result.push_back(Polynomial(-a2*l2*l2*v3+a2*v3+l2*l2*v3*a3+l2*l2*l3*d2+l2*v3*l3*a3+3*l2*d2,a2*l2*l2-a2-l2*l2*v3*l3*d2+l2*l2*a3+3*l2*v3*d2-l2*l3*a3));
    result.push_back(Polynomial(-a2*l2*l2+a2+l2*l2*v3*l3*d2-l2*l2*a3-3*l2*v3*d2+l2*l3*a3,-a2*l2*l2*v3+a2*v3+l2*l2*v3*a3+l2*l2*l3*d2+l2*v3*l3*a3+3*l2*d2));
    result.push_back(Polynomial(a2*l2*l2*l3-a2*l3+l2*l2*l3*a3+l2*l2*v3*d2+3*l2*l3*v3*d2+l2*a3,a2*l2*l2*l3*v3-a2*l3*v3-l2*l2*l3*v3*a3+l2*l2*d2-3*l2*l3*d2+l2*v3*a3));
    result.push_back(Polynomial(2*l2*l2*v3+2*l2*v3*l3,2*l2*l2-2*l2*l3));
    result.push_back(Polynomial(2*l2*l2*l3*v3-2*l2*v3,2*l2*l2*l3+2*l2));
    result.push_back(Polynomial(-2*l2*l2*l3-2*l2,2*l2*l2*l3*v3-2*l2*v3));
    result.push_back(Polynomial(-2*l2*l2+2*l2*l3,2*l2*l2*v3+2*l2*v3*l3));
  }
  else
  {
    result.push_back(Polynomial(-a1*l2*l2*l3*v3+a1*l2*v3+l2*l2*l3*v3*a2-l2*l2*l3*v3*a3-l2*l2*l3*l1*d2+l2*l2*v3*l1*a2-l2*l2*v3*l1*a3+l2*l2*d2-l2*l3*v3*l1*a3-l2*l3*d2+l2*v3*a3-l2*l1*d2-l3*v3*a2-v3*l1*a2,a1*l2*l2*l3+a1*l2-l2*l2*l3*v3*l1*d2-l2*l2*l3*a2-l2*l2*l3*a3-l2*l2*v3*d2+l2*l2*l1*a2+l2*l2*l1*a3-l2*l3*v3*d2-l2*l3*l1*a3+l2*v3*l1*d2-l2*a3+l3*a2-l1*a2));
    result.push_back(Polynomial(a1*l2*l2*v3+a1*l2*v3*l3+l2*l2*v3*l3*l1*a2-l2*l2*v3*l3*l1*a3-l2*l2*v3*a2+l2*l2*v3*a3+l2*l2*l3*d2+l2*l2*d2*l1+l2*v3*l3*a3+l2*v3*l1*a3-l2*l3*d2*l1+l2*d2-v3*l3*l1*a2+v3*a2,-a1*l2*l2+a1*l2*l3-l2*l2*v3*l3*d2+l2*l2*v3*d2*l1+l2*l2*l3*l1*a2+l2*l2*l3*l1*a3+l2*l2*a2+l2*l2*a3+l2*v3*l3*d2*l1+l2*v3*d2-l2*l3*a3+l2*l1*a3-l3*l1*a2-a2));
    result.push_back(Polynomial(a1*l2*l2-a1*l2*l3+l2*l2*v3*l3*d2-l2*l2*v3*l1*d2-l2*l2*l3*l1*a2-l2*l2*l3*l1*a3-l2*l2*a2-l2*l2*a3-l2*v3*l3*l1*d2-l2*v3*d2+l2*l3*a3-l2*l1*a3+l3*l1*a2+a2,a1*l2*l2*v3+a1*l2*v3*l3+l2*l2*v3*l3*l1*a2-l2*l2*v3*l3*l1*a3-l2*l2*v3*a2+l2*l2*v3*a3+l2*l2*l3*d2+l2*l2*l1*d2+l2*v3*l3*a3+l2*v3*l1*a3-l2*l3*l1*d2+l2*d2-v3*l3*l1*a2+v3*a2));
    result.push_back(Polynomial(-a1*l2*l2*l3-a1*l2+l2*l2*l3*v3*l1*d2+l2*l2*l3*a2+l2*l2*l3*a3+l2*l2*v3*d2-l2*l2*a2*l1-l2*l2*a3*l1+l2*l3*v3*d2+l2*l3*a3*l1-l2*v3*l1*d2+l2*a3-l3*a2+a2*l1,-a1*l2*l2*l3*v3+a1*l2*v3+l2*l2*l3*v3*a2-l2*l2*l3*v3*a3-l2*l2*l3*l1*d2+l2*l2*v3*a2*l1-l2*l2*v3*a3*l1+l2*l2*d2-l2*l3*v3*a3*l1-l2*l3*d2+l2*v3*a3-l2*l1*d2-l3*v3*a2-v3*a2*l1));
    result.push_back(Polynomial(-2*l2*l2*l3*l1*v3+2*l2*l2*v3+2*l2*l3*v3+2*l2*l1*v3,2*l2*l2*l3*l1+2*l2*l2-2*l2*l3+2*l2*l1));
    result.push_back(Polynomial(2*l2*l2*l3*v3+2*l2*l2*v3*l1+2*l2*l3*v3*l1-2*l2*v3,2*l2*l2*l3-2*l2*l2*l1+2*l2*l3*l1+2*l2));
    result.push_back(Polynomial(-2*l2*l2*l3+2*l2*l2*l1-2*l2*l3*l1-2*l2,2*l2*l2*l3*v3+2*l2*l2*v3*l1+2*l2*l3*v3*l1-2*l2*v3));
    result.push_back(Polynomial(-2*l2*l2*l3*l1-2*l2*l2+2*l2*l3-2*l2*l1,-2*l2*l2*l3*v3*l1+2*l2*l2*v3+2*l2*l3*v3+2*l2*v3*l1));
  }
  return result;
};

// jcapco Hyperplanes Tv2: 

std::vector<Polynomial> h1_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-l1*l3*l2*v3-l1*v3+l3*v3-l2*v3,l1*l3*l2-l1+l3+l2));
result.push_back(Polynomial(-l1*l3*v3+l1*v3*l2-l3*v3*l2-v3,-l1*l3-l1*l2+l3*l2-1));
result.push_back(Polynomial(l1*l3+l1*l2+l3*l2-1,-l1*l3*v3+l1*v3*l2+l3*v3*l2+v3));
result.push_back(Polynomial(-l1*l3*l2+l1+l3+l2,-l1*l3*v3*l2-l1*v3-l3*v3+v3*l2));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h2_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double v3 = a.v[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(l1*l3*l2-l1-l3-l2,l1*l3*v3*l2+l1*v3+l3*v3-v3*l2));
result.push_back(Polynomial(-l1*l3-l1*l2-l3*l2+1,l1*l3*v3-l1*v3*l2-l3*v3*l2-v3));
result.push_back(Polynomial(-l1*l3*v3+l1*v3*l2-l3*v3*l2-v3,-l1*l3-l1*l2+l3*l2-1));
result.push_back(Polynomial(-l1*l3*l2*v3-l1*v3+l3*v3-l2*v3,l1*l3*l2-l1+l3+l2));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h3_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(4*a2*l1*l1*l3*v3-4*a2*l1*l3*v3*l2-l1*l1*l1*l3*v3*l2*a3-2*l1*l1*l1*l3*d2-l1*l1*l1*v3*a3-l1*l1*l3*v3*a3-2*l1*l1*l3*l2*d2+l1*l1*v3*l2*a3+4*l1*l1*d2+l1*l3*v3*l2*a3+2*l1*l3*d2+l1*v3*a3+4*l1*l2*d2+l3*v3*a3+2*l3*l2*d2-v3*l2*a3,4*a2*l1*l1*l2-4*a2*l1-2*l1*l1*l1*l3*v3*d2-l1*l1*l1*l3*l2*a3+l1*l1*l1*a3-2*l1*l1*l3*v3*l2*d2+l1*l1*l3*a3+l1*l1*l2*a3-2*l1*l3*v3*d2+l1*l3*l2*a3-l1*a3-2*l3*v3*l2*d2-l3*a3-l2*a3));
result.push_back(Polynomial(-4*a2*l1*l1*v3+4*a2*l1*l2*v3-l1*l1*l1*l3*v3*a3+l1*l1*l1*l2*v3*a3+2*l1*l1*l1*d2+l1*l1*l3*l2*v3*a3+4*l1*l1*l3*d2+2*l1*l1*l2*d2+l1*l1*v3*a3+4*l1*l3*l2*d2+l1*l3*v3*a3-l1*l2*v3*a3-2*l1*d2-l3*l2*v3*a3-2*l2*d2-v3*a3,4*a2*l1*l1*l3*l2-4*a2*l1*l3+l1*l1*l1*l3*a3+l1*l1*l1*l2*a3+2*l1*l1*l1*v3*d2+l1*l1*l3*l2*a3+2*l1*l1*l2*v3*d2-l1*l1*a3-l1*l3*a3-l1*l2*a3+2*l1*v3*d2-l3*l2*a3+2*l2*v3*d2+a3));
result.push_back(Polynomial(-4*a2*l1*l1+4*a2*l1*l2-l1*l1*l1*l3*a3-2*l1*l1*l1*v3*d2-l1*l1*l1*l2*a3+4*l1*l1*l3*v3*d2+l1*l1*l3*l2*a3-2*l1*l1*v3*l2*d2-l1*l1*a3+4*l1*l3*v3*l2*d2+l1*l3*a3+2*l1*v3*d2+l1*l2*a3-l3*l2*a3+2*v3*l2*d2+a3,4*a2*l1*l1*l3*v3*l2-4*a2*l1*l3*v3-l1*l1*l1*l3*v3*a3+l1*l1*l1*v3*l2*a3+2*l1*l1*l1*d2-l1*l1*l3*v3*l2*a3-l1*l1*v3*a3+2*l1*l1*l2*d2+l1*l3*v3*a3-l1*v3*l2*a3+2*l1*d2+l3*v3*l2*a3+v3*a3+2*l2*d2));
result.push_back(Polynomial(4*a2*l1*l1*l3-4*a2*l1*l3*l2+2*l1*l1*l1*l3*v3*d2+l1*l1*l1*l3*l2*a3-l1*l1*l1*a3+2*l1*l1*l3*v3*l2*d2+l1*l1*l3*a3+4*l1*l1*v3*d2+l1*l1*l2*a3-2*l1*l3*v3*d2-l1*l3*l2*a3+4*l1*v3*l2*d2+l1*a3-2*l3*v3*l2*d2-l3*a3-l2*a3,4*a2*l1*l1*v3*l2-4*a2*l1*v3-l1*l1*l1*l3*v3*l2*a3-2*l1*l1*l1*l3*d2-l1*l1*l1*v3*a3+l1*l1*l3*v3*a3-2*l1*l1*l3*l2*d2-l1*l1*v3*l2*a3+l1*l3*v3*l2*a3-2*l1*l3*d2+l1*v3*a3-l3*v3*a3-2*l3*l2*d2+v3*l2*a3));
result.push_back(Polynomial(-2*l1*l1*l1*l3*v3+2*l1*l1*l1*v3*l2+2*l1*l1*l3*v3*l2+2*l1*l1*v3+2*l1*l3*v3-2*l1*v3*l2-2*l3*v3*l2-2*v3,2*l1*l1*l1*l3+2*l1*l1*l1*l2+2*l1*l1*l3*l2-2*l1*l1-2*l1*l3-2*l1*l2-2*l3*l2+2));
result.push_back(Polynomial(2*l1*l1*l1*l3*v3*l2+2*l1*l1*l1*v3+2*l1*l1*l3*v3-2*l1*l1*v3*l2-2*l1*l3*v3*l2-2*l1*v3-2*l3*v3+2*v3*l2,2*l1*l1*l1*l3*l2-2*l1*l1*l1-2*l1*l1*l3-2*l1*l1*l2-2*l1*l3*l2+2*l1+2*l3+2*l2));
result.push_back(Polynomial(-2*l1*l1*l1*l3*l2+2*l1*l1*l1-2*l1*l1*l3-2*l1*l1*l2+2*l1*l3*l2-2*l1+2*l3+2*l2,2*l1*l1*l1*l3*v3*l2+2*l1*l1*l1*v3-2*l1*l1*l3*v3+2*l1*l1*v3*l2-2*l1*l3*v3*l2-2*l1*v3+2*l3*v3-2*v3*l2));
result.push_back(Polynomial(-2*l1*l1*l1*l3-2*l1*l1*l1*l2+2*l1*l1*l3*l2-2*l1*l1+2*l1*l3+2*l1*l2-2*l3*l2+2,-2*l1*l1*l1*l3*v3+2*l1*l1*l1*v3*l2-2*l1*l1*l3*v3*l2-2*l1*l1*v3+2*l1*l3*v3-2*l1*v3*l2+2*l3*v3*l2+2*v3));

 return result;
 };

std::vector<Polynomial> h4_tc_v2 (Input& a)
{
    std::vector<Polynomial> result;
    double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double d2 = a.d[1], v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(4*a2*l1*l1*l2+4*a2*l1+l1*l1*l1*l3*l2*a3-2*l1*l1*l1*v3*l2*d2-l1*l1*l1*a3+4*l1*l1*l3*v3*l2*d2+l1*l1*l3*a3+2*l1*l1*v3*d2+l1*l1*l2*a3-4*l1*l3*v3*d2-l1*l3*l2*a3+2*l1*v3*l2*d2+l1*a3-l3*a3-2*v3*d2-l2*a3,-4*a2*l1*l1*l3*v3-4*a2*l1*l3*v3*l2-l1*l1*l1*l3*v3*l2*a3-l1*l1*l1*v3*a3-2*l1*l1*l1*l2*d2+l1*l1*l3*v3*a3-l1*l1*v3*l2*a3+2*l1*l1*d2+l1*l3*v3*l2*a3+l1*v3*a3-2*l1*l2*d2-l3*v3*a3+v3*l2*a3+2*d2));
result.push_back(Polynomial(4*a2*l1*l1*l3*l2+4*a2*l1*l3-2*l1*l1*l1*l3*l2*v3*d2-l1*l1*l1*l3*a3-l1*l1*l1*l2*a3+l1*l1*l3*l2*a3+2*l1*l1*l3*v3*d2-4*l1*l1*l2*v3*d2-l1*l1*a3+2*l1*l3*l2*v3*d2+l1*l3*a3+l1*l2*a3+4*l1*v3*d2-l3*l2*a3-2*l3*v3*d2+a3,4*a2*l1*l1*v3+4*a2*l1*l2*v3-2*l1*l1*l1*l3*l2*d2-l1*l1*l1*l3*v3*a3+l1*l1*l1*l2*v3*a3-l1*l1*l3*l2*v3*a3+2*l1*l1*l3*d2-l1*l1*v3*a3-2*l1*l3*l2*d2+l1*l3*v3*a3-l1*l2*v3*a3+l3*l2*v3*a3+2*l3*d2+v3*a3));
result.push_back(Polynomial(4*a2*l1*l1*l3*v3*l2+4*a2*l1*l3*v3+l1*l1*l1*l3*v3*a3+2*l1*l1*l1*l3*l2*d2-l1*l1*l1*v3*l2*a3-l1*l1*l3*v3*l2*a3-2*l1*l1*l3*d2-l1*l1*v3*a3-4*l1*l1*l2*d2-l1*l3*v3*a3-2*l1*l3*l2*d2+l1*v3*l2*a3+4*l1*d2+l3*v3*l2*a3+2*l3*d2+v3*a3,4*a2*l1*l1+4*a2*l1*l2-2*l1*l1*l1*l3*v3*l2*d2-l1*l1*l1*l3*a3-l1*l1*l1*l2*a3+2*l1*l1*l3*v3*d2-l1*l1*l3*l2*a3+l1*l1*a3-2*l1*l3*v3*l2*d2+l1*l3*a3+l1*l2*a3+2*l3*v3*d2+l3*l2*a3-a3));
result.push_back(Polynomial(4*a2*l1*l1*v3*l2+4*a2*l1*v3+l1*l1*l1*l3*v3*l2*a3+l1*l1*l1*v3*a3+2*l1*l1*l1*l2*d2+l1*l1*l3*v3*a3+4*l1*l1*l3*l2*d2-l1*l1*v3*l2*a3-2*l1*l1*d2-l1*l3*v3*l2*a3-4*l1*l3*d2-l1*v3*a3-2*l1*l2*d2-l3*v3*a3+v3*l2*a3+2*d2,-4*a2*l1*l1*l3-4*a2*l1*l3*l2+l1*l1*l1*l3*l2*a3-2*l1*l1*l1*v3*l2*d2-l1*l1*l1*a3-l1*l1*l3*a3+2*l1*l1*v3*d2-l1*l1*l2*a3-l1*l3*l2*a3-2*l1*v3*l2*d2+l1*a3+l3*a3+2*v3*d2+l2*a3));
result.push_back(Polynomial(-2*l1*l1*l1*l3-2*l1*l1*l1*l2+2*l1*l1*l3*l2-2*l1*l1+2*l1*l3+2*l1*l2-2*l3*l2+2,-2*l1*l1*l1*l3*v3+2*l1*l1*l1*v3*l2-2*l1*l1*l3*v3*l2-2*l1*l1*v3+2*l1*l3*v3-2*l1*v3*l2+2*l3*v3*l2+2*v3));
result.push_back(Polynomial(-2*l1*l1*l1*l3*l2+2*l1*l1*l1-2*l1*l1*l3-2*l1*l1*l2+2*l1*l3*l2-2*l1+2*l3+2*l2,2*l1*l1*l1*l3*v3*l2+2*l1*l1*l1*v3-2*l1*l1*l3*v3+2*l1*l1*v3*l2-2*l1*l3*v3*l2-2*l1*v3+2*l3*v3-2*v3*l2));
result.push_back(Polynomial(-2*l1*l1*l1*l3*l2*v3-2*l1*l1*l1*v3-2*l1*l1*l3*v3+2*l1*l1*l2*v3+2*l1*l3*l2*v3+2*l1*v3+2*l3*v3-2*l2*v3,-2*l1*l1*l1*l3*l2+2*l1*l1*l1+2*l1*l1*l3+2*l1*l1*l2+2*l1*l3*l2-2*l1-2*l3-2*l2));
result.push_back(Polynomial(2*l1*l1*l1*l3*v3-2*l1*l1*l1*v3*l2-2*l1*l1*l3*v3*l2-2*l1*l1*v3-2*l1*l3*v3+2*l1*v3*l2+2*l3*v3*l2+2*v3,-2*l1*l1*l1*l3-2*l1*l1*l1*l2-2*l1*l1*l3*l2+2*l1*l1+2*l1*l3+2*l1*l2+2*l3*l2-2));

 return result;
 };

// jcapco Hyperplanes RRP special case (planar): 

std::vector<Polynomial> h1_tsp (Input& a)
{
    std::vector<Polynomial> result;
    //double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(-l2-l3));
result.push_back(Polynomial(-l2*l3+1));
result.push_back(Polynomial(-l2*l3*v3-v3));
result.push_back(Polynomial(-l2*v3+v3*l3));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h2_tsp (Input& a)
{
    std::vector<Polynomial> result;
    //double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial(l2*v3-v3*l3));
result.push_back(Polynomial(l2*l3*v3+v3));
result.push_back(Polynomial(-l2*l3+1));
result.push_back(Polynomial(-l2-l3));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));
result.push_back(Polynomial(0));

 return result;
 };

std::vector<Polynomial> h3_tsp (Input& a)
{
    std::vector<Polynomial> result;
    //double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial((-a2*l2-a2*l3-l2*a3-l3*a3)/2));
result.push_back(Polynomial((-a2*l2*l3+a2-l2*l3*a3+a3)/2));
result.push_back(Polynomial((-a2*l2*l3*v3-a2*v3+l2*l3*v3*a3+v3*a3)/2));
result.push_back(Polynomial((-a2*l2*v3+a2*v3*l3+l2*v3*a3-v3*l3*a3)/2));
result.push_back(Polynomial(-l2*l3+1));
result.push_back(Polynomial(l2+l3));
result.push_back(Polynomial(-l2*v3+v3*l3));
result.push_back(Polynomial(l2*l3*v3+v3));

 return result;
 };

std::vector<Polynomial> h4_tsp (Input& a)
{
    std::vector<Polynomial> result;
    //double l1 = a.al[0], l2 = a.al[1], l3 = a.al[2];
    double l2 = a.al[1], l3 = a.al[2];
    //double d1 = a.d[0], d2 = a.d[1], v3 = a.v[2];
    double v3 = a.v[2];
    //double a1 = a.a[0], a2 = a.a[1], a3 = a.a[2];
    double a2 = a.a[1], a3 = a.a[2];
result.push_back(Polynomial((a2*l2*v3-a2*v3*l3-l2*v3*a3+v3*l3*a3)/2));
result.push_back(Polynomial((a2*l2*l3*v3+a2*v3-l2*l3*v3*a3-v3*a3)/2));
result.push_back(Polynomial((-a2*l2*l3+a2-l2*l3*a3+a3)/2));
result.push_back(Polynomial((-a2*l2-a2*l3-l2*a3-l3*a3)/2));
result.push_back(Polynomial(-l2*l3*v3-v3));
result.push_back(Polynomial(l2*v3-v3*l3));
result.push_back(Polynomial(l2+l3));
result.push_back(Polynomial(-l2*l3+1));

 return result;
 };

}
}
