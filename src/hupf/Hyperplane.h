/**
 * This file contains the Hyperplane equations
 * */
#pragma once

//PI defined in Input
#define PI 3.1415926535897932384626433

#include <hupf/Input.h>
#include <hupf/rrr/t2.h>
#include <hupf/rrr/left_3r.h>
#include <hupf/rrr/right_3r.h>
#include <hupf/rrp/left_rrp.h>

namespace LibHUPF
{

static double EPSILON=1e-6;
using namespace std; //jcapco todo: temporary

enum ThreeChain
{
  //RRP=1+2+0=3, RPR=1+0+4=5, PRR=0+2+4=6, PPR=4, RRR=1+2+4=7, 
  //PRP=0+2+0=2..etc. 8 possibilities!
  RRP_ = 3,
  RPR_ = 5,
  PRR_ = 6,
  PPR_ = 4,
  RRR_ = 7,
  PRP_ = 2
};

class Hyperplane
{

public:
  std::vector<std::vector<Polynomial> > h;
  //{1,2,3}x{4,5,6} corresponding to Tv/di x Tv/dj i=1,2,3; j=4,5,6
  int manifoldsUsed[2]; 

  Input a;
  size_t left_3r, right_3r;

  Hyperplane() {};

  Hyperplane(Input inputParms)
  {
    a = inputParms;
    manifoldsUsed[0] = 0;
    manifoldsUsed[1] = 0;

    //jcapco todo: this could be 9 (tcp) or 10 (tcw)
    left_3r = 0; right_3r =0;
    for (size_t i=0; i<3;++i)
    {
      left_3r += (inputParms.rots[i] << i);
      right_3r += (inputParms.rots[i+3] << i);
    }
    
    //RRP=1+2+0=3, RPR=1+0+4=5, PRR=0+2+4=6, PPR=4, RRR=1+2+4=7, 
    //PRP=0+2+0=2..etc. 8 possibilities!
    if (left_3r == RRR_) //RRR
    {
      //jcapco todo: convert everything to simplified form following the sequence 
      if (fabs(a.alpha[0])<EPSILON && fabs(a.alpha[1])<EPSILON) //RRR spherical special case
      {
    		manifoldsUsed[0]=0;    
      }
      else if (!(fabs(a.alpha[0])<EPSILON || fabs(a.a[0])<EPSILON)) //RRR Tv3 possible
      {
    		manifoldsUsed[0]=3;
        h.push_back(h1_tc_v3(a));
        h.push_back(h2_tc_v3(a));
        h.push_back(h3_tc_v3(a));
        h.push_back(h4_tc_v3(a));
      }
      else if (!(fabs(a.alpha[1])<EPSILON || fabs(a.a[1])<EPSILON)) //RRR Tv1 possible
      {
    		manifoldsUsed[0]=1;
        h.push_back(h1_tc_v1(a));
        h.push_back(h2_tc_v1(a));
        h.push_back(h3_tc_v1(a));
        h.push_back(h4_tc_v1(a));
      }
      else if (!((fabs(a.alpha[0] - Math::pi/2)<EPSILON || fabs(a.alpha[0]+Math::pi/2)<EPSILON) && (fabs(a.alpha[1]-Math::pi/2)<EPSILON || 
        fabs(a.alpha[1]+Math::pi/2)<EPSILON) && fabs(a.d[1])<EPSILON)) //RRR Tv2 possible
      {
		manifoldsUsed[0]=2;
        h.push_back(h1_tc_v2(a));
        h.push_back(h2_tc_v2(a));
        h.push_back(h3_tc_v2(a));
        h.push_back(h4_tc_v2(a));
      }    
    }
    else if (left_3r == RRP_) //RRP
    {
      if (!(fabs(a.alpha[0])<EPSILON || fabs(a.a[0])<EPSILON))
      {
        manifoldsUsed[0]=3;
		    //T(d3)
        h.push_back(RRP::h1_tc_d3(a)); //todo use same hyperplanes with switching!
        h.push_back(RRP::h2_tc_d3(a));
        h.push_back(RRP::h3_tc_d3(a));
        h.push_back(RRP::h4_tc_d3(a));
      }
      else if (!(fabs(a.alpha[1]+Math::pi/2)<EPSILON || fabs(a.alpha[1]-Math::pi/2)<EPSILON))
      {
    		manifoldsUsed[0]=1;
        //T(v1)
        h.push_back(RRP::h1_tc_v1(a));
        h.push_back(RRP::h2_tc_v1(a));
        h.push_back(RRP::h3_tc_v1(a));
        h.push_back(RRP::h4_tc_v1(a));        
      }
      else if (!(fabs(a.alpha[0])<EPSILON)) //Tv2, considering Td3 and Tv1 are not valid
      {
    		manifoldsUsed[0]=2;
        h.push_back(RRP::h1_tc_v2(a));
        h.push_back(RRP::h2_tc_v2(a));
        h.push_back(RRP::h3_tc_v2(a));
        h.push_back(RRP::h4_tc_v2(a));
      }
      else
      {
        manifoldsUsed[0]=0;
		    h.push_back(RRP::h1_tsp(a));
        h.push_back(RRP::h2_tsp(a));
        h.push_back(RRP::h3_tsp(a));
        h.push_back(RRP::h4_tsp(a));        
      }            
    }
    
    //todo right_3r, RRP        
    if (right_3r == RRR_)
    {
      // T(v6) usable when a[3]!=0 AND alpha[3]!=0
      // (Algorithm 3: T(v1)_left needs a2_rev=-a[3]!=0 and l2_rev=-tan(alpha[3]/2)!=0)
      if (!(fabs(a.alpha[3])<EPSILON || fabs(a.a[3])<EPSILON))
      {
        manifoldsUsed[1]=6;
        h.push_back(h1_v6q(a));
        h.push_back(h2_v6q(a));
        h.push_back(h3_v6q(a));
        h.push_back(h4_v6q(a));
      }
      // T(v4) usable when a[4]!=0 AND alpha[4]!=0
      // (Algorithm 3: T(v3)_left needs a1_rev=-a[4]!=0 and l1_rev=-tan(alpha[4]/2)!=0)
      else if (!(fabs(a.alpha[4])<EPSILON || fabs(a.a[4])<EPSILON))
	    {
		    manifoldsUsed[1]=4;
        h.push_back(h1_v4q(a));
        h.push_back(h2_v4q(a));
        h.push_back(h3_v4q(a));
        h.push_back(h4_v4q(a));
      }
      else if (!((fabs(a.alpha[3]-PI/2)<EPSILON || fabs(a.alpha[3]+PI/2)<EPSILON) && (fabs(a.alpha[4]-PI/2)<EPSILON || fabs(a.alpha[4]+PI/2)<EPSILON) && fabs(a.d[4])<EPSILON))
      {
        manifoldsUsed[1]=5;
		    //std::cout << "T(v5)\n";
        h.push_back(h1_v5q(a));
        h.push_back(h2_v5q(a));
        h.push_back(h3_v5q(a));
        h.push_back(h4_v5q(a));
      }    
    }
  };

};

}