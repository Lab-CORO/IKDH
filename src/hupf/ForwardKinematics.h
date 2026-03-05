#pragma once

#include <hupf/Hyperplane.h>

namespace LibHUPF
{

class ForwardKinematics
{

public:
  double a[6];
  double d[6];
  double alpha[6];
  double theta[6];
  Matrix eePose;
  Matrix transformationMatrix[6];
  ForwardKinematics(Input &input, vector<double> th)
  {
    eePose = Matrix(4,4);
    for(int i=0; i<6; i++)   //length param
    {
      a[i] = input.a[i];
      if (input.rots[i])
      {
        d[i] = input.d[i];
        theta[i] = th[i]*Math::deg2rad;
      }
      else
      {
        theta[i]= 0; //todo use theta given in input!
        d[i] = th[i];
      }
      
      alpha[i] = input.alpha[i];
    }
  };
  double maxElement(Hyperplane &h)
  {
    double result=0;
    Matrix eePoseDiff = calculateFwdKin() - h.a.eePose; //hyperplane variable
    if(d[0] != 0)
    {
      eePoseDiff.set(3,0,eePoseDiff.get(3,0)-d[0]);
    }
    for(int i=1; i<4; i++)
    {
      for(int j=0; j<4; j++)
      {
        if(fabs(eePoseDiff.get(i,j)) > result)
        {
          result = fabs(eePoseDiff.get(i,j));
        }
      }
    }
    return result;
  };
  Matrix calculateFwdKin()
  {
    for(int i=0; i<6; i++)
    {
      transformationMatrix[i] = Matrix(4, 4);
      transformationMatrix[i].set(0,0,1);
      transformationMatrix[i].set(0,1,0);
      transformationMatrix[i].set(0,2,0);
      transformationMatrix[i].set(0,3,0);
      transformationMatrix[i].set(1,1,cos(theta[i]));
      transformationMatrix[i].set(1,2,-cos(alpha[i]) * sin(theta[i]));
      transformationMatrix[i].set(1,3,sin(alpha[i]) * sin(theta[i]));
      transformationMatrix[i].set(1,0,a[i] * cos(theta[i]));
      transformationMatrix[i].set(2,1,sin(theta[i]));
      transformationMatrix[i].set(2,2,cos(alpha[i]) * cos(theta[i]));
      transformationMatrix[i].set(2,3,-sin(alpha[i]) * cos(theta[i]));
      transformationMatrix[i].set(2,0,a[i] * sin(theta[i]));
      transformationMatrix[i].set(3,1,0);
      transformationMatrix[i].set(3,2,sin(alpha[i]));
      transformationMatrix[i].set(3,3,cos(alpha[i]));
      transformationMatrix[i].set(3,0,d[i]);
    }
    eePose = transformationMatrix[0]*transformationMatrix[1]*transformationMatrix[2]*transformationMatrix[3]*transformationMatrix[4]*transformationMatrix[5];
    return eePose;
  };
};

}//namespace