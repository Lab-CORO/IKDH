/**
 * Provides IO functionality for the input parameters and the generated solution
 * */
#pragma once

#include <sstream>
#include <fstream>
#include <iomanip>
#include <hupf/Matrix.h>

namespace Math
{  
  const float pi = 3.141592653589793238462643383279502884f;
  const float deg2rad = float(pi/180.0f);
}

namespace LibHUPF
{
//#define PI 3.1415926535897932384626433
class Input
{
public:
  double a[6]; //!< link length
  double d[6]; //!< link offsets
  double v[6]; //!< link half-angle tangents (for P case), todo compute!
  double alpha[6]; //!< twist angles in radians
  double al[6]; //!< half tangent substitution value for twist angles
  Matrix eePose; //!< end effector position matrix
  std::string path;  //!< input file path
  bool rots[6]; //if a joint is prismatice it is false, otherwise true

  Input() : eePose(4,4)
  {

  };

  Matrix studyParams()
  {
    Matrix eePose_Study(8,1);
    if(eePose.numCols==4 && eePose.numRows==4)
    {
      eePose_Study.set(0,0,1+eePose.get(1,1)+eePose.get(2,2)+eePose.get(3,3));
      eePose_Study.set(1,0,eePose.get(3,2)-eePose.get(2,3));
      eePose_Study.set(2,0,eePose.get(1,3)-eePose.get(3,1));
      eePose_Study.set(3,0,eePose.get(2,1)-eePose.get(1,2));

      if(eePose_Study.get(0,0)==0 && eePose_Study.get(1,0)==0 &&
          eePose_Study.get(2,0)==0 && eePose_Study.get(3,0)==0)
      {
        eePose_Study.set(0,0,eePose.get(3,2)-eePose.get(2,3));
        eePose_Study.set(1,0,1+eePose.get(1,1)-eePose.get(2,2)-eePose.get(3,3));
        eePose_Study.set(2,0,eePose.get(2,1)+eePose.get(1,2));
        eePose_Study.set(4,0,eePose.get(1,3)+eePose.get(3,1));

        if(eePose_Study.get(0,0)==0 && eePose_Study.get(1,0)==0 &&
            eePose_Study.get(2,0)==0 && eePose_Study.get(3,0)==0)
        {
          eePose_Study.set(0,0,eePose.get(1,3)-eePose.get(3,1));
          eePose_Study.set(1,0,eePose.get(2,1)+eePose.get(1,2));
          eePose_Study.set(2,0,1-eePose.get(1,1)+eePose.get(2,2)-eePose.get(3,3));
          eePose_Study.set(3,0,eePose.get(2,3)+eePose.get(3,2));

          if(eePose_Study.get(0,0)==0 && eePose_Study.get(1,0)==0 &&
              eePose_Study.get(2,0)==0 && eePose_Study.get(3,0)==0)
          {
            eePose_Study.set(0,0,eePose.get(2,1)-eePose.get(1,2));
            eePose_Study.set(1,0,eePose.get(1,3)+eePose.get(3,1));
            eePose_Study.set(2,0,eePose.get(2,3)+eePose.get(3,2));
            eePose_Study.set(3,0,1-eePose.get(1,1)-eePose.get(2,2)+eePose.get(3,3));
          }
        }
      }
      eePose_Study.set(4,0,(eePose.get(1,0)*eePose_Study.get(1,0)+eePose.get(2,0)*eePose_Study.get(2,0)+eePose.get(3,0)*eePose_Study.get(3,0))/2);
      eePose_Study.set(5,0,(-eePose.get(1,0)*eePose_Study.get(0,0)+eePose.get(3,0)*eePose_Study.get(2,0)-eePose.get(2,0)*eePose_Study.get(3,0))/2);
      eePose_Study.set(6,0,(-eePose.get(2,0)*eePose_Study.get(0,0)-eePose.get(3,0)*eePose_Study.get(1,0)+eePose.get(1,0)*eePose_Study.get(3,0))/2);
      eePose_Study.set(7,0,(-eePose.get(3,0)*eePose_Study.get(0,0)+eePose.get(2,0)*eePose_Study.get(1,0)-eePose.get(1,0)*eePose_Study.get(2,0))/2);
    }
    //jcapco now: normalize
    return eePose_Study*(0.5/Matrix::FrobeniusNorm(eePose_Study));
    //return eePose_Study;
  };

  /**
   * tangent of half angle substitution.
   * */
  void halfAngleSubs()
  {
    for(int i=0; i<=5; i++)
    {
      al[i] = tan(alpha[i]/2);
    }
  };

  /**
   * inverse action of half angle substitution.
   * */
  static double inverseHalfAngleSubs(const double& v)
  {
    double result=0.0;
    result = 2*atan(v);
    return result;
  };

  Input rightSubs()
  {
    Input result;
    result.a[0] = result.a[4];
    result.a[1] = result.a[3];
    result.a[2] = 0;
    result.d[1] = result.d[4];
    result.d[2] = result.d[3];
    result.al[0] = result.al[4];
    result.al[1] = result.al[3];
    result.al[2] = 0;
    return result;
  };

  /**
   * parse input file. read out all the parameters.
   * @param path path to input file
   * @return pointer to Input object, note that this pointer should be deleted after usage. contains robot (and NOT endeffector informations)
   * */
  static Input* parseFile(const std::string& path)
  {
    using namespace std;

    //Input res;
    string line;
    Input* tmp = new Input();

    ifstream myfile(path.c_str());
    if(myfile.is_open())
    {

      //search [DH] block in file

      while(line.compare(0,4,"[DH]") == 0)
        if(!getline (myfile,line))
        {
          perror("File format not correct, or error while Reading0...");
          myfile.close();
          exit(1);
        }

      bool alphaIsRad=false;

      //read all 4 dh parameters
      for(int k=0; k<4; k++)
      {
        do
        {
          if(!getline (myfile,line))
          {
            perror("File format not correct, or error while Reading1...");
            myfile.close();
            exit(1);
          }
        }
        while((line.compare(0,3,"a =")!=0)&&(line.compare(0,3,"d =")!=0)&&(line.compare(0,9,"measure =")!=0)&&(line.compare(0,7,"alpha =")!=0)&&
              (line.compare(0,2,"a=")!=0)&&(line.compare(0,2,"d=")!=0)&&(line.compare(0,8,"measure=")!=0)&&(line.compare(0,6,"alpha=")!=0));

        if (line.compare(0,3,"a =")==0 || line.compare(0,2,"a=")==0)
        {
          string subline = line.substr(3,line.length());
          istringstream linestream(subline);
          string item;
          for(int i=0; i<6; i++)
          {
            if(getline(linestream, item, ','))
              tmp->a[i]=atof(item.c_str());
            else
            {
              perror("File format not correct!");
            }
          }

        }
        else if (line.compare(0,3,"d =")==0 || line.compare(0,2,"d=")==0)
        {
          string subline = line.substr(3,line.length());
          istringstream linestream(subline);
          string item;
          for(int i=0; i<6; i++)
          {
            if(getline(linestream, item, ','))
              tmp->d[i]=atof(item.c_str());
            else
            {
              perror("File format not correct!");
            }
          }
        }
        else if (line.compare(0,9,"measure =")==0 || line.compare(0,8,"measure=")==0)
        {
          if(line.compare(10,1,"d")==0)
            alphaIsRad=false;
          else if (line.compare(10,1,"r")==0)
            alphaIsRad=true;
          else
          {
            perror("Format of angle not correct!");
          }

        }
        else if (line.compare(0,7,"alpha =")==0 || line.compare(0,6,"alpha=")==0)
        {
          string subline = line.substr(7,line.length());
          istringstream linestream(subline);
          string item;
          for(int i=0; i<6; i++)
          {
            if(getline(linestream, item, ','))
              tmp->alpha[i]=atof(item.c_str());
            else
            {
              perror("File format not correct!");
            }
          }

          if(!alphaIsRad)
          {
            for(int i=0; i<6; i++)
            {
              tmp->alpha[i]=tmp->alpha[i]*Math::deg2rad;
              tmp->al[i]=tan(tmp->alpha[i]/2);
            }
          }
        }
      }
      myfile.seekg(0, ios_base::beg);  //reset filePointer

      //while(getline(myfile,line))
      getline(myfile,line);
      
      if (line.compare(0,4,"[EE]")==0) 		//new [EE] block detected
      {
        /*
        Input resInput;

        //populate new Input with DHs
        for(int k=0; k<6; k++)
        {
          resInput.a[k]=tmp.a[k];
          resInput.d[k]=tmp.d[k];
          resInput.alpha[k]=tmp.alpha[k];
          resInput.al[k]=tmp.al[k];
        }
        */
        
        for(int k=0; k<4; k++)
        {
          if(!getline(myfile,line))
          {
            perror("File format not correct!");
            myfile.close();
            exit(1);
          }

          istringstream linestream(line);
          string item;
          for(int i=0; i<4; i++)
          {
            if(getline(linestream, item, ','))
              //resInput.eePose.m_matrix[k][i]=atof(item.c_str());
              (tmp->eePose).m_matrix[k][i]=atof(item.c_str());
            else
            {
              perror("File format not correct!");
            }
          }
        }
        //res.push_back(resInput);
        //while ends
        
      }
    }
    else
    {
      perror("Error opening file!");
      myfile.close();
      exit(1);
    }
    myfile.close();
    return tmp;
  };
  
};

inline void clean_half_angle_tan(double& a)
{
  double eps = 1e-3;
  if (fabs(a)< eps) a=0;
  else if (fabs(a-1)<eps) a=1;
  else if (fabs(a+1)<eps) a=-1;
}

 /**
 * parse input file. read out all the parameters.
 * @param dh parameters: a is an array for length, d is an array for offset, alpha is the twist angle
 * @return pointer to Input object, note that this pointer should be deleted after usage. contains robot (and NOT endeffector informations)
 * */
static Input* getInput(double* a, double* d, double* theta, double* alpha, bool* rots)
{
  Input* tmp = new Input();

  bool alphaIsRad=true;

  //read all 4 dh parameters
  std::copy(a,a+6, tmp->a);
  std::copy(d,d+6, tmp->d);
  std::copy(rots, rots+6, tmp->rots);
  std::copy(alpha,alpha+6, tmp->alpha);

  if(!alphaIsRad)
    for(int i=0; i<6; i++)
    {
      tmp->alpha[i]=tmp->alpha[i]*Math::deg2rad;
    }

  for(int i=0; i<6; i++)
  {
    if (!(tmp->rots[i]))
      tmp->v[i]=tan(theta[i]/2);
    tmp->al[i]=tan(tmp->alpha[i]/2);
  }

  return tmp;
};

}//namespace 