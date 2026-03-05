/**
* Represents the surface in Kinematic Image Space after solving 7 pencils of hyperplanes
**/
#pragma once

#include <hupf/Hyperplane.h>
// omp.h removed for portability - timing stubs below
static inline double omp_get_wtime() { return 0.0; }

namespace LibHUPF
{

//adds p_add to polyCoefficient[0] of p and modifies the degree of p if p_add is 0 degree
inline void addPolyCoeff(Polynomial& p, Polynomial p_add)
{
  if (p_add.degree < 0) p.degree = -1;
  else
  {
    p.degree = 0;
    p.polyCoefficient[0]=p_add; 
  }
}

class KinematicSurface
{
public:
  // Equations for kinematic image surface. Polynomials are equations with polynomial coeff.
  Polynomial p0,p1,p2,p3,p4,p5,p6,p7;
  Hyperplane hyperplane;
  vector<Polynomial> r;
  KinematicSurface() {};
  
  /**
  * Work to do for Thread 1, calculates surfaces 0,1,4,5 and writes results into kinematicSurface kin
  * @param kin the kinematicSurface for writing the solutions
  **/
  static void *compute1(void* kin)
  {
    KinematicSurface* kinem=(KinematicSurface*)kin;    
    //preparing a processing the matrices parallel;
    Matrix mat[4];
    for (size_t i=0; i<4; ++i)
      mat[i] = Matrix(7,7,true);

    //solving SLE with projective space, assume x0=1, unknowns: x1,x2,x3,y0,y1,y2,y3
    //uses Cramer's Rule!
    //coeff if x0 will be the negative of the answer-vector of the SLE
    //for easier explanation assume T(u1) and T(u4) is chosen, so we choose 4 H-planes from T(u1) and 3 from T(u4)
    //constant coeff correspond to u1
        
    for(int i=0; i<=6; i++)
    {
      for(int j=0; j<=6; j++)
      {
        if(i>3)
        {
          //(R[u1])[u4] , but the coeff. of u4 are constants so it is R[u4]
          mat[0].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1]; //main coefficient matrix, j+1 because x0 is set to 1

          if(j==0)
            mat[1].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]); 
          else
            mat[1].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];

          if(j==3)
            mat[2].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]);
          else
            mat[2].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];

          if(j==4)
            mat[3].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]);
          else
            mat[3].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];
        }
        else
        {
          //constant coeff corresponding to u1 since A=(R[u1])[u4] but since there is no u[4] its the constant term of the polynomial A
          Polynomial p(0,true);          
          addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
          //p.polyCoefficient[0]=(kinem->hyperplane.h[i])[j+1]; 
          mat[0].p_matrix[i][j] = p;

          if(j==0)
          {
            //Polynomial p(0,true);
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            //p.polyCoefficient[0]=Polynomial::negPoly((kinem->hyperplane.h[i])[0]); 
            mat[1].p_matrix[i][j] = p;
          }
          else
          {
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(kinem->hyperplane.h[i])[j+1];
            mat[1].p_matrix[i][j] = p;
          }

          if(j==3)
          {
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            mat[2].p_matrix[i][j] = p;
          }
          else
          {
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=((kinem->hyperplane.h[i])[j+1]);
            mat[2].p_matrix[i][j] = p;
          }

          if(j==4)
          {
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            mat[3].p_matrix[i][j] = p;
          }
          else
          {
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=((kinem->hyperplane.h[i])[j+1]);
            mat[3].p_matrix[i][j] = p;
          }
        }
      }
    }    

    for (int i=0; i<4; ++i)
    {
      Polynomial p = Matrix::myDeterminant(mat[i]);
      for (size_t j=0; j<p.polyCoefficient.size();++j)
      {
        p.polyCoefficient[j].clean();
      }
      kinem->r[(i%2)+4*(i/2)]=p;
    }

    return 0;
  };

  /**
  * Work to do for Thread 2, calculates surfaces 2,3,6,7 and writes results into kinematicSurface kin
  * @param kin the kinematicSurface for writing the solutions
  **/
  static void *compute2(void* kin)
  {
    KinematicSurface* kinem=(KinematicSurface*)kin;
    Matrix mat[4];
    for (size_t i=0; i<4; ++i)
      mat[i] = Matrix(7,7,true);

    //solving SLE with projective space, assume x0=1, unknowns: x1,x2,x3,y0,y2,y3
    for(int i=0; i<=6; i++)
    {
      for(int j=0; j<=6; j++)
      {
        if(i>3)
        {
          if(j==1)
            mat[0].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]);
          else
            mat[0].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];

          if(j==2)
            mat[1].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]);
          else
            mat[1].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];

          if(j==5)
            mat[2].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]);
          else
            mat[2].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];

          if(j==6)
            mat[3].p_matrix[i][j] = Polynomial::negPoly((kinem->hyperplane.h[i])[0]);
          else
            mat[3].p_matrix[i][j] = (kinem->hyperplane.h[i])[j+1];

        }
        else
        {
          Polynomial p(0,true);
          if(j==1)
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            mat[0].p_matrix[i][j] = p;
          }
          else
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=((kinem->hyperplane.h[i])[j+1]);
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            mat[0].p_matrix[i][j] = p;
          }

          if(j==2)
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            mat[1].p_matrix[i][j] = p;
          }
          else
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=((kinem->hyperplane.h[i])[j+1]);
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            mat[1].p_matrix[i][j] = p;
          }

          if(j==5)
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            mat[2].p_matrix[i][j] = p;
          }
          else
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=((kinem->hyperplane.h[i])[j+1]);
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            mat[2].p_matrix[i][j] = p;
          }

          if(j==6)
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=(Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            addPolyCoeff(p, Polynomial::negPoly((kinem->hyperplane.h[i])[0]));
            mat[3].p_matrix[i][j] = p;
          }
          else
          {
            //Polynomial p(0,true);
            //p.polyCoefficient[0]=((kinem->hyperplane.h[i])[j+1]);
            addPolyCoeff(p, (kinem->hyperplane.h[i])[j+1]);
            mat[3].p_matrix[i][j] = p;
          }
        }
      }
    }

    for (int i=0; i<4; ++i)
    {
      Polynomial p = Matrix::myDeterminant(mat[i]);
      for (size_t j=0; j<p.polyCoefficient.size();++j)
      {
        p.polyCoefficient[j].clean();
      }
      kinem->r[(i%2)+2+4*(i/2)]=p;
    }

    return 0;
  };

  /**
  * constructor to calculate kinematicSurfaces, either parallel on two threads or sequential
  * @param hyperplane1 the hyperplane input
  **/
  KinematicSurface(Hyperplane &hyperplane1, vector<double> &times)
  {
    hyperplane=hyperplane1;
    r.resize(8);

    //sequential calculation    
    double end = 0;
    double start = omp_get_wtime();    
    compute1(this);
    //jcapco todo: debug, check rank by taking random u1 and u4
    end = omp_get_wtime();
    times.push_back(end-start);
    start = end;
    compute2(this);    
    times.push_back(omp_get_wtime()-start);
  }

  /**
  * Calculate the value of all r's by substituting the 1st and 2nd root.
  * @param va First Root
  * @param vb Second Root
  * @return All r's after substitution and evaluation
  **/
  vector<double> substituteValues(double va, double vb)
  {
    vector<double> result;
    for(int i=0; i<8; i++)
    {
      result.push_back(r[i].subs(va).eval(vb));
    }
    return result;
  };

  /**
  * Calculate the value of all r's by substituting the root for wrist partiioned manipulators
  * @param va root to be substitued
  * @return All r's after substitution and evaluation
  **/
  vector<double> substituteValues(double va)
  {
    vector<double> result;
    for(int i=0; i<8; i++)
    {
      result.push_back(r[i].subs(va).coefficient[0]);
    }
    return result;
  };
};

} //namespace