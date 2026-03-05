#pragma once
/**
*	Provides functionality to calculate angles after having eqns for Kinematic Imgage Space.
**/
#include <hupf/ForwardKinematics.h>
#include <hupf/KinematicSurface.h>
#include <hupf/rpoly.h>
#include <hupf/rrr/wrist_sol.h>
#include <hupf/poly_capture.h>

//debug
//#include <iostream>

namespace LibHUPF
{

class SolveForAngles
{
public:
  /**
  * Main Method in this class, checks the type of input (hyperplane) and calls approprite method (general case or special case)
  * After execution sorts the results
  * @param hyperplane The hyperplane input
  * @param kin The KinematicSurface input
  * @return the results of the main program
  **/
  static std::vector<std::vector<double> > solveKinSurface(Hyperplane &hyperplane, KinematicSurface &kin)
  {
    using namespace std;
    vector<vector<double> > result;

    // Need exactly 8 hyperplanes (4 left + 4 right); if not, robot geometry is unsupported
    if(hyperplane.h.size() < 8)
      return result;

    if(hyperplane.manifoldsUsed[0]<1 || hyperplane.manifoldsUsed[1]<1)
	  {
      //cout << "first manifold: " << hyperplane.manifoldsUsed[0] << endl;
      //cout << "second manifold: " << hyperplane.manifoldsUsed[1] << endl;
      result = solveForSpecialCase(hyperplane, kin); //hopefully our RP case does not end up here
    }
    else
    {
      result = solveForGeneralCase(hyperplane, kin);
    }
    SolveForAngles::sortResults(&result);
    return result;
  };

  /**
  * Calculation for wrist partitioned manipulators
  * @param hyperplane The hyperplane input
  * @param kin The KinematicSurface input
  * @return the results of the problem
  **/
  static vector<vector<double> > solveForSpecialCase(Hyperplane &hyperplane, KinematicSurface &kin)
  {
    vector<vector<double> > finalSolution;
    double complexPrecision = 1e-3;
    Matrix hypMatForWristPart(8,8,true);

    for(int i=0; i<hypMatForWristPart.numRows; i++)
    {
      for(int j=0; j<hypMatForWristPart.numCols; j++)
      {
        hypMatForWristPart.p_matrix[i][j] = hyperplane.h[i][j];
      }
    }
    //Jenkins Traub root finder --> no speedup
    //JenkinsTraubRoot jtr;
    //vector<Complex> wristRoots = jtr.calcRoots(Matrix::myDeterminant(hypMatForWristPart));

    //calculates all roots of the determinant of hypMatForWristPart
    vector<Complex> wristRoots = Matrix::qrCompanion(Matrix::balanceCompanionMatrix(Matrix::companionPoly(Matrix::myDeterminant(hypMatForWristPart))));


    vector<double> realRoots;
    //remove complex solutions
    for(size_t i=0; i<wristRoots.size(); i++)
    {
      if(fabs(wristRoots[i].imaginary) < complexPrecision)
        realRoots.push_back(wristRoots[i].real);
    }

    for(size_t i=0; i<realRoots.size(); i++)
    {
      vector<double> theta;
      theta.resize(6);
      vector<double> rValue = kin.substituteValues(realRoots[i]);
      double v1,v2,v3,v4,v5,v6;

      if(hyperplane.manifoldsUsed[0]==1)
      {
        v1=realRoots[i];
      }
      else
      {
        vector<Polynomial> p = h2_tc_v1(hyperplane.a);
        Polynomial tmpres(1,true);
        for(size_t j=0; j<p.size(); j++)
        {
          tmpres += p[j] * rValue[j];
        }
        v1=-tmpres.polyCoefficient[0].coefficient[0]/tmpres.polyCoefficient[1].coefficient[0];
      }

      if(hyperplane.manifoldsUsed[0]==2)
      {
        v2=realRoots[i];
      }
      else
      {
        vector<Polynomial> p = h2_tc_v2(hyperplane.a);
        Polynomial tmpres(1,true);
        for(size_t j=0; j<p.size(); j++)
        {
          tmpres += p[j] * rValue[j];
        }
        v2=-tmpres.polyCoefficient[0].coefficient[0]/tmpres.polyCoefficient[1].coefficient[0];
      }

      if(hyperplane.manifoldsUsed[0]==3)
      {
        v3=realRoots[i];
      }
      else
      {
        vector<Polynomial> p = h2_tc_v3(hyperplane.a);
        Polynomial tmpres(1,true);
        for(size_t j=0; j<p.size(); j++)
        {
          tmpres += p[j] * rValue[j];
        }
        v3=-tmpres.polyCoefficient[0].coefficient[0]/tmpres.polyCoefficient[1].coefficient[0];
      }

      vector<double> v6_value = new_wrist_sol1_v6(rValue, hyperplane.a);
      v6=v6_value[0];
      v5=new_wrist_sol1_v5(rValue,v6,hyperplane.a);
      v4=new_wrist_sol1_v4(rValue,v6,hyperplane.a);

      theta[0] = v1;
      theta[1] = v2;
      theta[2] = v3;
      theta[3] = v4;
      theta[4] = v5;
      theta[5] = v6;

      for(size_t j=0; j<6; j++)
      {
        theta[j] = Input::inverseHalfAngleSubs(theta[j])/Math::deg2rad;
      }
      finalSolution.push_back(theta);

      v6=v6_value[1];
      v5=new_wrist_sol1_v5(rValue,v6,hyperplane.a);
      v4=new_wrist_sol1_v4(rValue,v6,hyperplane.a);

      theta[0] = v1;
      theta[1] = v2;
      theta[2] = v3;
      theta[3] = v4;
      theta[4] = v5;
      theta[5] = v6;

      for(size_t j=0; j<6; j++)
      {
        theta[j] = Input::inverseHalfAngleSubs(theta[j])/Math::deg2rad;
      }
      finalSolution.push_back(theta);
    }
    return finalSolution;
  };

  /**
  * for general 6R manipulator
  * @param hyperplane The hyperplane input
  * @param kin The KinematicSurface input
  * @return the results of the problem, in degrees (for rotational angles)!
  **/
  static vector<vector<double> > solveForGeneralCase(Hyperplane &hyperplane, KinematicSurface &kin)
  {
    vector<vector<double> > finalSolution;
    //study quadric equation
    //instead try to use hyperplane rather than kin.. (this will allow solving for situation when we are in a kinsing.)
    Polynomial e1 = kin.r[0] * kin.r[4]
                    + kin.r[1] * kin.r[5]
                    + kin.r[2] * kin.r[6]
                    + kin.r[3] * kin.r[7];

    //last hyperplane h8
    Polynomial e2 = ((hyperplane.h[7])[0]) * kin.r[0] +
                    ((hyperplane.h[7])[1]) * kin.r[1] +
                    ((hyperplane.h[7])[2]) * kin.r[2] +
                    ((hyperplane.h[7])[3]) * kin.r[3] +
                    ((hyperplane.h[7])[4]) * kin.r[4] +
                    ((hyperplane.h[7])[5]) * kin.r[5] +
                    ((hyperplane.h[7])[6]) * kin.r[6] +
                    ((hyperplane.h[7])[7]) * kin.r[7];

    //jcapco todo: we can make some nice shortcuts if the desired value is prismatic value (not half-tangent from rotational joint)
    vector<pair<double,double> > realRoots = calculateRoots(e1,e2,kin.r[0],kin.r[1]);
    int numberofSolutions = int(realRoots.size());
    //cout << "Number of dSolutions: " << numberofSolutions << endl;

    vector<vector<double> > final_v = findAllAngles(realRoots, hyperplane, kin);
    for(int i=0; i<numberofSolutions; i++)
    {
      for(int j=0; j<6; j++)
      {
        if (hyperplane.a.rots[j])
          (final_v[i])[j] = Input::inverseHalfAngleSubs((final_v[i])[j])/Math::deg2rad;
      }
    }

    //checkSolutions with forwardKinematics
    double fwdPrecision = 1e-3;
    for(int i=0; i<numberofSolutions; i++)
    {
      vector<double> theta;
      theta.push_back((final_v[i])[0]);
      theta.push_back((final_v[i])[1]);
      theta.push_back((final_v[i])[2]);
      theta.push_back((final_v[i])[3]);
      theta.push_back((final_v[i])[4]);
      theta.push_back((final_v[i])[5]);

      ForwardKinematics fwd(hyperplane.a, theta);
      double maxElement = fwd.maxElement(hyperplane);
      if(maxElement<fwdPrecision)
        finalSolution.push_back(theta);
    }
    return finalSolution;
  };

  /**
  * calculate remaining 4 angles
  * @param va Vector of known roots
  * @param hyperplane The hyperplane input
  * @param kin The KinematicSurface input
  * @return A double array of all the 6 angles for all solutions
  **/
  static vector<vector<double> > findAllAngles(vector<pair<double,double> > va, Hyperplane &hyperplane, KinematicSurface &kin)
  {
    int numberOfSolutions = int(va.size());
    vector<vector<double> > final_v;
    int lo[3] = {0,0,0}; //l_order
    
    //order of solving depending on the manifold used
    switch (hyperplane.manifoldsUsed[0])
    {
    case 1: 
      { lo[0]=1-1; lo[1]=3-1; lo[2]=2-1; break; }
    case 2:
      { lo[0]=2-1; lo[1]=3-1; lo[2]=1-1; break; }
    case 3:
      { lo[0]=3-1; lo[1]=1-1; lo[2]=2-1; break; }
    }
  
    std::vector<LibHUPF::Polynomial> (*hplanes[3])(LibHUPF::Input&); //={0,0,0};
    if (hyperplane.left_3r == RRP_)
    {      
      hplanes[0]=RRP::h1_tc_v1;
      hplanes[1]=RRP::h1_tc_v2;
      hplanes[2]=RRP::h1_tc_d3;
    }
    else if (hyperplane.left_3r == RRR_)
    {
      //todo: add namespace to RRR
      hplanes[0]=LibHUPF::h1_tc_v1;
      hplanes[1]=LibHUPF::h1_tc_v2;
      hplanes[2]=LibHUPF::h1_tc_v3;
    }
    
    //jcapco now: correct hyperplanes for RRP
    //no parallelisation potential --> doesnt worth
    for(int i=0; i<numberOfSolutions; i++)
    {
      vector<double> tmp;
      final_v.push_back(tmp);
      vector<double> rValue = kin.substituteValues(va[i].first, va[i].second);

      //we can use the RRR hyperplanes in RRP computations!
      if(hyperplane.left_3r==RRP_ && hyperplane.manifoldsUsed[0]==3) 
      {
        hyperplane.a.d[2]=va[i].first;
        //cout << "Testing solved d3: " << va[i].first << endl;
      }
      
      double vsol[3] ={0,0,0};
      
      vsol[lo[0]]=va[i].first;
      
      for (size_t k=1; k<3;k++)
      {
        vector<Polynomial> p = hplanes[lo[k]](hyperplane.a);
        Polynomial tempresult(1,true);
        for(size_t j=0; j<p.size(); j++)
          tempresult = tempresult+p[j]*rValue[j];
        vsol[lo[k]] = -tempresult.polyCoefficient[0].coefficient[0]/tempresult.polyCoefficient[1].coefficient[0];
      }
      
      for (size_t j=0;j<3;++j)
        final_v[i].push_back(vsol[j]);

      //for v4, in thesis this is 1/v4
      if(hyperplane.manifoldsUsed[1]==4)
        final_v[i].push_back(va[i].second);
      else
      {
        vector<Polynomial> p = h2_v4q(hyperplane.a);
        Polynomial tempresult(1,true);
        for(size_t j=0; j<p.size(); j++)
          tempresult = tempresult+p[j]*rValue[j];
        (final_v[i]).push_back(-tempresult.polyCoefficient[0].coefficient[0] / tempresult.polyCoefficient[1].coefficient[0]);
      }
      (final_v[i])[3] = 1/(final_v[i])[3];

      //for v5, this is v5' in HUPF thesis and is replaced by v5'=-v5
      if(hyperplane.manifoldsUsed[1]==5)
        final_v[i].push_back(va[i].second);
      else
      {
        vector<Polynomial> p = h2_v5q(hyperplane.a);
        Polynomial tempresult(1,true);
        for(size_t j=0; j<p.size(); j++)
          tempresult = tempresult+p[j]*rValue[j];
        (final_v[i]).push_back(-tempresult.polyCoefficient[0].coefficient[0] / tempresult.polyCoefficient[1].coefficient[0]);
      }
      (final_v[i])[4] = -(final_v[i])[4];

      //for v6, in thesis this is 1/v6
      if(hyperplane.manifoldsUsed[1]==6)
        final_v[i].push_back(va[i].second);
      else
      {
        vector<Polynomial> p = h2_v6q(hyperplane.a);
        Polynomial tempresult(1,true);
        for(size_t j=0; j<p.size(); j++)
          tempresult = tempresult+p[j]*rValue[j];
        (final_v[i]).push_back(-tempresult.polyCoefficient[0].coefficient[0] / tempresult.polyCoefficient[1].coefficient[0]);
      }
      (final_v[i])[5] = 1/(final_v[i])[5];
    }
    return final_v;
  };

  /**
  * Returns the real roots for polynomials
  * @param e1 1st Polynomial for real Roots
  * @param e2 2nd Polynomial for real Roots
  * @param r0 1st Polynomial for wrong Roots
  * @param r1 2nd Polynomial for wrong Roots
  * @return the roots of the polynomials
  **/
  static vector<pair<double,double> > calculateRoots(Polynomial &e1, Polynomial &e2, Polynomial &r0, Polynomial &r1)
  {
    double comparisonPrecision=1e-3;
    double imaginaryPrecision=1e-3;
    JenkinsTraubRoot jtr;


    //doesnt worth parallelizing (1 thread 90%, 2nd thread 10%)

    vector<double> realRoots;
    vector<double> wrongRealRoots;

    //use companion matrix QR eigenvalue method to calculate right roots
    //(more robust than Jenkins-Traub for clustered or nearly-degenerate roots)
    Polynomial realPol=resultantMethod(e1,e2);
    poly_capture_ref() = realPol.coefficient; // cache degree-56 resultant for external access
    vector<Complex> roots = Matrix::qrCompanion(Matrix::balanceCompanionMatrix(Matrix::companionPoly(realPol)));

    for(size_t i=0; i<roots.size(); i++)
    {
      if(fabs(roots[i].imaginary) < imaginaryPrecision)
        realRoots.push_back(roots[i].real);
    }
    //use jenkins traub to calculate wrong roots
    Polynomial falsePol=resultantMethod(r0,r1);
    vector<Complex> wrongRoots = jtr.calcRoots(falsePol);

    for(size_t i=0; i<wrongRoots.size(); i++)
    {
      if(fabs(wrongRoots[i].imaginary) < imaginaryPrecision)
        wrongRealRoots.push_back(wrongRoots[i].real);
    }

    //remove wrong roots
    for(size_t i=0; i<wrongRealRoots.size(); i++)
    {
      size_t j=0;
      while(j<realRoots.size())
      {
        if(fabs(wrongRealRoots[i] - realRoots[j]) < comparisonPrecision)
          realRoots.erase(realRoots.begin()+j);
        else
          j++;
      }
    }

    vector<pair<double,double> > finalRoots;

    for(size_t i=0; i<realRoots.size(); i++)
    {
      Polynomial newE1 = e1.subs(realRoots[i]);
      vector<Complex> rootsE1 = jtr.calcRoots(newE1);//Matrix::qrCompanion(Matrix::balanceCompanionMatrix(Matrix::companionPoly(newE1)));
      //printf("%d ", rootsE1.size());

      Polynomial newE2 = e2.subs(realRoots[i]);
      vector<Complex> rootsE2 = jtr.calcRoots(newE2); //Matrix::qrCompanion(Matrix::balanceCompanionMatrix(Matrix::companionPoly(newE2)));
      //printf("%d ", rootsE2.size());

      for(size_t j=0; j<rootsE1.size(); j++)
        if(fabs(rootsE1[j].imaginary) < imaginaryPrecision)
          for(size_t k=0; k<rootsE2.size(); k++)
            if(fabs(rootsE2[k].imaginary) < imaginaryPrecision)
              if(fabs(rootsE1[j].real - rootsE2[k].real) < comparisonPrecision)
                finalRoots.push_back(pair<double,double>(realRoots[i],rootsE2[k].real));
    }

    size_t i2;
    bool allChecked;

    if(finalRoots.size()>0)
    {
      for(size_t i1=0; i1<finalRoots.size()-1; i1++)
      {
        i2=i1+1;
        allChecked=false;
        while(!allChecked)
        {
          if(fabs(finalRoots[i1].first-finalRoots[i2].first) < comparisonPrecision && fabs(finalRoots[i1].second-finalRoots[i2].second) < comparisonPrecision)
          {
            finalRoots.erase(finalRoots.begin()+i2);
            i2=i1;
          }
          i2++;
          if(i2 >= finalRoots.size())
            allChecked=true;
        }
      }
    }

    return finalRoots;
  };



  /**
  * Returns the determinant of Resultant Matrix created using two input polynomials.
  * @param e1 First bivariate Polynomial
  * @param e2 Second bivariate Polynomial
  * @return resultant Polynomial
  **/
  static Polynomial resultantMethod(Polynomial &e1, Polynomial &e2)
  {
    Matrix resultantMatrix(e1.degree+e2.degree, e1.degree+e2.degree, true);
    for(int i=0; i<e2.degree; i++)
      for(int j=0; j<=e1.degree; j++)
      {
        resultantMatrix.p_matrix[j+i][i] = e1.polyCoefficient[j];
      }

    for(int i=0; i<e1.degree; i++)
      for(int j=0; j<=e2.degree; j++)
      {
        resultantMatrix.p_matrix[j+i][e2.degree+i] = e2.polyCoefficient[j];
      }
    return Matrix::myDeterminant(resultantMatrix);
  };
  /**
   * Calculates degree of polynomial determinant.
   * @param mat Matrix to get determinant
   * @return the degree
   **/
  static double getMaxDegree(Matrix &mat)
  {
    double subtract_me;
    double res;
    if(mat.numRows != mat.numCols)
      perror("Only for square matrices.");
    if(mat.numRows==1 && mat.numCols==1)
      return mat.m_matrix[0][0];
    if(mat.numRows==2 && mat.numCols==2)
      return max(mat.m_matrix[0][0]+mat.m_matrix[1][1], mat.m_matrix[0][1]+mat.m_matrix[1][0]);
    for(int i = 1; i<mat.numRows; i++)
    {
      for(int j = 1; j<mat.numRows; j++)
      {
        subtract_me = mat.m_matrix[i][0] + mat.m_matrix[0][j];
        mat.m_matrix[i][j] = max(mat.m_matrix[i][j] + mat.m_matrix[0][0], subtract_me);
      }
    }
    Matrix subMat(mat.numRows-1,mat.numCols-1);
    for(int i=0; i<subMat.numRows; i++)
      for(int j=0; j<subMat.numCols; j++)
        subMat.m_matrix[i][j] = mat.m_matrix[i+1][j+1];
    res = getMaxDegree(subMat);
    for(int i=1; i<=mat.numRows-2; i++)
    {
      res = res - mat.m_matrix[0][0];
    }
    return res;
  };

  /**
  * Method of quickSort to part solutions in two parts
  * @param vals the vector of values, each value contains multiple angles
  * @param left the index of the most left element
  * @param right the index of the most right element
  * @param pivot the index of the pivot element
  * @param sortNr the index of the angle to use for sorting
  * @return the index of the new pivotelement
  **/
  static int partQuick(vector<vector<double> > *vals,int left, int right, int pivot, int sortNr)
  {
    double pivotVal=((*vals)[pivot])[sortNr];
    vector<double> h=(*vals)[pivot];
    (*vals)[pivot]=(*vals)[right];
    (*vals)[right]=h;

    int storeIndex=left;

    for(int i=left; i<right; i++)
      if(((*vals)[i])[sortNr]<pivotVal)
      {
        h=(*vals)[i];
        (*vals)[i]=(*vals)[storeIndex];
        (*vals)[storeIndex]=h;
        storeIndex++;
      }
    h=(*vals)[storeIndex];
    (*vals)[storeIndex]=(*vals)[right];
    (*vals)[right]=h;
    return storeIndex;
  };

  /**
  * QuickSort Method for the solutions
  * @param array vector of solutions, each solution contains 6 angles
  * @param left the index of the most left element
  * @param right the index of the most right element
  * @param sortNr the index of th eangle to use for sorting
  **/
  static void quickSort(vector<vector<double> >* array, int left, int right, int sortNr)
  {
    if(left<right)
    {
      int pivotIndex=partQuick(array,left,right,(left+right)/2,sortNr);
      quickSort(array,left,pivotIndex-1,sortNr);
      quickSort(array,pivotIndex+1,right,sortNr);
    }
  };

  /**
  * Sorting all results ascending
  * @param results the results to sort
  **/
  static void sortResults(vector<vector<double> >* results)
  {
    quickSort(results,0,int((*results).size())-1,0);		//sort by first angle

    for (int k=1; k<6; k++)
    {
      bool distGT0=false;
      double beforeVal=-10.0;
      int dist=0;
      for (int i=0; i<int((*results).size()); i++)
      {
        if (fabs(((*results)[i])[k-1]-beforeVal)<0.0000001)
          dist++;
        else
        {
          if(dist>0)
          {
            quickSort(results,i-dist-1,i-1,k);
            distGT0=true;
          }
          beforeVal=((*results)[i])[k-1];
          dist=0;
        }
      }
      if (!distGT0)		//in case of dist never gt0
        break;
    }
  };
};

} //namespace