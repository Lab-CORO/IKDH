/**
 *     Provides functionality to handle polymonial calculations.
 *     A Polynomial can have real or polynomial coefficients.
 **/
#pragma once

#include <vector>
//#include <iostream>
//#if defined(_WIN32) || defined(_WIN64)
#include<math.h>
//#endif

namespace LibHUPF{

class Polynomial
{

public:
  int degree; //!< polynomial degree
  std::vector<Polynomial> polyCoefficient; //nonempty -> bivariate polynomial
  std::vector<double> coefficient; //noempty -> univariate polynomial
  bool realCoeffType; //!< true if real coefficients polynomial, false if polynomial coefficients polynomial

  /**
   * Polynomial with degree -1 and polynomial coefficients, i.e. the zero polynomial
   * */
  Polynomial()
  {
    degree=-1;
    realCoeffType=false;
  };

  /**
   * Polynomial with degree n and real coefficients.
   * @param deg Degree of the polynomial
   * */
  Polynomial(int deg) : coefficient(deg+1)
  {
    degree=deg;
    realCoeffType=true;
  };

  /**
   * Creates a Polynomial with real coefficients from a given vector<double>.
   * @param coeff vector<double> filled with coefficients in ascending order.
   * */
  Polynomial(std::vector<double> coeff) : coefficient(coeff)
  {
    degree = int(coeff.size())-1;
    realCoeffType = true;
  };

  /**
   * Creates a polynomial with real coefficients from a given double pointer.
   * The degree of the polynomial must be given.
   * @param coeff double pointer
   * @param size degree of the polynomial
   * */
  Polynomial(double* coeff, int size)
  {
    degree = size;
    for(int i=0; i<=size; i++)
    {
      coefficient.push_back(coeff[i]);
    }
    realCoeffType = true;
  };

  /**
   * Creates a polynomial of degree 1 with real coefficients. This constructor is
   * used in the emplace_back blocks in the Hyperplane class.
   * @param v1 first coefficient
   * @param v2 second coefficient
   * */
  Polynomial(double v1, double v2)
  {
    //make sure we don't have unecessary zero polynomials!, this is used in hyperplane initialization!
    degree = 1;
    if (fabs(v2) < 1e-40)
    {
      degree--;
      if (fabs(v1) < 1e-40) degree--;
    }
    
    if (degree > -1)
      coefficient.push_back(v1);
    if (degree > 0)
      coefficient.push_back(v2);
    realCoeffType = true;
  };

  /**
   * Creates a polynomial of degree 0 with real coefficients. This constructor is
   * used in the emplace_back blocks in the Hyperplane class.
   * @param v1 first coefficient
   * */
  Polynomial(double v1)
  {
    degree=0;
    realCoeffType=true;
    coefficient.push_back(v1);
  };

  /**
   * Creates a polynomial of degree n with polynomial coefficients.
   * @param deg degree of the polynomial
   * @param poly true --> polynomial coefficient type
   * */
  Polynomial(int deg, bool poly) : polyCoefficient(deg+1)
  {
    if(poly)
    {
      degree=deg;
      realCoeffType=false;
    }
  };

  /**
   * Creates a polynomial with polynomial coefficients. The coefficients
   * are given in form of a vector<Polynomial>.
   * @param coeff polynomial coefficient vector
   * */
  Polynomial(std::vector<Polynomial> coeff) : polyCoefficient(coeff)
  {
    degree = int(coeff.size())-1;
    realCoeffType=false;
  }

  /**
   * Copy constructor
   * */
  Polynomial(const Polynomial &p) : coefficient(p.coefficient)
  {
    degree = p.degree;
    if(!p.realCoeffType)
    {
      polyCoefficient = p.polyCoefficient;
      realCoeffType = false;
    }
    else
    {
      //coefficient(p.coefficient);
      realCoeffType = true;
    }
  };

  /**
   * Get real coefficient value at position n
   * @param n position
   * @return real coefficient at position n
   * */
  double get(int n)
  {
    return coefficient[n];
  };

  /**
   * Set real coefficient value at position n
   * @param n position
   * @param v value
   * */
  void set(int n, double v)
  {
    coefficient[n] = v;
  };

  /**
   * Get polynomial coefficient value at position n
   * @param n position
   * @return polynomial coefficient at position n
   * */
  Polynomial getElement(int n)
  {
    return polyCoefficient[n];
  };

  /**
   * Set polynomial coefficient value at position n
   * @param n position
   * @param p polynomial coefficient value
   * */
  void setElement(int n, Polynomial p)
  {
    polyCoefficient[n] = p;
  };

  /**
   * negate a given Polynomial
   * @param p polynomial
   * @return negated polynomial
   * */
  static Polynomial negPoly(Polynomial p)
  {
    if (p.realCoeffType)
    {
      for (int i=0; i<=p.degree; i++)
        p.coefficient[i]=-p.coefficient[i];
    }
    else
    {
      for (int i=0; i<=p.degree; i++)
        p.polyCoefficient[i]=negPoly(p.polyCoefficient[i]);
    }
    return p;
  };

  /**
   * substitute a value in a polynomial with polynomial coefficients.
   * The polynomial coefficients are evaluated and the resulting real coefficient
   * polynomial is returned.
   * @param d number for which coefficient the polynomial need to be evaluated
   * @return substituted polynomial
   * */
  Polynomial subs(double d)
  {
    if (realCoeffType)
    {
      perror("Subs possible only for polynomials having poly type coeffecients");
      exit(EXIT_FAILURE);
    }
    Polynomial res(degree);
    for (int i = 0; i <= degree; i++)
      res.coefficient[i] = polyCoefficient[i].eval(d);
    return res;
  };

  /**
   * Evaluate the value of the polynomial by substituting with d.
   * @param d value
   * @return resulting value
   * */
  double eval(double d)
  {
    if (!realCoeffType)
    {
      printf("Only polynomial with Real coefficients allowed");
      exit(EXIT_FAILURE);
    }
    double total = 0.0;
    for (int i = 0; i <= degree; i++)
      total += coefficient[i] * pow(d, i);
    return total;
  };

//todo: remove the operators from class method
  /**
   * minus operator (Polynomial-Polynomial)
   * */
  Polynomial& operator-(const Polynomial &b)
  {
    if (degree==-1 && b.degree != -1)
      return (*(new Polynomial(negPoly(b))));
    if (degree!=-1 && b.degree == -1)
      return *this;

    Polynomial *res;
    int resDegree=0;

    if (degree >= b.degree)
      resDegree = degree;
    else
      resDegree = b.degree;

    if (realCoeffType && b.realCoeffType)
      res = new Polynomial(resDegree);
    else
      res = new Polynomial(resDegree,true);


    if(realCoeffType && b.realCoeffType)
    {
      double A=0;
      double B=0;
      for(int i=resDegree; i>=0; i--)
      {
        if(i<=degree)
          A=coefficient[i];

        if(i<=b.degree)
          B=b.coefficient[i];
        res->coefficient[i]=A-B;
      }
    }


    if (!realCoeffType && !b.realCoeffType)
    {
      Polynomial polA;
      Polynomial polB;

      for (int i = 0; i <= resDegree; i++)
      {
        if (i > degree)
          polA.degree=-1;
        else
          polA = polyCoefficient[i];

        if (i > b.degree)
          polB.degree=-1;
        else
          polB = b.polyCoefficient[i];

        if (polA.degree==-1 && polB.degree!=-1)
          res->polyCoefficient[i]=negPoly(polB);
        if (degree!=-1 && b.degree == -1)
          res->polyCoefficient[i] = polA;
        if (degree!=-1 && b.degree!=-1)
          res->polyCoefficient[i]=polA-polB;
      }
    }

    if (realCoeffType && !b.realCoeffType)
    {
      double helpA = 0;
      Polynomial helpB;
      for (int i = 0; i <= resDegree; i++)
      {
        if (i > degree)
          helpA = 0;
        else
          helpA = coefficient[i];

        if (i > b.degree)
          helpB.degree=-1;
        else
          helpB = b.polyCoefficient[i];
        res->polyCoefficient[i]=negPoly(helpB)+helpA;
      }

    }

    if (!realCoeffType && b.realCoeffType)
    {
      Polynomial polA;
      double helpB = 0;

      for (int i = 0; i <= resDegree; i++)
      {
        if (i > degree)
          polA.degree=-1;
        else
          polA = polyCoefficient[i];

        if (i > b.degree)
          helpB = 0;
        else
          helpB = b.coefficient[i];
        res->polyCoefficient[i]=polA-helpB;
      }

    }
    return *res;
  };

  /**
   * minus equals operator (Polynomial-=Polynomial)
   * */
  Polynomial& operator-=(const Polynomial &b)
  {
    if (degree==-1 && b.degree != -1)
    {
      Polynomial c(negPoly(b));
      degree=c.degree;
      realCoeffType=c.realCoeffType;
      if(realCoeffType)
        coefficient=c.coefficient;
      else
        polyCoefficient=c.polyCoefficient;
      return (*this);
    }
    if (degree!=-1 && b.degree == -1)
      return *this;
    int resDegree=0;

    if (degree >= b.degree)
      resDegree = degree;
    else
      resDegree = b.degree;

    if (realCoeffType && b.realCoeffType)
    {
      coefficient.resize(resDegree+1);
      //delete higher degrees in this
      if(b.degree > degree)
      {
        for(int i=degree+1; i<=b.degree; i++)
          coefficient[i]=0.0;
      }
      degree = resDegree;
    }
    else
    {
      degree = resDegree;
      polyCoefficient.resize(resDegree+1);
    }

    if (realCoeffType && b.realCoeffType)
    {
      for(int i=0; i<=resDegree; i++)
      {
        if(i>b.degree)
          break;
        if(i>degree)
          coefficient[i]=-b.coefficient[i];
        else
          coefficient[i]-=b.coefficient[i];
      }
    }

    if (!realCoeffType && !b.realCoeffType)
    {
      for(int i=0; i<=resDegree; i++)
      {
        if(i>b.degree)
          break;
        if(i>degree)
          polyCoefficient[i]=negPoly(b.polyCoefficient[i]);
        else
          polyCoefficient[i]-=b.polyCoefficient[i];
      }
    }

    if (realCoeffType && !b.realCoeffType)
    {
      double helpA = 0;
      Polynomial helpB;
      for (int i = 0; i <= resDegree; i++)
      {
        if (i > degree)
          helpA = 0;
        else
          helpA = coefficient[i];

        if (i > b.degree)
          helpB.degree=-1;
        else
          helpB = b.polyCoefficient[i];
        polyCoefficient[i]=negPoly(helpB)+helpA;
      }

    }

    if (!realCoeffType && b.realCoeffType)
    {
      Polynomial polA;
      double helpB = 0;

      for (int i = 0; i <= resDegree; i++)
      {
        if (i > degree)
          polA.degree=-1;
        else
          polA = polyCoefficient[i];

        if (i > b.degree)
          helpB = 0;
        else
          helpB = b.coefficient[i];
        polyCoefficient[i]=polA-helpB;
      }

    }
    return *this;
  };

  /**
   * plus equals operator (Polynomial+=Polynomial)
   * */
  Polynomial& operator+=(const Polynomial &b)
  {
    if(degree==-1 && b.degree!=-1)
    {
      degree=b.degree;
      realCoeffType=b.realCoeffType;
      if(b.realCoeffType)
        coefficient=b.coefficient;
      else
        polyCoefficient=b.polyCoefficient;
      return (*this);
    }
    if(degree!=-1 && b.degree==-1)
      return (*this);
    int rDegree =0;
    if(degree >= b.degree)
      rDegree=degree;
    else
      rDegree=b.degree;

    if(realCoeffType && b.realCoeffType)
    {
      coefficient.resize(rDegree+1);
      //delete higher degrees in this
      if(b.degree > degree)
      {
        for(int i=degree+1; i<=b.degree; i++)
          coefficient[i]=0.0;
      }
      degree = rDegree;
    }
    else
    {
      degree = rDegree;
      polyCoefficient.resize(rDegree+1);
    }

    if(realCoeffType && b.realCoeffType)
    {
      for(int i=0; i<=rDegree; i++)
      {
        if(i>b.degree)
          break;
        if(i>degree)
          coefficient[i] = b.coefficient[i];
        else
          coefficient[i] += b.coefficient[i];
      }
    }

    if(!realCoeffType && !b.realCoeffType)
    {
      for(int i=0; i<=rDegree; i++)
      {
        if(i>b.degree)
          break;
        if(i>degree)
          polyCoefficient[i]=b.polyCoefficient[i];
        else
          polyCoefficient[i]+=b.polyCoefficient[i];
      }
    }

    if(realCoeffType && !b.realCoeffType)
    {
      realCoeffType=false;
      double A=0;
      Polynomial B;
      for(int i=0; i<=rDegree; i++)
      {
        if(i>degree)
          A=0;
        else
          A=coefficient[i];
        if(i>b.degree)
          B.degree=-1;
        else
          B=b.polyCoefficient[i];
        polyCoefficient[i]=B+A;
      }
    }

    if(!realCoeffType && b.realCoeffType)
    {
      double B=0;
      for(int i=rDegree; i>=0; i--)
      {
        if(i<=b.degree)
          B=b.coefficient[i];
        if(i>degree)
          polyCoefficient[i]=B;
        else
          polyCoefficient[i]+=B;
      }
    }
    return *this;
  };

  /**
   * plus operator (Polynomial+Polynomial)
   * */
  Polynomial operator+(const Polynomial &b)
  {
    Polynomial result;
    if(degree==-1 && b.degree!=-1)
      return Polynomial(b);
    if(degree!=-1 && b.degree==-1)
      return (*this);
    int rDegree =0;
    if(degree >= b.degree)
      rDegree=degree;
    else
      rDegree=b.degree;

    if(realCoeffType && b.realCoeffType)
      result= Polynomial((int)rDegree);
    else
      result= Polynomial(rDegree, true);

    if(realCoeffType && b.realCoeffType)
    {
      double A=0;
      double B=0;
      for(int i=rDegree; i>=0; i--)
      {
        if(i<=degree)
          A=coefficient[i];

        if(i<=b.degree)
          B=b.coefficient[i];
        result.coefficient[i]=A+B;
      }
    }

    if(!realCoeffType && !b.realCoeffType)
    {
      Polynomial *A=NULL;
      Polynomial B;
      for(int i=rDegree; i>=0; i--)
      {
        if(i<=degree)
          A=&polyCoefficient[i];

        if(i<=b.degree)
          B=Polynomial(b.polyCoefficient[i]);

        if(A==NULL)
          result.polyCoefficient[i]=b.polyCoefficient[i];
        if(B.degree==-1)
          result.polyCoefficient[i]=polyCoefficient[i];
        if(A!=NULL && B.degree!=-1)
          result.polyCoefficient[i]=(*A)+(B);
      }
    }

    if(realCoeffType && !b.realCoeffType)
    {
      double A=0;
      Polynomial B;
      for(int i=0; i<=rDegree; i++)
      {
        if(i>degree)
          A=0;
        else
          A=coefficient[i];
        if(i>b.degree)
          B.degree=-1;
        else
          B=b.polyCoefficient[i];
        result.polyCoefficient[i]=B+A;
      }
    }

    if(!realCoeffType && b.realCoeffType)
    {
      double B=0;
      for(int i=0; i<=rDegree; i++)
      {
        if(i>b.degree)
          B=0;
        else
          B=b.coefficient[i];
        if(i>degree)
          result.polyCoefficient[i]=B;
        else
          result.polyCoefficient[i]=polyCoefficient[i]+B;
      }
    }
    return result;
  };

  /**
   * plus operator (Polynomial+double)
   * */
  Polynomial operator+(double a)
  {
    using namespace std;
    Polynomial res;
    if(degree==-1)
    {
      vector<double> help;
      help.push_back(a);
      return Polynomial(help);
    }
    else
    {
      if(realCoeffType)
      {
        res=*this;
        res.coefficient[0]+=a;
      }
      if(!realCoeffType)
      {
        res=*this;
        res.polyCoefficient[0] = polyCoefficient[0]+a;
      }
    }
    return res;
  };


  /**
   * minus operator (Polynomial-double)
   * */
  Polynomial& operator-(double b)
  {
    Polynomial *res;

    if (degree==-1)
    {
      res=new Polynomial(-b);
    }
    else
    {
      if (realCoeffType)
      {
        res=this;
        res->coefficient[0]-=b;
      }
      if (!realCoeffType)
      {
        res=this;
        res->polyCoefficient[0]=polyCoefficient[0]-b;
      }

    }
    return *res;

  };

  /**
   * multiplication operator (Polynomial*Polynomial)
   * */
  Polynomial operator*(const Polynomial &b)
  {
    Polynomial res;

    if (degree==-1 || b.degree==-1)
    {
      return res;
    }
    else
    {
      if (realCoeffType && b.realCoeffType)
      {
        res = Polynomial(degree+b.degree);
        for (int i=degree; i>=0; i--)
          for (int j=b.degree; j>=0; j--)
            res.coefficient[i+j]+=coefficient[i]*b.coefficient[j];
        return res;
      }
      else if (!realCoeffType && !b.realCoeffType)
      {
        res = Polynomial(degree+b.degree,true);
        for (int i=degree; i>=0; i--)
          for (int j=b.degree; j>=0; j--)
            res.polyCoefficient[i+j]+=polyCoefficient[i]*b.polyCoefficient[j];
        return res;
      }
      else if (!realCoeffType && b.realCoeffType)
      {
        res = Polynomial(degree+b.degree,true);
        for (int i=degree; i>=0; i--)
          for (int j=b.degree; j>=0; j--)
          {
            res.polyCoefficient[i+j]+=(polyCoefficient[i]*b.coefficient[j]);
          }
        return res;
      }
      else if (!b.realCoeffType && realCoeffType)
      {
        res = Polynomial(degree+b.degree,true);
        for (int i=degree; i>=0; i--)
          for (int j=b.degree; j>=0; j--)
          {
            Polynomial tmp1(b.polyCoefficient[j]);
            res.polyCoefficient[i+j]+=(tmp1*coefficient[i]);
          }
        return res;
      }
      else return res;
    }
  };

  /**
   * multiplication equals operator (Polynomial*=Polynomial)
   * */
  Polynomial& operator*=(const Polynomial &b)
  {
    if (realCoeffType && b.realCoeffType)
    {
      for (int i=degree; i>=0; i--)
        for (int j=b.degree; j>=0; j--)
          coefficient[i+j]+=coefficient[i]*b.coefficient[j];
    }
    else    if (!realCoeffType && !b.realCoeffType)
    {
      for (int i=degree; i>=0; i--)
        for (int j=b.degree; j>=0; j--)
          polyCoefficient[i+j]+=polyCoefficient[i]*b.polyCoefficient[j];
    }
    else    if (!realCoeffType && b.realCoeffType)
    {
      for (int i=degree; i>=0; i--)
        for (int j=b.degree; j>=0; j--)
        {
          polyCoefficient[i+j]+=(polyCoefficient[i]*b.coefficient[j]);
        }
    }
    else if (!b.realCoeffType && realCoeffType)
    {
      for (int i=degree; i>=0; i--)
        for (int j=b.degree; j>=0; j--)
        {
          Polynomial tmp1(b.polyCoefficient[j]);
          polyCoefficient[i+j]+=(tmp1*coefficient[i]);
        }
    }
    return *this;
  };

  /**
   * multiplication operator (Polynomial*double)
   * */
  Polynomial operator*(const double b)
  {
    Polynomial res;
    if (degree!=-1)
    {
      if (realCoeffType)
      {
        res=*this;
        for (int i=0; i<=res.degree; i++)
          res.coefficient[i]*=b;
      }
      if (!realCoeffType)
      {
        res=*this;
        for (int i=0; i<=res.degree; i++)
          res.polyCoefficient[i]*=b;
      }
    }
    return res;
  };

  /**
   * division operator (Polynomial/d)
   * */
  Polynomial operator/(double d)
  {
    Polynomial res(*this);

    if (degree==-1)
      res.degree=-1;
    else
    {
      if (realCoeffType)
      {
        for (int i=0; i<=res.degree; i++)
          res.coefficient[i]=res.coefficient[i]/d;
      }
      if (!realCoeffType)
      {
        for (int i=0; i<=res.degree; i++)
          res.polyCoefficient[i]=res.polyCoefficient[i]/d;
      }
    }

    return res;
  };

  /**
   * divison operator (Polynomial/Polynomial)
   * */
  //Note: This works only for polynomials with real coefficients
  Polynomial operator/(Polynomial b)
  {
    int degA=degree;
    int degB=b.degree;
    if (degB > degA) return Polynomial(0); //zero polynomial!

    int wasNormalized=0;
    Polynomial a1(*this);
    Polynomial b1=b;

    if (realCoeffType)
    {
      if (b1.coefficient[degB] != 1)
      {
        b1=b1/b1.coefficient[degB];
        wasNormalized=1;
      }

      b1=negPoly(b1);

      Polynomial result(degA-degB);

      for (int i = 0; i <= degA - degB; i++) // rows of division step.
      {
        result.coefficient[degA - degB - i] = a1.coefficient[degA - i];
        int index = degB - 1;

        for (int j = degA - i - 1; j >= degA - i - degB; j--) // columns of division step
        {
          a1.coefficient[j] = a1.coefficient[j] + b1.coefficient[index] * result.coefficient[degA - degB - i];
          index--;
        }
      }

      /* //we want division also with remainder.. we take only polynomial quotient
      for (int i = 0; i < degB; i++) // All these are remainder elements.
      {
        //cout << "i: " << i << endl;
        if (a1.coefficient[i] != 0)  //TODO: Replace with tolerance checking.
        {        
          //cout << "a1.coeff[i]: " << a1.coefficient[i] << endl;
          perror("Not Completely Divisible");
          exit(EXIT_FAILURE);
        }
      }
      */

      if (wasNormalized)
        result = result/b.coefficient[degB];

      return result;
    }
    else
    {
      b1=negPoly(b1);
      Polynomial result(degA-degB, true);

      for (int i = 0; i <= degA - degB; i++) // rows of division step.
      {
        result.polyCoefficient[degA - degB - i] = a1.polyCoefficient[degA - i]/b1.polyCoefficient[degB]; //assuming the division becomes exact!
        int index = degB - 1;

        for (int j = degA - i - 1; j >= degA - i - degB; j--) // columns of division step
        {
          a1.polyCoefficient[j] = a1.polyCoefficient[j] + b1.polyCoefficient[index] * result.polyCoefficient[degA - degB - i];
          index--;
        }
      }
      
      return result;
    }

  };

  /**
   * print polynomial coefficients.
   * only for testing and debugging purposes.
   * @param p polynomial
   * @param offset start at coefficient n
   * */
  static void printPolynomial(Polynomial p, int offset)
  {
    for (int i=0; i<offset; i++)
      printf("\t");

    printf("Polynomial { degree: %d, realCoeffType: %d\n",p.degree,p.realCoeffType);

    for (int i=0; i<offset; i++)
      printf("\t");

    if (p.realCoeffType)
    {
      printf("Coefficients: ");

      for (int i=0; i<=p.degree; i++)
      {
        printf("%f\n ",p.coefficient[i]);

        for (int k=0; k<offset; k++)
          printf("\t");
      }

      printf("}\n");

    }
    else
    {
      printf("PolyCoeffs:\n");

      for (int i=0; i<=p.degree; i++)
        printPolynomial(p.polyCoefficient[i],offset+1);

    }
  }

  void clean()
  {
    for (int i=degree; i>=0; --i)
    {
      if (fabs(coefficient[i]) < 1e-20)
      {        
        coefficient.pop_back();
        degree--;
      }
    }
  }
  /**
   * equals operator (Polynomial==Polynomial)
   * only for testing and debugging purposes
   * @param b polynomial
   * @return true if this polynomial equals polynomial b, otherwise false is returned
   * */
  bool operator==(const Polynomial &b) const
  {
    if(degree!=b.degree || realCoeffType!=b.realCoeffType)
      return false;
    if(degree==b.degree && realCoeffType==b.realCoeffType)
    {
      if(realCoeffType)
      {
        for(size_t i=0; i<coefficient.size(); i++)
          if((coefficient[i] - b.coefficient[i]) > 0.0001)
            return false;
      }
      else
      {
        for(size_t i=0; i<polyCoefficient.size(); i++)
          for(size_t j=0; j<polyCoefficient[i].coefficient.size(); j++)
            if(!(polyCoefficient[i]==b.polyCoefficient[i]))
              return false;
      }
    }
    return true;
  }

};

}
