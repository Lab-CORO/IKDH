/**
 * 	Provides functionality to handle complex number calculations.
 **/
#pragma once
#include <math.h>

namespace LibHUPF
{

class Complex
{
public:
  double real;  //!< real part
  double imaginary;  //!< imaginary part

  /**
   * Creates an empty complex number.
   * */
  Complex()
  {
    real=0.0;
    imaginary=0.0;
  };

  /**
   * Returns norm of a complex number.
   * @param a complex number a
   * @return norm of complex number
   * */
  static double abs(Complex a)
  {
    return sqrt(a.real*a.real+a.imaginary*a.imaginary);
  }

  /**
   * Creates a complex number an sets real and imaginary part to given values
   * @param rea real part
   * @param img imaginary part
   * */
  Complex(double rea,double img)
  {
    real=rea;
    imaginary=img;
  };

  /**
   * divison operator (Complex/Complex)
   * */
  Complex operator/(const Complex &b)
  {
    Complex res;
    res.real=(real*b.real+imaginary*b.imaginary)/(b.real*b.real+b.imaginary*b.imaginary);
    res.imaginary=(imaginary*b.real-real*b.imaginary)/(b.real*b.real+b.imaginary*b.imaginary);
    return res;
  }

  /**
   * multiplication operator (Complex*Complex)
   * */
  Complex operator*(const Complex &b)
  {
    Complex res;
    res.real=(real*b.real-imaginary*b.imaginary);
    res.imaginary=real*b.imaginary+imaginary*b.real;
    return res;
  }

  /**
   * multiplication operator (Complex*double)
   * */
  Complex operator*(const double &b)
  {
    Complex res;
    res.real=real*b;
    res.imaginary=imaginary*b;
    return res;
  }

  /**
   * multiplication operator (Complex*int)
   * */
  Complex operator*(const int &b)
  {
    Complex res;
    res.real=real*b;
    res.imaginary=imaginary*b;
    return res;
  }

  /**
   * plus operator (Complex+Complex)
   * */
  Complex operator+(const Complex &b)
  {
    Complex res;
    res.real=real+b.real;
    res.imaginary=(imaginary+b.imaginary);
    return res;
  }

  /**
   * plus operator (Complex+double)
   * */
  Complex operator+(const double &b)
  {
    Complex res;
    res.real=real+b;
    res.imaginary=(imaginary);
    return res;
  }

  /**
   * minus operator (Complex-Complex)
   * */
  Complex operator-(const Complex &b)
  {
    Complex res;
    res.real=real-b.real;
    res.imaginary=(imaginary-b.imaginary);
    return res;
  }

};

}