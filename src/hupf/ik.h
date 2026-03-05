#pragma once

#include <hupf/libdef.h>
#include <vector>
#include <string>
//a class is constructed by Input
//use below to populate input vector
//vector<Input> IP = Input::parseFile("robot info");

//#include <Input.h>
//if IP.empty() then do nothing, otherwise
//for every end effector (12 doubles) compute the angles in radians
//output angles (double)

#ifndef _MSC_VER
#include <stdint.h>
typedef int8_t __int8;
#endif

namespace LibHUPF
{

class LIBHUPF_LIBRARY_INTERFACE ik_solver
{
private:
  void* _input;
public:
  ik_solver(const std::string& ik_file);
  //with DH-parameters as 6 arrays for distance, offset and twist respectively.
  ik_solver(double* a, double* d, double* theta, double* alpha, bool* rots);
  ~ik_solver();

  //end effector as an 16-valued double array, assuming row-major!
  std::vector< std::vector<double> > solve(double* ee);
};

}

#ifdef __cplusplus
extern "C"
{ 
#endif

LIBHUPF_LIBRARY_INTERFACE void* create_ik_solver(double* a, double* d, double* theta, double* alpha, __int8* rots);

LIBHUPF_LIBRARY_INTERFACE void destroy_ik_solver(void* iks);

LIBHUPF_LIBRARY_INTERFACE int solve_ik(void* iks, double* ee, double* ret);

#ifdef __cplusplus
}
#endif