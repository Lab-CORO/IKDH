#include <hupf/ik.h>

#include <hupf/Input.h>
#include <hupf/Calculate.h>

#include <iostream>
namespace LibHUPF
{
  using namespace std;

  ik_solver::ik_solver(const string& ik_file)
  {
    //jcapco todo: refactor: create input without file
    _input = Input::parseFile(ik_file);
  }

  ik_solver::ik_solver(double* a, double* d, double* theta, double* alpha, bool* rots)
  {
    _input = getInput(a,d,theta,alpha,rots);
  }

  ik_solver::~ik_solver()
  {
    delete static_cast<Input*>(_input);
  }

  /*
    ee is a double array rep. end effector transformation
    with at most 12 entries (last row 0,0,0,1 is not necessary)
    output is the solution set (vector of many solutions in one given ee)
    it has the following 3x4 form (flat matrix with rows written first)
    R t , were R is the orthogonal 3x3 rotational matrix and t is the
    translation column vector
  */
  vector< vector<double> > ik_solver::solve(double* ee)
  {
    Input* inp = static_cast<Input*>(_input);

    //hupf has the form (1,0;t R)
    (inp->eePose).m_matrix[0][0]=1;
    (inp->eePose).m_matrix[0][1]=0;
    (inp->eePose).m_matrix[0][2]=0;
    (inp->eePose).m_matrix[0][3]=0;

    for (size_t i=0; i<3; ++i) //translation (//rows)
    {
      (inp->eePose).m_matrix[i+1][0]=ee[3+4*i];
      for (size_t j=0; j<3; ++j) //columns
        (inp->eePose).m_matrix[i+1][j+1]=ee[j+4*i];
    }

    vector<double> times; //profiling

    return Calculate::InverseKin(*inp, times);
  }
}

using namespace LibHUPF;

LIBHUPF_LIBRARY_INTERFACE void* create_ik_solver(double* a, double* d, double* theta, double* alpha, __int8* rots)
{
  bool brots[6] = {0,0,0,0,0,0};
  for (size_t i=0; i<6; ++i)
    if (rots[i]) brots[i] = true;
  return new ik_solver(a, d, theta, alpha, brots);
}

LIBHUPF_LIBRARY_INTERFACE void destroy_ik_solver(void* iks1)
{
  ik_solver* iks = static_cast<ik_solver*>(iks1);
  delete iks;
}

LIBHUPF_LIBRARY_INTERFACE int solve_ik(void* iks1, double* ee, double* ret)
{
  ik_solver* iks = static_cast<ik_solver*>(iks1);
  std::vector<std::vector<double>> sol = iks->solve(ee);
  size_t k=0;
  for (size_t i=0; i<sol.size(); ++i)
  {
    for (size_t j=0; j<sol[i].size();++j)
    {
      ret[k++] = sol[i][j];
    }
  }
  return int(sol.size());
}
