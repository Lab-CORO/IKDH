#define PARALLEL 0

#include <hupf/Calculate.h>
#include <hupf/SolveForAngles.h>

using namespace LibHUPF;
using namespace std;

vector<vector<double> > Calculate::InverseKin(Input inputParms, vector<double> &times)
{
  double start1 = omp_get_wtime();

  if(inputParms.d[0]!=0)
    inputParms.eePose.m_matrix[3][0] -= inputParms.d[0];
  inputParms.halfAngleSubs();

  Hyperplane h(inputParms);
  times.push_back(omp_get_wtime()-start1);

  // Guard: need 8 hyperplanes (4 left + 4 right); spherical wrists / unsupported geometries fall here
  if(h.h.size() < 8)
    return vector<vector<double> >();

  KinematicSurface surf(h, times);
  vector<vector<double> > result;

  start1 = omp_get_wtime();
  result=SolveForAngles::solveKinSurface(h,surf);
  times.push_back(omp_get_wtime()-start1);

  return result;
}
