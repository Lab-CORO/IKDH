//PI Defined in Input
//#define PI 3.1415926535897932384626433

//#define DEBUG 0
#define PARALLEL 0

#include <hupf/Calculate.h>
#include <hupf/SolveForAngles.h>

using namespace LibHUPF;
using namespace std;

//not to be exported
vector<vector<double> > Calculate::InverseKin(Input inputParms, vector<double> &times)
{

  double start1 = omp_get_wtime();
  //double end = 0;

  if(inputParms.d[0]!=0)
    inputParms.eePose.m_matrix[3][0] -= inputParms.d[0];
  inputParms.halfAngleSubs();

  Hyperplane h(inputParms);
  //adding hyperplane time
  times.push_back(omp_get_wtime()-start1);

  // Guard: need 8 hyperplanes (4 left + 4 right); spherical wrists / unsupported geometries fall here
  if(h.h.size() < 8)
    return vector<vector<double> >();

  //by parallelizing: time reduced to 30% from 45ms to 30ms
  KinematicSurface surf(h, times); //50%
  //adding kinematic surface time, this take the MOST time!
  vector<vector<double> > result;     
  
  start1 = omp_get_wtime();
  result=SolveForAngles::solveKinSurface(h,surf); //50%
  //adding solve kin surface time
  //free(h.manifoldsUsed);
  times.push_back(omp_get_wtime()-start1);

  return result;
}
