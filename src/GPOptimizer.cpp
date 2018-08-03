#include "GPOptimizer.h"
#include <NewtonMinimizerGradHessian.h>
#include <vector>
#include <Eigen/Core>
#include "GausProc.h"
#include "GPFitNewton.h"
#include "FunctionGradHessian.h"

using namespace std;
using namespace Eigen;
using namespace FitNewton;

void GPOptimizer::GPoptimize(int n, int m) {
  FitNewton::NewtonMinimizerGradHessian _minimizer;
  GPFitNewton fct(a, n, m);
  _minimizer.setFunction(&fct);
  double st[10] = {6, 12, 0, 0, 0, 0, 0, 0, 0, 0};
  VectorXd start_point = VectorXd::Zero(n + m);  // input initial guess
  for (int i = 0; i < n; i++) {
    if (par[0] != 0)
      start_point(i) = par[i];
    else
      start_point(i) = st[i];
    cout << "GPOp start par:" << i << " " << start_point(i) << endl;
  }
  // output storage
  VectorXd min_point =
      VectorXd::Zero(n + m);  // output vector for minimize method below

  // minimize
  _minimizer.minimize(start_point, min_point);

  // store output vertex spatial point
  for (int i = 0; i < n; i++) {
    par[i] = min_point(i);
  }
}
