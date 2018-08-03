#ifndef GPFITNEWTON_H
#define GPFITNEWTON_H

#include "FunctionGradHessian.h"

class GausProc;
class GPFitNewton : public FitNewton::FunctionGradHessian {
 public:
  GPFitNewton(GausProc *in, int n, int m);
  FitNewton::FunctionGradHessian *Clone() const { return 0; }

 private:
  bool calcValGradHessian(const Eigen::VectorXd &x, double &val,
                          Eigen::VectorXd &grad, Eigen::MatrixXd &hessian);
  GausProc *a;
  int npar;
  int nparFixed;
};

#endif
