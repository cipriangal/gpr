#include "GPFitNewton.h"
#include "GausProc.h"

using namespace std;
using namespace Eigen;
using namespace FitNewton;

GPFitNewton::GPFitNewton(GausProc *in, int n, int m)
    : FunctionGradHessian(n, m) {
  npar = n;
  nparFixed = m;
  a = in;
}

bool GPFitNewton::calcValGradHessian(const VectorXd &x, double &val,
                                     VectorXd &grad, MatrixXd &hessian) {
  for (int i = 0; i < npar; i++) {
    a->SetPar(i, x(i));
  }
  int prc = a->process();
  if (prc != 0) cout << "GausProc process says:" << prc << endl;

  val = -a->prob();
  // cout<<"GPO:val:"<<val<<endl;

  for (int i = 0; i < npar; i++) {
    double hes = a->probDer(i);
    grad(i) = hes;
    // cout<<"GPO:gr"<<i<<" "<<hes<<endl;
  }

  for (int i = 0; i < npar; i++)
    for (int j = 0; j < npar; j++) {
      double hes = a->probDer(i, j);
      // cout<<"GPO:hes"<<i<<" "<<j<<" "<<hes<<endl;
      hessian(i, j) = hes;
    }
  return true;
}
