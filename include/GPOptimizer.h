#ifndef GPOPTIMIZER_H
#define GPOPTIMIZER_H

#include <algorithm>

class GausProc;
class GPFitNewton;

class GPOptimizer {
 public:
  GPOptimizer(GausProc *in) : a(in) {
    std::fill_n(par, sizeof(par) / sizeof(par[0]), 0);
  }
  GPOptimizer(GausProc *in, double p1, double p2) : a(in) {
    std::fill_n(par, sizeof(par) / sizeof(par[0]), 0);
    par[0] = p1;
    par[1] = p2;
  }
  void GPoptimize(int n, int m);
  double getPar(int i) { return par[i]; }

 private:
  GausProc *a;
  double par[10];
};

#endif
