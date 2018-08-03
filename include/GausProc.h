#ifndef GAUSPROC_H
#define GAUSPROC_H

#include <iostream>
#include <string>
#include <vector>
#include <TMatrixDSym.h>

#ifndef __CINT__
#include <cassert>
#include <Eigen/Dense>
#endif

class GausProc {
 public:
  enum Kernel {
    RBF = 0,
    MATERN = 1
  };

  GausProc(const std::vector<double> &x, const std::vector<double> &y,
           const std::vector<double> &sigma_y, double xmin, double xmax,
           int n_predictions, const std::string &oname);

  static void Read(const std::string &fname, std::vector<double> &x,
                   std::vector<double> &y, std::vector<double> &dy,
                   double &xmin, double &xmax, unsigned &nPredictions);

  void SetKernel(Kernel mod);
  void SetRangeXlow(double x) { xlow = x; }
  void SetRangeXhigh(double x) { xhigh = x; }
  void SetPar(int n, double p) { par[n] = p; }
  void SetWarpPar(int n, double p) { warpPar[n] = p; }
  void SetVerbosity(int ver) { verbosity = ver; }

  TMatrixDSym* GetCovarianceMatrix();

  int process();
  void Write(int first, const char *add = "") const;
  double prob() const;
  double probDer(int n) const;
  double probDer(int n, int m) const;
  void warp(int wp = 0);
  void unwarp(int wp = 0);
  void Integral(const double _xmin, const double _xmax, 
		double &Integral, double &dIntegral);

#ifndef __CINT__
#define CHECKED_ACCESS 0
  class vector1 {
    // a simple dynamic one-dimensional array on top of std::vector
   public:
    vector1(int n) : data(n, 0) {}
    vector1(const std::vector<double> &d) : data(d) {}
    std::vector<double>::iterator begin() { return data.begin(); }
    std::vector<double>::iterator end() { return data.end(); }
    std::size_t rows() const { return data.size(); }

    double &operator()(std::size_t i) {
#ifndef EIGEN_NO_DEBUG
      return data.at(i);
#else
      return data[i];
#endif
    }
    const double &operator()(std::size_t i) const {
#ifndef EIGEN_NO_DEBUG
      return data.at(i);
#else
      return data[i];
#endif
    }

   private:
    std::vector<double> data;
  };
  typedef Eigen::MatrixXd vector2;
#endif

 private:
  double getUncert(double x1, double x2) const;
  double calcK(double x1, double x2) const;
  double calcKd1(double x1, double x2, int n) const;
  double calcKd2(double x1, double x2, int n, int m) const;

  double detK;
#ifndef __CINT__
  vector2 k;

  vector1 _y;
  vector1 _dy;
  vector2 y;
  vector1 yy;
  vector1 dy;
  vector1 x;

  vector2 yOut;
  vector2 dyOut;
  vector1 xOut;
#endif

  Kernel model;
  int npar;
  int nWarpPar;
  std::vector<double> par;
  std::vector<double> warpPar;
  int stepsOpt;
  int nrD;          // nr data points
  int nrE;          // nr of extrapolations
  double xlow;      // lower limit of extrapolations
  double xhigh;     // upper limit
  std::string fnm;  // output filename
  int verbosity;
};

#endif /* __GAUSPROC_H__ */
