#include "GausProc.h"
#include <algorithm>
#include <string>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <cstdlib>
#include "math.h"
#include "gsl/gsl_sf_bessel.h"

#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TDirectory.h>
#include <TMath.h>

using std::cout;
using std::endl;

namespace {
void printMatrix(const TMatrixD &b) {
  for (int i = 0; i < b.GetNrows(); i++) {
    for (int j = 0; j < b.GetNcols(); j++) {
      std::cout << b(i, j) << " ";
    }
    std::cout << std::endl;
  }
}
}

GausProc::GausProc(const std::vector<double> &X, const std::vector<double> &Y,
                   const std::vector<double> &sigmaY, double xmin, double xmax,
                   int n_predictions, const std::string &oname)
    : detK(0),
      k(vector2(X.size(),X.size())),
      _y(Y),
      _dy(sigmaY),
      y(vector2(X.size(), X.size())),  // FIXME, can be 1D?
      yy(Y),
      dy(sigmaY),
      x(X),
      yOut(vector2(n_predictions, n_predictions)),
      dyOut(vector2(n_predictions, n_predictions)),
      xOut(n_predictions),
      model(RBF),
      nWarpPar(0),
      par(std::vector<double>(10, 0)),
      warpPar(std::vector<double>(10, 0)),
      nrD(X.size()),
      nrE(n_predictions),
      xlow(xmin),
      xhigh(xmax) {
  assert(X.size() == Y.size() and "inconsistent number of points");
  assert(Y.size() == sigmaY.size() and "inconsistent number of points");
  fnm = oname;
  verbosity = 0;
  npar = 2;
  stepsOpt = 0;

  const double dx = (xmax - xmin) / nrE;
  int nrEE=0;
  for (int i = 0; i < nrE; i++) {
    int skip=0;
    double _xx=xmin + i * dx + dx / 2;
    for(int j = 0; j < nrD; j++ ){
      if(fabs( (_xx - x(j))/x(j) ) < 0.001)	
	skip++;
    }
    if(skip!=0)  continue;

    xOut(nrEE) = _xx;
    nrEE++;
  }
  nrE=nrEE;

  for (unsigned i = 0; i < X.size(); ++i) {
    y(i, 0) = Y[i];
  }
}

TMatrixDSym*
GausProc::GetCovarianceMatrix() {
  // For whatever reason using the initialization in construction
  // suffers from rounding errors, so load element by element...
  TMatrixDSym *covMatrix = new TMatrixDSym(nrE);

  for( int i = 0; i < nrE; i++ )
    for( int j = 0; j < nrE; j++ )
      (*covMatrix)[i][j] = dyOut(i,j);
  
  return covMatrix;
}

int GausProc::process() {
  TMatrixDSym mk(nrD);
  // TMatrixD mk(nrD,nrD);

  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) {
      k(i, j) = calcK(x(i), x(j));
      mk(i, j) = k(i, j);
    }

  vector2 ks(nrE, nrD);
  for (int i = 0; i < nrE; i++)
    for (int j = 0; j < nrD; j++) {
      ks(i, j) = calcK(xOut(i), x(j));
    }
  vector2 kss(nrE, nrE);
  for (int i = 0; i < nrE; i++)
    for (int j = 0; j < nrE; j++) {
      kss(i, j) = calcK(xOut(i), xOut(j));
    }

  if (verbosity == 3) printMatrix(mk);

  mk.InvertFast(&detK);
  if (verbosity == 1)
    cout << "Inverted matrix with determinant:" << detK << endl;
  if (verbosity == 3) {
    printMatrix(mk);
  }

  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) k(i, j) = mk(i, j);

  yOut = ks * k * y;

  if (verbosity == 2){
    for (int i = 0; i < nrE; i++)
      for (int j = 0; j < nrE; j++) {
	if(i==j)
	  cout << "i/j/dyOut/kss"
	       << " " << i << " " << j << " " << dyOut(i, j) << " "
	       << calcK(xOut(i), xOut(j)) << endl;
      }
  }

  dyOut = kss - ks * k * ks.transpose();

  if (verbosity == 2){
    for (int i = 0; i < nrE; i++)
      for (int j = 0; j < nrE; j++) {
	if(i==j)
	  cout << "i/j/dyOut/kss"
	       << " " << i << " " << j << " " << dyOut(i, j) << " "
	       << calcK(xOut(i), xOut(j)) << endl;
      }
  }

  return 0;
}

void GausProc::Write(int first, const char *add) const {
  TFile fout(fnm.c_str(), "UPDATE");
  TH1 *ho = new TH1F("ho", "GPR Output", nrE, xlow, xhigh);
  std::vector<double> yyo(nrE, 0);
  std::vector<double> dyyo(nrE, 0);
  for (int i = 0; i < nrE; i++) {
    yyo[i] = yOut(i, 0);
    dyyo[i] = sqrt(dyOut(i, i));//1 sigma, use 1.96* for 95% confidence level
    
    ho->SetBinContent(i + 1, yyo[i]);
    if (dyyo[i] != dyyo[i]) {
      cout << "Problem:" << i << " " << xOut(i) << " " << dyyo[i] << " "
           << yyo[i] << " "
           << " " << dyOut(i, i) << " " << calcK(xOut(i), xOut(i)) << endl;
      ho->SetBinError(i + 1, yyo[i]);
    } else {
      ho->SetBinError(i + 1, dyyo[i]);
    }
  }

  TGraphErrors *hi = new TGraphErrors(nrD, &x(0), &yy(0), 0, &dy(0));
  hi->SetName(Form("hi%s", add));
  hi->SetTitle("Gaussian Process Input");
  hi->SetLineColor(2);
  hi->SetMarkerColor(2);
  hi->SetMarkerStyle(20);
  if (first == 0) {
    hi->Write();
    fout.mkdir("steps");
  } else if (first == -1)
    hi->Write();

  ho->SetName(Form("ho%s", add));
  ho->SetTitle(Form("GPR Output l/sf: %4.2f %4.2f", par[1], par[0]));
  ho->SetLineColor(4);
  ho->SetMarkerColor(1);
  ho->SetMarkerStyle(31);
  ho->SetMarkerSize(0.5);
  if (first >= 0) {
    ho->SetName(Form("ho_%d", first));
    fout.cd("steps");
  }
  ho->Write();

  fout.Close();
}

double GausProc::prob() const {
  const double pi = TMath::Pi();

  const double dt =
      vector2(y.transpose() * k * y)(0, 0);

  if (verbosity) {
    for (int i = 0; i < npar; i++) cout << par[i] << " ";
    cout << endl;
    cout << "prob:dt/det " << dt << " " << detK << " " << log(detK) << endl;
  }
  double logp = -0.5 * (dt + log(detK) +
                        nrD * log(2 * pi));  // this is the log likelihood
  return logp;
}

double GausProc::probDer(int n) const {
  vector2 kd(nrD, nrD);
  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) {
      kd(i, j) = calcKd1(x(i), x(j), n);
    }

  const vector2 dm1 = k * kd;
  double term2 = -dm1.trace();
  const vector2 dn = dm1 * k;

  const double dt =
      vector2(y.transpose() * dn * y)(0, 0);

  const double logp = -0.5 * (dt + term2);  // this is the - log likelihood
  return logp;
}

double GausProc::probDer(int n, int m) const {
  vector2 kd(nrD, nrD);
  vector2 ks(nrD, nrD);
  vector2 ksT(nrD, nrD);  
  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) {
      kd(i, j) = calcKd2(x(i), x(j), n, m);
      ks(i, j) = calcKd1(x(i), x(j), n);
      ksT(i, j) = calcKd1(x(i), x(j), m);
    }

  double t1 = -(k * ks * k * ksT).trace();
  t1 += vector2(y.transpose()  * k * ks * k * ksT * k * y)(0, 0);
  t1 += vector2(y.transpose()  * k * ksT * k * ks * k * y)(0, 0);

  const vector2 dm4 = k * kd;
  t1 += dm4.trace();
  t1 -= vector2(y.transpose() * dm4 * k * y)(0, 0);

  double logp = 0.5 * (t1);  // this is the -log likelihood
  return logp;
}

double GausProc::calcK(double x1, double x2) const {
  const double uncert = getUncert(x1, x2);
  switch (model) {
    case RBF:
      return std::pow(par[0], 2) * exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) +
             std::pow(uncert, 2);
    case MATERN: {
      double v = sqrt(2.0 * par[0]) * std::pow(x1 - x2, 2) / par[1];
      if (v == 0 || v > 600) return std::pow(uncert, 2);
      if (verbosity == 2) {
        cout << "1/2/v/p0/bessel:" << x1 << " " << x2 << " " << v << " "
             << par[0] << " ";
        cout << gsl_sf_bessel_Knu(par[0], v) << endl;
      }
      return std::pow(2, 1 - par[0]) / (TMath::Gamma(par[0])) * std::pow(v, par[0]) *
                 (gsl_sf_bessel_Knu(par[0], v)) +
             std::pow(uncert, 2);
    }
    default:
      return std::pow((x1 + x2) / 2, par[2]) * std::pow(par[0], 2) *
                 exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) +
             std::pow(uncert, 2);
  }
}

// calculate the first derivative
double GausProc::calcKd1(double x1, double x2, int n) const {
  switch (model) {
    case RBF:
      if (n == 0)
        return 2 * par[0] * exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2)));
      else
        return std::pow(par[0], 2) * std::pow(x1 - x2, 2) *
               exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) / std::pow(par[1], 3);
    default:
      return 1;
  }
}
// calculate the second derivatives
double GausProc::calcKd2(double x1, double x2, int n, int m) const {
  switch (model) {
    case RBF:
      if (n == 0 && m == 0)  // sigmaf second derivative
        return 2 * exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2)));
      else if (n == 1 && m == 1)  // length second derivative
        return std::pow(par[0], 2) *
               (-3. * std::pow(x1 - x2, 2) *
                    exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) /
                    std::pow(par[1], 4) +
                std::pow(x1 - x2, 4) * exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) /
                    std::pow(par[1], 6));
      else  // cross terms are the same
        return par[0] * 2 * std::pow(x1 - x2, 2) *
               exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) / std::pow(par[1], 3);
    default:
      return 1;
  }
}

double GausProc::getUncert(double x1, double x2) const {
  double precision = 0.0001;
  if (x1 != 0) {
    if (fabs((x1 - x2) / x1) > precision) return 0;
  } else {
    if (fabs(x1 - x2) > precision) return 0;
  }

  for (int i = 0; i < nrD; i++)
    if (x1 != 0) {
      if (fabs((x(i) - x1) / x1) <= precision) return dy(i);
    } else {
      if (fabs(x(i) - x1) <= precision) return dy(i);
    }

  return 0;
}

void GausProc::SetKernel(GausProc::Kernel mod) {
  model = mod;
  switch (model) {
    case RBF:
      std::cout << "You have chosen an RBF kernel and the number of parameters "
                   "was set to two" << std::endl;
      npar = 2;
      break;
    case MATERN:
      std::cout
          << "You have chosen a Matern kernel and the number of parameters "
             "was set to two" << std::endl;
      npar = 2;
      break;
    default:
      std::cout << "You have chosen a silly RBF kernel and the number of "
                   "parameters was set to three" << std::endl;
      npar = 3;
      break;
  }
}

void GausProc::warp(int wp) {
  if (wp == 0) {
    for (int i = 0; i < nrD; i++) {
      if (_y(i) == 0) {
        yy(i) = 0;
        y(i, 0) = yy(i);
        dy(i) = 1;
        continue;
      }
      double dmy = log(_y(i));
      double dmdy = 0;
      double dmy1 = log(fabs(_y(i) + _dy(i)));
      double dmy2 = log(fabs(_y(i) - _dy(i)));
      if (fabs(dmy1 - dmy) > fabs(dmy2 - dmy))
        dmdy = fabs(dmy1 - dmy);
      else
        dmdy = fabs(dmy2 - dmy);

      yy(i) = dmy;
      y(i, 0) = yy(i);
      if (dmy == 0)
        dy(i) = 1;
      else
        dy(i) = dmdy;
      if (verbosity == 1) {
        cout <<"warp:"<< i << " " << yy(i) << " " << dy(i) << endl;
        cout <<"warp:"<< i << "~" << _y(i) << "~" << _dy(i) << endl;
      }
    }
  } else if(wp==1){
    for (int i = 0; i < nrD; i++) {
      double dmy = TMath::TanH(_y(i));
      double dmdy = 0;
      double dmy1 = TMath::TanH(_y(i) + _dy(i));
      double dmy2 = TMath::TanH(_y(i) - _dy(i));
      if (fabs(dmy1 - dmy) > fabs(dmy2 - dmy))
        dmdy = fabs(dmy1 - dmy);
      else
        dmdy = fabs(dmy2 - dmy);

      yy(i) = dmy;
      y(i, 0) = yy(i);
      if (dmy == 0)
        dy(i) = 1;
      else
        dy(i) = dmdy;
      if (verbosity == 1) {
        cout <<"warp:" << i << " " << yy(i) << " " << dy(i) << endl;
        cout <<"warp:" << i << "~" << _y(i) << "~" << _dy(i) << endl;
      }
    }    
  } else {
    for (int i = 0; i < nrD; i++) {
      double dmy1 = _y(i) + _dy(i);
      double dmy2 = _y(i) - _dy(i);
      int niter = nWarpPar / 3;
      for (int j = 0; j < niter; j++) {
        dmy1 += warpPar[j] *
                TMath::TanH(warpPar[j + 1] * (_y(i) + _dy(i) + warpPar[j + 2]));
        dmy2 += warpPar[j] *
                TMath::TanH(warpPar[j + 1] * (_y(i) - _dy(i) + warpPar[j + 2]));
      }
      yy(i) = (dmy1 + dmy2) / 2;
      y(i, 0) = yy(i);
      dy(i) = (dmy1 - dmy2) / 2;
    }
  }
}

void GausProc::unwarp(int wp) {
  if (wp == 0) { 
    for (int i = 0; i < nrD; i++) {
      yy(i) = _y(i);
      y(i, 0) = yy(i);
      dy(i) = _dy(i);
      
      if (verbosity == 2) {
        cout <<"unwarp:"<< i << " " << _y(i) << " " << _dy(i) << endl;
        cout <<"unwarp:"<< i << "~" << yy(i) << "~" << dy(i) << endl;
      }
    }
    vector2 covUnwarp(vector2(nrE,nrE));
    vector1 expY(nrE);
    
    for (int i = 0; i < nrE; i++) {
      for(int j=0;j<nrE;j++){
    	double cov_unwarp=exp(yOut(i,0)) * dyOut(i,j) * exp(yOut(j,0));
	covUnwarp(i,j)=cov_unwarp;
      }     
      expY(i) = exp(yOut(i,0));
    }
    for(int i=0;i<nrE;i++)
      {
	yOut(i,0)=expY(i);
	for(int j=0;j<nrE;j++) dyOut(i,j)=covUnwarp(i,j);
      }
  }else if(wp==1){
    for (int i = 0; i < nrD; i++) {
      yy(i) = _y(i);
      y(i, 0) = yy(i);
      dy(i) = _dy(i);
      
      if (verbosity == 2) {
        cout <<"unwarp:"<< i << " " << _y(i) << " " << _dy(i) << endl;
        cout <<"unwarp:"<< i << "~" << yy(i) << "~" << dy(i) << endl;
      }
    }
    for (int i = 0; i < nrE; i++) {

      double dmy = TMath::ATanH(yOut(i, 0));
      double dmdy = 0;
      double uncert = sqrt(dyOut(i, i));
      if (uncert != uncert) continue;
      double dmy1 = TMath::ATanH(yOut(i, 0) + uncert);
      double dmy2 = TMath::ATanH(yOut(i, 0) - uncert);
      if (fabs(dmy1 - dmy) > fabs(dmy2 - dmy))
        dmdy = fabs(dmy1 - dmy);
      else
        dmdy = fabs(dmy2 - dmy);

      if (verbosity == 2) {
        cout << i << " " << yOut(i, i) << " " << sqrt(dyOut(i, i)) << endl;
        cout << i << "~" << dmy << "~" << dmdy << endl;
      }
      yOut(i, 0) = dmy;
      dyOut(i, i) = std::pow(dmdy, 2);
    }
    
  } else {
    cout << "It would be great if you implement something here (say a Newton "
            "Raphson algorithm)! :)" << endl;
  }
  Write(-1, "_unwarp");
}

void GausProc::Integral(const double _xmin, const double _xmax, 
		       double &Integral, double &dIntegral)
{
  Integral=0;
  dIntegral=0;
  double width=0;

  for(int i=0;i<nrE;i++)
    if(xOut(i)>=_xmin && xOut(i)<_xmax)
      {
	width=fabs(xOut(i)-xOut(i+1));
	Integral+=width*yOut(i,0);
	double dm=0;
	for(int j=0;j<nrE;j++)
	  if(xOut(j)>=_xmin && xOut(j)<_xmax)
	    {
	      dm+=dyOut(i,j)*pow(width,2);
	    }
	dIntegral+=dm;
      }
  dIntegral=sqrt(dIntegral);
}

void GausProc::Read(const std::string &fname, std::vector<double> &X,
                    std::vector<double> &Y, std::vector<double> &dY,
                    double &xmin, double &xmax, unsigned &nPredictions) {
  double nrD;
  std::ifstream fin(fname.c_str());
  fin >> nrD >> nPredictions >> xmin >> xmax;
  X.resize(nrD);
  Y.resize(nrD);
  dY.resize(nrD);

  for (unsigned i = 0; i < nrD; i++) {
    fin >> X[i] >> Y[i] >> dY[i];
  }
  cout << "Number of data points:" << nrD << endl;
}
