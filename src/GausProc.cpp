#include "GausProc.h"
#include <algorithm>
#include <string>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cstdlib>
#include "math.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_errno.h"

#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TDecompChol.h>
#include <TRandom3.h>
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
  void printMatrix(const GausProc::vector2 &b) {
    for (int i = 0; i < b.rows(); i++) {
      for (int j = 0; j < b.cols(); j++) {
	std::cout<<std::setprecision(20) << b(i, j) << " ";
      }
      std::cout << std::endl;
    }
  }
}

GausProc::GausProc(const std::vector<double> &X, const std::vector<double> &Y,
                   const std::vector<double> &sigmaY, double xmin, double xmax,
                   int n_predictions, const std::string &oname)
  : detK(0),//determinant of the input covariance matrix
    k(vector2(X.size(),X.size())), //input covariance matrix
    _y(Y), //original y values of the input data 
    _dy(sigmaY), //original uncertainty of the y value of the input data 
    yy(Y), //y values that will be used in the GPR (they may get warped)
    dy(sigmaY), //uncertainties of the y values above (they may also get warped)
    x(X), //input x values
      //y(vector2(X.size(), X.size())),  
    y(vector2(X.size(), 1)),  //collumn vector that will be used in the matrix multiplication
    yOut(vector2(n_predictions, n_predictions)), //output y values
    dyOut(vector2(n_predictions, n_predictions)), //uncertainties of the output y values
    cholDecomp(vector2(n_predictions, n_predictions)), //the cholenski decomposed matrix
    xOut(n_predictions), //x positions for the GPR predictions
    model(RBF), //this is the variable that knows the kernel you have (for now only RBF)
    nWarpPar(0), //if you wish to warp by tanh functions you should say here how many you want (multiples of 3)
    par(std::vector<double>(10, 0)), //vector for the kernel hyperparameters
    warpPar(std::vector<double>(10, 0)),//the warping parameters
    nrD(X.size()), //number of data points
    nrE(n_predictions), //number of predictions
    xlow(xmin), //min value of x for the prediction
    xhigh(xmax) { //max values of x for the prediction
  assert(X.size() == Y.size() and "inconsistent number of points");
  assert(Y.size() == sigmaY.size() and "inconsistent number of points");
  fnm = oname;
  verbosity = 0;
  npar = 2;
  stepsOpt = 0;

  //cout<<"xOut "<<endl;
  //sometimes a prediction point falls on a input data point --> this skips it
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
    //cout<<xOut(nrEE)<<" ";
    nrEE++;
  }
  //cout<<endl;
  nrE=nrEE;

  cout<<"xIn "<<endl;
  for (unsigned i = 0; i < X.size(); ++i) {
    y(i, 0) = Y[i];
    cout<<x(i)<<" ";
  }
  cout<<endl;
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

///run the GPR for one set of hyperparameters (that you can set from SetPar
int GausProc::process() {
  TMatrixDSym mk(nrD);
  // TMatrixD mk(nrD,nrD);

  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) {
      k(i, j) = calcK(x(i), x(j));//calculate the covariance matrix for the input
      mk(i, j) = k(i, j);
    }

  vector2 ks(nrE, nrD);
  for (int i = 0; i < nrE; i++)
    for (int j = 0; j < nrD; j++) {
      ks(i, j) = calcK(xOut(i), x(j));//calculate the covariance between input and output points
    }
  vector2 kss(nrE, nrE);
  for (int i = 0; i < nrE; i++)
    for (int j = 0; j < nrE; j++) {
      kss(i, j) = calcK(xOut(i), xOut(j));//calculate the covariance between output points
    }

  if (verbosity >= 3) {
    cout<<"covariance matrix:"<<endl;
    printMatrix(mk);
    cout<<"ks matrix:"<<endl;
    printMatrix(ks);
    cout<<"kss matrix:"<<endl;
    printMatrix(kss);
    cout<<"y matrix:"<<endl;
    printMatrix(y);
  }

  mk.InvertFast(&detK); //Root invert function
  if (verbosity >= 1)
    cout << "Inverted matrix with determinant:" << detK << endl;
  if (verbosity >= 3) {
    cout<<"Inverted covariance matrix:"<<endl;
    printMatrix(mk);
  }

  //copy over the inverted matrix into the original covariance matrix
  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) k(i, j) = mk(i, j); 

  yOut = ks * k * y; //calculate the y out values

  if (verbosity >= 4){
    for (int i = 0; i < nrE; i++)
      for (int j = 0; j < nrE; j++) {
	if(i==j)
	  cout << "i/j/dyOut/kss"
	       << " " << i << " " << j << " " << dyOut(i, j) << " "
	       << kss(i,j) << endl;
      }
  }

  dyOut = kss - ks * k * ks.transpose(); //calculate the uncertainties of the y values

  if (verbosity >= 4){
    for (int i = 0; i < nrE; i++)
      for (int j = 0; j < nrE; j++) {
	if(i==j)
	  cout << "i/j/dyOut/kss"
	       << " " << i << " " << j << " " << dyOut(i, j) << " "
	       << kss(i,j) << endl;
      }
  }

  if (verbosity >= 3) {
    cout<<"Output y:"<<endl;
    printMatrix(yOut);
    cout<<"Output dy:"<<endl;
    printMatrix(dyOut);
  }


  return 0;
}

///write function that creates root file with Tgraph for the input and TH1 for the output
//it still has an old feature that can be used to save histos for ouput for different sets of hyperparameters
//default is add=""; but can save other histos using the name ho(output)_"*add" or hi(input)_"*add"  
void GausProc::Write(int first, const char *add, float xmin, float xmax){
  TFile fout(fnm.c_str(), "UPDATE");
  TH1 *ho = new TH1F("ho", "GPR Output", nrE, xlow, xhigh);
  std::vector<double> yyo(nrE, 0);
  std::vector<double> dyyo(nrE, 0);
  double integ=-5; //integral value
  double dinteg=0; //uncertainty of integral value
  if(xmin!=xmax)
    Integral(xmin,xmax,integ,dinteg);

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
  if(xmin!=xmax)
    ho->SetTitle(Form("GPR Output l/sf (int/dint): %4.2f %4.2f (%4.2f %4.2f)", 
		      par[1], par[0],integ,dinteg));
  else
    ho->SetTitle(Form("GPR Output l/sf: %4.2f %4.2f",par[1],par[0])); 
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

///calculate the log likelihood for a certain set of hyperparameters 
//this is the equation that needs to be minimized (feeds into GPOptimizer)
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

///first derivative for the log likelihood eq
//feeds into GPOptimizer
double GausProc::probDer(int n) const {
  vector2 kd(nrD, nrD);
  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) {
      kd(i, j) = calcKd1(x(i), x(j), n); //calculate the first derivative for the RBF kernel
    }

  const vector2 dm1 = k * kd;
  double term2 = -dm1.trace();
  const vector2 dn = dm1 * k;

  const double dt =
      vector2(y.transpose() * dn * y)(0, 0);

  const double logp = -0.5 * (dt + term2);  // this is the - log likelihood
  return logp;
}

///second derivatives for the log likelihood eq
//feeds into GPOptimizer
double GausProc::probDer(int n, int m) const {
  vector2 kd(nrD, nrD);
  vector2 ks(nrD, nrD);
  vector2 ksT(nrD, nrD);  
  for (int i = 0; i < nrD; i++)
    for (int j = 0; j < nrD; j++) {//calculate the second derivatives for the RBF kernel
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

///calculates the covariance between two points
//here we can implement more kernels
double GausProc::calcK(double x1, double x2) const {
  const double uncert = getUncert(x1, x2);
  switch (model) {
    case RBF:
      return std::pow(par[0], 2) * exp(-std::pow(x1 - x2, 2) / (2 * std::pow(par[1], 2))) +
             std::pow(uncert, 2);
    case MATERN: {
      double v = sqrt(2.0 * par[0]) * std::pow(x1 - x2, 2) / par[1];
      if (v == 0 || v > 600) return std::pow(uncert, 2);
      if (verbosity >= 2) {
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

/// calculate the first derivative of the RBF kernel
//for different kernels the derivatives should be implemented to get the optimezer to work
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
// calculate the second derivatives of the RBF kernel
//for different kernels the derivatives should be implemented to get the optimezer to work
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

///get uncertainty when calculating the covariance between x1 and x2
// it returns something only if x1==x2 and x1 is from the input data
double GausProc::getUncert(double x1, double x2) const {
  double precision = 0.0001;
  if(verbosity>=4) cout<<"getUncert: "<<x1<<" "<<x2;    
  if (x1 != 0) {
    if (fabs((x1 - x2) / x1) > precision){
      if(verbosity>=4) cout<<" return (a) 0"<<endl;
      return 0;
    }
  } else {
    if (fabs(x1 - x2) > precision){
      if(verbosity>=4) cout<<" return (b) 0"<<endl;
      return 0;
    }
  }

  for (int i = 0; i < nrD; i++)
    if (x1 != 0) {
      if (fabs((x(i) - x1) / x1) <= precision){
	if(verbosity>=4) cout<<" return (c) "<<dy(i)<<endl;
	return dy(i);
      }
    } else {
      if (fabs(x(i) - x1) <= precision){
	if(verbosity>=4) cout<<" return (d) "<<dy(i)<<endl;
	return dy(i);
      }
    }

  if(verbosity>=4) cout<<" return (e) 0 "<<endl;
  return 0;
}

///set kernel method
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

///warping function
//for now only the log is properly implemented (so it should be used only for data that is positive only)
void GausProc::warp(int wp) {
  if (wp == 0) {
    for (int i = 0; i < nrD; i++) {
      double dmy;
      if (_y(i) == 0) 
	dmy = -10;
      else
	dmy = log(_y(i));

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
      if (verbosity >= 1) {
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
      if (verbosity >= 1) {
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

///un-warping function (writes the unwarped histos at the end)
void GausProc::unwarp(int wp) {
  if (wp == 0) { 
    for (int i = 0; i < nrD; i++) {
      yy(i) = _y(i);
      y(i, 0) = yy(i);
      dy(i) = _dy(i);
      
      if (verbosity >= 1) {
        cout <<"unwarp:"<< i << " " << _y(i) << " " << _dy(i) << endl;
        cout <<"unwarp:"<< i << "~" << yy(i) << "~" << dy(i) << endl;
      }
    }
    vector2 covUnwarp(vector2(nrE,nrE));
    vector1 expY(nrE);
    
    for (int i = 0; i < nrE; i++) {
      for(int j=0;j<nrE;j++){
	double cov_unwarp= exp(yOut(i,0)) * dyOut(i,j) * exp(yOut(j,0));
	covUnwarp(i,j)=cov_unwarp;
      }     
      expY(i) = exp(yOut(i,0));
      if(verbosity>=1)
	cout<<"i:og:unwp"<<i<<" "<<yOut(i,0)<<" "<<expY(i)<<endl;
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
      
      if (verbosity >= 1) {
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

      if (verbosity >= 2) {
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
  Write(-1, "_unwarp",30,50);
  //Write(-1, "_unwarp"); <---- normal version shoul fix this
}

///cholenski decomposition function
void GausProc::DecomposeChol()
{
  cout<<"matrix dimension "<<nrE<<endl;
  if(nrE>5000){
    cout<<"DecomposeChol: you want to do a Cholesky decomposition for "<<nrE
	<<" points! You better get some popcorn (this will take a while)!"<<endl;
    return;
  }

  //TMatrixDSym cholMat(nrE);
  gsl_matrix * m = gsl_matrix_alloc (nrE, nrE);
  for(int i=0;i<nrE;i++){
    //for(int j=0;j<=i;j++){
    //for(int j=i;j<nrE;j++){
    for(int j=0;j<nrE;j++){
      double val=dyOut(i,j);
      if(i==j)
	val+=pow(10,-9);
      //14 > 50 <100
      //15 < 17

      // if(val>=0)
      // 	val=(floor(val*pow(10,12))/pow(10,12));
      // else
      // 	val=-(floor(-val*pow(10,12))/pow(10,12));
      
      //cout<<std::setprecision(20)<<dyOut(i,j)<<" "<<val<<endl;
      gsl_matrix_set (m, i, j, val ); 
      //gsl_matrix_set (m, j, i, val ); 
      //cholMat(i,j)=val;
      //if(i<=j && dyOut(i,j)!=dyOut(j,i))
      //cout<<i<<" "<<j<<" input "<<std::setprecision(20)<<gsl_matrix_get(m,i,j)<<endl;
	//<<" "<<dyOut(i,j)<<" "<<dyOut(j,i)<<" "<<dyOut(j,i)-dyOut(i,j)<<endl;
    }
  }

  // TDecompChol decompRoot(cholMat);
  // decompRoot.Decompose();
  // cout<<"Decomposed root!"<<endl;
  //TMatrixD Uchol=decompRoot.GetU();
  // cout<<"Uchol"<<endl;
  // printMatrix(Uchol);

  // vector2 lt2(nrE,nrE);
  // for(int i=0;i<nrE;i++)
  //   for(int j=0;j<nrE;j++){
  // 	cholDecomp(i,j)=Uchol(j,i);
  // 	lt2(i,j)=Uchol(i,j);
  //     }
  // printMatrix(cholDecomp);
  // vector2 recomp2(nrE,nrE);
  // recomp2=cholDecomp * lt2;
  // cout<<"dyOut re-composed"<<endl;
  // printMatrix(recomp2);
  // vector2 diff2(nrE,nrE);
  // diff2=recomp2-dyOut;
  // cout<<"recomposed - dyOut"<<endl;
  // printMatrix(diff2);
  // double max2=-500;
  // double min2=500;
  // for(int i=0;i<nrE;i++)
  //   for(int j=0;j<nrE;j++)
  //     if(diff2(i,j)>max2) max2=diff2(i,j);
  //     else if(diff2(i,j)<min2) min2=diff2(i,j);
  // cout<<"max:min "<<max2<<" "<<min2<<endl;

  // gsl_error_handler_t* old_handler;
  // old_handler = gsl_set_error_handler_off();
  // if(verbosity>=1) cout<<"You are decomposing without the gsl handler!"<<endl;

  int testI = nrE/2;
  cout<<"\t"<<testI<<" input "<<std::setprecision(20)<<gsl_matrix_get(m,testI,testI)
      <<" <> "<<dyOut(testI,testI)<<endl;

  // try{
    gsl_linalg_cholesky_decomp (m);
  // }catch(const std::exception& exception){
  //   cout<<"Problem with Cholensky decomposition: "<<exception.what()<<endl;
  //   int testI = nrE/2;
  //   cout<<"\t"<<testI<<" input "<<std::setprecision(20)<<gsl_matrix_get(m,testI,testI)
  // 	<<" <> "<<dyOut(testI,testI)<<endl;
  //   exit(2);
  // }

  cout<<"Decomposed gsl!"<<endl;
  
  vector2 lt(nrE,nrE);
  for(int i=0;i<nrE;i++)
    for(int j=0;j<nrE;j++)
      if(j<i){
	cholDecomp(i,j)=gsl_matrix_get(m,i,j);
	lt(i,j)=0;
      }
      else if(j>i){
	cholDecomp(i,j)=0;
	lt(i,j)=gsl_matrix_get(m,i,j);
      }
      else{
	cholDecomp(i,j)=gsl_matrix_get(m,i,j);
	lt(i,j)=gsl_matrix_get(m,i,j);
      }

  vector2 recomp(nrE,nrE);
  recomp=cholDecomp * lt;
  vector2 diff(nrE,nrE);
  diff=recomp-dyOut;
  if(verbosity>4){
    cout<<"gsl decomp"<<endl;
    printMatrix(cholDecomp);
    cout<<"dyOut re-composed"<<endl;
    printMatrix(recomp);
    cout<<"recomposed - dyOut"<<endl;
    printMatrix(diff);
    double max=-500;
    double min=500;
    for(int i=0;i<nrE;i++)
      for(int j=0;j<nrE;j++)
	if(diff(i,j)>max) max=diff(i,j);
	else if(diff(i,j)<min) min=diff(i,j);
    cout<<"max:min "<<max<<" "<<min<<endl;
  }
  gsl_matrix_free (m);
}

///if you have the cholenski decomposed matrix you can sample a function with this
void GausProc::GetOneFunction(int seed,std::vector<double> &Xf, std::vector<double> &Yf, int &props)
{
  Xf.resize(nrE);
  Yf.resize(nrE);
  vector2 mean(vector2(nrE,1));
  
  for(int i=0;i<nrE;i++){
    Xf[i]=xOut(i);
    mean(i,0)=yOut(i,0);
    //cout<<i<<" x mean "<<Xf[i]<<" "<<mean(i,0)<<endl;
  }

  vector2 covFct(vector2(nrE,nrE));
  for(int i=0;i<nrE;i++)
    for(int j=0;j<nrE;j++)
      {
  	covFct(i,j)=cholDecomp(i,j);
	// if(cholDecomp(i,j)!=cholDecomp(j,i) && j>=i)
	//   cout<<"cholDecomp "<<i<<" "<<j<<" "<<cholDecomp(i,j)<<" "<<cholDecomp(j,i)<<endl;
	//cout<<i<<" "<<j<<" cov "<<covFct(i,j)<<" "<<cholDecomp(i,j)<<" "<<endl;
      }
  
  vector2 randX(vector2(nrE,1));
  vector2 _yf(vector2(nrE,1));
  TRandom3 sd(seed);
  for(int i=0;i<nrE;i++){
    randX(i,0)=sd.Gaus(0,1);
    //cout<<i<<" "<<randX(i,0)<<endl;
  }
  
  _yf=covFct*randX + mean;      
  for(int i=0;i<nrE;i++){
    Yf[i]=_yf(i,0);
    //cout<<i<<" y "<<Yf[i]<<endl;
  }
  
  props=FunctionProp(Xf,Yf);
}

///implementation of function properties (like monotonicity, convexity .. etc)
int GausProc::FunctionProp(std::vector<double> xf, std::vector<double> yf)
{
  int props=0;
  int nr=xf.size();

  std::vector<double> derSn;
  derSn.resize(nr);
  for(int i=0;i<nr-1;i++)
    derSn[i]=(yf[i+1]-yf[i])/(xf[i+1]-xf[i]);
  
  double sgn=derSn[0];
  int consistent=0;
  int soft=0;
  for(int i=0;i<nr-1;i++)
    {
      if(sgn*derSn[i]<0) consistent++;
      if(derSn[i]==0) soft++;
    }

  if(consistent==0){
    if(sgn>0){
      if(soft) props+=1;
      else props+=2;
    }
    if(sgn<0){
      if(soft) props+=4;
      else props+=8;
    }
  }

  //add convex/convave 

  return props;
}

///proper integral for a range in the GPR prediction that takes into account the correlations between those points
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

///one way to read points in
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
