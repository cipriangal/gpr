void example2()
{
  gSystem->Load("../install/lib/libGausProc.so");

  int verbosity=0;

  string outfile="o_gp_ex2.root";
  const unsigned nData=30;
  double _x[nData];
  double _y[nData];
  double _dy[nData];
  std::vector<double> x,y,sigma_y;

  TF1 *fS=new TF1("fS","gaus(0)",-6,6);
  fS->SetParameter(0,200);
  fS->SetParameter(1,1);
  fS->SetParameter(2,1);
  TF1 *fB = new TF1("fB","pol3(0)",-6,6);
  fB->SetParameters(1000,-0.2,0.2,-1);

  for(int i=0;i<nData;i++){
    double xTemp = -6 + i* 12./nData;
    double bkg = fB->Eval(xTemp);
    double sig = fS->Eval(xTemp);

    if(xTemp>-1 && xTemp<3) continue;

    x.push_back(xTemp);
    y.push_back(bkg + sig);
    sigma_y.push_back(sqrt(bkg+sig));
  }
  double xmin=-7;
  double xmax=7;
  unsigned nPredictions=1000;
  
  GausProc a(x, y, sigma_y, xmin, xmax, nPredictions, outfile.c_str());
  gSystem->Exec(Form("rm -f %s",outfile.c_str()));
  
  cout<<"Input data:"<<endl;
  for(int i=0;i<x.size();i++)
    cout<<x[i]<<" "<<y[i]<<" "<<sigma_y[i]<<endl;

  a.SetVerbosity(verbosity);
  a.SetKernel(GausProc::RBF);
  
  GPOptimizer c(&a,500,2);
  c.GPoptimize(2,0);
  cout<<c.getPar(0)<<" "<<c.getPar(1)<<endl;  
  
  GausProc d(x, y, sigma_y, xmin, xmax, nPredictions, outfile.c_str());
  d.SetPar(0,c.getPar(0));
  d.SetPar(1,c.getPar(1));
  d.process();
  d.Write(-1);

  double integ(0), dinteg(0);
  d.Integral(-1,3,integ,dinteg);
  cout<<"Integral between -1 and 3:\t"<<integ<<"\t"<<dinteg<<endl;

}

void drawExample(){
  TFile *fin=TFile::Open("o_gp_ex2.root","READ");
  TF1 *fS=new TF1("fS","gaus(0)",-6,6);
  fS->SetParameter(0,200);
  fS->SetParameter(1,1);
  fS->SetParameter(2,1);
  TF1 *fB = new TF1("fB","pol3(0)",-6,6);
  fB->SetParameters(1000,-0.2,0.2,-1);

  double x[30],y[30],dx[30],dy[30];
  for(int i=0;i<30;i++){
    double xTemp = -6 + i* 12./30;
    double bkg = fB->Eval(xTemp);
    double sig = fS->Eval(xTemp);

    x[i]=xTemp;
    y[i]=bkg + sig;
    dy[i]=sqrt(bkg+sig);
    dx[i]=0;
  }
  TGraphErrors *allPts = new TGraphErrors(30,x,y,dx,dy);
  allPts->SetMarkerStyle(20);
  allPts->SetMarkerColor(4);
  allPts->SetLineColor(4);
  allPts->GetXaxis()->SetLimits(-8,8);
  allPts->GetYaxis()->SetRangeUser(0,1300);

  TGraphErrors *hi=(TGraphErrors*)fin->Get("hi");
  hi->GetXaxis()->SetLimits(-8,8);

  TH1F *ho=(TH1F*)fin->Get("ho");

  TCanvas *c1=new TCanvas("c1","c1");
  allPts->Draw("AP");
  fB->SetLineColor(2);
  fB->SetLineWidth(3);
  fB->Draw("L&&same");
  fS->SetLineColor(3);
  fS->SetLineWidth(3);
  fS->Draw("L&&same");
  cout<<"background integral: "<<fB->Integral(-1,3)<<endl;

  TCanvas *c2=new TCanvas("c2","c2");
  hi->Draw("AP");
  ho->DrawCopy("same");
  hi->Draw("P");

  TCanvas *c1=new TCanvas("c3","c3");
  ho->DrawCopy();
  fB->Draw("L&&same");
  fS->Draw("L&&same");

  fin->Close();
  delete fin;
}

