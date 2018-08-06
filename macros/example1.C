void example1()
{
  gSystem->Load("../install/lib/libGausProc.so");

  int verbosity=0;

  char *outfile="o_gp_ex1.root";
  const unsigned nData=6;
  double _x[nData] = {1,2,3,4,5,6};
  double _y[nData] = {-1.8, -1.3, -0.5, 0.2, 0.8, 1.3};
  double _dy[nData] = {0.1,0.1,0.1,0.1,0.1,0.1};
  std::vector<double> x,y,sigma_y;
  for(int i=0;i<nData;i++){
    x.push_back(_x[i]);
    y.push_back(_y[i]);
    sigma_y.push_back(_dy[i]);
  }
  double xmin=7;
  double xmax=7;
  unsigned nPredictions=1;
  
  GausProc a(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  gSystem->Exec(Form("rm -f %s",outfile));
  
  cout<<"Input data:"<<endl;
  for(int i=0;i<x.size();i++)
    cout<<x[i]<<" "<<y[i]<<" "<<sigma_y[i]<<endl;

  a.SetVerbosity(verbosity);
  a.SetKernel(GausProc::RBF);
  
  GPOptimizer c(&a,2,10);
  c.GPoptimize(2,0);
  cout<<c.getPar(0)<<" "<<c.getPar(1)<<endl;  
  
  GausProc d(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  d.SetPar(0,c.getPar(0));
  d.SetPar(1,c.getPar(1));
  d.process();
  d.Write(-1);
}

void drawExample(){
  TFile *fin=TFile::Open("o_gp_ex1.root","READ");
  TGraphErrors *hi = (TGraphErrors*) fin->Get("hi");
  TH1F *ho= (TH1F*) fin->Get("ho");
  hi->GetXaxis()->SetLimits(0,8);
  hi->GetYaxis()->SetRangeUser(-2.5,2);
  hi->Draw("AP");
  ho->SetMarkerStyle(21);
  ho->DrawCopy("same");
  fin->Close();
  delete fin;
}

