void example3()
{
  gSystem->Load("../install/lib/libGausProc.so");

  int verbosity=0;

  char *outfile="o_gp_ex3.root";
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
  double xmin=0;
  double xmax=7;
  unsigned nPredictions=400;
  
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

  cout<<"Cholensky decompostion starting:"<<endl;
  d.DecomposeChol();
  cout<<"Cholensky decompostion finished"<<endl;

  TFile *fout=TFile::Open(outfile,"UPDATE");
  const int nFct=300;
  TGraph *gr[nFct];
  for(int i=0;i<nFct;i++){
    //cout<<"function "<<i<<endl;
    gr[i]=new TGraph();
    std::vector<double> x,y;
    int props;
    d.GetOneFunction(i+300,x,y,props);
    for(int j=0;j<x.size();j++)
      gr[i]->SetPoint(j,x[j],y[j]);
    gr[i]->SetName(Form("oneFct_%d",i));
    fout->cd();
    gr[i]->Write();
  }
  fout->Close();

}

void drawExample(){
  TFile *fin=TFile::Open("o_gp_ex3.root","READ");
  TGraphErrors *hi = (TGraphErrors*) fin->Get("hi");
  TH1F *ho= (TH1F*) fin->Get("ho");
  hi->GetXaxis()->SetLimits(0,8);
  hi->GetYaxis()->SetRangeUser(-2.5,2);
  hi->Draw("AP");
  ho->SetMarkerStyle(21);
  ho->DrawCopy("same");
  hi->Draw("P");

  TCanvas *c2=new TCanvas("c2","c2");
  hi->Draw("AP");
  TGraph *gr;
  for(int i=0;i<5;i++){
    gr=(TGraph*)fin->Get(Form("oneFct_%d",i));
    gr->SetLineColor(i);
    gr->SetLineWidth(2);
    gr->Draw("L");
  }
  hi->Draw("P");
  gStyle->SetOptStat(0);

  TCanvas *c3=new TCanvas("c3","c3");
  hi->Draw("AP");
  TH2D *h2d = new TH2D("h2d","heat map histo",200,0,8,200,-2.5,2);
  TGraph *gr;
  for(int i=0;i<300;i++){
    gr=(TGraph*)fin->Get(Form("oneFct_%d",i));
    int npts = gr->GetN();
    double x(0),y(0);
    for(int j=0;j<npts;j++){
      gr->GetPoint(j,x,y);
      h2d->Fill(x,y);
    }
  }
  h2d->DrawCopy("colz");
  hi->Draw("P");

  TCanvas *c4=new TCanvas("c4","c4");
  h2d->DrawCopy("colz");
  ho->DrawCopy("same");

  fin->Close();
  delete fin;
}




