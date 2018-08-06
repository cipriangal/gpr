void gp()
{
  gSystem->Load("../install/lib/libGausProc.so");

  int verbosity=0;

  char *outfile="o_gp.root";
  const unsigned nData=30;

  double _x[30]={0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,6.1,6.3,6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1,8.3,8.5,8.7,8.9,9.1,9.3,9.5,9.7,9.9};
  double _y[30]={351,371,331,352,330,323,322,356,333,325,205,204,203,233,220,225,217,238,223,247,249,258,284,285,323,301,360,366,407,402};
  double _dy[30]={18.735,19.2614,18.1934,18.7617,18.1659,17.9722,17.9444,18.868,18.2483,18.0278,14.3178,14.2829,14.2478,15.2643,14.8324,15,14.7309,15.4272,14.9332,15.7162,15.7797,16.0624,16.8523,16.8819,17.9722,17.3494,18.9737,19.1311,20.1742,20.0499};
  std::vector<double> x,y,sigma_y;
  x.resize(nData);
  y.resize(nData);
  sigma_y.resize(nData);
  for(int i=0;i<nData;i++){
    x[i]=_x[i];
    y[i]=_y[i];
    sigma_y[i]=_dy[i];
  }
  double xmin=-0.5;
  double xmax=11;
  unsigned nPredictions=1000;
  
  GausProc a(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  gSystem->Exec(Form("rm -f %s",outfile));
  
  for(int i=0;i<x.size();i++)
    cout<<x[i]<<" "<<y[i]<<" "<<sigma_y[i]<<endl;

  a.SetVerbosity(verbosity);
  a.SetKernel(GausProc::RBF);
  //  a.warp(0);//should be used only if you have data spanning several orders of magnitude
  
  GPOptimizer c(&a,2,10);
  c.GPoptimize(2,0);
  cout<<c.getPar(0)<<" "<<c.getPar(1)<<endl;  
  
  GausProc d(x, y, sigma_y, xmin, xmax, nPredictions, outfile);
  //d.warp(0);//see comment above
  d.SetPar(0,c.getPar(0));
  d.SetPar(1,c.getPar(1));
  d.process();
  d.Write(-1);
  //d.unwarp(0);//see comment above
}

