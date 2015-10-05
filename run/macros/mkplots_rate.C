void do_nhits();
void do_ereco();

void mkplots_rate(){
  //do_nhits();
  plot_var("nhit", 12, 0.0, 12.0, "1");
  plot_var("nhit", 50, 0.0, 100.0, "1");
  //plot_var("mcore", 50, 0.0, 100.0, "1");
  //plot_var("llr", 50, 0.0, 1.0, "1");
  //plot_var("ereco", 50, 18.0, 22.0, "(llr>0.07)");
  plot_var("eshower", 50, 0.0, 25.0, "(nhit>7)");
}

void plot_var(const char * var, double nbin, double xmin, double xmax, const char * cuts){
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas();

  //double emin = 18.0;
  //double emax = 22.0;
  //int nbin = 50;
  TH1F * hframe = new TH1F("hframe", "", nbin, xmin, xmax);
  TH1F * h0 = new TH1F("h0", "", nbin, xmin, xmax);
  TH1F * h1 = new TH1F("h1", "", nbin, xmin, xmax);
  TH1F * h2 = new TH1F("h2", "", nbin, xmin, xmax);
  TH1F * h3 = new TH1F("h3", "", nbin, xmin, xmax);
  h0->Sumw2();
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();


  c->SetLogy();
  
  _file0->ls();
  TTree * rate = _file0->Get("rate");
  rate->Print();


  char cmda[100];
  char cmdb[100];

  for (int i=0; i<4; i++){
    sprintf(cmda, "%s>>h%d", var, i);
    sprintf(cmdb, "wgt*((sample==%d)&&%s)", i, cuts);
    cout << cmda << "\n";
    cout << cmdb << "\n";
    rate->Draw(cmda,cmdb);
  }

  h0->SetLineColor(1);
  h1->SetLineColor(2);
  h2->SetLineColor(4);
  h3->SetLineColor(3);

  hframe->SetMaximum(100.0);
  hframe->SetMinimum(1E-20);
  hframe->Draw();
  hframe->SetXTitle(var);


  h0->Draw("HSAME");
  h1->Draw("epSAME");
  h2->Draw("epSAME");
  h3->Draw("epSAME");

  char name[100];
  sprintf(name, "%s.pdf", var);
  c->SaveAs(name);

  double r0 = h0->GetSumOfWeights();
  double r1 = h1->GetSumOfWeights();
  double r2 = h2->GetSumOfWeights();
  double r3 = h3->GetSumOfWeights();

  cout << "integrated total rate sample 0 (Hz/km2):  " << r0 << "\n";
  cout << "integrated total rate sample 1 (Hz/km2):  " << r1 << "\n";
  cout << "integrated total rate sample 2 (Hz/km2):  " << r2 << "\n";
  cout << "integrated total rate sample 3 (Hz/km2):  " << r3 << "\n";

  cout << "integrated total rate sample 0 (Hz/m2):  " << r0*1E-6 << "\n";
  cout << "integrated total rate sample 1 (Hz/m2):  " << r1*1E-6 << "\n";
  cout << "integrated total rate sample 2 (Hz/m2):  " << r2*1E-6 << "\n";
  cout << "integrated total rate sample 3 (Hz/m2):  " << r3*1E-6 << "\n";

  cout << "integrated total rate sample 0 (Hz/m2/sr):  " << r0*1E-6/(2*3.1415) << "\n";
  cout << "integrated total rate sample 1 (Hz/m2/sr):  " << r1*1E-6/(2*3.1415) << "\n";
  cout << "integrated total rate sample 2 (Hz/m2/sr):  " << r2*1E-6/(2*3.1415) << "\n";
  cout << "integrated total rate sample 3 (Hz/m2/sr):  " << r3*1E-6/(2*3.1415) << "\n";


}
