void do_nhits();
void do_ereco();

void mkplots_rate(){
  //do_nhits();
  do_ereco();
}


void do_nhits(){
  gStyle->SetOptStat(0);
  int MAXHITS = 50;
  //int MAXHITS = 10;
  TCanvas * c = new TCanvas();
  TH1F * hframe = new TH1F("hframe", "", MAXHITS, 0., MAXHITS);

  TH1F * h1 = new TH1F("h1", "", MAXHITS, 0., MAXHITS);
  TH1F * h2 = new TH1F("h2", "", MAXHITS, 0., MAXHITS);
  TH1F * h3 = new TH1F("h3", "", MAXHITS, 0., MAXHITS);
  TH1F * h4 = new TH1F("h4", "", MAXHITS, 0., MAXHITS);
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h4->Sumw2();


  c->SetLogy();
  



  _file0->ls();
  TTree * rate = _file0->Get("rate");
  rate->Print();

  rate->Draw("nhit>>h1","wgt*(sample==0)");
  rate->Draw("nhit>>h2","wgt*(sample==1)");
  rate->Draw("nhit>>h3","wgt*(sample==2)");
  rate->Draw("nhit>>h4","wgt*(sample==3)");

  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(4);
  h4->SetLineColor(3);

  hframe->SetMaximum(10.0);
  //hframe->SetMinimum(1E-13);
  hframe->SetMinimum(1E-25);
  hframe->Draw();
  //h1->Draw("SAME");
  //h2->Draw("SAME");
  //h3->Draw("SAME");
  h1->Draw("HSAME");
  h2->Draw("epSAME");
  h3->Draw("epSAME");
  h4->Draw("epSAME");

  double r1 = h1->GetSumOfWeights();
  double r2 = h2->GetSumOfWeights();
  double r3 = h3->GetSumOfWeights();
  double r4 = h4->GetSumOfWeights();

  cout << "integrated total rate sample 1 (Hz/km2):  " << r1 << "\n";
  cout << "integrated total rate sample 2 (Hz/km2):  " << r2 << "\n";
  cout << "integrated total rate sample 3 (Hz/km2):  " << r3 << "\n";
  cout << "integrated total rate sample 4 (Hz/km2):  " << r4 << "\n";

  cout << "integrated total rate sample 1 (Hz/m2):  " << r1*1E-6 << "\n";
  cout << "integrated total rate sample 2 (Hz/m2):  " << r2*1E-6 << "\n";
  cout << "integrated total rate sample 3 (Hz/m2):  " << r3*1E-6 << "\n";
  cout << "integrated total rate sample 4 (Hz/m2):  " << r4*1E-6 << "\n";

  cout << "integrated total rate sample 1 (Hz/m2/sr):  " << r1*1E-6/(2*3.1415) << "\n";
  cout << "integrated total rate sample 2 (Hz/m2/sr):  " << r2*1E-6/(2*3.1415) << "\n";
  cout << "integrated total rate sample 3 (Hz/m2/sr):  " << r3*1E-6/(2*3.1415) << "\n";
  cout << "integrated total rate sample 4 (Hz/m2/sr):  " << r4*1E-6/(2*3.1415) << "\n";
}


void do_ereco(){
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas();

  double emin = 18.0;
  double emax = 22.0;
  int nbin = 50;
  TH1F * hframe = new TH1F("hframe", "", nbin, emin, emax);
  TH1F * h1 = new TH1F("h1", "", nbin, emin, emax);
  TH1F * h2 = new TH1F("h2", "", nbin, emin, emax);
  TH1F * h3 = new TH1F("h3", "", nbin, emin, emax);
  TH1F * h4 = new TH1F("h4", "", nbin, emin, emax);
  h1->Sumw2();
  h2->Sumw2();
  h3->Sumw2();
  h4->Sumw2();


  c->SetLogy();
  
  _file0->ls();
  TTree * rate = _file0->Get("rate");
  rate->Print();

  rate->Draw("ereco>>h1","wgt*(sample==0)");
  rate->Draw("ereco>>h2","wgt*(sample==1)");
  rate->Draw("ereco>>h3","wgt*(sample==2)");
  rate->Draw("ereco>>h4","wgt*(sample==3)");

  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(4);
  h4->SetLineColor(3);

  hframe->SetMaximum(10.0);
  hframe->SetMinimum(1E-25);
  hframe->Draw();

  h1->Draw("HSAME");
  h2->Draw("epSAME");
  h3->Draw("epSAME");
  h4->Draw("epSAME");
}
