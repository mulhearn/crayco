


void do_nhits();
void do_ereco();

void mkplots_rate(){
  //do_nhits();
  
  
  plot_var("nhit", "Number of Hits", 10, 0.0, 10.0, 1E-15, 100.0, "1",0);
  //plot_var("nhit", "Number of Hits", 20, 0.0, 50.0, 1E-15, 100.0, "1",0);
  //plot_var("mcore", 100, 0.0, 100.0, 1E-15, 100.0, "1",0);
  plot_var("llr", "Log Likelihood Ratio", 20, 0.0, 1.0, 1E-15, 10E-5, "(llr>0.0)",0);
  plot_var("ereco", "Reconstructed Energy [eV]", 10, 19.5, 21.0, 1E-17, 1E-8.0, "(llr>0.07)",1);
  //plot_var("ereco", 10, 0.0, 21.0, 1E-17, 1E2.0, "1");
  //plot_var("eshower", 50, 0.0, 25.0, 1E-15, 100.0, "(nhit>7)");
}


TGraphErrors * get_log10x_graph(TH1F * h){
  int max = h->GetNbinsX();
  double * x  = new double[max];
  double * y  = new double[max];
  double * dx = new double[max];
  double * dy = new double[max];
  int cnt = 0;
  
  for (int i=1; i<=h->GetNbinsX(); i++){
    double log10_x = h->GetBinLowEdge(i) + (0.5 * (h->GetBinWidth(i)));
    x[cnt]  = pow(10.0, log10_x);
    dx[cnt] = 0.0;
    y[cnt]  = h->GetBinContent(i);
    dy[cnt] = h->GetBinError(i);
    if (y[cnt] > 0.0){
      cout << "x:  " << x[cnt] << "y:  " << y[cnt] << " dy:  " << dy[cnt] << "\n";
      cnt++;
    }
  }
  cout << "INFO: graph will have " << cnt << " entries.\n";
  if (cnt > 0){
    return new TGraphErrors(cnt, x, y, dx, dy);
  } else {
    return NULL;
  }
}

void plot_var(const char * var, const char * xtitle, double nbin, double xmin, double xmax, double ymin, double ymax, const char * cuts, int mode_logx){
  enum {COLOR_A=1, COLOR_B=2, COLOR_C=4};
  enum {STYLE_A=20, STYLE_B=21, STYLE_C=22};
  enum {WIDTH=2};
  enum {LEGWIDTH=3};
  const double MSIZE=1.5;

  TH1F* dummy1 = new TH1F("dummy","",100, 0, 10);
  //TH1F* dummy2 = new TH1F("dummy","",100, 0, 10);
  //TH1F* dummy3 = new TH1F("dummy","",100, 0, 10);
  //TH1F* dummy4 = new TH1F("dummy","",100, 0, 10);
  //TH1F* dummy5 = new TH1F("dummy","",100, 0, 10);
  //TH1F* dummy6 = new TH1F("dummy","",100, 0, 10);
  
  //dummy1->SetLineWidth(LEGWIDTH);
  //dummy1->SetMarkerSize(MSIZE);
  //dummy1->SetLineColor(COLOR_A);
  //leg1->AddEntry(dummy1, "500 devices/km^{2}", "l");
  //leg1->AddEntry(dummy6, "#epsilon A = 1e-9 (#gamma), 1e-5 (#mu)", "p");
  //leg1->AddEntry(dummy5, "#epsilon A = 5e-9 (#gamma), 5e-5 (#mu)", "p");
  //leg1->AddEntry(dummy4, "#epsilon A = 1e-7");    

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas();

  //double emin = 18.0;
  //double emax = 22.0;
  //int nbin = 50;
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


  h0->SetLineColor(COLOR_A);
  h0->SetMarkerColor(COLOR_A);
  h0->SetLineWidth(WIDTH);
  h0->SetMarkerSize(MSIZE);
  h0->SetMarkerStyle(STYLE_A);

  h1->SetLineColor(COLOR_B);
  h1->SetMarkerColor(COLOR_B);
  h1->SetLineWidth(WIDTH);
  h1->SetMarkerSize(MSIZE);
  h1->SetMarkerStyle(STYLE_B);

  h2->SetLineColor(COLOR_C);
  h2->SetMarkerColor(COLOR_C);
  h2->SetLineWidth(WIDTH);
  h2->SetMarkerSize(MSIZE);
  h2->SetMarkerStyle(STYLE_C);

  TLegend *leg = new TLegend(0.60,0.62,0.90,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h0, "Comb Bkg", "p");
  leg->AddEntry(h1, "E > 5 #times 10^{19} eV",      "p");
  leg->AddEntry(h2, "E > 1 #times 10^{20} eV",      "p");

  if (! mode_logx){
    TH1F * hframe = new TH1F("hframe", "", nbin, xmin, xmax);
    hframe->SetMaximum(ymax);
    hframe->SetMinimum(ymin);
    hframe->Draw();
    c->SetBottomMargin(0.12);
    hframe->GetXaxis()->SetTitleSize(0.06); 
    hframe->GetXaxis()->SetTitleOffset(0.95);
    hframe->GetYaxis()->SetTitleOffset(0.75);  
    hframe->GetYaxis()->SetTitleSize(0.06);


    hframe->SetXTitle(xtitle);
    hframe->SetYTitle("Rate [Hz/km^{2}]");

    h0->Draw("HPSAME");
    h1->Draw("HPSAME");
    h2->Draw("HPSAME");
    //h3->Draw("HPSAME");

    leg->Draw();

  } else { 
    //cout << "INFO: logx mode not yet prepared\n" << "\n";
    TGraphErrors * g0 = get_log10x_graph(h0);
    TGraphErrors * g1 = get_log10x_graph(h1);
    TGraphErrors * g2 = get_log10x_graph(h2);
    TGraphErrors * g3 = get_log10x_graph(h3);
    c->SetLogx();


    if (g0){
      g0->SetLineColor(COLOR_A);
      g0->SetMarkerColor(COLOR_A);
      g0->SetLineWidth(WIDTH);
      g0->SetMarkerSize(MSIZE);
      g0->SetMarkerStyle(STYLE_A);
    }
    if (g1){
      g1->SetLineColor(COLOR_B);
      g1->SetMarkerColor(COLOR_B);
      g1->SetLineWidth(WIDTH);
      g1->SetMarkerSize(MSIZE);
      g1->SetMarkerStyle(STYLE_B);
    }
    if (g2){
      g2->SetLineColor(COLOR_C);
      g2->SetMarkerColor(COLOR_C);
      g2->SetLineWidth(WIDTH);
      g2->SetMarkerSize(MSIZE);
      g2->SetMarkerStyle(STYLE_C);
    }

    TH1F * hframe = new TH1F("hframe", "", nbin, pow(10.0, xmin), pow(10.0, xmax));
    hframe->SetMaximum(ymax);
    hframe->SetMinimum(ymin);


    c->SetBottomMargin(0.12);
    hframe->GetXaxis()->SetTitleSize(0.06); 
    hframe->GetXaxis()->SetTitleOffset(0.95);
    hframe->GetYaxis()->SetTitleOffset(0.75);  
    hframe->GetYaxis()->SetTitleSize(0.06);

    hframe->Draw();
    hframe->SetYTitle("Rate [Hz/km^{2}]");
    hframe->SetXTitle(xtitle);

    hframe->Draw();
    if (g0){ g0->Draw("PSAME"); }
    if (g1){ g1->Draw("LPESAME"); }
    if (g2){ g2->Draw("LPSAME"); }
    if (g3){ g3->Draw("LPSAME"); }
    
    leg->Draw();

  }
  
  char name[100];
  sprintf(name, "%s.pdf", var);
  c->SaveAs(name);
  


}
