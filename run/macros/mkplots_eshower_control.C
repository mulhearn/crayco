


void scale_log(TH1F * h){
  static const double log10 = log(10.0);
  h->Print();
  int nbin = h->GetNbinsX();
  cout << nbin << "\n";
  for (int i=1; i<=nbin; i++){
    double a = h->GetBinLowEdge(i);
    double b = a + h->GetBinWidth(i);
    double Ea = exp(a * log10);
    double Eb = exp(b * log10);
    double dE = Eb - Ea;
    double x  = h->GetBinContent(i);
    double dx = h->GetBinError(i);
    if (dE > 0.0){
      h->SetBinContent(i, (1E9)*(x/dE)/1E6/(2*3.1415));
      h->SetBinError(i, (1E9)*(dx/dE)/1E6/(2*3.1415));
    }
  }
}



mkplots_eshower_control(){
  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas();
  c->SetLogy();
  TFile *_file0 = TFile::Open("eshower_control.root");
  TH1F * h = (TH1F *) _file0->Get("h_eshower_control_eshower_wgt");
  //h->Draw("H");

  // the cumulative version will be automatically normalized correctly, since each event is weighted correctly.
  TH1F * hc = (TH1F *) h->GetCumulative(kFALSE); // kFALSE to set direction of integral
  hc->Scale(1E-6);
  hc->Draw("H");

  //TH1F * hd = (TH1F *) h->Clone("hd");
  //scale_log(hd);
  //hd->Draw("H");
}
