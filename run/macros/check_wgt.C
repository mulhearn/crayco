{

  // input is from e.g.:
  // ./bin/rate_study 1 10000 1000 0.0 0.0 0.2


  gStyle->SetErrorX(0.0001);
  gStyle->SetOptStat(0);

  TFile *_file0 = TFile::Open("noise_control.root");
  TCanvas * c = new TCanvas();
  c->SetLogy();
  TH1F * hframe = new TH1F("hframe", "", 20, 0.0, 20.0);
  hframe->SetMinimum(1E-8);
  hframe->SetMaximum(1E5);

  h_noise_control_nhits_unwgt->SetLineColor(1);
  h_noise_control_nhits_wgt->SetLineColor(2);  
  g_noise_control_nhits->SetLineColor(4);

  h_noise_control_nhits_unwgt->SetMarkerColor(1);
  h_noise_control_nhits_wgt->SetMarkerColor(2);  

  h_noise_control_nhits_unwgt->SetMarkerStyle(21);
  h_noise_control_nhits_wgt->SetMarkerStyle(22);  

  hframe->Draw();
  h_noise_control_nhits_unwgt->Draw("epSAME");
  h_noise_control_nhits_wgt->Draw("epSAME");  
  g_noise_control_nhits->Draw("SAME");
}
