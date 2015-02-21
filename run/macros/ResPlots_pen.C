void ResPlots_pen(){
   enum {COLOR_A=3, COLOR_B=4, COLOR_C=1};
   enum {STYLE_A=1, STYLE_B=2, STYLE_C=3};
   enum {WIDTH=2};

   gStyle->SetOptStat(0);

  TFile *f_1000  = new TFile("plots_1000_final.root");
  TFile *f_5000  = new TFile("plots_5000_final.root");
  TFile *f_10000 = new TFile("plots_10000_final.root");


  TGraphErrors *geres_1000_8    = (TGraphErrors*)f_1000->Get("eres8");
  TGraphErrors *gtres_1000_8    = (TGraphErrors*)f_1000->Get("tres8");
  TGraphErrors *gpres_1000_8    = (TGraphErrors*)f_1000->Get("pres8");
  TGraphErrors *geres_1000_9    = (TGraphErrors*)f_1000->Get("eres9");
  TGraphErrors *gtres_1000_9    = (TGraphErrors*)f_1000->Get("tres9");
  TGraphErrors *gpres_1000_9    = (TGraphErrors*)f_1000->Get("pres9"); 

  TGraphErrors *geres_5000_8    = (TGraphErrors*)f_5000->Get("eres8");
  TGraphErrors *gtres_5000_8    = (TGraphErrors*)f_5000->Get("tres8");
  TGraphErrors *gpres_5000_8    = (TGraphErrors*)f_5000->Get("pres8");
  TGraphErrors *geres_5000_9    = (TGraphErrors*)f_5000->Get("eres9");
  TGraphErrors *gtres_5000_9    = (TGraphErrors*)f_5000->Get("tres9");
  TGraphErrors *gpres_5000_9    = (TGraphErrors*)f_5000->Get("pres9");

  TGraphErrors *geres_10000_8    = (TGraphErrors*)f_10000->Get("eres8");
  TGraphErrors *gtres_10000_8    = (TGraphErrors*)f_10000->Get("tres8");
  TGraphErrors *gpres_10000_8    = (TGraphErrors*)f_10000->Get("pres8");
  TGraphErrors *geres_10000_9    = (TGraphErrors*)f_10000->Get("eres9");
  TGraphErrors *gtres_10000_9    = (TGraphErrors*)f_10000->Get("tres9");
  TGraphErrors *gpres_10000_9    = (TGraphErrors*)f_10000->Get("pres9");

  geres_1000_8->SetLineColor(COLOR_A);
  gtres_1000_8->SetLineColor(COLOR_A);
  gpres_1000_8->SetLineColor(COLOR_A);
  geres_1000_9->SetLineColor(COLOR_A);
  gtres_1000_9->SetLineColor(COLOR_A);
  gpres_1000_9->SetLineColor(COLOR_A);
  //geres_5000_8->SetLineColor(COLOR_B);
  //gtres_5000_8->SetLineColor(COLOR_B);
  //gpres_5000_8->SetLineColor(COLOR_B);
  //geres_5000_9->SetLineColor(COLOR_B);
  //gtres_5000_9->SetLineColor(COLOR_B);
  //gpres_5000_9->SetLineColor(COLOR_B);
  geres_10000_8->SetLineColor(COLOR_C);
  gtres_10000_8->SetLineColor(COLOR_C);
  gpres_10000_8->SetLineColor(COLOR_C);
  geres_10000_9->SetLineColor(COLOR_C);
  gtres_10000_9->SetLineColor(COLOR_C);
  gpres_10000_9->SetLineColor(COLOR_C);

  geres_1000_8->SetLineWidth(WIDTH);
  gtres_1000_8->SetLineWidth(WIDTH);
  gpres_1000_8->SetLineWidth(WIDTH);
  geres_1000_9->SetLineWidth(WIDTH);
  gtres_1000_9->SetLineWidth(WIDTH);
  gpres_1000_9->SetLineWidth(WIDTH);

  //geres_5000_8->SetLineWidth(WIDTH);
  //gtres_5000_8->SetLineWidth(WIDTH);
  //gpres_5000_8->SetLineWidth(WIDTH);
  //geres_5000_9->SetLineWidth(WIDTH);
  //gtres_5000_9->SetLineWidth(WIDTH);
  //gpres_5000_9->SetLineWidth(WIDTH);

  geres_10000_8->SetLineWidth(WIDTH);
  gtres_10000_8->SetLineWidth(WIDTH);
  gpres_10000_8->SetLineWidth(WIDTH);
  geres_10000_9->SetLineWidth(WIDTH);
  gtres_10000_9->SetLineWidth(WIDTH);
  gpres_10000_9->SetLineWidth(WIDTH);

  ghits_1000_8->SetLineWidth(WIDTH);
  ghits_1000_9->SetLineWidth(WIDTH);
  ghits_5000_8->SetLineWidth(WIDTH);
  ghits_5000_9->SetLineWidth(WIDTH);
  ghits_10000_8->SetLineWidth(WIDTH);
  ghits_10000_9->SetLineWidth(WIDTH);

  geres_1000_8->SetLineStyle(STYLE_B);
  gtres_1000_8->SetLineStyle(STYLE_B); 
  gpres_1000_8->SetLineStyle(STYLE_B);
  geres_1000_9->SetLineStyle(STYLE_C);
  gtres_1000_9->SetLineStyle(STYLE_C);
  gpres_1000_9->SetLineStyle(STYLE_C);
  geres_5000_8->SetLineStyle(STYLE_B);
  gtres_5000_8->SetLineStyle(STYLE_B); 
  gpres_5000_8->SetLineStyle(STYLE_B);
  geres_5000_9->SetLineStyle(STYLE_C);
  gtres_5000_9->SetLineStyle(STYLE_C);
  gpres_5000_9->SetLineStyle(STYLE_C);
  geres_10000_8->SetLineStyle(STYLE_B);
  gtres_10000_8->SetLineStyle(STYLE_B); 
  gpres_10000_8->SetLineStyle(STYLE_B);
  geres_10000_9->SetLineStyle(STYLE_C);
  gtres_10000_9->SetLineStyle(STYLE_C);
  gpres_10000_9->SetLineStyle(STYLE_C);


  TGraphErrors *geres2_1000_8    = (TGraphErrors*)f_1000->Get("eres8");
  TGraphErrors *gtres2_1000_8    = (TGraphErrors*)f_1000->Get("tres8");
  TGraphErrors *gpres2_1000_8    = (TGraphErrors*)f_1000->Get("pres8");
  TGraphErrors *geres2_1000_9    = (TGraphErrors*)f_1000->Get("eres9");
  TGraphErrors *gtres2_1000_9    = (TGraphErrors*)f_1000->Get("tres9");
  TGraphErrors *gpres2_1000_9    = (TGraphErrors*)f_1000->Get("pres9"); 

  TGraphErrors *geres2_5000_8    = (TGraphErrors*)f_5000->Get("eres8");
  TGraphErrors *gtres2_5000_8    = (TGraphErrors*)f_5000->Get("tres8");
  TGraphErrors *gpres2_5000_8    = (TGraphErrors*)f_5000->Get("pres8");
  TGraphErrors *geres2_5000_9    = (TGraphErrors*)f_5000->Get("eres9");
  TGraphErrors *gtres2_5000_9    = (TGraphErrors*)f_5000->Get("tres9");
  TGraphErrors *gpres2_5000_9    = (TGraphErrors*)f_5000->Get("pres9");

  TGraphErrors *geres2_10000_8    = (TGraphErrors*)f_10000->Get("eres8");
  TGraphErrors *gtres2_10000_8    = (TGraphErrors*)f_10000->Get("tres8");
  TGraphErrors *gpres2_10000_8    = (TGraphErrors*)f_10000->Get("pres8");
  TGraphErrors *geres2_10000_9    = (TGraphErrors*)f_10000->Get("eres9");
  TGraphErrors *gtres2_10000_9    = (TGraphErrors*)f_10000->Get("tres9");
  TGraphErrors *gpres2_10000_9    = (TGraphErrors*)f_10000->Get("pres9");

  geres2_1000_8->SetLineColor(COLOR_A);
  gtres2_1000_8->SetLineColor(COLOR_A);
  gpres2_1000_8->SetLineColor(COLOR_A);
  geres2_1000_9->SetLineColor(COLOR_A);
  gtres2_1000_9->SetLineColor(COLOR_A);
  gpres2_1000_9->SetLineColor(COLOR_A);
  geres2_5000_8->SetLineColor(COLOR_B);
  gtres2_5000_8->SetLineColor(COLOR_B);
  gpres2_5000_8->SetLineColor(COLOR_B);
  geres2_5000_9->SetLineColor(COLOR_B);
  gtres2_5000_9->SetLineColor(COLOR_B);
  gpres2_5000_9->SetLineColor(COLOR_B);
  geres2_10000_8->SetLineColor(COLOR_C);
  gtres2_10000_8->SetLineColor(COLOR_C);
  gpres2_10000_8->SetLineColor(COLOR_C);
  geres2_10000_9->SetLineColor(COLOR_C);
  gtres2_10000_9->SetLineColor(COLOR_C);
  gpres2_10000_9->SetLineColor(COLOR_C);

  geres2_1000_8->SetLineWidth(WIDTH);
  gtres2_1000_8->SetLineWidth(WIDTH);
  gpres2_1000_8->SetLineWidth(WIDTH);
  geres2_1000_9->SetLineWidth(WIDTH);
  gtres2_1000_9->SetLineWidth(WIDTH);
  gpres2_1000_9->SetLineWidth(WIDTH);
  geres2_5000_8->SetLineWidth(WIDTH);
  gtres2_5000_8->SetLineWidth(WIDTH);
  gpres2_5000_8->SetLineWidth(WIDTH);
  geres2_5000_9->SetLineWidth(WIDTH);
  gtres2_5000_9->SetLineWidth(WIDTH);
  gpres2_5000_9->SetLineWidth(WIDTH);
  geres2_10000_8->SetLineWidth(WIDTH);
  gtres2_10000_8->SetLineWidth(WIDTH);
  gpres2_10000_8->SetLineWidth(WIDTH);
  geres2_10000_9->SetLineWidth(WIDTH);
  gtres2_10000_9->SetLineWidth(WIDTH);
  gpres2_10000_9->SetLineWidth(WIDTH);

  geres2_1000_8->SetLineStyle(1);
  gtres2_1000_8->SetLineStyle(1);
  gpres2_1000_8->SetLineStyle(1);
  geres2_1000_9->SetLineStyle(1);
  gtres2_1000_9->SetLineStyle(1);
  gpres2_1000_9->SetLineStyle(1);
  geres2_5000_8->SetLineStyle(1);
  gtres2_5000_8->SetLineStyle(1);
  gpres2_5000_8->SetLineStyle(1);
  geres2_5000_9->SetLineStyle(1);
  gtres2_5000_9->SetLineStyle(1);
  gpres2_5000_9->SetLineStyle(1);
  geres2_10000_8->SetLineStyle(1);
  gtres2_10000_8->SetLineStyle(1);
  gpres2_10000_8->SetLineStyle(1);
  geres2_10000_9->SetLineStyle(1);
  gtres2_10000_9->SetLineStyle(1);
  gpres2_10000_9->SetLineStyle(1);

  TLegend *leg1 = new TLegend(0.12,0.12,0.50,0.40);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);

  TH1F* dummy1 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy2 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy3 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy4 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy5 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy6 = new TH1F("dummy","",100, 0, 10);
  
  dummy1->SetLineWidth(WIDTH);
  dummy2->SetLineWidth(WIDTH);
  dummy3->SetLineWidth(WIDTH);
  dummy4->SetLineWidth(WIDTH);
  dummy5->SetLineWidth(WIDTH);
  dummy6->SetLineWidth(WIDTH);

  dummy1->SetLineColor(COLOR_A);
  leg1->AddEntry(dummy1, "1000/km^{2}", "l");
  dummy2->SetLineColor(COLOR_B);
  leg1->AddEntry(dummy2, "5000/km^{2}", "l");
  dummy3->SetLineColor(COLOR_C);
  leg1->AddEntry(dummy3, "10000/km^{2}", "l");

  dummy4->SetLineStyle(STYLE_A);
  dummy4->SetLineColor(1);
  //leg1->AddEntry(dummy4, "#epsilon A = 1e-7");

  dummy5->SetLineStyle(STYLE_B);
  dummy5->SetLineColor(1);
  leg1->AddEntry(dummy5, "#epsilon A = 1e-8");
  dummy6->SetLineStyle(STYLE_C);
  dummy6->SetLineColor(1);
  leg1->AddEntry(dummy6, "#epsilon A = 1e-9");
  
  TCanvas *ce = new TCanvas("ce");
  TCanvas *ct = new TCanvas("ct");
  TCanvas *cp = new TCanvas("cp");
  TCanvas *ch = new TCanvas("ch");
  
  ce->cd();
  ce->SetLogx();
  TH2F * hframe = new TH2F("hframe", "", 10, 2E18, 1E21, 10, 0.0, 0.75);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]");
  hframe->GetXaxis()->SetTitleOffset(1.3);
  hframe->GetYaxis()->SetTitle("Fractional Energy Resolution");
  hframe->GetYaxis()->SetTitleOffset(1.2);
  hframe->Draw();
  geres2_1000_8 ->Draw("P");
  geres2_1000_9 ->Draw("P");
  geres2_5000_8 ->Draw("P");
  geres2_5000_9 ->Draw("P");
  geres2_10000_8->Draw("P");
  geres2_10000_9->Draw("P");
  geres_1000_8 ->Draw("LP");
  geres_1000_9 ->Draw("LP");
  geres_5000_8 ->Draw("LP");
  geres_5000_9 ->Draw("LP");
  geres_10000_8->Draw("LP");
  geres_10000_9->Draw("LP");
  leg1->Draw();
  ce->SaveAs("energy_res.png");

  return 0;
   
  ct->cd();
  ct->SetLogx();
  TH2F * hframe = new TH2F("hframe", "", 10, 2E18, 1E21, 10, 0.0, 0.23);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]");
  hframe->GetXaxis()->SetTitleOffset(1.3);
  hframe->GetYaxis()->SetTitle("Theta Resolution [Radians]");
  hframe->GetYaxis()->SetTitleOffset(1.2);
  hframe->Draw();
  gtres_1000_8 ->Draw("LP");
  gtres_1000_9 ->Draw("LP");
  gtres_5000_8 ->Draw("LP");
  gtres_5000_9 ->Draw("LP");
  gtres_10000_8->Draw("LP");
  gtres_10000_9->Draw("LP");
  gtres2_1000_8 ->Draw("P");
  gtres2_1000_9 ->Draw("P");
  gtres2_5000_8 ->Draw("P");
  gtres2_5000_9 ->Draw("P");
  gtres2_10000_8->Draw("P");
  gtres2_10000_9->Draw("P");
  leg1->Draw();
  ct->SaveAs("theta_res.png");


  cp->cd();
  cp->SetLogx();
  TH2F * hframe = new TH2F("hframe", "", 10, 2E18, 1E21, 10, 0.0, 0.9);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]");
  hframe->GetXaxis()->SetTitleOffset(1.3);
  hframe->GetYaxis()->SetTitle("Phi Resolution [Radians]");
  hframe->GetYaxis()->SetTitleOffset(1.2);
  hframe->Draw();
  gpres_1000_8 ->Draw("LP");
  gpres_1000_9 ->Draw("LP");
  gpres_5000_8 ->Draw("LP");
  gpres_5000_9 ->Draw("LP");
  gpres_10000_8->Draw("LP");
  gpres_10000_9->Draw("LP");
  gpres2_1000_8 ->Draw("P");
  gpres2_1000_9 ->Draw("P");
  gpres2_5000_8 ->Draw("P");
  gpres2_5000_9 ->Draw("P");
  gpres2_10000_8->Draw("P");
  gpres2_10000_9->Draw("P");
  leg1->Draw();
  cp->SaveAs("phi_res.png");

}
