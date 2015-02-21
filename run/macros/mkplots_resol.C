void mkplots_resol(){
   enum {COLOR_A=2, COLOR_B=4, COLOR_C=6};
   enum {STYLE_A=2, STYLE_B=24, STYLE_C=20};
   enum {WIDTH=2};
   enum {LEGWIDTH=3};
   const double MSIZE=1.5;

   gStyle->SetOptStat(0);

  TFile *f_1000  = new TFile("resplots_500.root");
  TFile *f_5000  = new TFile("resplots_1000.root");
  TFile *f_10000 = new TFile("resplots_5000.root");

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
  geres_5000_8->SetLineColor(COLOR_B);
  gtres_5000_8->SetLineColor(COLOR_B);
  gpres_5000_8->SetLineColor(COLOR_B);
  geres_5000_9->SetLineColor(COLOR_B);
  gtres_5000_9->SetLineColor(COLOR_B);
  gpres_5000_9->SetLineColor(COLOR_B);
  geres_10000_8->SetLineColor(COLOR_C);
  gtres_10000_8->SetLineColor(COLOR_C);
  gpres_10000_8->SetLineColor(COLOR_C);
  geres_10000_9->SetLineColor(COLOR_C);
  gtres_10000_9->SetLineColor(COLOR_C);
  gpres_10000_9->SetLineColor(COLOR_C);

  geres_1000_8->SetMarkerColor(COLOR_A);
  gtres_1000_8->SetMarkerColor(COLOR_A);
  gpres_1000_8->SetMarkerColor(COLOR_A);
  geres_1000_9->SetMarkerColor(COLOR_A);
  gtres_1000_9->SetMarkerColor(COLOR_A);
  gpres_1000_9->SetMarkerColor(COLOR_A);
  geres_5000_8->SetMarkerColor(COLOR_B);
  gtres_5000_8->SetMarkerColor(COLOR_B);
  gpres_5000_8->SetMarkerColor(COLOR_B);
  geres_5000_9->SetMarkerColor(COLOR_B);
  gtres_5000_9->SetMarkerColor(COLOR_B);
  gpres_5000_9->SetMarkerColor(COLOR_B);
  geres_10000_8->SetMarkerColor(COLOR_C);
  gtres_10000_8->SetMarkerColor(COLOR_C);
  gpres_10000_8->SetMarkerColor(COLOR_C);
  geres_10000_9->SetMarkerColor(COLOR_C);
  gtres_10000_9->SetMarkerColor(COLOR_C);
  gpres_10000_9->SetMarkerColor(COLOR_C);

  geres_1000_8->SetLineWidth(WIDTH);
  gtres_1000_8->SetLineWidth(WIDTH);
  gpres_1000_8->SetLineWidth(WIDTH);
  geres_1000_9->SetLineWidth(WIDTH);
  gtres_1000_9->SetLineWidth(WIDTH);
  gpres_1000_9->SetLineWidth(WIDTH);
  geres_5000_8->SetLineWidth(WIDTH);
  gtres_5000_8->SetLineWidth(WIDTH);
  gpres_5000_8->SetLineWidth(WIDTH);
  geres_5000_9->SetLineWidth(WIDTH);
  gtres_5000_9->SetLineWidth(WIDTH);
  gpres_5000_9->SetLineWidth(WIDTH);
  geres_10000_8->SetLineWidth(WIDTH);
  gtres_10000_8->SetLineWidth(WIDTH);
  gpres_10000_8->SetLineWidth(WIDTH);
  geres_10000_9->SetLineWidth(WIDTH);
  gtres_10000_9->SetLineWidth(WIDTH);
  gpres_10000_9->SetLineWidth(WIDTH);


  geres_1000_8->SetMarkerSize(MSIZE);
  gtres_1000_8->SetMarkerSize(MSIZE);
  gpres_1000_8->SetMarkerSize(MSIZE);
  geres_1000_9->SetMarkerSize(MSIZE);
  gtres_1000_9->SetMarkerSize(MSIZE);
  gpres_1000_9->SetMarkerSize(MSIZE);
  geres_5000_8->SetMarkerSize(MSIZE);
  gtres_5000_8->SetMarkerSize(MSIZE);
  gpres_5000_8->SetMarkerSize(MSIZE);
  geres_5000_9->SetMarkerSize(MSIZE);
  gtres_5000_9->SetMarkerSize(MSIZE);
  gpres_5000_9->SetMarkerSize(MSIZE);
  geres_10000_8->SetMarkerSize(MSIZE);
  gtres_10000_8->SetMarkerSize(MSIZE);
  gpres_10000_8->SetMarkerSize(MSIZE);
  geres_10000_9->SetMarkerSize(MSIZE);
  gtres_10000_9->SetMarkerSize(MSIZE);
  gpres_10000_9->SetMarkerSize(MSIZE);

  geres_1000_8->SetMarkerStyle(STYLE_B);
  gtres_1000_8->SetMarkerStyle(STYLE_B); 
  gpres_1000_8->SetMarkerStyle(STYLE_B);
  geres_1000_9->SetMarkerStyle(STYLE_C);
  gtres_1000_9->SetMarkerStyle(STYLE_C);
  gpres_1000_9->SetMarkerStyle(STYLE_C);
  geres_5000_8->SetMarkerStyle(STYLE_B);
  gtres_5000_8->SetMarkerStyle(STYLE_B); 
  gpres_5000_8->SetMarkerStyle(STYLE_B);
  geres_5000_9->SetMarkerStyle(STYLE_C);
  gtres_5000_9->SetMarkerStyle(STYLE_C);
  gpres_5000_9->SetMarkerStyle(STYLE_C);
  geres_10000_8->SetMarkerStyle(STYLE_B);
  gtres_10000_8->SetMarkerStyle(STYLE_B); 
  gpres_10000_8->SetMarkerStyle(STYLE_B);
  geres_10000_9->SetMarkerStyle(STYLE_C);
  gtres_10000_9->SetMarkerStyle(STYLE_C);
  gpres_10000_9->SetMarkerStyle(STYLE_C);

  TLegend *leg1 = new TLegend(0.12,0.14,0.50,0.42);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);

  TH1F* dummy1 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy2 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy3 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy4 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy5 = new TH1F("dummy","",100, 0, 10);
  TH1F* dummy6 = new TH1F("dummy","",100, 0, 10);
  
  dummy1->SetLineWidth(LEGWIDTH);
  dummy2->SetLineWidth(LEGWIDTH);
  dummy3->SetLineWidth(LEGWIDTH);
  dummy4->SetLineWidth(LEGWIDTH);
  dummy5->SetLineWidth(LEGWIDTH);
  dummy6->SetLineWidth(LEGWIDTH);

  dummy1->SetMarkerSize(MSIZE);
  dummy2->SetMarkerSize(MSIZE);
  dummy3->SetMarkerSize(MSIZE);
  dummy4->SetMarkerSize(MSIZE);
  dummy5->SetMarkerSize(MSIZE);
  dummy6->SetMarkerSize(MSIZE);



  dummy1->SetLineColor(COLOR_A);
  leg1->AddEntry(dummy1, "500 devices/km^{2}", "l");
  dummy2->SetLineColor(COLOR_B);
  leg1->AddEntry(dummy2, "1000 devices/km^{2}", "l");
  dummy3->SetLineColor(COLOR_C);
  leg1->AddEntry(dummy3, "5000 devices/km^{2}", "l");

  dummy4->SetMarkerStyle(STYLE_A);
  dummy4->SetLineColor(1);
  dummy5->SetMarkerStyle(STYLE_B);
  dummy5->SetLineColor(1);
  dummy6->SetMarkerStyle(STYLE_C);
  dummy6->SetLineColor(1);

  leg1->AddEntry(dummy6, "#epsilon A = 1e-9 (#gamma), 1e-5 (#mu)", "p");
  leg1->AddEntry(dummy5, "#epsilon A = 5e-9 (#gamma), 5e-5 (#mu)", "p");
  //leg1->AddEntry(dummy4, "#epsilon A = 1e-7");    
  TCanvas *ct = new TCanvas("ct");
  TCanvas *cp = new TCanvas("cp");
  TCanvas *ce = new TCanvas("ce");
  //TCanvas *ch = new TCanvas("ch");
  
  ce->cd();
  ce->SetLogx();
  ce->SetBottomMargin(0.12);
  TH2F * hframe = new TH2F("hframe", "", 10, 2E18, 1E21, 10, 0.0, 0.6);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]"); 
  hframe->GetXaxis()->SetTitleSize(0.06); 
  hframe->GetXaxis()->SetTitleOffset(0.95);
  hframe->GetYaxis()->SetTitle("Fractional Energy Resolution");
  hframe->GetYaxis()->SetTitleOffset(0.75);  
  hframe->GetYaxis()->SetTitleSize(0.06);


  hframe->Draw();
  geres_1000_8 ->Draw("LP");
  geres_1000_9 ->Draw("LP");
  geres_5000_8 ->Draw("LP");
  geres_5000_9 ->Draw("LP");
  geres_10000_8->Draw("LP");
  geres_10000_9->Draw("LP");
  leg1->Draw();
  ce->SaveAs("energy_res.png");
  ce->SaveAs("energy_res.eps");



  ct->cd();
  ct->SetLogx();
  ct->SetBottomMargin(0.12);
  TH2F * hframe = new TH2F("hframe", "", 10, 2E18, 1E21, 10, 0.0, 0.275);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]");
  hframe->GetXaxis()->SetTitleSize(0.06); 
  hframe->GetXaxis()->SetTitleOffset(0.95);
  hframe->GetYaxis()->SetTitle("Theta Resolution [Radians]");
  hframe->GetYaxis()->SetTitleOffset(0.75);  
  hframe->GetYaxis()->SetTitleSize(0.06);

  hframe->Draw();
  gtres_1000_8 ->Draw("LP");
  gtres_1000_9 ->Draw("LP");
  gtres_5000_8 ->Draw("LP");
  gtres_5000_9 ->Draw("LP");
  gtres_10000_8->Draw("LP");
  gtres_10000_9->Draw("LP");
  leg1->Draw();
  ct->SaveAs("theta_res.png");
  ct->SaveAs("theta_res.eps");

  cp->cd();
  cp->SetLogx();
  cp->SetBottomMargin(0.12);
  TH2F * hframe = new TH2F("hframe", "", 10, 2E18, 1E21, 10, 0.0, 0.8);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]");
  hframe->GetXaxis()->SetTitleSize(0.06); 
  hframe->GetXaxis()->SetTitleOffset(0.95);
  hframe->GetYaxis()->SetTitle("Phi Resolution [Radians]");
  hframe->GetYaxis()->SetTitleOffset(0.75);  
  hframe->GetYaxis()->SetTitleSize(0.06);

  hframe->Draw();
  gpres_1000_8 ->Draw("LP");
  gpres_1000_9 ->Draw("LP");
  gpres_5000_8 ->Draw("LP");
  gpres_5000_9 ->Draw("LP");
  gpres_10000_8->Draw("LP");
  gpres_10000_9->Draw("LP");
  leg1->Draw();
  cp->SaveAs("phi_res.png");
  cp->SaveAs("phi_res.eps");
}
