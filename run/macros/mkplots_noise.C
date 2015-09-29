
void nhits_plots();
void ereco_plots();
void plots(const char * prefix, const char * var, const char * title, double min, double max);

mkplots_noise(){
  //plots("h_04", "llr", "LLR ", 0.0, 1.0);
  //plots("h_04", "sigf", "LogE Significance / sqrt(N)", 0.0, 5.0);
  //nhits_plots();
  ereco_plots();
}


TGraph * add_hist(TFile * f, const char * tagn, const char * tagsa, const char * tagsb, int color){
   TH1F * hsa = f->Get(tagsa);
   TH1F * hsb = f->Get(tagsb);
   TH1F * hn = f->Get(tagn);

   int nsa = hsa->GetNbinsX();
   int nsb = hsb->GetNbinsX();
   int nn = hn->GetNbinsX();
   
   double * xsa = new double[2*nsa];
   double * ysa = new double[2*nsa];
   double * xsb = new double[2*nsb];
   double * ysb = new double[2*nsb];
   double * xn  = new double[nn];
   double * yn  = new double[nn];
   double * dyn = new double[nn];

   for (int i=0; i<nsa; i++){
     ysa[2*i] = hsa->GetBinContent(i);
     xsa[2*i] = pow(10.0, hsa->GetBinLowEdge(i)); 
     ysa[2*i+1] = hsa->GetBinContent(i);
     xsa[2*i+1] = pow(10.0, hsa->GetBinLowEdge(i)+ hsa->GetBinWidth(i)); 
   }
   for (int i=0; i<nsb; i++){
     ysb[2*i] = hsb->GetBinContent(i);
     xsb[2*i] = pow(10.0, hsb->GetBinLowEdge(i)); 
     ysb[2*i+1] = hsb->GetBinContent(i);
     xsb[2*i+1] = pow(10.0, hsb->GetBinLowEdge(i)+ hsb->GetBinWidth(i)); 
   }
   for (int i=0; i<nn; i++){
     xn[i] = pow(10.0, hn->GetBinLowEdge(i)+0.5*hn->GetBinWidth(i)); 
     yn[i] = hn->GetBinContent(i);
     dyn[i] = hn->GetBinError(i);
     cout << hn->GetBinContent(i) << " +/- " << hn->GetBinError(i) << "\n";

   }

   gsa = new TGraph(2*nsa, xsa, ysa);
   gsb = new TGraph(2*nsb, xsb, ysb);
   gn = new TGraphErrors(nn, xn, yn, 0, dyn);

   gsa->SetLineStyle(2);
   gsb->SetLineStyle(2);
   gn->SetLineStyle(1);
   gsa->SetLineColor(4);
   gsb->SetLineColor(2);
   gn->SetLineColor(1);
   gn->SetMarkerStyle(7);
   gn->SetMarkerColor(1);
   gsa->Draw("L");
   gsb->Draw("LSAME");
   gn->Draw("EP");
   return gn;
}

void ereco_plots(){
   gStyle->SetOptStat(0);

   TFile * f = new TFile("noise.root");
   f->ls();

   TLegend *leg = new TLegend(0.5,0.7,0.85, 0.90);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);

   
   TCanvas * c1 = new TCanvas;
   c1->SetLogy();
   c1->SetLogx();
   c1->SetBottomMargin(0.13);
   TH1F * hframe = new TH1F("hframe", "", 10, 1E18, 1E22);

   hframe->SetMaximum(1E1);  
   hframe->SetMinimum(1E-17);  
   hframe->SetXTitle("Reconstructed Primary Energy [eV]");
   hframe->SetYTitle("Rate per bin [Hz / km^2]");
   //hframe->SetYTitle("Rate");
   hframe->GetXaxis()->SetTitleSize(0.06); 
   hframe->GetYaxis()->SetTitleSize(0.06);

   hframe->GetXaxis()->SetTitleOffset(0.90);
   hframe->GetYaxis()->SetTitleOffset(0.75);  

   hframe->Draw();

   add_hist(f, "h_noise_ereco", "h_signala_ereco", "h_signalb_ereco", 4);


   hframe->Draw("SAME");

   TH1F * hdummy1 = new TH1F("hdummy1", "", 100, 0.0, 1.0);
   TH1F * hdummy2 = new TH1F("hdummy2", "", 100, 0.0, 1.0);
   TH1F * hdummy3 = new TH1F("hdummy3", "", 100, 0.0, 1.0);
   hdummy1->SetLineStyle(1);
   hdummy1->SetMarkerStyle(7);
   hdummy2->SetLineStyle(2);
   hdummy2->SetLineColor(4);
   hdummy3->SetLineStyle(2);
   hdummy3->SetLineColor(2);

   leg->AddEntry(hdummy1,  "noise only", "p");
   leg->AddEntry(hdummy2,  "noise and low-E simulated shower", "l");
   leg->AddEntry(hdummy3,  "noise and high-E simulated shower", "l");
   leg->Draw();

   c1->SaveAs("plots/noise.pdf");
   c1->SaveAs("plots/noise.eps");
}

void nhits_plots(){
  gStyle->SetOptStat(0);
  
  TFile * f = new TFile("noise.root");
  f->ls();
  
  TH1F * h = new TH1F("hframe", "", 100, 0., 100.0);
  h->SetMinimum(1E-11);
  h->SetMaximum(h_02_noise_nhits_unwgt->GetMaximum()*1.5);

  TCanvas * c = new TCanvas();
  c.SetLogy();
  
  h->SetXTitle("Number of Hits");
  h->Draw();

  TGraph * ga = f->Get("g_02_noise_nhits");
  TGraph * gb = f->Get("g_2_noise_nhits");
  ga->SetLineColor(2);
  gb->SetLineColor(4);
  ga->Draw("SAME");
  gb->Draw("SAME");

  TH1F * ha = f->Get("h_02_noise_nhits_unwgt");
  TH1F * hb = f->Get("h_2_noise_nhits_unwgt");
  TH1F * hc = f->Get("h_02_noise_nhits_wgt");
  TH1F * hd = f->Get("h_2_noise_nhits_wgt");
  ha->SetMarkerStyle(2);
  hb->SetMarkerStyle(2);
  hc->SetMarkerStyle(7);
  hd->SetMarkerStyle(7);
  ha->Draw("EPSAME");  
  hb->Draw("EPSAME");
  hc->Draw("EPSAME");  
  hd->Draw("EPSAME");
  c->SaveAs("plots/nhits.pdf");
  c->SaveAs("plots/nhits.eps");
}


void plots(const char * prefix, const char * var, const char * title, double min, double max){
  gStyle->SetOptStat(0);
  
  TFile * f = new TFile("noise.root");
  f->ls();
  
  TH1F * h = new TH1F("hframe", "", 100, min, max);
  h->SetMinimum(1E-15);
  h->SetMaximum(10000.0);

  TCanvas * c = new TCanvas();
  c->SetLogy();
  
  h->SetXTitle(title);
  h->Draw();

  const char name[100];
  sprintf(name, "%s_noise_%s", prefix, var);
  TH1F * ha = f->Get(name);
  sprintf(name, "%s_signal_%s", prefix, var);
  TH1F * hb = f->Get(name);  
  sprintf(name, "%s_signalb_%s", prefix, var);
  TH1F * hc = f->Get(name);  

  ha->SetLineColor(2);
  hb->SetLineColor(4);
  hc->SetLineColor(3);
  ha->Draw("EPSAME");  
  hb->Draw("EPSAME");
  cb->Draw("EPSAME");

  sprintf(name, "plots/%s.pdf", var);
  c->SaveAs(name);
  sprintf(name, "plots/%s.pdf", var);
  c->SaveAs(name);
}
