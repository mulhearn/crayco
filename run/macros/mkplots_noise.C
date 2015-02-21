

TGraph * add_hist(TFile * f, const char * tagn, const char * tags, int color, double scale_s, double scale_n){
   TH1F * hs = f->Get(tags);
   TH1F * hn = f->Get(tagn);

   int ns = hs->GetNbinsX();
   int nn = hn->GetNbinsX();
   
   double * xs = new double[2*ns];
   double * ys = new double[2*ns];
   double * xn = new double[2*nn];
   double * yn = new double[2*nn];

   for (int i=0; i<ns; i++){
     ys[2*i] = scale_s * hs->GetBinContent(i);
     xs[2*i] = pow(10.0, 9.0 + hs->GetBinLowEdge(i)); 
     ys[2*i+1] = scale_s * hs->GetBinContent(i);
     xs[2*i+1] = pow(10.0, 9.0 + hs->GetBinLowEdge(i)+ hs->GetBinWidth(i)); 
   }
   for (int i=0; i<nn; i++){
     yn[2*i] = scale_n * hn->GetBinContent(i);
     xn[2*i] = pow(10.0, 9.0 + hn->GetBinLowEdge(i)); 
     yn[2*i+1] = scale_n * hn->GetBinContent(i);
     xn[2*i+1] = pow(10.0, 9.0 + hn->GetBinLowEdge(i)+ hs->GetBinWidth(i)); 
   }

   gs = new TGraph(2*ns, xs, ys);
   gn = new TGraph(2*nn, xn, yn);

   gs->SetLineStyle(2);
   gn->SetLineStyle(1);
   gs->SetLineColor(color);
   gn->SetLineColor(color);
   gs->Draw("L");
   gn->Draw("L");
   return gn;
}


void mkplots_noise(){
   gStyle->SetOptStat(0);

   TFile * f = new TFile("noise.root");
   f->ls();

   TLegend *leg = new TLegend(0.1,0.6,0.60, 0.85);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);

   
   TCanvas * c1 = new TCanvas;
   c1->SetLogy();
   c1->SetLogx();
   c1->SetBottomMargin(0.13);
   TH1F * hframe = new TH1F("hframe", "", 10, 1E18, 1E22);

   hframe->SetMaximum(10.0);  
   hframe->SetMinimum(0.0000000001);  
   hframe->SetXTitle("Reconstructed Energy [eV]");
   hframe->SetYTitle("Rate [Hz / km^2]");
   hframe->GetXaxis()->SetTitleSize(0.06); 
   hframe->GetYaxis()->SetTitleSize(0.06);

   hframe->GetXaxis()->SetTitleOffset(0.90);
   hframe->GetYaxis()->SetTitleOffset(0.75);  

   hframe->Draw();
   double peryrHz = 3.17E-8;
   double N_n = 1000;
   double N_s = 100;
   double scale_n = 0.5 / N_n;
   double scale_s = peryrHz / N_s; // 10^19 eV
   //double scale_s = 1E-3 *  peryrHz / N; // 10^20 eV

   leg->AddEntry(add_hist(f, "h_02_noise", "h_02_signal", 1, scale_s, scale_n), "rate=0.02", "l");
   leg->AddEntry(add_hist(f, "h_2_noise", "h_2_signal", 2, scale_s, scale_n), "rate=0.2", "l");
   hframe->Draw("SAME");

   TH1F * hdummy1 = new TH1F("hdummy1", "", 100, 0.0, 1.0);
   TH1F * hdummy2 = new TH1F("hdummy2", "", 100, 0.0, 1.0);
   hdummy1->SetLineStyle(1);
   hdummy2->SetLineStyle(2);

   leg->AddEntry(hdummy1,  "noise only", "l");
   leg->AddEntry(hdummy2,  "noise and simulated shower", "l");
   leg->Draw();

   c1->SaveAs("plots/noise.png");
   c1->SaveAs("plots/noise.eps");
}
