

TH1F * add_hist(TFile * f, const char * tagn, const char * tags, int color){
   TH1F * hs = f->Get(tags);
   TH1F * hn = f->Get(tagn);
   hs->SetLineStyle(2);
   hn->SetLineStyle(1);
   hs->SetLineColor(color);
   hn->SetLineColor(color);
   hs->Draw("HSAME");
   hn->Draw("HSAME");
   return hn;
}


void mkplots_noise(){
   gStyle->SetOptStat(0);

   TFile * f = new TFile("noise.root");
   f->ls();

   TLegend *leg = new TLegend(0.2,0.6,0.70, 0.85);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);

   
   TCanvas * c1 = new TCanvas;
   c1->SetLogy();
   c1->SetBottomMargin(0.13);
   TH1F * hframe = new TH1F("hframe", "", 10, 9.0, 13.0);

   hframe->SetMaximum(1000.0);  
   hframe->SetMinimum(0.1);  
   hframe->SetXTitle("Fit Log(Energy)");
   hframe->SetYTitle("Events");
   hframe->GetXaxis()->SetTitleSize(0.06); 
   hframe->GetYaxis()->SetTitleSize(0.06);

   hframe->GetXaxis()->SetTitleOffset(0.90);
   hframe->GetYaxis()->SetTitleOffset(0.75);  


   hframe->Draw();
   leg->AddEntry(add_hist(f, "h_01_noise", "h_01_signal", 1), "rate=0.01", "l");
   leg->AddEntry(add_hist(f, "h_05_noise", "h_05_signal", 2), "rate=0.05", "l");
   leg->AddEntry(add_hist(f, "h_1_noise",  "h_1_signal", 3),  "rate=0.10", "l");
   hframe->Draw("SAME");

   TH1F * hdummy1 = new TH1F("hdummy1", "", 100, 0.0, 1.0);
   TH1F * hdummy2 = new TH1F("hdummy2", "", 100, 0.0, 1.0);
   hdummy1->SetLineStyle(1);
   hdummy2->SetLineStyle(2);

   leg->AddEntry(hdummy1,  "no shower", "l");
   leg->AddEntry(hdummy2,  "with shower", "l");
   leg->Draw();

   c1->SaveAs("noise.png");
   c1->SaveAs("noise.eps");

}
