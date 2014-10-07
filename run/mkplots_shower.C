TH1F * do_plots(const char * finename, const char * hname, const char * pname);

void mkplots_shower(){
   hrata = do_plots("shower_7_1k.root", "hrata", "shower_7_1k.png"); 
   hratb = do_plots("shower_7_10k.root", "hratb", "shower_7_10k.png");  
   hratc = do_plots("shower_8_10k.root", "hratc", "shower_8_10k.png");    
   
   hrata->SetLineColor(1);
   hratb->SetLineColor(2);
   hratc->SetLineColor(3);

   TCanvas * c2 = new TCanvas;
   TH1F * hframe2 = new TH1F("hframe2", "", 100, 10.0, 12.0);
   hframe->SetMaximum(2.0);
   hframe->SetMinimum(0.0);
   hframe->SetXTitle("Fit Log(Energy)");
   hframe->GetXaxis()->SetTitleOffset(1.3);
   hframe->SetYTitle("Fraction from Lower Energy");
   hframe->GetYaxis()->SetTitleOffset(1.2);

   TLegend *leg = new TLegend(0.2,0.7,0.70, 0.85);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);

   leg->AddEntry(hrata, "#epsilon A =10^-7, N=1000", "l");
   leg->AddEntry(hratb, "#epsilon A =10^-7, N=10,000", "l");
   leg->AddEntry(hratc, "#epsilon A =10^-8, N=10,000", "l");

   hframe->Draw();

   hrata->Draw("epSAME");
   hratb->Draw("epSAME");
   hratc->Draw("epSAME");

   leg->Draw();
   hframe->Draw("SAME");
   c2->SaveAs("misreco.png");
   c2->SaveAs("misreco.eps");
   
}


TH1F * do_plots(const char * filename, const char * hname, const char * pname){
   gStyle->SetOptStat(0);

   TFile * f = new TFile(filename);
   f->ls();
   TH1F * hframe = new TH1F("hframe", "", 100, 10.0, 12.0);

   TH1F * htrue = (TH1F *) f->Get("hshower_true");
   TH1F * hreco = (TH1F *) f->Get("hshower_reco");
   TH1F * hnear = (TH1F *) f->Get("hshower_near");

   TH1F * hrat  = hnear->Clone(hname);
   hrat->Divide(hnear, hreco, 1.0, 1.0, "B");

   TCanvas * c1 = new TCanvas;
   c1->SetLogy();
   hframe->SetMaximum(1E1);
   hframe->SetMinimum(1E-7);

   hframe->SetXTitle("Fit Log(Energy)");
   hframe->GetXaxis()->SetTitleOffset(1.3);
   hframe->SetYTitle("Showers");
   hframe->GetYaxis()->SetTitleOffset(1.2);

   hframe->Draw();
   htrue->Draw("HSAME");
   hnear->SetLineColor(2);
   hnear->Draw("HSAME");
   hreco->SetLineColor(3);
   hreco->Draw("HSAME");
   c1->SaveAs(pname);

   return hrat;

}
