

mkplots_ereco(){
   gStyle->SetOptStat(0);

   TFile * f = new TFile("edist.root");
   f->ls();

   TH1F * hereco = (TH1F *) f->Get("hereco");
   TCanvas * c1 = new TCanvas;
   c1->SetLogy();
   hereco->SetXTitle("Fit Log(Energy)");
   hereco->GetXaxis()->SetTitleOffset(1.3);
   hereco->SetYTitle("Events");
   hereco->GetYaxis()->SetTitleOffset(1.2);
   hereco->Draw("ep");
   c1->SaveAs("hereco.png");

}
