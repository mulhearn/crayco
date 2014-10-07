

mkplots(){
   gStyle->SetOptStat(1);

   TFile * f = new TFile("plots.root");
   f->ls();

   TH1F * hpn = (TH1F *) f->Get("hpn");
   TH1F * hpt = (TH1F *) f->Get("hpt");
   TH1F * hpp = (TH1F *) f->Get("hpp");
   TH1F * hde = (TH1F *) f->Get("hde");
   TH1F * hdp = (TH1F *) f->Get("hdp");
   TH1F * hdt = (TH1F *) f->Get("hdt");

   TCanvas * c1 = new TCanvas;
   hpn->SetXTitle("pull shower log(N)");
   hpn->Draw("ep");
   c1->SaveAs("plots/pn.png");

   TCanvas * c2 = new TCanvas;
   hpt->SetXTitle("pull shower eta)");
   hpt->Draw("ep");
   c2->SaveAs("plots/pt.png");

   TCanvas * c3 = new TCanvas;
   hpp->SetXTitle("pull shower phi");
   hpp->Draw("ep");
   c3->SaveAs("plots/pp.png");

   TCanvas * c4 = new TCanvas;
   hde->SetXTitle("shower E fractional difference");
   hde->Draw("ep");
   c4->SaveAs("plots/dn.png");

   TCanvas * c5 = new TCanvas;
   hdt->SetXTitle("shower theta fractional difference");
   hdt->Draw("ep");
   c5->SaveAs("plots/dt.png");

   TCanvas * c6 = new TCanvas;
   hdp->SetXTitle("shower phi fractional difference");
   hdp->Draw("ep");
   c6->SaveAs("plots/dp.png");


}
