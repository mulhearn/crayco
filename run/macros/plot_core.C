


{  
  gStyle->SetOptStat(0);
  TCanvas c;
  TH1F h1("h1","",100,18.5,20.0);
  TH1F h2("h2","",100,18.5,20.0);
  TH1F h3("h3","",100,18.5,20.0);
  TH1F h4("h4","",100,18.5,20.0);
  TH1F h5("h5","",100,18.5,20.0);
  TH1F h6("h6","",100,18.5,20.0);
  h1->SetLineColor(1);
  h2->SetLineColor(2);
  h3->SetLineColor(3);
  h4->SetLineColor(4);
  h5->SetLineColor(5);
  h6->SetLineColor(6);

  TLegend *leg = new TLegend(0.5,0.7,0.85, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);


  //leg->AddEntry(&h1, "E near 1st peak", "l");
  //leg->AddEntry(&h2, "E near 2nd peak", "l");

  noise->Draw("ereco>>h1","wgt*((eshower==0)&&(ncore==0))");
  noise->Draw("ereco>>h2","wgt*((eshower==0)&&(ncore==1))");
  noise->Draw("ereco>>h3","wgt*((eshower==0)&&(ncore==2))");
  noise->Draw("ereco>>h4","wgt*((eshower==0)&&(ncore==3))");
  noise->Draw("ereco>>h5","wgt*((eshower==0)&&(ncore==4))");
  noise->Draw("ereco>>h6","wgt*((eshower==0)&&(ncore==5))");
  
  double s1 = h1->GetSumOfWeights();
  double s2 = h2->GetSumOfWeights();
  double s3 = h3->GetSumOfWeights();
  double s4 = h4->GetSumOfWeights();
  double s5 = h5->GetSumOfWeights();
  double s6 = h6->GetSumOfWeights();

  
  h1->Scale(100.0/s1);
  h2->Scale(100.0/s2);
  h3->Scale(100.0/s3);
  h4->Scale(100.0/s4);
  h5->Scale(100.0/s5);
  h6->Scale(100.0/s6);

  h1->SetMaximum(60.0);

  h1->SetXTitle("num hits in core");
  h1->Draw();
  h2->Draw("SAME");
  h3->Draw("SAME");
  //h4->Draw("SAME");
  //h5->Draw("SAME");
  //h6->Draw("SAME");
  //leg->Draw();
  c->SaveAs("binomial.pdf");
}
