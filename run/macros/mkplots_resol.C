
//#define OPTA "short/resolution_res_short_500_opt.root"
//#define OPTB "short/resolution_res_short_1000_opt.root"
//#define OPTC "short/resolution_res_short_5000_opt.root"

//#define PSMA "short/resolution_res_short_500_psm.root"
//#define PSMB "short/resolution_res_short_1000_psm.root"
//#define PSMC "short/resolution_res_short_5000_psm.root"

//#define PSMA "short/resolution_res_short_500_mdm.root"
//#define PSMB "short/resolution_res_short_1000_mdm.root"
//#define PSMC "short/resolution_res_short_5000_mdm.root"

#define OPTA "root/full/resolution_full_500_opt.root"
#define OPTB "root/full/resolution_full_1000_opt.root"
#define OPTC "root/full/resolution_full_5000_opt.root"

//#define PSMA "root/full/resolution_full_500_psm.root"
//#define PSMB "root/full/resolution_full_1000_psm.root"
//#define PSMC "root/full/resolution_full_5000_psm.root"

#define PSMA "root/full/resolution_full_500_mdm.root"
#define PSMB "root/full/resolution_full_1000_mdm.root"
#define PSMC "root/full/resolution_full_5000_mdm.root"

void mkplots_var(const char * var);

void mkplots_resol(){
  mkplots_var("eres", "Fractional Energy Resolution", 0.52);
  mkplots_var("pres", "Phi Resolution [rad]]", 0.9);
  mkplots_var("tres", "Theta Resolution [rad]", 0.4);
}

void mkplots_var(const char * var, const char * ytitle, double ymax){
   enum {COLOR_A=2, COLOR_B=4, COLOR_C=6};
   enum {STYLE_A=2, STYLE_B=24, STYLE_C=20};
   enum {WIDTH=2};
   enum {LEGWIDTH=3};
   const double MSIZE=1.5;

   gStyle->SetOptStat(0);

  TFile *fopta  = new TFile(OPTA);
  TFile *foptb  = new TFile(OPTB);
  TFile *foptc  = new TFile(OPTC);
  TFile *fpsma  = new TFile(PSMA);
  TFile *fpsmb  = new TFile(PSMB);
  TFile *fpsmc  = new TFile(PSMC);

  TGraphErrors *gopta    = (TGraphErrors*)fopta->Get(var);
  TGraphErrors *goptb    = (TGraphErrors*)foptb->Get(var);
  TGraphErrors *goptc    = (TGraphErrors*)foptc->Get(var);

  TGraphErrors *gpsma    = (TGraphErrors*)fpsma->Get(var);
  TGraphErrors *gpsmb    = (TGraphErrors*)fpsmb->Get(var);
  TGraphErrors *gpsmc    = (TGraphErrors*)fpsmc->Get(var);

  gopta->SetLineColor(COLOR_A);
  gopta->SetMarkerColor(COLOR_A);
  gopta->SetLineWidth(WIDTH);
  gopta->SetMarkerSize(MSIZE);
  gopta->SetMarkerStyle(STYLE_B);
  
  goptb->SetLineColor(COLOR_B);
  goptb->SetMarkerColor(COLOR_B);
  goptb->SetLineWidth(WIDTH);
  goptb->SetMarkerSize(MSIZE);
  goptb->SetMarkerStyle(STYLE_B);

  goptc->SetLineColor(COLOR_C);
  goptc->SetMarkerColor(COLOR_C);
  goptc->SetLineWidth(WIDTH);
  goptc->SetMarkerSize(MSIZE);
  goptc->SetMarkerStyle(STYLE_B);



  gpsma->SetLineColor(COLOR_A);
  gpsma->SetMarkerColor(COLOR_A);
  gpsma->SetLineWidth(WIDTH);
  gpsma->SetMarkerSize(MSIZE);
  gpsma->SetMarkerStyle(STYLE_C);

  gpsmb->SetLineColor(COLOR_B);
  gpsmb->SetMarkerColor(COLOR_B);
  gpsmb->SetLineWidth(WIDTH);
  gpsmb->SetMarkerSize(MSIZE);
  gpsmb->SetMarkerStyle(STYLE_C);

  gpsmc->SetLineColor(COLOR_C);
  gpsmc->SetMarkerColor(COLOR_C);
  gpsmc->SetLineWidth(WIDTH);
  gpsmc->SetMarkerSize(MSIZE);
  gpsmc->SetMarkerStyle(STYLE_C);

  TLegend *leg1 = new TLegend(0.12,0.14,0.55,0.40);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);

  TH1F* dummy1 = new TH1F("dummy1","",100, 0, 10);
  TH1F* dummy2 = new TH1F("dummy2","",100, 0, 10);
  TH1F* dummy3 = new TH1F("dummy3","",100, 0, 10);
  TH1F* dummy4 = new TH1F("dummy4","",100, 0, 10);
  TH1F* dummy5 = new TH1F("dummy5","",100, 0, 10);
  
  
  dummy1->SetLineWidth(LEGWIDTH);
  dummy1->SetMarkerSize(MSIZE);
  dummy1->SetLineColor(COLOR_A);
  dummy2->SetLineWidth(LEGWIDTH);
  dummy2->SetMarkerSize(MSIZE);
  dummy2->SetLineColor(COLOR_B);
  dummy3->SetLineWidth(LEGWIDTH);
  dummy3->SetMarkerSize(MSIZE);
  dummy3->SetLineColor(COLOR_C);

  dummy4->SetLineWidth(LEGWIDTH);
  dummy4->SetMarkerSize(MSIZE);
  dummy4->SetMarkerStyle(STYLE_B);

  dummy5->SetLineWidth(LEGWIDTH);
  dummy5->SetMarkerSize(MSIZE);
  dummy5->SetMarkerStyle(STYLE_C);


  leg1->AddEntry(dummy1, "500 devices/km^{2}", "l");
  leg1->AddEntry(dummy2, "1000 devices/km^{2}", "l");
  leg1->AddEntry(dummy3, "5000 devices/km^{2}", "l");
  leg1->AddEntry(dummy4, "high efficiency", "p");
  leg1->AddEntry(dummy5, "low efficiency", "p");

  TCanvas *c = new TCanvas();
  c->cd();
  c->SetLogx();
  c->SetBottomMargin(0.12);
  TH2F * hframe = new TH2F("hframe", "", 10, 8E18, 2E21, 10, 0.0, ymax);
  hframe->GetXaxis()->SetTitle("Primary Energy [eV]"); 
  hframe->GetXaxis()->SetTitleSize(0.06); 
  hframe->GetXaxis()->SetTitleOffset(0.95);
  hframe->GetYaxis()->SetTitle(ytitle);
  hframe->GetYaxis()->SetTitleOffset(0.75);  
  hframe->GetYaxis()->SetTitleSize(0.06);


  hframe->Draw();
  gopta->Draw("LPE");
  goptb->Draw("LPE");
  goptc->Draw("LPE");
  gpsma->Draw("LPE");
  gpsmb->Draw("LPE");
  gpsmc->Draw("LPE");

  leg1->Draw();
  char name[100];
  sprintf(name,"%s.pdf", var);
  c->SaveAs(name);
}
