#include <iostream>
#include <fstream>
#include <set>
#include<TFile.h>
#include<TTree.h>
#include<TH1.h>
#include<TH2.h>
#include<THStack.h>
#include<TGraphErrors.h>
#include<TGraphAsymmErrors.h>
#include<TCanvas.h>
#include<TFrame.h>
#include<TLegend.h>
#include<vector>
#include<iostream>
#include<TMath.h>
#include<TROOT.h>
#include<TInterpreter.h>
#include<TStyle.h>
#include<TChain.h>
#include<TString.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TLatex.h>
#include "TEfficiency.h"
#include "TAxis.h"
//#include "TOTTOTerror.h"
//#include "JERerror.h"
//#include "JESerror.h"
#include "staterror.h"
#include "TRandom3.h"
#include "TBuffer.h"
#include "TRandom2.h"
#include "TUUID.h"


void bdt_com_v1_stat_2016_mergebin(){
  Bool_t syst_on = false;
  gStyle->SetLineWidth(2);
  TCanvas *c1 =new TCanvas("c1", " ", 0, 0,700,800);

  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
//  c1->SetGridx();
//  c1->SetGridy();
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->Draw();

//  c1->SetLogy();

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
//  pad1->SetLogy();
//  pad1->SetGridx();
//  pad1->SetGridy();
  pad1->SetTickx(1);
  pad1->SetTicky(1);

//  pad1->cd();



  THStack *hs = new THStack("hs","");
  THStack *hs_mod = new THStack("hs_mod","");

  Double_t error;
  ////////////////////Data rootfiles/////////////////////


  TFile *fdata1 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/Data/BDT_UL18_MC_hists.root");
  TH1F  *h1 = (TH1F*)fdata1->Get("BDT_UL18_MC_test_bkg_mvaHist");
  h1->Rebin(2);
  h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX())+h1->GetBinContent(h1->GetNbinsX()+1));
  h1->SetBinContent(1, h1->GetBinContent(1)+h1->GetBinContent(0));


  double xmax = h1->GetXaxis()->GetXmax();
  double xmin = h1->GetXaxis()->GetXmin();
  double xmin1 = 0.0;
  double xmax1 = 1.0;

  //////////////////////MC files//////////////////////


  TFile *f1 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/TTGJets/BDT_UL18_MC_hists.root");
  TH1F *hmc1 = (TH1F*)f1->Get("BDT_UL18_MC_test_bkg_mvaHist");
  hmc1->Rebin(2);
  hmc1->SetBinContent(hmc1->GetNbinsX(), hmc1->GetBinContent(hmc1->GetNbinsX())+hmc1->GetBinContent(hmc1->GetNbinsX()+1));
  hmc1->SetBinContent(1, hmc1->GetBinContent(1)+hmc1->GetBinContent(0));

  TFile *f2 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/TTSL/BDT_UL18_MC_hists.root");
  TH1F *hmc2 = (TH1F*)f2->Get("BDT_UL18_MC_test_bkg_mvaHist");
  hmc2->Rebin(2);
  hmc2->SetBinContent(hmc2->GetNbinsX(), hmc2->GetBinContent(hmc2->GetNbinsX())+hmc2->GetBinContent(hmc2->GetNbinsX()+1));
  hmc2->SetBinContent(1, hmc2->GetBinContent(1)+hmc2->GetBinContent(0));

  TFile *f3 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/TT2L2Nu/BDT_UL18_MC_hists.root");
  TH1F *hmc3 = (TH1F*)f3->Get("BDT_UL18_MC_test_bkg_mvaHist");
  hmc3->Rebin(2);
  hmc3->SetBinContent(hmc3->GetNbinsX(), hmc3->GetBinContent(hmc3->GetNbinsX())+hmc3->GetBinContent(hmc3->GetNbinsX()+1));
  hmc3->SetBinContent(1, hmc3->GetBinContent(1)+hmc3->GetBinContent(0));

  TFile *f4 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/DY/BDT_UL18_MC_hists.root");
  TH1F *hmc4 = (TH1F*)f4->Get("BDT_UL18_MC_test_bkg_mvaHist");
  hmc4->Rebin(2);
  hmc4->SetBinContent(hmc4->GetNbinsX(), hmc4->GetBinContent(hmc4->GetNbinsX())+hmc4->GetBinContent(hmc4->GetNbinsX()+1));
  hmc4->SetBinContent(1, hmc4->GetBinContent(1)+hmc4->GetBinContent(0));

  TFile *f5 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/WH20/BDT_UL18_MC_hists.root");
  TH1F *hsig1 = (TH1F*)f5->Get("BDT_UL18_MC_test_sig_mvaHist");
  hsig1->Rebin(2);

  TFile *f6 = TFile::Open("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2016_merg/WH55/BDT_UL18_MC_hists.root");
  TH1F *hsig2 = (TH1F*)f6->Get("BDT_UL18_MC_test_sig_mvaHist");
  hsig2->Rebin(2);


  // === Merge bins in [0.9, 1.0] into a single bin for data, background and signal ===
  auto mergeHighBDT = [](TH1F* h) -> TH1F* {
    int nOrig = h->GetNbinsX();
    // Find the lower bin edge nearest to 0.9 (robust against floating-point binning)
    int binMergeStart = 1;
    double minDist = TMath::Abs(h->GetXaxis()->GetBinLowEdge(1) - 0.9);
    for (int i = 2; i <= nOrig; i++) {
      double dist = TMath::Abs(h->GetXaxis()->GetBinLowEdge(i) - 0.9);
      if (dist < minDist) { minDist = dist; binMergeStart = i; }
    }
    cout << "mergeHighBDT: " << h->GetName()
         << "  nearest edge to 0.9 → bin " << binMergeStart
         << "  (edge = " << h->GetXaxis()->GetBinLowEdge(binMergeStart) << ")" << endl;
    // New bin count: (binMergeStart-1) unchanged bins + 1 merged bin
    int nNew = binMergeStart;
    std::vector<double> edges(nNew + 1);
    for (int i = 1; i <= binMergeStart; i++)
      edges[i - 1] = h->GetXaxis()->GetBinLowEdge(i);
    edges[nNew] = h->GetXaxis()->GetXmax();   // right edge = 1.0
    return (TH1F*)h->Rebin(nNew, Form("%s_rb", h->GetName()), edges.data());
  };

  h1   = mergeHighBDT(h1);    // data
  hmc1 = mergeHighBDT(hmc1);  // TTGJets
  hmc2 = mergeHighBDT(hmc2);  // TTSL
  hmc3 = mergeHighBDT(hmc3);  // TT2L2Nu
  hmc4 = mergeHighBDT(hmc4);  // DY
  hsig1 = mergeHighBDT(hsig1); // Signal WH20
  hsig2 = mergeHighBDT(hsig2); // Signal WH55
  // =========================================================================

  TH1F *totalMC = (TH1F*)hmc1->Clone();
  totalMC->Add(hmc2);
  totalMC->Add(hmc3);
  totalMC->Add(hmc4);


  totalMC->SetFillColor(14);
  gStyle->SetHatchesSpacing(0.7);
  gStyle->SetHatchesLineWidth(1.);
  totalMC->SetFillStyle(3354);

  hmc1->SetFillColor(TColor::GetColor("#5790fc"));
  hmc1->SetFillStyle(1001);
  hmc1->SetLineColor(TColor::GetColor("#5790fc"));

  hmc2->SetFillColor(TColor::GetColor("#f89c20"));
  hmc2->SetFillStyle(1001);
  hmc2->SetLineColor(TColor::GetColor("#f89c20"));

  hmc3->SetFillColor(TColor::GetColor("#e42536"));
  hmc3->SetFillStyle(1001);
  hmc3->SetLineColor(TColor::GetColor("#e42536"));

  hmc4->SetFillColor(TColor::GetColor("#964a8b"));
  hmc4->SetFillStyle(1001);
  hmc4->SetLineColor(TColor::GetColor("#964a8b"));



  TH1F *totalMC_mod = (TH1F*)totalMC->Clone();

  TH1F *h_data = (TH1F*)h1->Clone();
  TH1F *hmc1_mod = (TH1F*)hmc1->Clone();
  TH1F *hmc2_mod = (TH1F*)hmc2->Clone();
  TH1F *hmc3_mod = (TH1F*)hmc3->Clone();
  TH1F *hmc4_mod =  (TH1F*)hmc4->Clone();




  h1->SetBinErrorOption(TH1::kPoisson);



int n = h1->GetNbinsX();
  float q = (1.-0.6827)/2.;
  float N, dN, yscale;
  Double_t x[n], y[n], EYlow[n], EYhigh[n], EXlow[n], EXhigh[n];
  for(int i = 0; i < n; i++){
    N = h1->GetBinContent(i+1);
    dN = h1->GetBinError(i+1);
    if(N > 0. && dN > 0. && abs((dN*dN)/N-1) > 0.0001){
      yscale = ((dN*dN)/N);
      N = (N/dN)*(N/dN);
    }
    else yscale = 1;
    x[i] = h1->GetXaxis()->GetBinCenter(i+1);
    y[i] = yscale*N;
    if(N > 0) EYlow[i] = yscale*(N - ROOT::Math::chisquared_quantile_c(1-q,2*N)/2.);
    else EYlow[i] = 0.;
    EYhigh[i] = yscale*(ROOT::Math::chisquared_quantile_c(q,2*(N+1))/2.- N);
    EXlow[i] = h1->GetXaxis()->GetBinUpEdge(i+1) - x[i];
    EXhigh[i] = x[i] - h1->GetXaxis()->GetBinLowEdge(i+1);
    EXlow[i] = 0.0;
    EXhigh[i] = 0.0;
  }
  TGraphAsymmErrors *hdata1 = new TGraphAsymmErrors(n,x,y,EXlow,EXhigh,EYlow,EYhigh);



  hs->Add(hmc4,"hist");
  hs->Add(hmc1,"hist");
  hs->Add(hmc3,"hist");
  hs->Add(hmc2,"hist");

  hs_mod->Add(hmc4_mod,"hist");
  hs_mod->Add(hmc1_mod,"hist");
  hs_mod->Add(hmc3_mod,"hist");
  hs_mod->Add(hmc2_mod,"hist");

  double a = h1->GetMaximum();
  hs->Draw("histo");


    TBox *shadedBox = new TBox(0.0, 0.0, 0.78, 200.0);
  shadedBox->SetFillColorAlpha(kGray, 0.60);  // gray, 40% opacity
  shadedBox->SetFillStyle(1001);
  shadedBox->SetLineWidth(0);                  // no border
//  shadedBox->Draw("same");




TH1F *totalMC_band = (TH1F*)totalMC->Clone("totalMC_band");
  totalMC_band->SetFillColor(TColor::GetColor("#9c9ca1"));
  totalMC_band->SetFillStyle(3001);
  totalMC_band->SetLineColor(TColor::GetColor("#9c9ca1"));
  totalMC_band->SetLineWidth(0);
  totalMC_band->SetMarkerSize(0);
  totalMC_band->Draw("E2 same");

  hs->Draw("histo same");



// Dashed line at BDT = 0.8
TLine *vline1 = new TLine(0.78, 0.0, 0.78, 200.0);
vline1->SetLineColor(kBlack);
vline1->SetLineWidth(2);
vline1->SetLineStyle(2);
//vline1->Draw("same");

// Dashed line at BDT = 0.9
TLine *vline2 = new TLine(0.9, 0.0, 0.9, 200.0);
vline2->SetLineColor(kBlack);
vline2->SetLineWidth(2);
vline2->SetLineStyle(2);
//vline2->Draw("same");



  hdata1->SetLineColor(1);
hdata1->SetMarkerColor(1);
hdata1->SetMarkerSize(1.0);
hdata1->SetMarkerStyle(20);
hdata1->Draw("pz same");


  hsig1->SetLineColor(TColor::GetColor("#832db6"));
  hsig1->SetLineWidth(3);
  hsig2->SetLineColor(TColor::GetColor("#92dadd"));
  hsig2->SetLineWidth(3);
  hsig1->Draw("hist same");
  hsig2->Draw("hist same");

  hs->GetYaxis()->SetLabelSize(0.045);

  hs->GetYaxis()->SetTitle("Events / 0.04");
  hs->GetYaxis()->SetTitleOffset(1.0);
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->SetMinimum(0.);
  hs->SetMaximum(200.);
  hs->GetXaxis()->SetLabelSize(0);
  hs->GetXaxis()->SetRangeUser(xmin1,xmax1);


  TLegend *legend1 = new TLegend(0.5, 0.5, 0.7, 0.85);
  legend1->SetTextFont(42);
  legend1->SetLineColor(0);
  legend1->SetTextSize(0.04);
  legend1->SetFillColor(0);
  legend1->AddEntry(hdata1, "Data", "lep");
  legend1->AddEntry(hmc4, "Drell-Yan", "f");
  legend1->AddEntry(hmc1, "t#bar{t}+#gamma+jets", "f");
  legend1->AddEntry(hmc3, "t#bar{t} (dilepton)", "f");
  legend1->AddEntry(hmc2, "t#bar{t} (semilepton)", "f");
  legend1->AddEntry(hsig1, "M_{A} 20 x 100", "l");
  legend1->AddEntry(hsig2, "M_{A} 55 x 100", "l");
  legend1->Draw();

  // "CMS" bold — inside plot, upper left
  TLatex *t_cms = new TLatex(0.14, 0.83, "#bf{CMS}");
  t_cms->SetNDC();
  t_cms->SetTextFont(42);
  t_cms->SetTextSize(0.058);
  t_cms->SetTextAlign(11);
  t_cms->Draw("same");

  // "Preliminary" italic — just to the right of CMS, same height
  TLatex *t_prelim = new TLatex(0.26, 0.83, "#it{Preliminary}");
  t_prelim->SetNDC();
  t_prelim->SetTextFont(42);
  t_prelim->SetTextSize(0.046);
  t_prelim->SetTextAlign(11);
  t_prelim->Draw("same");

  // Lumi + energy — upper right, above the frame
  TLatex *t_lumi = new TLatex(0.9, 0.92, "35.9 fb^{-1} (13 TeV)");
  t_lumi->SetNDC();
  t_lumi->SetTextFont(42);
  t_lumi->SetTextSize(0.046);
  t_lumi->SetTextAlign(31);
  t_lumi->Draw("same");

pad1->RedrawAxis();


  c1->cd();


  
//  TH1F *hist_data = (TH1F*)hdata1->Clone();
  Double_t ratioy[n], ratioEYlow[n], ratioEYhigh[n], ratioEXlow[n], ratioEXhigh[n];

  for(int i = 0; i < n; i++){
      double stackcontent = ((TH1F*)(hs->GetStack()->Last()))->GetBinContent(i+1);
      double stackerror = ((TH1F*)(hs->GetStack()->Last()))->GetBinError(i+1);

      double datacontent = y[i];
      double dataerrorYup = EYhigh[i];
      double dataerrorYdn = EYlow[i];

    if( (stackcontent!=0) && (datacontent !=0) ) {
      ratioy[i] = (datacontent / stackcontent) ;
      ratioEYhigh[i] = ratioy[i]*sqrt(pow((dataerrorYup/datacontent),2) + pow((stackerror/stackcontent),2));
      ratioEYlow[i] = ratioy[i]*sqrt(pow((dataerrorYdn/datacontent),2) + pow((stackerror/stackcontent),2));

      ratioEXlow[i] = 0.0;
      ratioEXhigh[i] = 0.0;
    }
  }

  TGraphAsymmErrors *hist_data = new TGraphAsymmErrors(n,x,ratioy,EXlow,EXhigh,ratioEYlow,ratioEYhigh);


  // Ratio Plot


  TPad *pad2 = new TPad("pad2", "newpad",0,0,1,0.3);
  pad2->Draw();
  pad2->cd();
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.25);
  pad2->SetRightMargin(0.9);
  pad2->SetFillStyle(0);
  pad2->SetTickx(1);
  pad2->SetTicky(1);

//  pad2->cd();


  hist_data->SetLineWidth(1.);
  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(1.0);
  hist_data->SetFillStyle(0);

  TH1F *band = (TH1F*)totalMC->Clone();
  band->Divide(totalMC);
  band->SetFillColor(TColor::GetColor("#9c9ca1"));
  band->SetLineColor(TColor::GetColor("#9c9ca1"));
  band->SetFillStyle(1001);

//cout<< totalMC->GetNbinsX() << "  " << h1->GetBinContent(1) << endl;


  TH1F *h_JetEnUp; TH1F *h_JetEnDown; TH1F *h_JetResUp; TH1F *h_JetResDown; TH1F *h_UnclusteredEnUp;TH1F *h_UnclusteredEnDown;
  TGraphAsymmErrors *tottoterree; TGraphAsymmErrors *jererree; TGraphAsymmErrors *erree;

  TGraphAsymmErrors *staterree = (TGraphAsymmErrors*)staterror(h1,totalMC);
  staterree->SetFillColor(kRed-7);
  staterree->SetLineColor(kRed-7);
  staterree->SetFillStyle(1001);
  band->GetXaxis()->SetTitle("BDT score");
  band->GetXaxis()->SetLabelSize(0.10);
  band->GetXaxis()->SetTitleSize(0.12);
  band->GetXaxis()->SetTitleOffset(0.9);
  band->GetXaxis()->SetRangeUser(xmin1,xmax1);
  band->GetXaxis()->SetTickLength(0.075);
  band->GetYaxis()->SetTitle("Data / MC");
  band->GetYaxis()->SetLabelSize(0.1);
  band->GetYaxis()->SetTitleSize(0.12);
  band->GetYaxis()->SetNdivisions(505);
  band->GetYaxis()->SetTitleOffset(0.4);
  band->GetYaxis()->SetRangeUser(0.0,2.0);
  band->SetTitle("");
  band->SetStats(0);

  band->Draw("AXIS");


  const int nb = totalMC->GetNbinsX();

TGraphAsymmErrors *gr = new TGraphAsymmErrors(nb);
int ip=0;
for (int i=1;i<=nb;++i){
  const double mc  = totalMC->GetBinContent(i);
  const double emc = totalMC->GetBinError(i);
  const double x   = totalMC->GetXaxis()->GetBinCenter(i);
  const double bw  = totalMC->GetXaxis()->GetBinWidth(i);
  const double rel = (mc>0.0)? emc/mc : 0.0;

  gr->SetPoint(ip, x, 1.0);
  gr->SetPointEXlow(ip,  0.5*bw);    // <-- fills continuously
  gr->SetPointEXhigh(ip, 0.5*bw);
  gr->SetPointEYlow(ip,  rel);
  gr->SetPointEYhigh(ip, rel);
  ++ip;
}
gr->SetFillColor(TColor::GetColor("#9c9ca1"));
gr->SetFillStyle(3001);
gr->SetLineColor(TColor::GetColor("#9c9ca1"));
gr->SetLineWidth(0);
gr->Draw("E2 SAME");


  hist_data->Draw("same p");


  TLegend* leg_ratio = new TLegend(0.20, 0.80, 0.38, 0.92);
leg_ratio->SetBorderSize(0);
leg_ratio->SetTextSize(0.07);
leg_ratio->AddEntry(gr, "Stat", "f");
leg_ratio->Draw();

    // Unity line
    TLine* line = new TLine(xmin1, 1.0, xmax1, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->SetLineWidth(1.);
    line->Draw();

    pad2->RedrawAxis();


    c1->SaveAs("bdt_com_stat_pas_2016_mergebin.png");
    c1->SaveAs("bdt_com_stat_pas_2016_mergebin.pdf");

}
