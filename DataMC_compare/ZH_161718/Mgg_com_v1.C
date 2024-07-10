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


void Mgg_com_v1(){
  Bool_t syst_on = false;
  gStyle->SetLineWidth(2);
  TCanvas *c1 =new TCanvas("c1", " ", 0, 0,700,800);

  c1->Range(0,0,1,1);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameBorderMode(0);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->Draw();

//  c1->SetLogy();

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();
  pad1->SetGridx();
  pad1->SetGridy();
  pad1->SetTickx(1);
  pad1->SetTicky(1);


  double lumi_16_APV = 19500;
  double lumi_16 = 16500;
  double lumi_17 = 41500;
  double lumi_18 = 59730;

  THStack *hs = new THStack("hs","");
  THStack *hs_mod = new THStack("hs_mod","");

  Double_t error;
  ////////////////////Data rootfiles/////////////////////

  TFile *fdata1 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016/full_reg/V1/DoubleEG_2016v2_isEleMu.root");
  TH1F  *h1 = (TH1F*)fdata1->Get("h_dipho_invmass");
  TFile *fdata2 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016/full_reg/V1/DoubleMuon_2016v2_isEleMu.root");
  TH1F  *h2 = (TH1F*)fdata2->Get("h_dipho_invmass");
  TFile *fdata3 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016APV/full_reg/V1/DoubleEG_2016PreAPVv2_isEleMu.root");                           // 2016
  TH1F  *h3 = (TH1F*)fdata3->Get("h_dipho_invmass");
  TFile *fdata4 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016APV/full_reg/V1/DoubleMuon_2016PreAPVv2_isEleMu.root");
  TH1F  *h4 = (TH1F*)fdata4->Get("h_dipho_invmass");
  h1->Add(h2);
  h1->Add(h3);
  h1->Add(h4);

  TFile *fdata5 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2017/full_reg/V1/DoubleEG_2017v2_isEleMu.root");
  TH1F  *h5 = (TH1F*)fdata5->Get("h_dipho_invmass");
  TFile *fdata6 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2017/full_reg/V1/DoubMu_2017V2_isEleMu.root");                             // 2017
  TH1F  *h6 = (TH1F*)fdata6->Get("h_dipho_invmass");
  h1->Add(h5);
  h1->Add(h6);

  TFile *fdata7 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2018/full_reg/V1/DoubMu_2018v3_isEleMu.root");
  TH1F  *h7 = (TH1F*)fdata7->Get("h_dipho_invmass");
  TFile *fdata8 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2018/full_reg/V1/egamma_2018v3_isEleMu.root");                                // 2018
  TH1F  *h8 = (TH1F*)fdata8->Get("h_dipho_invmass");
  h1->Add(h7);
  h1->Add(h8);

  h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX())+h1->GetBinContent(h1->GetNbinsX()+1));
  h1->SetBinContent(1, h1->GetBinContent(1)+h1->GetBinContent(0));


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
  }
  TGraphAsymmErrors *hdata1 = new TGraphAsymmErrors(n,x,y,EXlow,EXhigh,EYlow,EYhigh);

  double xmax = h1->GetXaxis()->GetXmax();
  double xmin = h1->GetXaxis()->GetXmin();
  double xmin1 = 0.0;
  double xmax1 = 200.0;

  //////////////////////MC files//////////////////////


  /////////////////////////////////////////////////////////////////    2016 MC     ////////////////////////////////////////////////////////////////////////


  double w1 = (lumi_16*4.078)/(5.0593e+06);
  TFile *f1 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016/full_reg/V1/TTGJets_2016v2_isEleMu.root");
  TH1F *hmc1 = (TH1F*)f1->Get("h_dipho_invmass");
  hmc1->SetBinContent(hmc1->GetNbinsX(), hmc1->GetBinContent(hmc1->GetNbinsX())+hmc1->GetBinContent(hmc1->GetNbinsX()+1));
  hmc1->SetBinContent(1, hmc1->GetBinContent(1)+hmc1->GetBinContent(0));
  hmc1->Scale(w1);

  double w2 = (lumi_16*365.34)/(4.28189e+10);
  TFile *f2 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016/full_reg/V1/TTSL_2016v2_isEleMu.root");
  TH1F *hmc2 = (TH1F*)f2->Get("h_dipho_invmass");
  hmc2->SetBinContent(hmc2->GetNbinsX(), hmc2->GetBinContent(hmc2->GetNbinsX())+hmc2->GetBinContent(hmc2->GetNbinsX()+1));
  hmc2->SetBinContent(1, hmc2->GetBinContent(1)+hmc2->GetBinContent(0));
  hmc2->Scale(w2);

  double w3 = (lumi_16*88.29)/(2.00013e+09);
  TFile *f3 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016/full_reg/V1/TTTo2L2Nu_2016v2_isEleMu.root");
  TH1F *hmc3 = (TH1F*)f3->Get("h_dipho_invmass");
  hmc3->SetBinContent(hmc3->GetNbinsX(), hmc3->GetBinContent(hmc3->GetNbinsX())+hmc3->GetBinContent(hmc3->GetNbinsX()+1));
  hmc3->SetBinContent(1, hmc3->GetBinContent(1)+hmc3->GetBinContent(0));
  hmc3->Scale(w3);

  double w4 = (lumi_16*5343)/(8.08308e+07);
  TFile *f4 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016/full_reg/V1/DYJetsToLL_2016v2_isEleMu.root");
  TH1F *hmc4 = (TH1F*)f4->Get("h_dipho_invmass");
  hmc4->SetBinContent(hmc4->GetNbinsX(), hmc4->GetBinContent(hmc4->GetNbinsX())+hmc4->GetBinContent(hmc4->GetNbinsX()+1));
  hmc4->SetBinContent(1, hmc4->GetBinContent(1)+hmc4->GetBinContent(0));
  hmc4->Scale(w4);

  double w7 = (lumi_16_APV*4.078)/(8.50926e+06);
  TFile *f7 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016APV/full_reg/V1/TTGJets_2016APVv2_isEleMu.root");
  TH1F *hmc7 = (TH1F*)f7->Get("h_dipho_invmass");
  hmc7->SetBinContent(hmc7->GetNbinsX(), hmc7->GetBinContent(hmc7->GetNbinsX())+hmc7->GetBinContent(hmc7->GetNbinsX()+1));
  hmc7->SetBinContent(1, hmc7->GetBinContent(1)+hmc7->GetBinContent(0));
  hmc7->Scale(w7);

  double w8 = (lumi_16_APV*365.34)/(3.97048e+10);
  TFile *f8 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016APV/full_reg/V1/TTSL_2016APVv2_isEleMu.root");
  TH1F *hmc8 = (TH1F*)f8->Get("h_dipho_invmass");
  hmc8->SetBinContent(hmc8->GetNbinsX(), hmc8->GetBinContent(hmc8->GetNbinsX())+hmc8->GetBinContent(hmc8->GetNbinsX()+1));
  hmc8->SetBinContent(1, hmc8->GetBinContent(1)+hmc8->GetBinContent(0));
  hmc8->Scale(w8);

  double w9 = (lumi_16_APV*88.29)/(1.85203e+09);
  TFile *f9 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016APV/full_reg/V1/TTTo2L2Nu_2016APVv2_isEleMu.root");
  TH1F *hmc9 = (TH1F*)f9->Get("h_dipho_invmass");
  hmc9->SetBinContent(hmc9->GetNbinsX(), hmc9->GetBinContent(hmc9->GetNbinsX())+hmc9->GetBinContent(hmc9->GetNbinsX()+1));
  hmc9->SetBinContent(1, hmc9->GetBinContent(1)+hmc9->GetBinContent(0));
  hmc9->Scale(w9);

  double w10 = (lumi_16_APV*5343)/(7.53999e+07);
  TFile *f10 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2016APV/full_reg/V1/DYJetsToLL_2016APVv2_isEleMu.root");
  TH1F *hmc10 = (TH1F*)f10->Get("h_dipho_invmass");
  hmc10->SetBinContent(hmc10->GetNbinsX(), hmc10->GetBinContent(hmc10->GetNbinsX())+hmc10->GetBinContent(hmc10->GetNbinsX()+1));
  hmc10->SetBinContent(1, hmc10->GetBinContent(1)+hmc10->GetBinContent(0));
  hmc10->Scale(w10);

  /////////////////////////////////////////////////////////////////    2017 MC     ////////////////////////////////////////////////////////////////////////


  double w13 = (lumi_17*4.078)/(1.60904e+07);
  TFile *f13 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2017/full_reg/V1/TTGJets_2017v2_isEleMu.root");
  TH1F *hmc13 = (TH1F*)f13->Get("h_dipho_invmass");
  hmc13->SetBinContent(hmc13->GetNbinsX(), hmc13->GetBinContent(hmc13->GetNbinsX())+hmc13->GetBinContent(hmc13->GetNbinsX()+1));
  hmc13->SetBinContent(1, hmc13->GetBinContent(1)+hmc13->GetBinContent(0));
  hmc13->Scale(w13);

  double w14 = (lumi_17*365.34)/(1.00204e+11);
  TFile *f14 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2017/full_reg/V1/TTSL_2017v2_isEleMu.root");
  TH1F *hmc14 = (TH1F*)f14->Get("h_dipho_invmass");
  hmc14->SetBinContent(hmc14->GetNbinsX(), hmc14->GetBinContent(hmc14->GetNbinsX())+hmc14->GetBinContent(hmc14->GetNbinsX()+1));
  hmc14->SetBinContent(1, hmc14->GetBinContent(1)+hmc14->GetBinContent(0));
  hmc14->Scale(w14);

  double w15 = (lumi_17*88.29)/(7.49227e+09);
  TFile *f15 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2017/full_reg/V1/TTTo2L2Nu_2017v2_isEleMu.root");
  TH1F *hmc15 = (TH1F*)f15->Get("h_dipho_invmass");
  hmc15->SetBinContent(hmc15->GetNbinsX(), hmc15->GetBinContent(hmc15->GetNbinsX())+hmc15->GetBinContent(hmc15->GetNbinsX()+1));
  hmc15->SetBinContent(1, hmc15->GetBinContent(1)+hmc15->GetBinContent(0));
  hmc15->Scale(w15);

  double w16 = (lumi_17*5343)/(6.18062e+07);
  TFile *f16 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2017/full_reg/V1/DYJetsToLL_2017v2_isEleMu.root");
  TH1F *hmc16 = (TH1F*)f16->Get("h_dipho_invmass");
  hmc16->SetBinContent(hmc16->GetNbinsX(), hmc16->GetBinContent(hmc16->GetNbinsX())+hmc16->GetBinContent(hmc16->GetNbinsX()+1));
  hmc16->SetBinContent(1, hmc16->GetBinContent(1)+hmc16->GetBinContent(0));
  hmc16->Scale(w16);

  /////////////////////////////////////////////////////////////////    2018 MC     ////////////////////////////////////////////////////////////////////////

  double w19 = (lumi_18*4.078)/(2.66557e+07);
  TFile *f19 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2018/full_reg/V1/TTGJets_v2_isEleMu.root");
  TH1F *hmc19 = (TH1F*)f19->Get("h_dipho_invmass");
  hmc19->SetBinContent(hmc19->GetNbinsX(), hmc19->GetBinContent(hmc19->GetNbinsX())+hmc19->GetBinContent(hmc19->GetNbinsX()+1));
  hmc19->SetBinContent(1, hmc19->GetBinContent(1)+hmc19->GetBinContent(0));
  hmc19->Scale(w19);

  double w20 = (lumi_18*365.34)/(1.48113e+11);
  TFile *f20 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2018/full_reg/V1/TTSL_v2_isEleMu.root");
  TH1F *hmc20 = (TH1F*)f20->Get("h_dipho_invmass");
  hmc20->SetBinContent(hmc20->GetNbinsX(), hmc20->GetBinContent(hmc20->GetNbinsX())+hmc20->GetBinContent(hmc20->GetNbinsX()+1));
  hmc20->SetBinContent(1, hmc20->GetBinContent(1)+hmc20->GetBinContent(0));
  hmc20->Scale(w20);

  double w21 = (lumi_18*88.29)/(1.03533e+10);
  TFile *f21 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2018/full_reg/V1/TTTo2L2Nu_v2_isEleMu.root");
  TH1F *hmc21 = (TH1F*)f21->Get("h_dipho_invmass");
  hmc21->SetBinContent(hmc21->GetNbinsX(), hmc21->GetBinContent(hmc21->GetNbinsX())+hmc21->GetBinContent(hmc21->GetNbinsX()+1));
  hmc21->SetBinContent(1, hmc21->GetBinContent(1)+hmc21->GetBinContent(0));
  hmc21->Scale(w21);

  double w22 = (lumi_18*5343)/(9.14154e+07);
  TFile *f22 = TFile::Open("/home/anirban/data/hh/workspace/root_files/CR_studies/new_v1/MVA_modified/ZH_2018/full_reg/V1/DYJetsToLL_v2_isEleMu.root");
  TH1F *hmc22 = (TH1F*)f22->Get("h_dipho_invmass");
  hmc22->SetBinContent(hmc22->GetNbinsX(), hmc22->GetBinContent(hmc22->GetNbinsX())+hmc22->GetBinContent(hmc22->GetNbinsX()+1));
  hmc22->SetBinContent(1, hmc22->GetBinContent(1)+hmc22->GetBinContent(0));
  hmc22->Scale(w22);


  TH1F *totalMC = (TH1F*)hmc1->Clone();
  totalMC->Add(hmc2);
  totalMC->Add(hmc3);
  totalMC->Add(hmc4);
  totalMC->Add(hmc7);
  totalMC->Add(hmc8);
  totalMC->Add(hmc9);
  totalMC->Add(hmc10);
  totalMC->Add(hmc13);
  totalMC->Add(hmc14);
  totalMC->Add(hmc15);
  totalMC->Add(hmc16);
  totalMC->Add(hmc19);
  totalMC->Add(hmc20);
  totalMC->Add(hmc21);
  totalMC->Add(hmc22);


  double scaletodata = h1->Integral()/totalMC->Integral();

  hmc1->Add(hmc7);
  hmc1->Add(hmc13);
  hmc1->Add(hmc19);
  hmc2->Add(hmc8);
  hmc2->Add(hmc14);
  hmc2->Add(hmc20);
  hmc3->Add(hmc9);
  hmc3->Add(hmc15);
  hmc3->Add(hmc21);
  hmc4->Add(hmc10);
  hmc4->Add(hmc16);
  hmc4->Add(hmc22);
/*  hsig1->Add(hsig3);
  hsig1->Add(hsig5);
  hsig1->Add(hsig7);
  hsig2->Add(hsig4);
  hsig2->Add(hsig6);
  hsig2->Add(hsig8);
*/

  totalMC->SetFillColor(14);
  gStyle->SetHatchesSpacing(0.7);
  gStyle->SetHatchesLineWidth(1.);
  totalMC->SetFillStyle(3354);

  hmc1->Scale(scaletodata);
//  hmc1->Scale(0.8);
  hmc1->SetFillColor(kAzure+2);
  hmc1->SetFillStyle(1001);
  hmc1->SetLineColor(kBlack);

  hmc2->Scale(scaletodata);
  hmc2->Scale(1.1);
  hmc2->SetFillColor(kTeal-9);
  hmc2->SetFillStyle(1001);
  hmc2->SetLineColor(kBlack);

  hmc3->Scale(scaletodata);
//  hmc3->Scale(0.9);
  hmc3->SetFillColor(kYellow-9);
  hmc3->SetFillStyle(1001);
  hmc3->SetLineColor(kBlack);

  hmc4->Scale(scaletodata);
  hmc4->SetFillColor(kRed-7);
  hmc4->SetFillStyle(1001);
  hmc4->SetLineColor(kBlack);

  totalMC->Scale(scaletodata);

  cout<<scaletodata<<endl;

  TH1F *totalMC_mod = (TH1F*)totalMC->Clone();

  TH1F *h_data = (TH1F*)h1->Clone();
  TH1F *hmc1_mod = (TH1F*)hmc1->Clone();
  TH1F *hmc2_mod = (TH1F*)hmc2->Clone();
  TH1F *hmc3_mod = (TH1F*)hmc3->Clone();
  TH1F *hmc4_mod =  (TH1F*)hmc4->Clone();

//Divide the bin content with bin width

  for(int ij=1; ij<=hmc1->GetNbinsX(); ++ij){
        double nmc1 = 1.0*hmc1_mod->GetBinContent(ij)/hmc1_mod->GetBinWidth(ij);
        hmc1_mod->SetBinContent(ij,nmc1);

        double nmc2 = 1.0*hmc2_mod->GetBinContent(ij)/hmc2_mod->GetBinWidth(ij);
        hmc2_mod->SetBinContent(ij,nmc2);

        double nmc3 = 1.0*hmc3_mod->GetBinContent(ij)/hmc3_mod->GetBinWidth(ij);
        hmc3_mod->SetBinContent(ij,nmc3);

        double nmc4 = 1.0*hmc4_mod->GetBinContent(ij)/hmc4_mod->GetBinWidth(ij);
        hmc4_mod->SetBinContent(ij,nmc4);

        double nData = 1.0*h_data->GetBinContent(ij)/h_data->GetBinWidth(ij);
        double BinError = h_data->GetBinError(ij)/h_data->GetBinWidth(ij);
        h_data->SetBinContent(ij,nData);
        h_data->SetBinError(ij,BinError);

        double ntot = 1.0*totalMC_mod->GetBinContent(ij)/totalMC_mod->GetBinWidth(ij);
        double mcBinError = totalMC_mod->GetBinError(ij)/totalMC_mod->GetBinWidth(ij);
        totalMC_mod->SetBinContent(ij,ntot); totalMC_mod->SetBinError(ij,mcBinError);
  }


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

  h1->SetLineColor(1);
  h1->Draw("pez same");
  h1->SetMarkerSize(0.8);
  h1->SetMarkerStyle(20);
/*
  hsig1->SetLineColor(kMagenta+2);
  hsig1->SetLineWidth(3);
  hsig2->SetLineColor(kGreen+3);
  hsig2->SetLineWidth(3);
  hsig1->Draw("hist same");
  hsig2->Draw("hist same");
*/
  hs->GetYaxis()->SetLabelSize(0.045);

  hs->GetYaxis()->SetTitle("Events / 10 GeV");
  hs->GetYaxis()->SetTitleOffset(1.0);
  hs->GetYaxis()->SetTitleSize(0.05);
  hs->SetMinimum(1.);
  hs->SetMaximum(300.);
  hs->GetXaxis()->SetLabelSize(0);
  hs->GetXaxis()->SetRangeUser(xmin1,xmax1);


  TLegend *legend1 = new TLegend(0.7, 0.7, 0.87, 0.87);
  legend1->SetTextFont(42);
  legend1->SetLineColor(0);
  legend1->SetTextSize(0.03);
  legend1->SetFillColor(0);
  legend1->AddEntry(h1, "Data", "lep");
  legend1->AddEntry(hmc4, "DY", "f");
  legend1->AddEntry(hmc1, "TTGJets", "f");
  legend1->AddEntry(hmc3, "TTTo2L2Nu", "f");
  legend1->AddEntry(hmc2, "TTToSemiLepton", "f");
//  legend1->AddEntry(hsig1, "WH mA 20 X 20", "l");
//  legend1->AddEntry(hsig2, "WH mA 55 X 20", "l");
  legend1->Draw();

  TLatex *t2a = new TLatex(0.51,0.9," #bf{CMS} #it{Preliminary}             137 fb^{-1} (13 TeV, full Run II) ");
  t2a->SetNDC();
  t2a->SetTextFont(42);
  t2a->SetTextSize(0.05);
  t2a->SetTextAlign(20);
  t2a->Draw("same");

  c1->cd();

//  TH1F *hist_data = (TH1F*)hdata1->Clone();
  Double_t ratioy[n], ratioEYlow[n], ratioEYhigh[n];

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

    }
  }

  TGraphAsymmErrors *hist_data = new TGraphAsymmErrors(n,x,ratioy,EXlow,EXhigh,ratioEYlow,ratioEYhigh);

  // Ratio Plot

  TPad *pad2 = new TPad("pad2", "newpad",0,0,1,0.3);
  pad2->Draw();
  pad2->cd();
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.3);
  pad2->SetRightMargin(0.9);
  pad2->SetFillStyle(0);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->SetTickx(1);
  pad2->SetTicky(1);
  TLine *line = new TLine(xmin1, 1.,xmax1, 1.);
  line->SetLineColor(kBlack);
  line->SetLineWidth(1.);

  hist_data->SetLineWidth(1.);
  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.8);
  hist_data->SetFillStyle(0);

  TH1F *band = (TH1F*)totalMC->Clone();
  band->Divide(totalMC);
  band->SetFillColor(kGray);
  band->SetLineColor(kWhite);
  band->SetFillStyle(1001);

//cout<< totalMC->GetNbinsX() << "  " << h1->GetBinContent(1) << endl;


  TH1F *h_JetEnUp; TH1F *h_JetEnDown; TH1F *h_JetResUp; TH1F *h_JetResDown; TH1F *h_UnclusteredEnUp;TH1F *h_UnclusteredEnDown;
  TGraphAsymmErrors *tottoterree; TGraphAsymmErrors *jererree; TGraphAsymmErrors *erree;

  TGraphAsymmErrors *staterree = (TGraphAsymmErrors*)staterror(h1,totalMC);
  staterree->SetFillColor(kRed-10);
  staterree->SetLineColor(kRed-10);
  staterree->SetFillStyle(1001);
  band->GetXaxis()->SetTitle("M_{#gamma#gamma} (GeV)");
  band->GetXaxis()->SetLabelSize(0.1);
  band->GetXaxis()->SetTitleSize(0.13);
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

  band->Draw("");
  staterree->Draw ("e2 same");
  hist_data->Draw("same pez");

  line->Draw("same");
  pad2->RedrawAxis();

  TLegend* legratio(0);
  legratio = new TLegend(0.12,0.35,0.4,0.45);
  legratio->SetFillColor(0);
  legratio->SetBorderSize(0);
  legratio-> SetNColumns(4);
  legratio->SetTextSize(0.07);
  legratio->AddEntry(staterree, "Sim. stat.","f");
//  legratio->Draw();
  c1->SaveAs("Mgg_com_v1.png");
}

