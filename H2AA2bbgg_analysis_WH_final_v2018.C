#include <iostream>
#include<fstream>
#include<string>
#include "TObject.h"
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"



using namespace std;

double PhiInRange(const double& phi) {
  double phiout = phi;

  if( phiout > 2*M_PI || phiout < -2*M_PI) {
    phiout = fmod( phiout, 2*M_PI);
  }
  if (phiout <= -M_PI) phiout += 2*M_PI;
  else if (phiout >  M_PI) phiout -= 2*M_PI;

  return phiout;
}


double delta2R(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(pow(eta1 - eta2,2) +pow(PhiInRange(phi1 - phi2),2));
}

//------------------implementation of pdf weight-----------------------//

float getTheoryEsystematics_PDF(int nLHEPDFWeights_b, float* LHEPDFWeights_b, bool isHessian)
{
	float pdferr=0;

        if(isHessian){
                for(int ilhe=1; ilhe<nLHEPDFWeights_b; ilhe++){
                        pdferr +=  pow(abs(LHEPDFWeights_b[ilhe]-LHEPDFWeights_b[0]),2);
                }
        }
        else{

                float pdfmean = 0;
                for(int ilhe=1; ilhe<nLHEPDFWeights_b; ilhe++){
                        pdfmean += LHEPDFWeights_b[ilhe];
                }
                pdfmean *= 1./(nLHEPDFWeights_b-1);   // nLHEPDFWeights-1 is used since nLHEPDFWeights includes the central value

                for(int ilhe=1; ilhe<nLHEPDFWeights_b; ilhe++){
                        pdferr +=  pow(abs(LHEPDFWeights_b[ilhe]-pdfmean),2);
                }
                pdferr *= 1./(nLHEPDFWeights_b-2.);  // nLHEPDFWeights-2 is used since nLHEPDFWeights includes the central value

                }

        pdferr = sqrt(pdferr);

        return pdferr;
}


float* Muon_SF(TFile *file_mu_sf, string id, float pt, float eta){

	char name[100];

	sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt",id.c_str());
	TH2F *h_SF = (TH2F*)file_mu_sf->Get(name);
	sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt_stat",id.c_str());
	TH2F *h_SF_stat = (TH2F*)file_mu_sf->Get(name);
	sprintf(name,"NUM_%sID_DEN_TrackerMuons_abseta_pt_syst",id.c_str());
	TH2F *h_SF_sys = (TH2F*)file_mu_sf->Get(name);

	int eta_bin_id = h_SF->GetXaxis()->FindBin(fabs(eta));
	int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);

	float sf, sf_stat, sf_sys, sf_err;

	if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
		sf = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
		sf_err = h_SF->GetBinError(eta_bin_id,pt_bin_id);
		sf_stat = h_SF_stat->GetBinContent(eta_bin_id,pt_bin_id);
		sf_sys = h_SF_sys->GetBinContent(eta_bin_id,pt_bin_id);
	}else{
		sf = 1;
		sf_err = 0;
		sf_stat = sf_sys = 1;
		}

	static float sfvalues[5];
	sfvalues[0] = sf;
	sfvalues[1] = sf+sf_err;
	sfvalues[2] = sf-sf_err;
	sfvalues[3] = sf_stat;
	sfvalues[4] = sf_sys;

	return sfvalues;
}


float* Electron_SF(TFile *file_el_sf, float pt, float eta){

	char name[100];

	TH2F *h_SF = (TH2F*)file_el_sf->Get("EGamma_SF2D");
	TH2F *h_SF_statData = (TH2F*)file_el_sf->Get("statData");
	TH2F *h_SF_statMC = (TH2F*)file_el_sf->Get("statMC");
	TH2F *h_SF_altBkgModel = (TH2F*)file_el_sf->Get("altBkgModel");
	TH2F *h_SF_altSignalModel = (TH2F*)file_el_sf->Get("altSignalModel");
	TH2F *h_SF_altMCEff = (TH2F*)file_el_sf->Get("altMCEff");
	TH2F *h_SF_altTagSelection = (TH2F*)file_el_sf->Get("altTagSelection");

	int eta_bin_id = h_SF->GetXaxis()->FindBin(eta);
	int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);

	static float sfvalues[9] = {-100,-100,-100,-100,-100,-100,-100,-100,-100};

	if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
		sfvalues[0] = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[1] = sfvalues[0]+h_SF->GetBinError(eta_bin_id,pt_bin_id);
		sfvalues[2] = sfvalues[0]-h_SF->GetBinError(eta_bin_id,pt_bin_id);
		sfvalues[3] = h_SF_statData->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[4] = h_SF_statMC->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[5] = h_SF_altBkgModel->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[6] = h_SF_altSignalModel->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[7] = h_SF_altMCEff->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[8] = h_SF_altTagSelection->GetBinContent(eta_bin_id,pt_bin_id);
	}
	else{
		sfvalues[0] = sfvalues[1] = sfvalues[2] = 1;
		sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = sfvalues[7] = sfvalues[8] = 0;
		}

	return sfvalues;
}


float* Photon_SF(TFile *file_pho_sf, float pt, float eta){

        char name[100];

        TH2F *h_SF = (TH2F*)file_pho_sf->Get("EGamma_SF2D");
        TH2F *h_SF_statData = (TH2F*)file_pho_sf->Get("statData");
        TH2F *h_SF_statMC = (TH2F*)file_pho_sf->Get("statMC");
        TH2F *h_SF_altBkgModel = (TH2F*)file_pho_sf->Get("altBkgModel");
        TH2F *h_SF_altSignalModel = (TH2F*)file_pho_sf->Get("altSignalModel");
        TH2F *h_SF_altMCEff = (TH2F*)file_pho_sf->Get("altMCEff");
        TH2F *h_SF_altTagSelection = (TH2F*)file_pho_sf->Get("altTagSelection");

        int eta_bin_id = h_SF->GetXaxis()->FindBin(eta);
        int pt_bin_id = h_SF->GetYaxis()->FindBin(pt);

        static float sfvalues[9] = {-100,-100,-100,-100,-100,-100,-100,-100,-100};

        if(eta_bin_id>0 && eta_bin_id<=(h_SF->GetNbinsX()) && pt_bin_id>0 && pt_bin_id<=(h_SF->GetNbinsY())){
                sfvalues[0] = h_SF->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[1] = sfvalues[0]+h_SF->GetBinError(eta_bin_id,pt_bin_id);
                sfvalues[2] = sfvalues[0]-h_SF->GetBinError(eta_bin_id,pt_bin_id);
                sfvalues[3] = h_SF_statData->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[4] = h_SF_statMC->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[5] = h_SF_altBkgModel->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[6] = h_SF_altSignalModel->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[7] = h_SF_altMCEff->GetBinContent(eta_bin_id,pt_bin_id);
                sfvalues[8] = h_SF_altTagSelection->GetBinContent(eta_bin_id,pt_bin_id);
        }
        else{
                sfvalues[0] = sfvalues[1] = sfvalues[2] = 1;
                sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = sfvalues[7] = sfvalues[8] = 0;
                }

        return sfvalues;
}
 

float* Get_PU_Weights(TFile *file_pu_ratio, int npu){

	TH1F *h_data = (TH1F*)file_pu_ratio->Get("pileup_weight");
	TH1F *h_data_plus = (TH1F*)file_pu_ratio->Get("pileup_plus_weight");
	TH1F *h_data_minus = (TH1F*)file_pu_ratio->Get("pileup_minus_weight");

	int bin_id = h_data->FindBin(npu);

	static float puweight[3] = {0,0,0};
	if(bin_id>=0 && bin_id<100){
		puweight[0] = h_data->GetBinContent(bin_id);
		puweight[1] = h_data_plus->GetBinContent(bin_id);
		puweight[2] = h_data_minus->GetBinContent(bin_id);
	}
	return puweight;
}



  static const int njetmx = 100;
  static const int njetmxAK8 =100;
  static const int npartmx = 100;
  static const int nconsmax = 1000;
  static const int njetconsmax = 3;
  static const int ngenjetAK8mx =100;
  static const int nlhepdfmax = 103;

  BTagCalibration calib_deepcsv, calib_deepflav;
  BTagCalibrationReader reader_deepcsv, reader_deepflav;

  double Generator_weight;
  double weights[njetmx];

  double prefiringweight, prefiringweightup, prefiringweightdown;

  bool isMC;
  bool isFastSIM;
  bool isDATA;

  int nPFJetAK4, PFJetAK4_hadronflav[njetmx], PFJetAK4_partonflav[njetmx];
  float PFJetAK4_pt[njetmx], PFJetAK4_eta[njetmx], PFJetAK4_y[njetmx], PFJetAK4_phi[njetmx], PFJetAK4_mass[njetmx], PFJetAK4_energy[njetmx], PFJetAK4_PUID[njetmx];
  float PFJetAK4_btag_DeepCSV[njetmx], PFJetAK4_btag_DeepFlav[njetmx];
  bool PFJetAK4_jetID[njetmx], PFJetAK4_jetID_tightlepveto[njetmx];
  float PFJetAK4_btag_DeepCSV_SF[njetmx], PFJetAK4_btag_DeepCSV_SF_up[njetmx], PFJetAK4_btag_DeepCSV_SF_dn[njetmx];
  float PFJetAK4_btag_DeepFlav_SF[njetmx], PFJetAK4_btag_DeepFlav_SF_up[njetmx], PFJetAK4_btag_DeepFlav_SF_dn[njetmx];
  float PFJetAK4_reso[njetmx], PFJetAK4_resoup[njetmx], PFJetAK4_resodn[njetmx];
  float PFJetAK4_pt_reg[njetmx], PFJetAK4_mass_reg[njetmx], PFJetAK4_energy_reg[njetmx], PFJetAK4_bcorr[njetmx], PFJetAK4_breso[njetmx];
  float PFJetAK4_JEC[njetmx];
  float PFJetAK4_jesup_Total[njetmx], PFJetAK4_jesdn_Total[njetmx];

  int nGenPart;
  int GenPart_status[npartmx], GenPart_pdg[npartmx], GenPart_mompdg[npartmx], GenPart_momstatus[npartmx], GenPart_grmompdg[npartmx], GenPart_momid[npartmx], GenPart_daugno[npartmx];
  float GenPart_pt[npartmx], GenPart_eta[npartmx], GenPart_phi[npartmx], GenPart_mass[npartmx]; //GenPart_q[npartmx];
  bool GenPart_fromhard[npartmx], GenPart_fromhardbFSR[npartmx], GenPart_isPromptFinalState[npartmx], GenPart_isLastCopyBeforeFSR[npartmx], GenPart_isDirectPromptTauDecayProductFinalState[npartmx];

  double LHE_weight, Rho;
  float BTAG_SF, BTAG_jes_up, BTAG_jes_dn;

  int nMuon;
  float Muon_minchiso[njetmx], Muon_minnhiso[njetmx], Muon_minphiso[njetmx], Muon_minisoall[njetmx];
  float Muon_charge[njetmx], Muon_p[njetmx], Muon_pt[njetmx], Muon_eta[njetmx], Muon_phi[njetmx], Muon_e[njetmx], Muon_dz[njetmx], Muon_ip3d[njetmx], Muon_ptErr[njetmx], Muon_chi[njetmx], Muon_ecal[njetmx], Muon_hcal[njetmx];
  float Muon_posmatch[njetmx], Muon_trkink[njetmx], Muon_segcom[njetmx], Muon_pfiso[njetmx], Muon_dxy[njetmx], Muon_dxyErr[njetmx], Muon_hit[njetmx], Muon_pixhit[njetmx], Muon_mst[njetmx], Muon_trklay[njetmx], Muon_valfrac[njetmx],Muon_dxy_sv[njetmx];
  int Muon_ndf[njetmx];
  bool Muon_isPF[njetmx], Muon_isGL[njetmx], Muon_isTRK[njetmx];
  bool Muon_isGoodGL[njetmx], Muon_isTight[njetmx], Muon_isHighPt[njetmx], Muon_isHighPttrk[njetmx], Muon_isMed[njetmx], Muon_isMedPr[njetmx], Muon_isLoose[njetmx], Muon_TightID[njetmx];
  float Muon_corrected_pt[njetmx], Muon_correctedUp_pt[njetmx], Muon_correctedDown_pt[njetmx];

  int nElectron;
  bool Electron_mvaid_Fallv2WP90[njetmx], Electron_mvaid_Fallv2WP90_noIso[njetmx];
  bool Electron_mvaid_Fallv2WP80[njetmx], Electron_mvaid_Fallv2WP80_noIso[njetmx];
  float Electron_charge[njetmx], Electron_pt[njetmx], Electron_eta[njetmx], Electron_phi[njetmx], Electron_e[njetmx], Electron_e_ECAL[njetmx], Electron_p[njetmx];
  float Electron_dxy[njetmx],  Electron_dxyErr[njetmx], Electron_dxy_sv[njetmx], Electron_dz[njetmx], Electron_dzErr[njetmx], Electron_ip3d[njetmx];
  float Electron_hovere[njetmx], Electron_qovrper[njetmx], Electron_chi[njetmx]; //Electron_emiso03[njetmx], Electron_hadiso03[njetmx], Electron_emiso04[njetmx], Electron_hadiso04[njetmx];
  float Electron_eoverp[njetmx], Electron_ietaieta[njetmx], Electron_etain[njetmx], Electron_phiin[njetmx], Electron_fbrem[njetmx];
  float Electron_nohits[njetmx], Electron_misshits[njetmx];
  float Electron_pfiso_drcor[njetmx];
  float Electron_pfiso_eacor[njetmx];
  float Electron_pfiso04_eacor[njetmx];
  int Electron_ndf[njetmx];
  float Electron_eccalTrkEnergyPostCorr[njetmx];
  float Electron_energyScaleValue[njetmx];
  float Electron_energyScaleUp[njetmx];
  float Electron_energyScaleDown[njetmx];
  float Electron_energySigmaValue[njetmx];
  float Electron_energySigmaUp[njetmx];
  float Electron_energySigmaDown[njetmx];
  float Electron_supcl_eta[njetmx];
  float Electron_supcl_phi[njetmx];
  float Electron_supcl_e[njetmx];
  float Electron_supcl_rawE[njetmx];
  float Electron_sigmaieta[njetmx];
  float Electron_sigmaiphi[njetmx];
  float Electron_r9full[njetmx];
  float Electron_supcl_etaw[njetmx];
  float Electron_supcl_phiw[njetmx];
  float Electron_hcaloverecal[njetmx];
  float Electron_cloctftrkn[njetmx];
  float Electron_cloctftrkchi2[njetmx];
  float Electron_e1x5bye5x5[njetmx];
  float Electron_normchi2[njetmx];
  float Electron_hitsmiss[njetmx];
  float Electron_trkmeasure[njetmx];
  float Electron_convtxprob[njetmx];
  float Electron_ecloverpout[njetmx];
  float Electron_ecaletrkmomentum[njetmx];
  float Electron_deltaetacltrkcalo[njetmx];
  float Electron_supcl_preshvsrawe[njetmx];
  bool Electron_convVeto[njetmx];
  float Electron_pfisolsumphet[njetmx];
  float Electron_pfisolsumchhadpt[njetmx];
  float Electron_pfsiolsumneuhadet[njetmx];
  float Electron_minchiso[njetmx];
  float Electron_minnhiso[njetmx];
  float Electron_minphiso[njetmx];
  float Electron_minisoall[njetmx];
  float Electron_ecalEnergyPreCorr[njetmx];
  float Electron_ecalEnergyPostCorr[njetmx];

  int nPhoton;
  bool Photon_passEveto[njetmx];
  bool Photon_PixelSeed[njetmx];
  bool Photon_mvaid_Fall17V2_WP90[njetmx];
  bool Photon_mvaid_Fall17V2_WP80[njetmx];
  float Photon_mvaid_Fall17V2_raw[njetmx];
  bool Photon_mvaid_Spring16V1_WP90[njetmx];
  bool Photon_mvaid_Spring16V1_WP80[njetmx];
  float Photon_e[njetmx];
  float Photon_pt[njetmx];
  float Photon_eta[njetmx];
  float Photon_phi[njetmx];
  float Photon_e1by9[njetmx];
  float Photon_e9by25[njetmx];
  float Photon_hadbyem[njetmx];
  float Photon_trkiso[njetmx];
  float Photon_emiso[njetmx];
  float Photon_hadiso[njetmx];
  float Photon_chhadiso[njetmx];
  float Photon_neuhadiso[njetmx];
  float Photon_PUiso[njetmx];
  float Photon_phoiso[njetmx];
  float Photon_ietaieta[njetmx];
  float Photon_ecalEnergyPreCorr[njetmx];
  float Photon_ecalEnergyPostCorr[njetmx];

  int nGenJetAK4;
  float GenJetAK4_pt[njetmx];
  float GenJetAK4_eta[njetmx];
  float GenJetAK4_phi[njetmx];
  float GenJetAK4_mass[njetmx];
  int GenJetAK4_hadronflav[njetmx];

  float Generator_qscale, Generator_x1, Generator_x2, Generator_xpdf1, Generator_xpdf2, Generator_scalePDF;
  int Generator_id1, Generator_id2;
  int nLHEPDFWeights;
  float LHEPDFWeights[nlhepdfmax], LHEPDF_err;

  float miset , misphi , sumEt, misetsig, Met, SumEt, MetPhi ;
  int jet_no;

  int npu_vert;
  int npu_vert_true;
  
  int trig_value;


  float leppt, lepy, lepeta, lepphi, lepe;

  int b1_hadFlv, b2_hadFlv, b1_partFlv, b2_partFlv;
  float bjet_no, b1pt, b1y, b1eta, b1phi, b1e, b2pt, b2y, b2eta, b2phi, b2e, bb_inv_mass, b1_DeepFlv, b2_DeepFlv, b1_PUid, b2_PUid, b1pt_reg, b1e_reg, b2pt_reg, b2e_reg, b1pt_scl, b1e_scl, b2pt_scl, b2e_scl;
  float b1pt_jecup, b1e_jecup, b2pt_jecup, b2e_jecup, b1pt_jecdn, b1e_jecdn, b2pt_jecdn, b2e_jecdn, b1pt_resup, b1e_resup, b2pt_resup, b2e_resup, b1pt_resdn, b1e_resdn, b2pt_resdn, b2e_resdn;
  float b1pt_final, b1y_final, b1eta_final, b1phi_final, b1e_final, b2pt_final, b2y_final, b2eta_final, b2phi_final, b2e_final, bb_inv_mass_final, invmassbbgg_final;
  float b1pt_jesup, b1y_jesup, b1eta_jesup, b1phi_jesup, b1e_jesup, b2pt_jesup, b2y_jesup, b2eta_jesup, b2phi_jesup, b2e_jesup, bb_inv_mass_jesup;
  float b1pt_jesdn, b1y_jesdn, b1eta_jesdn, b1phi_jesdn, b1e_jesdn, b2pt_jesdn, b2y_jesdn, b2eta_jesdn, b2phi_jesdn, b2e_jesdn, bb_inv_mass_jesdn;
  float b1pt_resoup, b1y_resoup, b1eta_resoup, b1phi_resoup, b1e_resoup, b2pt_resoup, b2y_resoup, b2eta_resoup, b2phi_resoup, b2e_resoup, bb_inv_mass_resoup;
  float b1pt_resodn, b1y_resodn, b1eta_resodn, b1phi_resodn, b1e_resodn, b2pt_resodn, b2y_resodn, b2eta_resodn, b2phi_resodn, b2e_resodn, bb_inv_mass_resodn;

  float pho_no, pho1pt, pho1y, pho1eta, pho1phi, pho1e, pho2pt, pho2y, pho2eta, pho2phi, pho2e, pho1r9, pho2r9, dipho_invmass, pho1MVA, pho2MVA;

  float invmassbbgg, invmassbbgg_jesup, invmassbbgg_jesdn, invmassbbgg_resoup, invmassbbgg_resodn;
  float bb_gg_dphi, bb_gg_dphi_jesup, bb_gg_dphi_jesdn, bb_gg_dphi_resoup, bb_gg_dphi_resodn;

  double xsec_weight, weight_nom, weight_PU_up, weight_leptonsf_up, weight_photonsf_up, weight_prefiring_up, weight_pujid_up, weight_pujid_dn, weight_bSF_up, weight_trig_up, weight_PU_dn, weight_leptonsf_dn, weight_photonsf_dn, weight_prefiring_dn, weight_bSF_dn, weight_trig_dn, weight_haspixsf_up, weight_haspixsf_dn;
  double jecup_wt, jecdn_wt, jerup_wt, jerdn_wt;
  float Jet_pt_nom[njetmx], Jet_mass_nom[njetmx], Jet_pt_jesup[njetmx], Jet_mass_jesup[njetmx], Jet_pt_jesdn[njetmx], Jet_mass_jesdn[njetmx], Jet_pt_resoup[njetmx], Jet_mass_resoup[njetmx], Jet_pt_resodn[njetmx], Jet_mass_resodn[njetmx];
  float dR_leadpho_genpart[npartmx], dR_subleadpho_genpart[npartmx];
  int leadpho_corr_pdg, subleadpho_corr_pdg;

  float arr1[npartmx], arr2[npartmx], sval1, sval2;
  int pos1, pos2, leadpho_matched_genpdg, subleadpho_matched_genpdg, leadpho_matched_genmompdg, subleadpho_matched_genmompdg;

  float genb1_pt, genb1_eta, genb1_phi, genb1_mass, genb2_pt, genb2_eta, genb2_phi, genb2_mass, genb_invmass;
  double evt_rho, puwt, pujidwt;
  int puvtx;
  int aa, ab, bb;
  int bb_jecup, bb_jecdn, bb_jerup, bb_jerdn;
  bool isEle, isMu;
  bool no_CR_1, no_CR_2;


  int year = 2018; 
//  int year = 2017;
//  int year = 2016;

  string muon_id_name = "Tight";
  string electron_id_name = "wp90noiso";

  float DAK4_T = 0.71;
  float DAK4_M = 0.2783;
  float DAK4_L = 0.0490; 

//  2018       
  bool hlt_IsoMu24;
  bool hlt_Ele32_WPTight_Gsf;
  bool hlt_Mu17_Mu8;
  bool hlt_Ele23_Ele12;
  bool hlt_DoubleEle33;
  bool hlt_DoubleEle25_CaloIdL_MW; 


//  2017
/*  bool hlt_IsoMu24;
  bool hlt_Ele32_WPTight_Gsf;
  bool hlt_Mu17_Mu8;
  bool hlt_Ele23_Ele12; 
*/
//  2016
/*  bool hlt_IsoMu24;
  bool hlt_IsoTkMu24;
  bool hlt_Ele27_WPTight_Gsf;
  bool hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  bool hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  bool hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  bool hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  bool hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  bool hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  bool hlt_DoubleEle25_CaloIdL_MW; 
*/

int main(int argc, char *argv[])
{
	isMC = true;
//	isDATA = true;
        isFastSIM = false;
	char fOut[50];
        string inputFile=argv[3];
        string path="/home/abala/cms/CMSSW_10_5_0/src/condor_job/";

        calib_deepflav = BTagCalibration("DeepJet", "BtagRecommendation106XUL18/DeepJet_106XUL18SF_WPonly_V1p1.csv");
        reader_deepflav = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
        reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
        reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
        reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");

	TFile *f_Eff = new TFile("effcyPUID_81Xtraining.root");
	TH2D* Eff           = (TH2D*)f_Eff->Get("h2_eff_mc2018_M");
	TH2D* MisTag        = (TH2D*)f_Eff->Get("h2_mistag_mc2018_M");

	TFile *f_SF = new TFile("scalefactorsPUID_81Xtraining.root");
    	TH2D* SF_Eff = (TH2D * ) f_SF -> Get("h2_eff_sf2018_M");
    	TH2D* SF_MisTag = (TH2D * ) f_SF -> Get("h2_mistag_sf2018_M");
    	TH2D* SF_Eff_Err = (TH2D * ) f_SF -> Get("h2_eff_sf2018_M_Systuncty");
    	TH2D* SF_MisTag_Err = (TH2D * ) f_SF -> Get("h2_mistag_sf2018_M_Systuncty");

        char name[1000];

//	2018

        TFile *file_mu_sf;
        sprintf(name,"data/Efficiencies_muon_generalTracks_Z_Run%i_UL_ID.root",year);
        file_mu_sf = new TFile(name,"read");
   
        TFile *file_el_sf;
        sprintf(name,"data/egammaEffi.txt_Ele_%s_EGM2D_UL%i.root",electron_id_name.c_str(),year);
        file_el_sf = new TFile(name,"read");

	TFile *file_pho_sf;
        sprintf(name,"EGamma_SF/egammaEffi.txt_EGM2D_Pho_wp90.root_UL18.root");
        file_pho_sf = new TFile(name,"read");
   
        TFile *file_pu_ratio;
        sprintf(name,"data/pileup/RatioPileup-UL%i-100bins.root",year);
        file_pu_ratio = new TFile(name,"read");

	TFile *f_pho_hasPix = new TFile("EGamma_SF/HasPix_SummaryPlot_UL18.root");
        TH1F *h1 = (TH1F*) f_pho_hasPix->Get("MVAID/SF_HasPix_MVAID");
        TH1F *h2 = (TH1F*) f_pho_hasPix->Get("MVAID/Staunc_HasPix_MVAID");
        TH1F *h3 = (TH1F*) f_pho_hasPix->Get("MVAID/PUunc_HasPix_MVAID");
        TH1F *h4 = (TH1F*) f_pho_hasPix->Get("MVAID/Modelunc_HasPix_MVAID");

	float SF_HasPix_EB = h1->GetBinContent(1);
        float SF_HasPix_EE = h1->GetBinContent(4);
 	float Staunc_HasPix_EB = h2->GetBinContent(1);
        float Staunc_HasPix_EE = h2->GetBinContent(4);
	float PUunc_HasPix_EB = h3->GetBinContent(1);
        float PUunc_HasPix_EE = h3->GetBinContent(4);
	float Modelunc_HasPix_EB = h4->GetBinContent(1);
        float Modelunc_HasPix_EE = h4->GetBinContent(4);

	float SF_HasPix_EB_up = SF_HasPix_EB + sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB + Modelunc_HasPix_EB * Modelunc_HasPix_EB);
        float SF_HasPix_EB_dn = SF_HasPix_EB - sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB + Modelunc_HasPix_EB * Modelunc_HasPix_EB);
	float SF_HasPix_EE_up = SF_HasPix_EE + sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE + Modelunc_HasPix_EE * Modelunc_HasPix_EE);
        float SF_HasPix_EE_dn = SF_HasPix_EE - sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE + Modelunc_HasPix_EE * Modelunc_HasPix_EE);
	
	
//	2017
	
/*	TFile *file_mu_sf;
        sprintf(name,"data/Efficiencies_muon_generalTracks_Z_Run%i_UL_ID.root",year);
        file_mu_sf = new TFile(name,"read");

        TFile *file_el_sf;
        sprintf(name,"data/egammaEffi.txt_Ele_%s_EGM2D_UL%i.root",electron_id_name.c_str(),year);
        file_el_sf = new TFile(name,"read");

        TFile *file_pho_sf;
        sprintf(name,"EGamma_SF/egammaEffi.txt_EGM2D_PHO_MVA90_UL17.root");
        file_pho_sf = new TFile(name,"read");

        TFile *file_pu_ratio;
        sprintf(name,"data/pileup/RatioPileup-UL%i-100bins.root",year);
        file_pu_ratio = new TFile(name,"read");

        TFile *f_pho_hasPix = new TFile("EGamma_SF/HasPix_SummaryPlot_UL17.root");
        TH1F *h1 = (TH1F*) f_pho_hasPix->Get("MVAID/SF_HasPix_MVAID");
        TH1F *h2 = (TH1F*) f_pho_hasPix->Get("MVAID/Staunc_HasPix_MVAID");
        TH1F *h3 = (TH1F*) f_pho_hasPix->Get("MVAID/PUunc_HasPix_MVAID");
        TH1F *h4 = (TH1F*) f_pho_hasPix->Get("MVAID/Modelunc_HasPix_MVAID");

        float SF_HasPix_EB = h1->GetBinContent(1);
        float SF_HasPix_EE = h1->GetBinContent(4);
        float Staunc_HasPix_EB = h2->GetBinContent(1);
        float Staunc_HasPix_EE = h2->GetBinContent(4);
        float PUunc_HasPix_EB = h3->GetBinContent(1);
        float PUunc_HasPix_EE = h3->GetBinContent(4);
        float Modelunc_HasPix_EB = h4->GetBinContent(1);
        float Modelunc_HasPix_EE = h4->GetBinContent(4);

        float SF_HasPix_EB_up = SF_HasPix_EB + sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB + Modelunc_HasPix_EB * Modelunc_HasPix_EB);
        float SF_HasPix_EB_dn = SF_HasPix_EB - sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB + Modelunc_HasPix_EB * Modelunc_HasPix_EB);
        float SF_HasPix_EE_up = SF_HasPix_EE + sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE + Modelunc_HasPix_EE * Modelunc_HasPix_EE);
        float SF_HasPix_EE_dn = SF_HasPix_EE - sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE + Modelunc_HasPix_EE * Modelunc_HasPix_EE);
*/	
	
//	2016APV
	
/*	TFile *file_mu_sf;
        sprintf(name,"data/Efficiencies_muon_generalTracks_Z_Run%i_UL_HIPM_ID.root",year);
        file_mu_sf = new TFile(name,"read");

        TFile *file_el_sf;
        sprintf(name,"data/egammaEffi.txt_Ele_%s_preVFP_EGM2D_UL%i.root",electron_id_name.c_str(),year);
        file_el_sf = new TFile(name,"read");

        TFile *file_pho_sf;
        sprintf(name,"EGamma_SF/egammaEffi.txt_EGM2D_Pho_wp90_UL16.root");
        file_pho_sf = new TFile(name,"read");

        TFile *file_pu_ratio;
        sprintf(name,"data/pileup/RatioPileup-UL%i-100bins.root",year);
        file_pu_ratio = new TFile(name,"read");

        TFile *f_pho_hasPix = new TFile("EGamma_SF/HasPix_SummaryPlot_UL16_preVFP.root");
        TH1F *h1 = (TH1F*) f_pho_hasPix->Get("MVAID/SF_HasPix_MVAID");
        TH1F *h2 = (TH1F*) f_pho_hasPix->Get("MVAID/Staunc_HasPix_MVAID");
        TH1F *h3 = (TH1F*) f_pho_hasPix->Get("MVAID/PUunc_HasPix_MVAID");

        float SF_HasPix_EB = h1->GetBinContent(1);
        float SF_HasPix_EE = h1->GetBinContent(4);
        float Staunc_HasPix_EB = h2->GetBinContent(1);
        float Staunc_HasPix_EE = h2->GetBinContent(4);
        float PUunc_HasPix_EB = h3->GetBinContent(1);
        float PUunc_HasPix_EE = h3->GetBinContent(4);

        float SF_HasPix_EB_up = SF_HasPix_EB + sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB);
        float SF_HasPix_EB_dn = SF_HasPix_EB - sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB);
        float SF_HasPix_EE_up = SF_HasPix_EE + sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE);
        float SF_HasPix_EE_dn = SF_HasPix_EE - sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE);
*/	

//	2016
	
/*	TFile *file_mu_sf;
        sprintf(name,"data/Efficiencies_muon_generalTracks_Z_Run%i_UL_ID.root",year);
        file_mu_sf = new TFile(name,"read");

        TFile *file_el_sf;
        sprintf(name,"data/egammaEffi.txt_Ele_%s_postVFP_EGM2D_UL%i.root",electron_id_name.c_str(),year);
        file_el_sf = new TFile(name,"read");

        TFile *file_pho_sf;
        sprintf(name,"EGamma_SF/egammaEffi.txt_EGM2D_Pho_MVA90_UL16_postVFP.root");
        file_pho_sf = new TFile(name,"read");

        TFile *file_pu_ratio;
        sprintf(name,"data/pileup/RatioPileup-UL%i-100bins.root",year);
        file_pu_ratio = new TFile(name,"read");

        TFile *f_pho_hasPix = new TFile("EGamma_SF/HasPix_SummaryPlot_UL16_postVFP.root");
        TH1F *h1 = (TH1F*) f_pho_hasPix->Get("MVAID/SF_HasPix_MVAID");
        TH1F *h2 = (TH1F*) f_pho_hasPix->Get("MVAID/Staunc_HasPix_MVAID");
        TH1F *h3 = (TH1F*) f_pho_hasPix->Get("MVAID/PUunc_HasPix_MVAID");

        float SF_HasPix_EB = h1->GetBinContent(1);
        float SF_HasPix_EE = h1->GetBinContent(4);
        float Staunc_HasPix_EB = h2->GetBinContent(1);
        float Staunc_HasPix_EE = h2->GetBinContent(4);
        float PUunc_HasPix_EB = h3->GetBinContent(1);
        float PUunc_HasPix_EE = h3->GetBinContent(4);

        float SF_HasPix_EB_up = SF_HasPix_EB + sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB);
        float SF_HasPix_EB_dn = SF_HasPix_EB - sqrt(Staunc_HasPix_EB * Staunc_HasPix_EB + PUunc_HasPix_EB * PUunc_HasPix_EB);
        float SF_HasPix_EE_up = SF_HasPix_EE + sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE);
        float SF_HasPix_EE_dn = SF_HasPix_EE - sqrt(Staunc_HasPix_EE * Staunc_HasPix_EE + PUunc_HasPix_EE * PUunc_HasPix_EE);
*/	


	if(inputFile=="2018_vf2/WH_mA20_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_20_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/WH_mA25_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_25_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WH_mA30_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_30_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WH_mA35_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_35_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WH_mA35_new_v2t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_35_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WH_mA40_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_40_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/WH_mA45_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_45_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/WH_mA50_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_50_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/WH_mA55_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_55_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/WH_mA60_new_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2018/WH/WH_mA_60_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/TTGJets_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGJets/2018/WH/TTGJets_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/TTGG_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGG/2018/WH/TTGG_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/DYJetsToLL_M10to50_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2018/WH/DYJetsToLL_M10to50_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/DYJetsToLL_M50_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2018/WH/DYJetsToLL_M50_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/TTTo2L2Nu_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTTo2L2Nu/2018/WH/TTTo2L2Nu_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/TTToSemiLeptonic_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTSL/2018/WH/TTSL_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WW_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/WW/2018/WH/WW_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/ZZ_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/ZZ/2018/WH/ZZ_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WZ_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/WZ/2018/WH/WZ_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/WGToLNuG_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/WGToLNuG/2018/WH/WGToLNuG_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/ZGToLLG_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/ZGToLLG/2018/WH/ZGToLLG_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/DYJetsToLL_M10to50_nlo_2018.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2018/WH/DYJetsToLL_M10to50_nlo_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/DYJetsToLL_M50_nlo_2018.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2018/WH/DYJetsToLL_M50_nlo_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2018_vf2/Egamma_2018_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/egamma/2018/WH/Egamma_2018_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2018_vf2/SingleMuon_new_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinMu/2018/SinMu_2018_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA20_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_20_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA25_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_25_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA30_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_30_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA35_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_35_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA40_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_40_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA45_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_45_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA50_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_50_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA55_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_55_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/WH_mA60_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016APV/WH/WH_mA_60_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/TTGJets_16APV_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGJets/2016APV/WH/TTGJets_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2016_vf2_APV/DYJetsToLL_16APV_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2016APV/WH/DYJetsToLL_M50_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2016_vf2_APV/TTTo2L2Nu_16APV_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTTo2L2Nu/2016APV/WH/TTTo2L2Nu_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2016_vf2_APV/TTToSemiLeptonic_16APV_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTSL/2016APV/WH/TTSL_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/SingleElectron_16APV_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinEle/2016APV/SingleElectron_2016APV_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2_APV/SingleMuon_16APV_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinMu/2016APV/SinMu_2016APV_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA20_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_20_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA25_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_25_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA30_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_30_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA35_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_35_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA40_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_40_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA45_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_45_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA50_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_50_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA55_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_55_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/WH_mA60_16_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2016/WH/WH_mA_60_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/TTGJets_16_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGJets/2016/WH/TTGJets_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2016_vf2/DYJetsToLL_16_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2016/WH/DYJetsToLL_M50_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2016_vf2/TTTo2L2Nu_16_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTTo2L2Nu/2016/WH/TTTo2L2Nu_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2016_vf2/TTToSemiLeptonic_16_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTSL/2016/WH/TTSL_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/SingleElectron_16_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinEle/2016/SingleElectron_2016_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2016_vf2/SingleMuon_16_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinMu/2016/SinMu_2016_%s_%s.root",argv[1],argv[2]);
        }
	
	else if(inputFile=="2017_vf3/WH_mA20_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_20_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA25_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_25_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA30_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_30_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA35_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_35_%s_%s.root",argv[1],argv[2]);
        }
	
	else if(inputFile=="2017_vf3/WH_mA40_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_40_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA45_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_45_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA50_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_50_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA55_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_55_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/WH_mA60_17_v1t.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/signal/2017/WH/WH_mA_60_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/TTGJets_17_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTGJets/2017/WH/TTGJets_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/DYJetsToLL_17_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/DYJetsToLL/2017/WH/DYJetsToLL_M50_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2017_vf3/TTTo2L2Nu_17_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTTo2L2Nu/2017/WH/TTTo2L2Nu_%s_%s.root",argv[1],argv[2]);
        }

        else if(inputFile=="2017_vf3/TTToSemiLeptonic_17_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/TTSL/2017/WH/TTSL_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/SingleMuon_17_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinMu/2017/SinMu_2017_%s_%s.root",argv[1],argv[2]);
        }

	else if(inputFile=="2017_vf3/SingleElectron_17_v1.log"){
                sprintf(fOut,"/home/abala/t3store3/Higgs/SinEle/2017/SingleElectron_2017_%s_%s.root",argv[1],argv[2]);
        }

        else{
                cout<<"Input file does not exist"<<endl;
                exit(0);
        }



	TFile *fSF_singlemuon;
        TH2F *h_singlemuon; TH2F *h_singlemuon_stat;  TH2F *h_singlemuon_syst; 

	TFile *fSF_singleelectron;
	TH2F  *h_singleele;


//	2018
	
	if(!isDATA)
        {
        fSF_singlemuon = TFile::Open("Trigger_SF/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root");
        h_singlemuon      = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt");
        h_singlemuon_stat = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt_stat");
        h_singlemuon_syst = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt_syst");

        fSF_singleelectron  = TFile::Open("Trigger_SF/Efficiencies_electron_HLT_Ele32_WPTight_Gsf_Run2018_UL.root");
        h_singleele    = (TH2F*)fSF_singleelectron->Get("EGamma_SF2D");
        }

	
	//2017
	
/*	if(!isDATA)
        {
        fSF_singlemuon = TFile::Open("Trigger_SF/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root");
        h_singlemuon      = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt");
        h_singlemuon_stat = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt_stat");
        h_singlemuon_syst = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt_syst");

        fSF_singleelectron  = TFile::Open("Trigger_SF/Efficiencies_electron_HLT_Ele32_WPTight_Gsf_Run2017_UL.root");
        h_singleele    = (TH2F*)fSF_singleelectron->Get("EGamma_SF2D");
        }
*/
	//2016
	
/*	if(!isDATA)
        {
        fSF_singlemuon = TFile::Open("Trigger_SF/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root");
        h_singlemuon      = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt");
        h_singlemuon_stat = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt_stat");
        h_singlemuon_syst = (TH2F*)fSF_singlemuon->Get("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_abseta_pt_syst");

        fSF_singleelectron  = TFile::Open("Trigger_SF/Efficiencies_electron_HLT_Ele32_WPTight_Gsf_Run2016_UL.root");
        h_singleele    = (TH2F*)fSF_singleelectron->Get("EGamma_SF2D");
        }
*/

	TFile *fout = new TFile(fOut,"RECREATE");
        TTree *Tout = new TTree("Tout", "Info");


	Tout->Branch("leppt", &leppt, "leppt/F");
        Tout->Branch("lepy", &lepy, "lepy/F");
	Tout->Branch("lepeta", &lepeta, "lepeta/F");
	Tout->Branch("lepphi", &lepphi, "lepphi/F");
	Tout->Branch("lepe", &lepe, "lepe/F");

	Tout->Branch("Met", &Met, "Met/F");
        Tout->Branch("MetPhi", &MetPhi, "MetPhi/F");

        Tout->Branch("jet_no", &jet_no, "jet_no/I");
        Tout->Branch("b1pt", &b1pt, "b1pt/F");
        Tout->Branch("b1y", &b1y, "b1y/F");
        Tout->Branch("b1eta", &b1eta, "b1eta/F");
        Tout->Branch("b1phi", &b1phi, "b1phi/F");
        Tout->Branch("b1e", &b1e, "b1e/F");
        Tout->Branch("b2pt", &b2pt, "b2pt/F");
        Tout->Branch("b2y", &b2y, "b2y/F");
        Tout->Branch("b2eta", &b2eta, "b2eta/F");
        Tout->Branch("b2phi", &b2phi, "b2phi/F");
        Tout->Branch("b2e", &b2e, "b2e/F");
        Tout->Branch("bb_inv_mass", &bb_inv_mass, "bb_inv_mass/F");
	Tout->Branch("b1pt_jecup", &b1pt_jecup, "b1pt_jecup/F");
	Tout->Branch("b1e_jecup", &b1e_jecup, "b1e_jecup/F");
	Tout->Branch("b2pt_jecup", &b2pt_jecup, "b2pt_jecup/F");
        Tout->Branch("b2e_jecup", &b2e_jecup, "b2e_jecup/F");
	Tout->Branch("b1pt_jecdn", &b1pt_jecdn, "b1pt_jecdn/F");
        Tout->Branch("b1e_jecdn", &b1e_jecdn, "b1e_jecdn/F");
        Tout->Branch("b2pt_jecdn", &b2pt_jecdn, "b2pt_jecdn/F");
        Tout->Branch("b2e_jecdn", &b2e_jecdn, "b2e_jecdn/F");
	Tout->Branch("b1pt_resup", &b1pt_resup, "b1pt_resup/F");
        Tout->Branch("b1e_resup", &b1e_resup, "b1e_resup/F");
        Tout->Branch("b2pt_resup", &b2pt_resup, "b2pt_resup/F");
        Tout->Branch("b2e_resup", &b2e_resup, "b2e_resup/F");
	Tout->Branch("b1pt_resdn", &b1pt_resdn, "b1pt_resdn/F");
        Tout->Branch("b1e_resdn", &b1e_resdn, "b1e_resdn/F");
        Tout->Branch("b2pt_resdn", &b2pt_resdn, "b2pt_resdn/F");
        Tout->Branch("b2e_resdn", &b2e_resdn, "b2e_resdn/F");
	Tout->Branch("b1pt_reg", &b1pt_reg, "b1pt_reg/F");
	Tout->Branch("b1e_reg", &b1e_reg, "b1e_reg/F");
	Tout->Branch("b2pt_reg", &b2pt_reg, "b2pt_reg/F");
	Tout->Branch("b2e_reg", &b2e_reg, "b2e_reg/F");
	Tout->Branch("b1pt_scl", &b1pt_scl, "b1pt_scl/F");
        Tout->Branch("b1e_scl", &b1e_scl, "b1e_scl/F");
        Tout->Branch("b2pt_scl", &b2pt_scl, "b2pt_scl/F");
        Tout->Branch("b2e_scl", &b2e_scl, "b2e_scl/F");
	Tout->Branch("b1pt_final", &b1pt_final, "b1pt_final/F");
        Tout->Branch("b1y_final", &b1y_final, "b1y_final/F");
        Tout->Branch("b1eta_final", &b1eta_final, "b1eta_final/F");
        Tout->Branch("b1phi_final", &b1phi_final, "b1phi_final/F");
        Tout->Branch("b1e_final", &b1e_final, "b1e_final/F");
        Tout->Branch("b2pt_final", &b2pt_final, "b2pt_final/F");
        Tout->Branch("b2y_final", &b2y_final, "b2y_final/F");
        Tout->Branch("b2eta_final", &b2eta_final, "b2eta_final/F");
        Tout->Branch("b2phi_final", &b2phi_final, "b2phi_final/F");
        Tout->Branch("b2e_final", &b2e_final, "b2e_final/F");
	Tout->Branch("bb_inv_mass_final", &bb_inv_mass_final, "bb_inv_mass_final/F");
        Tout->Branch("b1_DeepFlv", &b1_DeepFlv, "b1_DeepFlv/F");
        Tout->Branch("b2_DeepFlv", &b2_DeepFlv, "b2_DeepFlv/F");
        Tout->Branch("b1_hadFlv", &b1_hadFlv, "b1_hadFlv/I");
	Tout->Branch("b2_hadFlv", &b2_hadFlv, "b2_hadFlv/I");
	Tout->Branch("b1_partFlv", &b1_partFlv, "b1_partFlv/I");
	Tout->Branch("b2_partFlv", &b2_partFlv, "b2_partFlv/I");
	Tout->Branch("b1_PUid", &b1_PUid, "b1_PUid/F");
	Tout->Branch("b2_PUid", &b2_PUid, "b2_PUid/F");

//        Tout->Branch("pho_no", &pho_no, "pho_no/I");
        Tout->Branch("pho1pt", &pho1pt, "pho1pt/F");
        Tout->Branch("pho2pt", &pho2pt, "pho2pt/F");
        Tout->Branch("pho1y", &pho1y, "pho1y/F");
        Tout->Branch("pho2y", &pho2y, "pho2y/F");
        Tout->Branch("pho1eta", &pho1eta, "pho1eta/F");
        Tout->Branch("pho1phi", &pho1phi, "pho1phi/F");
        Tout->Branch("pho1e", &pho1e, "pho1e/F");
        Tout->Branch("pho2eta", &pho2eta, "pho2eta/F");
        Tout->Branch("pho2phi", &pho2phi, "pho2phi/F");
        Tout->Branch("pho2e", &pho2e, "pho2e/F");
	Tout->Branch("pho1MVA", &pho1MVA, "pho1MVA/F");
	Tout->Branch("pho2MVA", &pho2MVA, "pho2MVA/F");
	Tout->Branch("pho1r9", &pho1r9, "pho1r9/F");
	Tout->Branch("pho2r9", &pho2r9, "pho2r9/F");
        Tout->Branch("dipho_invmass", &dipho_invmass, "dipho_invmass/F");
        Tout->Branch("invmassbbgg", &invmassbbgg, "invmassbbgg/F");
	Tout->Branch("invmassbbgg_final", &invmassbbgg_final, "invmassbbgg_final/F");
	Tout->Branch("evt_rho", &evt_rho, "evt_rho/D");
	Tout->Branch("puvtx", &puvtx, "puvtx/I");
	Tout->Branch("puwt", &puwt, "puwt/D");
	Tout->Branch("pujidwt", &pujidwt, "pujidwt/D");

	Tout->Branch("genb1_pt", &genb1_pt, "genb1_pt/F");
	Tout->Branch("genb1_eta", &genb1_eta, "genb1_eta/F");
	Tout->Branch("genb1_phi", &genb1_phi, "genb1_phi/F");
	Tout->Branch("genb1_mass", &genb1_mass, "genb1_mass/F");
	Tout->Branch("genb2_pt", &genb2_pt, "genb2_pt/F");
        Tout->Branch("genb2_eta", &genb2_eta, "genb2_eta/F");
        Tout->Branch("genb2_phi", &genb2_phi, "genb2_phi/F");
        Tout->Branch("genb2_mass", &genb2_mass, "genb2_mass/F");
	Tout->Branch("genb_invmass", &genb_invmass, "genb_invmass/F");

	Tout->Branch("xsec_weight", &xsec_weight, "xsec_weight/D");
	Tout->Branch("weight_nom", &weight_nom, "weight_nom/D");
	Tout->Branch("weight_PU_up", &weight_PU_up, "weight_PU_up/D");
	Tout->Branch("weight_leptonsf_up", &weight_leptonsf_up, "weight_leptonsf_up/D");
	Tout->Branch("weight_photonsf_up", &weight_photonsf_up, "weight_photonsf_up/D");
	Tout->Branch("weight_prefiring_up", &weight_prefiring_up, "weight_prefiring_up/D");
	Tout->Branch("weight_bSF_up", &weight_bSF_up, "weight_bSF_up/D");
	Tout->Branch("weight_trig_up", &weight_trig_up, "weight_trig_up/D");
	Tout->Branch("weight_haspixsf_up", &weight_haspixsf_up, "weight_haspixsf_up/D");
	Tout->Branch("weight_pujid_up", &weight_pujid_up, "weight_pujid_up/D");
	Tout->Branch("weight_PU_dn", &weight_PU_dn, "weight_PU_dn/D");
	Tout->Branch("weight_leptonsf_dn", &weight_leptonsf_dn, "weight_leptonsf_dn/D");
	Tout->Branch("weight_photonsf_dn", &weight_photonsf_dn, "weight_photonsf_dn/D");
	Tout->Branch("weight_prefiring_dn", &weight_prefiring_dn, "weight_prefiring_dn/D");
	Tout->Branch("weight_bSF_dn", &weight_bSF_dn, "weight_bSF_dn/D");
	Tout->Branch("weight_trig_dn", &weight_trig_dn, "weight_trig_dn/D");
	Tout->Branch("weight_haspixsf_dn", &weight_haspixsf_dn, "weight_haspixsf_dn/D");
	Tout->Branch("weight_pujid_dn", &weight_pujid_dn, "weight_pujid_dn/D");

        Tout->Branch("ab", &ab, "ab/I");
        Tout->Branch("bb", &bb, "bb/I");

	Tout->Branch("isEle", &isEle, "isEle/O");
	Tout->Branch("isMu", &isMu, "isMu/O");

	Tout->Branch("LHEPDF_err", &LHEPDF_err, "LHEPDF_err/F");

	TH1F *cut_flow_el = new TH1F("cut_flow_el","cut_flow_el",20,0,20);
        cut_flow_el->Sumw2();

        TH1F *cut_flow_mu = new TH1F("cut_flow_mu","cut_flow_mu",20,0,20);
        cut_flow_mu->Sumw2();

	TH1F *cut_flow_ph = new TH1F("cut_flow_ph","cut_flow_ph",20,0,20);
        cut_flow_ph->Sumw2();

        Double_t AK4PtBINS[] =  {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000, 7000};
        const int nAK4PtBINS = sizeof(AK4PtBINS)/sizeof(AK4PtBINS[0])-1; 
        Double_t EtaBINS[] = {0,0.6,1.2,2.5}; 
        const int nEtaBINS = sizeof(EtaBINS)/sizeof(EtaBINS[0])-1;

        TH2F* h_Ak4_b_flv = new TH2F("h_Ak4_b_flv", "h_Ak4_b_flv", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv = new TH2F("h_Ak4_c_flv", "h_Ak4_c_flv", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv = new TH2F("h_Ak4_l_flv", "h_Ak4_l_flv", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);

        TH2F* h_Ak4_b_flv_pass_L = new TH2F("h_Ak4_b_flv_pass_L", "h_Ak4_b_flv_pass_L", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv_pass_L = new TH2F("h_Ak4_c_flv_pass_L", "h_Ak4_c_flv_pass_L", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv_pass_L = new TH2F("h_Ak4_l_flv_pass_L", "h_Ak4_l_flv_pass_L", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
       
        TH2F* h_Ak4_b_flv_pass_M = new TH2F("h_Ak4_b_flv_pass_M", "h_Ak4_b_flv_pass_M", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv_pass_M = new TH2F("h_Ak4_c_flv_pass_M", "h_Ak4_c_flv_pass_M", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv_pass_M = new TH2F("h_Ak4_l_flv_pass_M", "h_Ak4_l_flv_pass_M", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);

        TH2F* h_Ak4_b_flv_pass_T = new TH2F("h_Ak4_b_flv_pass_T", "h_Ak4_b_flv_pass_T", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_c_flv_pass_T = new TH2F("h_Ak4_c_flv_pass_T", "h_Ak4_c_flv_pass_T", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);
        TH2F* h_Ak4_l_flv_pass_T = new TH2F("h_Ak4_l_flv_pass_T", "h_Ak4_l_flv_pass_T", nAK4PtBINS, AK4PtBINS, nEtaBINS, EtaBINS);



   int TotEvt=0,count=0;
   long int processedEvt=0;
   string fileName;
   ifstream infile;
   infile.open(path+argv[3]);
   while(!infile.eof()){
   count = count+1;
   getline(infile,fileName);

   int L_lim = stof(argv[1]);
   int H_lim = stof(argv[2]);
   if(count<=L_lim)continue;
   if(count>H_lim)continue;
   TFile *f = TFile::Open(fileName.data());
   if(f==0) continue;



   TTree *T1 = (TTree*)f->Get("Events");

   //2018
   T1->SetBranchAddress("trig_value",&trig_value);
   T1->SetBranchAddress("hlt_IsoMu24",&hlt_IsoMu24);
   T1->SetBranchAddress("hlt_Ele32_WPTight_Gsf",&hlt_Ele32_WPTight_Gsf);
   T1->SetBranchAddress("hlt_Mu17_Mu8",&hlt_Mu17_Mu8);
   T1->SetBranchAddress("hlt_Ele23_Ele12",&hlt_Ele23_Ele12);
   T1->SetBranchAddress("hlt_DoubleEle33",&hlt_DoubleEle33);
   T1->SetBranchAddress("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW); 


   //2017
/*   T1->SetBranchAddress("trig_value",&trig_value);
   T1->SetBranchAddress("hlt_IsoMu24",&hlt_IsoMu24);
   T1->SetBranchAddress("hlt_Ele32_WPTight_Gsf",&hlt_Ele32_WPTight_Gsf);
   T1->SetBranchAddress("hlt_Mu17_Mu8",&hlt_Mu17_Mu8);
   T1->SetBranchAddress("hlt_Ele23_Ele12",&hlt_Ele23_Ele12); 
*/
   //2016
/*   T1->SetBranchAddress("trig_value",&trig_value);
   T1->SetBranchAddress("hlt_IsoMu24",&hlt_IsoMu24);
   T1->SetBranchAddress("hlt_IsoTkMu24",&hlt_IsoTkMu24);
   T1->SetBranchAddress("hlt_Ele27_WPTight_Gsf",&hlt_Ele27_WPTight_Gsf);
   T1->SetBranchAddress("hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   T1->SetBranchAddress("hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   T1->SetBranchAddress("hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",&hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
   T1->SetBranchAddress("hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
   T1->SetBranchAddress("hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   T1->SetBranchAddress("hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   T1->SetBranchAddress("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW); 
*/

   T1->SetBranchAddress("prefiringweight",&prefiringweight);
   T1->SetBranchAddress("prefiringweightup",&prefiringweightup);
   T1->SetBranchAddress("prefiringweightdown",&prefiringweightdown);
   
   T1->SetBranchAddress("Rho",&Rho);

   T1->SetBranchAddress("CHSMET_pt",&miset);
   T1->SetBranchAddress("CHSMET_sumEt",&sumEt);
   T1->SetBranchAddress("CHSMET_phi",&misphi);

   T1->SetBranchAddress("nPFJetAK4",&nPFJetAK4);
   T1->SetBranchAddress("PFJetAK4_pt",PFJetAK4_pt);
   T1->SetBranchAddress("PFJetAK4_eta",PFJetAK4_eta);
   T1->SetBranchAddress("PFJetAK4_y",PFJetAK4_y);
   T1->SetBranchAddress("PFJetAK4_phi",PFJetAK4_phi);
   T1->SetBranchAddress("PFJetAK4_mass",PFJetAK4_mass);
   T1->SetBranchAddress("PFJetAK4_energy",PFJetAK4_energy);
   T1->SetBranchAddress("PFJetAK4_jetID",PFJetAK4_jetID);
   T1->SetBranchAddress("PFJetAK4_jetID_tightlepveto",PFJetAK4_jetID_tightlepveto);
   T1->SetBranchAddress("PFJetAK4_JEC",PFJetAK4_JEC);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV",PFJetAK4_btag_DeepCSV);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav",PFJetAK4_btag_DeepFlav);
   T1->SetBranchAddress("PFJetAK4_JER",PFJetAK4_reso);
   T1->SetBranchAddress("PFJetAK4_JERup",PFJetAK4_resoup);
   T1->SetBranchAddress("PFJetAK4_JERdn",PFJetAK4_resodn);
   T1->SetBranchAddress("PFJetAK4_jesup_Total",PFJetAK4_jesup_Total);
   T1->SetBranchAddress("PFJetAK4_jesdn_Total",PFJetAK4_jesdn_Total);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF",PFJetAK4_btag_DeepCSV_SF);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF_up",PFJetAK4_btag_DeepCSV_SF_up);
   T1->SetBranchAddress("PFJetAK4_btag_DeepCSV_SF_dn",PFJetAK4_btag_DeepCSV_SF_dn);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF",PFJetAK4_btag_DeepFlav_SF);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF_up",PFJetAK4_btag_DeepFlav_SF_up);
   T1->SetBranchAddress("PFJetAK4_btag_DeepFlav_SF_dn",PFJetAK4_btag_DeepFlav_SF_dn);
   T1->SetBranchAddress("PFJetAK4_hadronflav",PFJetAK4_hadronflav);
   T1->SetBranchAddress("PFJetAK4_partonflav",PFJetAK4_partonflav);
   T1->SetBranchAddress("PFJetAK4_PUID",PFJetAK4_PUID);
   T1->SetBranchAddress("PFJetAK4_pt_reg",PFJetAK4_pt_reg);
   T1->SetBranchAddress("PFJetAK4_mass_reg",PFJetAK4_mass_reg);
   T1->SetBranchAddress("PFJetAK4_energy_reg",PFJetAK4_energy_reg);
   T1->SetBranchAddress("PFJetAK4_bcorr",PFJetAK4_bcorr);
   T1->SetBranchAddress("PFJetAK4_breso",PFJetAK4_breso);

  T1->SetBranchAddress("nMuon",&nMuon);
  T1->SetBranchAddress("Muon_isPF",Muon_isPF);
  T1->SetBranchAddress("Muon_isGL",Muon_isGL);
  T1->SetBranchAddress("Muon_isTRK",Muon_isTRK);
  T1->SetBranchAddress("Muon_isLoose",Muon_isLoose);
  T1->SetBranchAddress("Muon_isGoodGL",Muon_isGoodGL);
  T1->SetBranchAddress("Muon_isMed",Muon_isMed);
  T1->SetBranchAddress("Muon_isMedPr",Muon_isMedPr);
  T1->SetBranchAddress("Muon_isTight",Muon_isTight);
  T1->SetBranchAddress("Muon_isHighPt",Muon_isHighPt);
  T1->SetBranchAddress("Muon_isHighPttrk",Muon_isHighPttrk);
  T1->SetBranchAddress("Muon_TightID",Muon_TightID);
  T1->SetBranchAddress("Muon_pt",Muon_pt);
  T1->SetBranchAddress("Muon_p",Muon_p);
  T1->SetBranchAddress("Muon_eta",Muon_eta);
  T1->SetBranchAddress("Muon_phi",Muon_phi);
  T1->SetBranchAddress("Muon_e",Muon_e);
  T1->SetBranchAddress("Muon_minisoch", Muon_minchiso);
  T1->SetBranchAddress("Muon_minisonh", Muon_minnhiso);
  T1->SetBranchAddress("Muon_minisoph", Muon_minphiso);
  T1->SetBranchAddress("Muon_minisoall", Muon_minisoall);
  T1->SetBranchAddress("Muon_dxy",Muon_dxy);
  T1->SetBranchAddress("Muon_dz",Muon_dz);
  T1->SetBranchAddress("Muon_dxyErr",Muon_dxyErr);
  T1->SetBranchAddress("Muon_ip3d",Muon_ip3d);
  T1->SetBranchAddress("Muon_ptErr",Muon_ptErr);
  T1->SetBranchAddress("Muon_chi",Muon_chi);
  T1->SetBranchAddress("Muon_ndf",Muon_ndf);
  T1->SetBranchAddress("Muon_ecal",Muon_ecal);
  T1->SetBranchAddress("Muon_hcal",Muon_hcal);
  T1->SetBranchAddress("Muon_pfiso",Muon_pfiso);
  T1->SetBranchAddress("Muon_posmatch",Muon_posmatch);
  T1->SetBranchAddress("Muon_trkink",Muon_trkink);
  T1->SetBranchAddress("Muon_segcom",Muon_segcom);
  T1->SetBranchAddress("Muon_hit",Muon_hit);
  T1->SetBranchAddress("Muon_pixhit",Muon_pixhit);
  T1->SetBranchAddress("Muon_mst",Muon_mst);
  T1->SetBranchAddress("Muon_trklay",Muon_trklay);
  T1->SetBranchAddress("Muon_valfrac",Muon_valfrac);
  T1->SetBranchAddress("Muon_dxy_sv",Muon_dxy_sv);
  T1->SetBranchAddress("Muon_corrected_pt",Muon_corrected_pt);
  T1->SetBranchAddress("Muon_correctedUp_pt",Muon_correctedUp_pt);
  T1->SetBranchAddress("Muon_correctedDown_pt",Muon_correctedDown_pt);

  T1->SetBranchAddress("nElectron",&nElectron);
  T1->SetBranchAddress("Electron_pt",Electron_pt);
  T1->SetBranchAddress("Electron_eta",Electron_eta);
  T1->SetBranchAddress("Electron_phi",Electron_phi);
  T1->SetBranchAddress("Electron_p",Electron_p);
  T1->SetBranchAddress("Electron_e",Electron_e);
  T1->SetBranchAddress("Electron_e_ECAL",Electron_e_ECAL);
  T1->SetBranchAddress("Electron_ecalEnergyPreCorr",Electron_ecalEnergyPreCorr);
  T1->SetBranchAddress("Electron_ecalEnergyPostCorr",Electron_ecalEnergyPostCorr);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP90",Electron_mvaid_Fallv2WP90);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP90_noIso",Electron_mvaid_Fallv2WP90_noIso);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP80",Electron_mvaid_Fallv2WP80);
  T1->SetBranchAddress("Electron_mvaid_Fallv2WP80_noIso",Electron_mvaid_Fallv2WP80_noIso);
  T1->SetBranchAddress("Electron_dxy",Electron_dxy);
  T1->SetBranchAddress("Electron_dxyErr",Electron_dxyErr);
  T1->SetBranchAddress("Electron_dz",Electron_dz);
  T1->SetBranchAddress("Electron_dzErr",Electron_dzErr);
  T1->SetBranchAddress("Electron_ip3d",Electron_ip3d);
  T1->SetBranchAddress("Electron_dxy_sv",Electron_dxy_sv);
  T1->SetBranchAddress("Electron_hovere",Electron_hovere);
  T1->SetBranchAddress("Electron_chi",Electron_chi);
  T1->SetBranchAddress("Electron_ndf",Electron_ndf);
  T1->SetBranchAddress("Electron_eoverp",Electron_eoverp);
  T1->SetBranchAddress("Electron_ietaieta",Electron_ietaieta);
  T1->SetBranchAddress("Electron_misshits",Electron_misshits);
  T1->SetBranchAddress("Electron_pfiso_drcor",Electron_pfiso_drcor);
  T1->SetBranchAddress("Electron_pfiso_eacor",Electron_pfiso_eacor);
  T1->SetBranchAddress("Electron_pfiso04_eacor",Electron_pfiso04_eacor);
  T1->SetBranchAddress("Electron_r9full", Electron_r9full);
  T1->SetBranchAddress("Electron_hcaloverecal", Electron_hcaloverecal);
  T1->SetBranchAddress("Electron_hitsmiss", Electron_hitsmiss);
  T1->SetBranchAddress("Electron_ecloverpout", Electron_ecloverpout);
  T1->SetBranchAddress("Electron_convVeto", Electron_convVeto);

  T1->SetBranchAddress("nPhoton",&nPhoton);
  T1->SetBranchAddress("Photon_passEveto",Photon_passEveto);
  T1->SetBranchAddress("Photon_PixelSeed",Photon_PixelSeed);
  T1->SetBranchAddress("Photon_e",Photon_e);
  T1->SetBranchAddress("Photon_pt",Photon_pt);
  T1->SetBranchAddress("Photon_eta",Photon_eta);
  T1->SetBranchAddress("Photon_phi",Photon_phi);
  T1->SetBranchAddress("Photon_ecalEnergyPreCorr",Photon_ecalEnergyPreCorr);
  T1->SetBranchAddress("Photon_ecalEnergyPostCorr",Photon_ecalEnergyPostCorr);
  T1->SetBranchAddress("Photon_mvaid_Fall17V2_raw",Photon_mvaid_Fall17V2_raw);
  T1->SetBranchAddress("Photon_mvaid_Fall17V2_WP90",Photon_mvaid_Fall17V2_WP90);
  T1->SetBranchAddress("Photon_mvaid_Fall17V2_WP80",Photon_mvaid_Fall17V2_WP80);
  T1->SetBranchAddress("Photon_mvaid_Spring16V1_WP90",Photon_mvaid_Spring16V1_WP90);
  T1->SetBranchAddress("Photon_mvaid_Spring16V1_WP80",Photon_mvaid_Spring16V1_WP80);
  T1->SetBranchAddress("Photon_e1by9",Photon_e1by9);
  T1->SetBranchAddress("Photon_e9by25",Photon_e9by25);
  T1->SetBranchAddress("Photon_trkiso",Photon_trkiso);
  T1->SetBranchAddress("Photon_emiso",Photon_emiso);
  T1->SetBranchAddress("Photon_hadiso",Photon_hadiso);
  T1->SetBranchAddress("Photon_chhadiso",Photon_chhadiso);
  T1->SetBranchAddress("Photon_neuhadiso",Photon_neuhadiso);
  T1->SetBranchAddress("Photon_phoiso",Photon_phoiso);
  T1->SetBranchAddress("Photon_PUiso",Photon_PUiso);
  T1->SetBranchAddress("Photon_hadbyem",Photon_hadbyem);
  T1->SetBranchAddress("Photon_ietaieta",Photon_ietaieta);

  if(isMC){

  // generator-related info //

  T1->SetBranchAddress("Generator_weight", &Generator_weight);
  T1->SetBranchAddress("Generator_qscale",&Generator_qscale);
  T1->SetBranchAddress("Generator_x1",&Generator_x1);
  T1->SetBranchAddress("Generator_x2",&Generator_x2);
  T1->SetBranchAddress("Generator_xpdf1",&Generator_xpdf1);
  T1->SetBranchAddress("Generator_xpdf2",&Generator_xpdf2);
  T1->SetBranchAddress("Generator_id1",&Generator_id1);
  T1->SetBranchAddress("Generator_id2",&Generator_id2);
  T1->SetBranchAddress("Generator_scalePDF",&Generator_scalePDF);

  T1->SetBranchAddress("nGenJetAK4", &nGenJetAK4);
  T1->SetBranchAddress("GenJetAK4_pt",GenJetAK4_pt);
  T1->SetBranchAddress("GenJetAK4_eta",GenJetAK4_eta);
  T1->SetBranchAddress("GenJetAK4_phi",GenJetAK4_phi);
  T1->SetBranchAddress("GenJetAK4_mass",GenJetAK4_mass);
  T1->SetBranchAddress("GenJetAK4_hadronflav",GenJetAK4_hadronflav);

  T1->SetBranchAddress("nGenPart",&nGenPart);
  T1->SetBranchAddress("GenPart_pt",GenPart_pt);
  T1->SetBranchAddress("GenPart_eta",GenPart_eta);
  T1->SetBranchAddress("GenPart_phi",GenPart_phi);
  T1->SetBranchAddress("GenPart_mass",GenPart_mass);
  T1->SetBranchAddress("GenPart_status",GenPart_status);
  T1->SetBranchAddress("GenPart_pdgId",GenPart_pdg);
  T1->SetBranchAddress("GenPart_mompdgId",GenPart_mompdg);
  T1->SetBranchAddress("GenPart_fromhard",GenPart_fromhard);

  T1->SetBranchAddress("npu_vert",&npu_vert);
  T1->SetBranchAddress("npu_vert_true",&npu_vert_true);

  T1->SetBranchAddress("LHE_weight",&LHE_weight);
  T1->SetBranchAddress("BTAG_SF",&BTAG_SF);
  T1->SetBranchAddress("BTAG_jes_up",&BTAG_jes_up);
  T1->SetBranchAddress("BTAG_jes_dn",&BTAG_jes_dn);

  T1->SetBranchAddress("nLHEPDFWeights",&nLHEPDFWeights);
  T1->SetBranchAddress("LHEPDFWeights",LHEPDFWeights);
  }


  int nevents;
  nevents=T1->GetEntries();

  for(int iev=0; iev<nevents; iev++)
  {
     
	  T1->GetEntry(iev);

          cut_flow_el->Fill(0);
	  cut_flow_mu->Fill(0);


        //-------------------------------------------event weight info-------------------------------------------//


	  double event_weight;

	      if(isMC){

                if (isFastSIM)
                {
                        event_weight = 1.0;
                }

		if(fabs(LHE_weight)>1.e-12) { event_weight = LHE_weight; }
		else { event_weight = Generator_weight;	}
        	      }
	      else{
		event_weight = 1.;
	          }


          //-----------------------------------------------------AK4 histogram filling for btag SF------------------------------------------------//
            
            
              if(isMC){
		for(unsigned ijet=0; ijet<nPFJetAK4; ijet++){
			double dR = 9999.9;
			for(unsigned gjet=0; gjet<nGenJetAK4; gjet++)
			{
				double temp_dR = delta2R(PFJetAK4_eta[ijet],PFJetAK4_phi[ijet],GenJetAK4_eta[gjet],GenJetAK4_phi[gjet]) ;
				if (temp_dR < dR )
				{
					dR = temp_dR;
				}
			}
			if(dR < 0.4)
			{
				if( abs(PFJetAK4_hadronflav[ijet]) == 5 )  {  
					h_Ak4_b_flv->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); 
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_T  )  { h_Ak4_b_flv_pass_T->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); } 
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_M ) 	{ h_Ak4_b_flv_pass_M->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); } 
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_L ) 	{ h_Ak4_b_flv_pass_L->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
				}    
				else if( abs(PFJetAK4_hadronflav[ijet]) == 4 )  {  
					h_Ak4_c_flv->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight);
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_T  )  { h_Ak4_c_flv_pass_T->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_M ) 	{ h_Ak4_c_flv_pass_M->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_L ) 	{ h_Ak4_c_flv_pass_L->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
				}
				else if( abs(PFJetAK4_hadronflav[ijet]) == 0 )  {  
					h_Ak4_l_flv->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight);
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_T  )  { h_Ak4_l_flv_pass_T->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_M ) 	{ h_Ak4_l_flv_pass_M->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
					if (PFJetAK4_btag_DeepFlav[ijet]   > DAK4_L ) 	{ h_Ak4_l_flv_pass_L->Fill(PFJetAK4_pt[ijet],fabs(PFJetAK4_eta[ijet]),event_weight); }
				}
			}		                                   
		}      
                }

       

	  //---------------------------------------trigger info---------------------------------------//


	  if (!(hlt_Ele32_WPTight_Gsf || hlt_IsoMu24)) continue;          
//	  if (!(hlt_Ele27_WPTight_Gsf || hlt_IsoMu24 || hlt_IsoTkMu24)) continue;          
          

          //------------------------------------------------lepton info--------------------------------------------//	      


	  int nele = 0;
          int nmu = 0;

          int EleIndex = -10;
          int MuIndex = -10;

          TLorentzVector leptonvec;

	  isEle = false;

	  if (hlt_Ele32_WPTight_Gsf) {
//	  if (hlt_Ele27_WPTight_Gsf) {
	  cut_flow_el->Fill(1);

	  for(int i=0; i<nElectron; i++)
	  {
		  if (Electron_hitsmiss[i]>1) continue;
                  if (fabs(Electron_dxy[i])>0.045) continue;
                  if (Electron_pfiso_drcor[i]>=0.15) continue;
                  if(fabs(Electron_eta[i]) > 2.5) continue;                                    //electron channel, for 2016 28, for all other years 34;
                  if(Electron_pt[i]<34) continue;
                  EleIndex = i;
                  nele++;
	  }

	  if (nele != 1) continue;

	  cut_flow_el->Fill(2);

	  int nlooseE_el = 0;
          for(int i=0; i<nElectron; i++)
          {
                  if (i==EleIndex) continue;
                  if (Electron_hitsmiss[i]>1) continue;
                  if (fabs(Electron_dxy[i])>0.045) continue;
                  if (Electron_pfiso_drcor[i]>=0.2) continue;
                  if(fabs(Electron_eta[i]) > 2.5) continue;
                  if(Electron_pt[i]<15) continue;
                  nlooseE_el++;
          }

          if (nlooseE_el > 0) continue;
          cut_flow_el->Fill(3);

          int nlooseMu_el = 0;
          for (int i=0; i<nMuon; i++)
          {
                  if (!Muon_isPF[i]) continue;
                  if (Muon_pfiso[i]>=0.2) continue;
                  if(fabs(Muon_eta[i]) > 2.5) continue;
                  if(Muon_pt[i]<15) continue;
                  nlooseMu_el++;
          }

          if (nlooseMu_el > 0) continue;
          cut_flow_el->Fill(4);

	  isEle = true;      }


	  isMu = false;

	  if (hlt_IsoMu24) {
//	  if (hlt_IsoMu24 || hlt_IsoTkMu24) {
	  cut_flow_mu->Fill(1);

	  for (int i=0; i<nMuon; i++)
          {
                  if (!Muon_isPF[i]) continue;
                  if (Muon_pfiso[i]>=0.15) continue;
                  if(fabs(Muon_eta[i]) > 2.5) continue;
                  if(Muon_pt[i]<26) continue;
                  MuIndex = i;                                                                  //muon channel
                  nmu++;
          }

	  if (nmu != 1) continue;

          cut_flow_mu->Fill(2);

	  int nlooseE_mu = 0;
          for(int i=0; i<nElectron; i++)
          {
                  if (Electron_hitsmiss[i]>1) continue;
                  if (fabs(Electron_dxy[i])>0.045) continue;
                  if (Electron_pfiso_drcor[i]>=0.2) continue;
                  if(fabs(Electron_eta[i]) > 2.5) continue;
                  if(Electron_pt[i]<15) continue;
                  nlooseE_mu++;
          }

          if (nlooseE_mu > 0) continue;
          cut_flow_mu->Fill(3);

          int nlooseMu_mu = 0;
          for (int i=0; i<nMuon; i++)
          {
                  if (i==MuIndex) continue;
                  if (!Muon_isPF[i]) continue;
                  if (Muon_pfiso[i]>=0.2) continue;
                  if(fabs(Muon_eta[i]) > 2.5) continue;
                  if(Muon_pt[i]<15) continue;
                  nlooseMu_mu++;
          }

          if (nlooseMu_mu > 0) continue;
          cut_flow_mu->Fill(4);

	  isMu = true;     }


	  if (isEle) {
          leptonvec.SetPtEtaPhiE(Electron_pt[EleIndex]*(Electron_ecalEnergyPostCorr[EleIndex]/Electron_ecalEnergyPreCorr[EleIndex]), Electron_eta[EleIndex], Electron_phi[EleIndex], Electron_ecalEnergyPostCorr[EleIndex]); }

	  if (isMu) {
          leptonvec.SetPtEtaPhiE(Muon_corrected_pt[MuIndex], Muon_eta[MuIndex], Muon_phi[MuIndex], Muon_e[MuIndex]*(Muon_corrected_pt[MuIndex]/Muon_pt[MuIndex])); }
     


	  //-------------------------------------------------------AK4 jet info-------------------------------------------------//


	  int indx = 0;

	  for(int i=0; i<nPFJetAK4; i++)
	  {
                  Jet_pt_nom[i] = PFJetAK4_pt[i];
                  Jet_mass_nom[i] = PFJetAK4_mass[i];

                  Jet_pt_jesup[i] = PFJetAK4_pt[i];
                  Jet_mass_jesup[i] = PFJetAK4_mass[i];

                  Jet_pt_jesdn[i] = PFJetAK4_pt[i];
                  Jet_mass_jesdn[i] = PFJetAK4_mass[i];

                  Jet_pt_resoup[i] = PFJetAK4_pt[i];
                  Jet_mass_resoup[i] = PFJetAK4_mass[i];

                  Jet_pt_resodn[i] = PFJetAK4_pt[i];
                  Jet_mass_resodn[i] = PFJetAK4_mass[i];


                  Jet_pt_nom[i] *= PFJetAK4_JEC[i];
                  Jet_mass_nom[i] *= PFJetAK4_JEC[i];
                  Jet_pt_nom[i] *= (1.+PFJetAK4_reso[i]);
                  Jet_mass_nom[i] *= (1.+PFJetAK4_reso[i]);

                  Jet_pt_jesup[i] *= PFJetAK4_jesup_Total[i];
                  Jet_mass_jesup[i] *= PFJetAK4_jesup_Total[i];
                  Jet_pt_jesup[i] *= (1.+PFJetAK4_reso[i]);
                  Jet_mass_jesup[i] *= (1.+PFJetAK4_reso[i]);

                  Jet_pt_jesdn[i] *= PFJetAK4_jesdn_Total[i];
                  Jet_mass_jesdn[i] *= PFJetAK4_jesdn_Total[i];
                  Jet_pt_jesdn[i] *= (1.+PFJetAK4_reso[i]);
                  Jet_mass_jesdn[i] *= (1.+PFJetAK4_reso[i]);

                  Jet_pt_resoup[i] *= PFJetAK4_JEC[i];
                  Jet_mass_resoup[i] *= PFJetAK4_JEC[i];
                  Jet_pt_resoup[i] *= (1.+PFJetAK4_resoup[i]);
                  Jet_mass_resoup[i] *= (1.+PFJetAK4_resoup[i]);

                  Jet_pt_resodn[i] *= PFJetAK4_JEC[i];
                  Jet_mass_resodn[i] *= PFJetAK4_JEC[i];
                  Jet_pt_resodn[i] *= (1.+PFJetAK4_resodn[i]);
                  Jet_mass_resodn[i] *= (1.+PFJetAK4_resodn[i]);


                  if(!PFJetAK4_jetID_tightlepveto[i]) continue;
                  if(fabs(PFJetAK4_y[i]) > 2.5) continue;

		  //2017 & 18
		  if (Jet_pt_nom[i] >= 10 && Jet_pt_nom[i] < 20)
                  {
                        if (PFJetAK4_PUID[i] < 0.26) continue;
                  }

		  if (Jet_pt_nom[i] >= 20 && Jet_pt_nom[i] < 30)
                  {
                        if (PFJetAK4_PUID[i] < 0.68) continue;
                  }

                  if (Jet_pt_nom[i] >= 30 && Jet_pt_nom[i] < 40)
                  {
                        if (PFJetAK4_PUID[i] < 0.90) continue;
                  }

                  if (Jet_pt_nom[i] >= 40 && Jet_pt_nom[i] < 50)
                  {
                        if (PFJetAK4_PUID[i] < 0.96) continue;
                  } 


		  //2016
/*		  if (Jet_pt_nom[i] >= 10 && Jet_pt_nom[i] < 20)
                  {
                        if (PFJetAK4_PUID[i] < 0.2) continue;
                  }

                  if (Jet_pt_nom[i] >= 20 && Jet_pt_nom[i] < 30)
                  {
                        if (PFJetAK4_PUID[i] < 0.62) continue;
                  }

                  if (Jet_pt_nom[i] >= 30 && Jet_pt_nom[i] < 40)
                  {
                        if (PFJetAK4_PUID[i] < 0.86) continue;
                  }

                  if (Jet_pt_nom[i] >= 40 && Jet_pt_nom[i] < 50)
                  {
                        if (PFJetAK4_PUID[i] < 0.93) continue;
                  } 
*/	

                  TLorentzVector b_vec;
                  b_vec.SetPtEtaPhiM(Jet_pt_nom[i], PFJetAK4_eta[i], PFJetAK4_phi[i], Jet_mass_nom[i]);
                        
                  if (leptonvec.DeltaR(b_vec) < 0.4) continue; 

                  if(isMC){
		
			BTagEntry::JetFlavor btv_flav;
			if(abs(PFJetAK4_hadronflav[i])==5){ btv_flav = BTagEntry::FLAV_B; }
			else if (abs(PFJetAK4_hadronflav[i])==4){ btv_flav = BTagEntry::FLAV_C; }
			else { btv_flav = BTagEntry::FLAV_UDSG; }
		
			PFJetAK4_btag_DeepFlav_SF[i] = reader_deepflav.eval_auto_bounds("central",btv_flav,fabs(PFJetAK4_y[i]),Jet_pt_nom[i]); 
			PFJetAK4_btag_DeepFlav_SF_up[i] = reader_deepflav.eval_auto_bounds("up",btv_flav,fabs(PFJetAK4_y[i]),Jet_pt_nom[i]);
			PFJetAK4_btag_DeepFlav_SF_dn[i] = reader_deepflav.eval_auto_bounds("down",btv_flav,fabs(PFJetAK4_y[i]),Jet_pt_nom[i]);
		
	 	  }
		
		  else{
			
			PFJetAK4_btag_DeepFlav_SF[i] = PFJetAK4_btag_DeepFlav_SF_up[i] = PFJetAK4_btag_DeepFlav_SF_dn[i] = 0;
			
			}

                  Jet_pt_nom[indx] = Jet_pt_nom[i];
                  Jet_mass_nom[indx] = Jet_mass_nom[i];
                  Jet_pt_jesup[indx] = Jet_pt_jesup[i];
                  Jet_mass_jesup[indx] = Jet_mass_jesup[i];
                  Jet_pt_jesdn[indx] = Jet_pt_jesdn[i];
                  Jet_mass_jesdn[indx] = Jet_mass_jesdn[i];
                  Jet_pt_resoup[indx] = Jet_pt_resoup[i];
                  Jet_mass_resoup[indx] = Jet_mass_resoup[i];
                  Jet_pt_resodn[indx] = Jet_pt_resodn[i];
                  Jet_mass_resodn[indx] = Jet_mass_resodn[i];

                  PFJetAK4_y[indx] = PFJetAK4_y[i];
                  PFJetAK4_eta[indx] = PFJetAK4_eta[i];
                  PFJetAK4_phi[indx] = PFJetAK4_phi[i];
                  PFJetAK4_btag_DeepFlav[indx] = PFJetAK4_btag_DeepFlav[i];
                  PFJetAK4_btag_DeepFlav_SF[indx] = PFJetAK4_btag_DeepFlav_SF[i];
                  PFJetAK4_btag_DeepFlav_SF_up[indx] = PFJetAK4_btag_DeepFlav_SF_up[i];
                  PFJetAK4_btag_DeepFlav_SF_dn[indx] = PFJetAK4_btag_DeepFlav_SF_dn[i];
		  PFJetAK4_hadronflav[indx] = PFJetAK4_hadronflav[i];
		  PFJetAK4_partonflav[indx] = PFJetAK4_partonflav[i];
		  PFJetAK4_PUID[indx] = PFJetAK4_PUID[i];
		  PFJetAK4_bcorr[indx] = PFJetAK4_bcorr[i];
		  PFJetAK4_breso[indx] = PFJetAK4_breso[i];

                  if (++indx >= njetmx) break;
	  }


	  nPFJetAK4 = indx;
	  int AK4coll = indx;

          if (nPFJetAK4 < 2) continue;

	  if (isEle)
	  {
		  cut_flow_el->Fill(5);
	  }
	  if (isMu)
          {
                  cut_flow_mu->Fill(5);
          }


	  for (int i=0; i<nPFJetAK4; i++)
          {
                  for (int j=i+1; j<nPFJetAK4; j++)
                  {
                           if(PFJetAK4_btag_DeepFlav[i]<PFJetAK4_btag_DeepFlav[j])
                           {
                                    float a = PFJetAK4_btag_DeepFlav[i];
                                    PFJetAK4_btag_DeepFlav[i] = PFJetAK4_btag_DeepFlav[j];         
                                    PFJetAK4_btag_DeepFlav[j] = a;

                                    a = Jet_pt_nom[i];
                                    Jet_pt_nom[i] = Jet_pt_nom[j];                   //nominal
                                    Jet_pt_nom[j] = a;

                                    a = PFJetAK4_y[i];
                                    PFJetAK4_y[i] = PFJetAK4_y[j];
                                    PFJetAK4_y[j] = a;

                                    a = PFJetAK4_eta[i];
                                    PFJetAK4_eta[i] = PFJetAK4_eta[j];
                                    PFJetAK4_eta[j] = a;

                                    a = PFJetAK4_phi[i];
                                    PFJetAK4_phi[i] = PFJetAK4_phi[j];
                                    PFJetAK4_phi[j] = a;

                                    a = Jet_mass_nom[i];
                                    Jet_mass_nom[i] = Jet_mass_nom[j];               //nominal
                                    Jet_mass_nom[j] = a;

                                    a = PFJetAK4_btag_DeepFlav_SF[i];
                                    PFJetAK4_btag_DeepFlav_SF[i] = PFJetAK4_btag_DeepFlav_SF[j];
                                    PFJetAK4_btag_DeepFlav_SF[j] = a;

                                    a = PFJetAK4_btag_DeepFlav_SF_up[i];
                                    PFJetAK4_btag_DeepFlav_SF_up[i] = PFJetAK4_btag_DeepFlav_SF_up[j];
                                    PFJetAK4_btag_DeepFlav_SF_up[j] = a;

                                    a = PFJetAK4_btag_DeepFlav_SF_dn[i];
                                    PFJetAK4_btag_DeepFlav_SF_dn[i] = PFJetAK4_btag_DeepFlav_SF_dn[j];
                                    PFJetAK4_btag_DeepFlav_SF_dn[j] = a;

				    a = PFJetAK4_hadronflav[i];
                                    PFJetAK4_hadronflav[i] = PFJetAK4_hadronflav[j];
                                    PFJetAK4_hadronflav[j] = a;

				    a = PFJetAK4_partonflav[i];
                                    PFJetAK4_partonflav[i] = PFJetAK4_partonflav[j];
                                    PFJetAK4_partonflav[j] = a;

				    a = PFJetAK4_PUID[i];
                                    PFJetAK4_PUID[i] = PFJetAK4_PUID[j];
                                    PFJetAK4_PUID[j] = a;

				    a = PFJetAK4_bcorr[i];
                                    PFJetAK4_bcorr[i] = PFJetAK4_bcorr[j];
                                    PFJetAK4_bcorr[j] = a;

				    a = PFJetAK4_breso[i];
                                    PFJetAK4_breso[i] = PFJetAK4_breso[j];
                                    PFJetAK4_breso[j] = a;

                                    a = Jet_pt_jesup[i];
                                    Jet_pt_jesup[i] = Jet_pt_jesup[j];                    //jesup
                                    Jet_pt_jesup[j] = a;

                                    a = Jet_mass_jesup[i];
                                    Jet_mass_jesup[i] = Jet_mass_jesup[j];                //jesup
                                    Jet_mass_jesup[j] = a;

                                    a = Jet_pt_jesdn[i];
                                    Jet_pt_jesdn[i] = Jet_pt_jesdn[j];                 //jesdn
                                    Jet_pt_jesdn[j] = a;

                                    a = Jet_mass_jesdn[i];
                                    Jet_mass_jesdn[i] = Jet_mass_jesdn[j];               //jesdn
                                    Jet_mass_jesdn[j] = a;

                                    a = Jet_pt_resoup[i];
                                    Jet_pt_resoup[i] = Jet_pt_resoup[j];                    //resoup
                                    Jet_pt_resoup[j] = a;

                                    a = Jet_mass_resoup[i];
                                    Jet_mass_resoup[i] = Jet_mass_resoup[j];                 //resoup
                                    Jet_mass_resoup[j] = a;

                                    a = Jet_pt_resodn[i];
                                    Jet_pt_resodn[i] = Jet_pt_resodn[j];                  //resodn
                                    Jet_pt_resodn[j] = a;

                                    a = Jet_mass_resodn[i];
                                    Jet_mass_resodn[i] = Jet_mass_resodn[j];                      //resodn
                                    Jet_mass_resodn[j] = a;
                           }
                   }
           }
 


	  TLorentzVector leadbvec_nom, subleadbvec_nom, leadbvec_jecup, subleadbvec_jecup, leadbvec_jecdn, subleadbvec_jecdn, leadbvec_resup, subleadbvec_resup, leadbvec_resdn, subleadbvec_resdn;
          TLorentzVector lbvec_reg, slbvec_reg;

          bb = 0;
          ab = 0;
          aa = 0;
    
          if (Jet_pt_nom[0] >= 10 && Jet_pt_nom[1] >= 10 && (PFJetAK4_btag_DeepFlav[0] >= 0.2783 && PFJetAK4_btag_DeepFlav[1] >= 0.2783)) {bb = 1;}
          if (Jet_pt_nom[0] >= 10 && Jet_pt_nom[1] >= 10 && (PFJetAK4_btag_DeepFlav[0] >= 0.2783 && PFJetAK4_btag_DeepFlav[1] > 0.0490 && PFJetAK4_btag_DeepFlav[1] < 0.2783)) {ab = 1;}
          if (Jet_pt_nom[0] >= 10 && Jet_pt_nom[1] >= 10 && (PFJetAK4_btag_DeepFlav[0] < 0.2783 && PFJetAK4_btag_DeepFlav[1] < 0.2783)) {aa = 1;}


	  if (isEle)
	  {
		  if (bb == 1) {cut_flow_el->Fill(6);}
          	  if (ab == 1) {cut_flow_el->Fill(7);}
	  }

	  if (isMu)
          {
                  if (bb == 1) {cut_flow_mu->Fill(6);}
                  if (ab == 1) {cut_flow_mu->Fill(7);}
          }


          if (Jet_pt_nom[0] >= 10 && Jet_pt_nom[1] >= 10 && (PFJetAK4_btag_DeepFlav[0] >= 0.2783 && PFJetAK4_btag_DeepFlav[1] >= 0.2783)) 
	  {  
          if (Jet_pt_nom[0] >= Jet_pt_nom[1])
          {
                   leadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
                   subleadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);

		   leadbvec_jecup.SetPtEtaPhiM(Jet_pt_jesup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesup[0]);
                   subleadbvec_jecup.SetPtEtaPhiM(Jet_pt_jesup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesup[1]);

		   leadbvec_jecdn.SetPtEtaPhiM(Jet_pt_jesdn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesdn[0]);
                   subleadbvec_jecdn.SetPtEtaPhiM(Jet_pt_jesdn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesdn[1]);

		   leadbvec_resup.SetPtEtaPhiM(Jet_pt_resoup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resoup[0]);
                   subleadbvec_resup.SetPtEtaPhiM(Jet_pt_resoup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resoup[1]);

                   leadbvec_resdn.SetPtEtaPhiM(Jet_pt_resodn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resodn[0]);
                   subleadbvec_resdn.SetPtEtaPhiM(Jet_pt_resodn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resodn[1]);

		   lbvec_reg.SetPtEtaPhiM(Jet_pt_nom[0]*PFJetAK4_bcorr[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]*PFJetAK4_bcorr[0]);
		   slbvec_reg.SetPtEtaPhiM(Jet_pt_nom[1]*PFJetAK4_bcorr[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]*PFJetAK4_bcorr[1]);
          }                                                             

          else
          {
                   subleadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
                   leadbvec_nom.SetPtEtaPhiM(Jet_pt_nom[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);

		   leadbvec_jecup.SetPtEtaPhiM(Jet_pt_jesup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesup[1]);
                   subleadbvec_jecup.SetPtEtaPhiM(Jet_pt_jesup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesup[0]);

                   leadbvec_jecdn.SetPtEtaPhiM(Jet_pt_jesdn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_jesdn[1]);
                   subleadbvec_jecdn.SetPtEtaPhiM(Jet_pt_jesdn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_jesdn[0]);

                   leadbvec_resup.SetPtEtaPhiM(Jet_pt_resoup[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resoup[1]);
                   subleadbvec_resup.SetPtEtaPhiM(Jet_pt_resoup[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resoup[0]);

                   leadbvec_resdn.SetPtEtaPhiM(Jet_pt_resodn[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_resodn[1]);
                   subleadbvec_resdn.SetPtEtaPhiM(Jet_pt_resodn[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_resodn[0]);

		   lbvec_reg.SetPtEtaPhiM(Jet_pt_nom[1]*PFJetAK4_bcorr[1], PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]*PFJetAK4_bcorr[1]);
                   slbvec_reg.SetPtEtaPhiM(Jet_pt_nom[0]*PFJetAK4_bcorr[0], PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]*PFJetAK4_bcorr[0]);
          }
	  }


	  //2018
	  float scl_nom = -0.018;
	  float scl_up = 0.001;
	  float scl_dn = -0.037;
	  float smr_nom = 0.05;
          float smr_up = 0.129;
          float smr_dn = -0.029; 


	  //2017
/*	  float scl_nom = 0.011;
          float scl_up = 0.033;
          float scl_dn = -0.011;
          float smr_nom = 0.051;
          float smr_up = 0.119;
          float smr_dn = -0.017; 
*/
	  //2016
/*	  float scl_nom = 0.004;
          float scl_up = 0.022;
          float scl_dn = -0.014;
          float smr_nom = -0.044;
          float smr_up = 0.017;
          float smr_dn = -0.105; 
*/


	  TLorentzVector reg_scale_nom_lead, reg_scale_nom_sublead, reg_scale_up_lead, reg_scale_up_sublead, reg_scale_dn_lead, reg_scale_dn_sublead;

	  reg_scale_nom_lead.SetPtEtaPhiE(lbvec_reg.Pt()*(1+scl_nom), lbvec_reg.Eta(), lbvec_reg.Phi(), lbvec_reg.E()*(1+scl_nom));
	  reg_scale_nom_sublead.SetPtEtaPhiE(slbvec_reg.Pt()*(1+scl_nom), slbvec_reg.Eta(), slbvec_reg.Phi(), slbvec_reg.E()*(1+scl_nom));

	  reg_scale_up_lead.SetPtEtaPhiE(lbvec_reg.Pt()*(1+scl_up), lbvec_reg.Eta(), lbvec_reg.Phi(), lbvec_reg.E()*(1+scl_up));
          reg_scale_up_sublead.SetPtEtaPhiE(slbvec_reg.Pt()*(1+scl_up), slbvec_reg.Eta(), slbvec_reg.Phi(), slbvec_reg.E()*(1+scl_up));

	  reg_scale_dn_lead.SetPtEtaPhiE(lbvec_reg.Pt()*(1+scl_dn), lbvec_reg.Eta(), lbvec_reg.Phi(), lbvec_reg.E()*(1+scl_dn));
          reg_scale_dn_sublead.SetPtEtaPhiE(slbvec_reg.Pt()*(1+scl_dn), slbvec_reg.Eta(), slbvec_reg.Phi(), slbvec_reg.E()*(1+scl_dn));


	  float mindr1, mindr2;
	  float dr1[njetmx], dr2[njetmx];
	  for(int i=0; i<nGenJetAK4; i++)
	  {
	  	TLorentzVector gj;
		gj.SetPtEtaPhiM(GenJetAK4_pt[i], GenJetAK4_eta[i], GenJetAK4_phi[i], GenJetAK4_mass[i]);
                dr1[i] = leadbvec_nom.DeltaR(gj);
                dr2[i] = subleadbvec_nom.DeltaR(gj);
	  }

	  int indx1 = 0;
	  int indx2 = 0;	  
	  mindr1 = dr1[0];
	  mindr2 = dr2[0];
	  for(int i=1; i<nGenJetAK4; i++)
	  {
	  	if(mindr1 > dr1[i])
		{
			mindr1 = dr1[i];
			indx1 = i;
		
		}
		if(mindr2 > dr2[i])
                {
                        mindr2 = dr2[i];
			indx2 = i;
                }
	  }

	  float pt_nom_lead, pt_up_lead, pt_dn_lead, pt_nom_sublead, pt_up_sublead, pt_dn_sublead;
	  if (mindr1 < 0.4 && leadbvec_nom.Pt()/GenJetAK4_pt[indx1] < 1.2 && leadbvec_nom.Pt()/GenJetAK4_pt[indx1] > 0.8)
	  {
		  pt_nom_lead = (GenJetAK4_pt[indx1] + (reg_scale_nom_lead.Pt()-GenJetAK4_pt[indx1])*(1 + smr_nom));
	   	  pt_up_lead = (GenJetAK4_pt[indx1] + (reg_scale_nom_lead.Pt()-GenJetAK4_pt[indx1])*(1 + smr_up));
	  	  pt_dn_lead = (GenJetAK4_pt[indx1] + (reg_scale_nom_lead.Pt()-GenJetAK4_pt[indx1])*(1 + smr_dn));
	  }
	  else
	  {
		  pt_nom_lead = reg_scale_nom_lead.Pt();
		  pt_up_lead = reg_scale_up_lead.Pt();
		  pt_dn_lead = reg_scale_dn_lead.Pt();
	  }

	  if (mindr2 < 0.4 && subleadbvec_nom.Pt()/GenJetAK4_pt[indx2] < 1.2 && subleadbvec_nom.Pt()/GenJetAK4_pt[indx2] > 0.8)
          {
                  pt_nom_sublead = (GenJetAK4_pt[indx2] + (reg_scale_nom_sublead.Pt()-GenJetAK4_pt[indx2])*(1 + smr_nom));
                  pt_up_sublead = (GenJetAK4_pt[indx2] + (reg_scale_nom_sublead.Pt()-GenJetAK4_pt[indx2])*(1 + smr_up));
                  pt_dn_sublead = (GenJetAK4_pt[indx2] + (reg_scale_nom_sublead.Pt()-GenJetAK4_pt[indx2])*(1 + smr_dn));
          }
          else
          {
                  pt_nom_sublead = reg_scale_nom_sublead.Pt();
                  pt_up_sublead = reg_scale_up_sublead.Pt();
                  pt_dn_sublead = reg_scale_dn_sublead.Pt();
          }

	  TLorentzVector final_nom_lead, final_nom_sublead, final_up_lead, final_up_sublead, final_dn_lead, final_dn_sublead;

	  final_nom_lead.SetPtEtaPhiM(pt_nom_lead, PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
          final_nom_sublead.SetPtEtaPhiM(pt_nom_sublead, PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);
          final_up_lead.SetPtEtaPhiM(pt_up_lead, PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
	  final_up_sublead.SetPtEtaPhiM(pt_up_sublead, PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);
	  final_dn_lead.SetPtEtaPhiM(pt_dn_lead, PFJetAK4_eta[0], PFJetAK4_phi[0], Jet_mass_nom[0]);
          final_dn_sublead.SetPtEtaPhiM(pt_dn_sublead, PFJetAK4_eta[1], PFJetAK4_phi[1], Jet_mass_nom[1]);


	  //-----------------------------------------------------photon info-----------------------------------------------------//
	
	
	  indx = 0;

	  for(int i=0; i<nPhoton; i++)
	  {
		  if (Photon_PixelSeed[i]) continue;
	//	  if (!Photon_passEveto[i]) continue;
		  if (Photon_pt[i] < 10) continue;
                  if (fabs(Photon_eta[i]) > 2.5) continue;
                  if (fabs(Photon_eta[i]) > 1.44 && fabs(Photon_eta[i]) < 1.57) continue;
      
		  if (fabs(Photon_eta[i]) <= 1.44)
		  {
		  	if (Photon_mvaid_Fall17V2_raw[i] < -0.02) continue;	//90% signal efficiency
		  }
		  if (fabs(Photon_eta[i]) >= 1.57)
                  {
                        if (Photon_mvaid_Fall17V2_raw[i] < -0.26) continue;	//90% signal efficiency 
                  }

		  TLorentzVector phovec;
                  phovec.SetPtEtaPhiE(Photon_pt[i]*(Photon_ecalEnergyPostCorr[i]/Photon_ecalEnergyPreCorr[i]), Photon_eta[i], Photon_phi[i], Photon_ecalEnergyPostCorr[i]);

		  if (leptonvec.DeltaR(phovec) < 0.4) continue;
                  if ((leadbvec_nom.DeltaR(phovec) < 0.4) || (subleadbvec_nom.DeltaR(phovec) < 0.4)) continue;   

                  Photon_pt[indx] = Photon_pt[i]*(Photon_ecalEnergyPostCorr[i]/Photon_ecalEnergyPreCorr[i]);
		  Photon_e[indx] = Photon_ecalEnergyPostCorr[i];
                  Photon_eta[indx] = Photon_eta[i];
                  Photon_phi[indx] = Photon_phi[i];
		  Photon_mvaid_Fall17V2_raw[indx] = Photon_mvaid_Fall17V2_raw[i];
		  Photon_e9by25[indx] = Photon_e9by25[i];

                  if (++indx >= njetmx) break;
	  }

	  nPhoton = indx;
          if (nPhoton < 2) continue;

	  if (isEle)
	  {
		  if (bb == 1) {cut_flow_el->Fill(8);}
                  if (ab == 1) {cut_flow_el->Fill(9);}
	  }

	  if (isMu)
          {
                  if (bb == 1) {cut_flow_mu->Fill(8);}
                  if (ab == 1) {cut_flow_mu->Fill(9);}
          }

	  for (int i=0; i<nPhoton; i++)
	  {
		  for (int j=i+1; j<nPhoton; j++)
		  {
			  if(Photon_pt[i]<Photon_pt[j])
			  {
				  float a = Photon_e[i];
                                  Photon_e[i] = Photon_e[j];
                                  Photon_e[j] = a;

                                  a = Photon_pt[i];
                                  Photon_pt[i] = Photon_pt[j];
                                  Photon_pt[j] = a;

                                  a = Photon_eta[i];
                                  Photon_eta[i] = Photon_eta[j];
                                  Photon_eta[j] = a;

                                  a = Photon_phi[i];
                                  Photon_phi[i] = Photon_phi[j];
                                  Photon_phi[j] = a;

				  a = Photon_mvaid_Fall17V2_raw[i];
                                  Photon_mvaid_Fall17V2_raw[i] = Photon_mvaid_Fall17V2_raw[j];
                                  Photon_mvaid_Fall17V2_raw[j] = a;

				  a = Photon_e9by25[i];
                                  Photon_e9by25[i] = Photon_e9by25[j];
                                  Photon_e9by25[j] = a;
			  }
		  }
	  }

	  TLorentzVector leadphovec, subleadphovec;

          leadphovec.SetPtEtaPhiE(Photon_pt[0], Photon_eta[0], Photon_phi[0], Photon_e[0]);
          subleadphovec.SetPtEtaPhiE(Photon_pt[1], Photon_eta[1], Photon_phi[1], Photon_e[1]);


	  //-----------------------------------------------gen level study-----------------------------------------------//


	
	  int indx_gen = 0;
	  for (int j=0; j<nGenPart; j++)
	  {
	  	if (fabs(GenPart_pdg[j]) == 5 && GenPart_mompdg[j] == 35)
		{
			GenPart_pt[indx_gen] = GenPart_pt[j];
			GenPart_eta[indx_gen] = GenPart_eta[j];
			GenPart_phi[indx_gen] = GenPart_phi[j];
			GenPart_mass[indx_gen] = GenPart_mass[j];
			indx_gen++;
		}
	  }

	  int bpart = indx_gen;

	  for (int j=0; j<bpart; j++)
	  {
	  	for (int k=j+1; k<bpart; k++)
		{
			if(GenPart_pt[j]<GenPart_pt[k])
			{
				float a = GenPart_pt[j];
                                GenPart_pt[j] = GenPart_pt[k];
                                GenPart_pt[k] = a;

				a = GenPart_eta[j];
                                GenPart_eta[j] = GenPart_eta[k];
                                GenPart_eta[k] = a;
				
				a = GenPart_phi[j];
                                GenPart_phi[j] = GenPart_phi[k];
                                GenPart_phi[k] = a;

				a = GenPart_mass[j];
                                GenPart_mass[j] = GenPart_mass[k];
                                GenPart_mass[k] = a;
			}
		}
	  }

	  TLorentzVector b1g, b2g;
	  b1g.SetPtEtaPhiM(GenPart_pt[0], GenPart_eta[0], GenPart_phi[0], GenPart_mass[0]);
	  b2g.SetPtEtaPhiM(GenPart_pt[1], GenPart_eta[1], GenPart_phi[1], GenPart_mass[1]);

	  float mindr1_b, mindr2_b;
          float dr1_b[njetmx], dr2_b[njetmx];
          int b1 = 0;
          int b2 = 0;

	  for (int j=0; j<nGenJetAK4; j++)
	  {
	  	TLorentzVector gj;
		gj.SetPtEtaPhiM(GenJetAK4_pt[j], GenJetAK4_eta[j], GenJetAK4_phi[j], GenJetAK4_mass[j]);
		dr1_b[j] = b1g.DeltaR(gj);
                dr2_b[j] = b2g.DeltaR(gj);
	  }

          mindr1_b = dr1_b[0];
          mindr2_b = dr2_b[0];

          for(int j=1; j<nGenJetAK4; j++)
          {
                if(mindr1_b > dr1_b[j])
                {
                        mindr1_b = dr1_b[j];
                        b1 = j;

                }
                if(mindr2_b > dr2_b[j])
                {
                        mindr2_b = dr2_b[j];
                        b2 = j;
                }
          }
	  
	  TLorentzVector gb1, gb2;
	  gb1.SetPtEtaPhiM(GenJetAK4_pt[b1], GenJetAK4_eta[b1], GenJetAK4_phi[b1], GenJetAK4_mass[b1]);
	  gb2.SetPtEtaPhiM(GenJetAK4_pt[b2], GenJetAK4_eta[b2], GenJetAK4_phi[b2], GenJetAK4_mass[b2]); 



	  //-----------------------------------------------all the SFs & etc---------------------------------------------//



          double leptonsf_weight = 1.0;
	  double leptonsf_weight_up = 1.0;
          double leptonsf_weight_dn = 1.0;
	  double leptonsf_weight_stat = 1.0;
	  double leptonsf_weight_syst = 1.0;

	  double photonsf_weight = 1.0;
          double photonsf_weight_up = 1.0;
          double photonsf_weight_dn = 1.0;
          double photonsf_weight_stat = 1.0;
          double photonsf_weight_syst = 1.0;

	  double haspix_wt = 1.0;
	  double haspix_wt_up = 1.0;
	  double haspix_wt_dn = 1.0;

	  double JetPuSF       =1;
   	  double PuJetIDSF     =1;
   	  double JetPuSF_Up    =1;
   	  double PuJetIDSF_Up  =1;
   	  double JetPuSF_Down  =1;
   	  double PuJetIDSF_Down=1;

	  if(isMC){
	  for(int i=0; i<nPFJetAK4; i++)
	  for(int i=0; i<AK4coll; i++)
	  {
	  	if ((PFJetAK4_pt[i] < 10) || (PFJetAK4_pt[i] > 50)) continue;
		if(!PFJetAK4_jetID_tightlepveto[i]) continue;
                if(fabs(PFJetAK4_eta[i]) > 2.5) continue;
	  
		bool matched =false;
		double eff_Pu       = 0;
     		double WJET_Pu      = 1;
     		double WJET_Pu_Up   = 1;
     		double WJET_Pu_Down = 1;

     		Int_t bin_pt  = SF_Eff->GetXaxis()->FindBin(PFJetAK4_pt[i]);
     		Int_t bin_eta = SF_Eff->GetYaxis()->FindBin(PFJetAK4_eta[i]);

		TLorentzVector vpfj;
		vpfj.SetPtEtaPhiM(PFJetAK4_pt[i], PFJetAK4_eta[i], PFJetAK4_phi[i], PFJetAK4_mass[i]);

		for(int kk=0; kk<nGenJetAK4; kk++)
      		{
			TLorentzVector vgenj;
                	vgenj.SetPtEtaPhiM(GenJetAK4_pt[kk], GenJetAK4_eta[kk], GenJetAK4_phi[kk], GenJetAK4_mass[kk]);
       			if(vpfj.DeltaR(vgenj) < 0.4)
        		{
         			matched=true;
         			break;
        		}
      		}

		if(matched)
        	{
         		JetPuSF      = SF_Eff->GetBinContent(bin_pt,bin_eta);
         		JetPuSF_Up   = JetPuSF+SF_Eff_Err->GetBinContent(bin_pt,bin_eta);
         		JetPuSF_Down = JetPuSF-SF_Eff_Err->GetBinContent(bin_pt,bin_eta);
         		eff_Pu       = Eff->GetBinContent(bin_pt,bin_eta);
        	}
      		if(!matched)
        	{
         		JetPuSF      = SF_MisTag->GetBinContent(bin_pt,bin_eta);
         		JetPuSF_Up   = JetPuSF+SF_MisTag_Err->GetBinContent(bin_pt,bin_eta);
         		JetPuSF_Down = JetPuSF-SF_MisTag_Err->GetBinContent(bin_pt,bin_eta);
         		eff_Pu       = MisTag->GetBinContent(bin_pt,bin_eta);
        	}

		WJET_Pu      = JetPuSF;
           	WJET_Pu_Up   = JetPuSF_Up;
           	WJET_Pu_Down = JetPuSF_Down;

		PuJetIDSF      *= WJET_Pu;
           	PuJetIDSF_Up   *= WJET_Pu_Up;
           	PuJetIDSF_Down *= WJET_Pu_Down;
	  }
	  }

	  double puWeight, puWeightup, puWeightdown;
	  puWeight = puWeightup = puWeightdown = 1.0;

	  if(isMC){

		if(npu_vert_true>=0 && npu_vert_true<100){
			float *puweights = Get_PU_Weights(file_pu_ratio, npu_vert_true);
			puWeight = puweights[0];
			puWeightup = puweights[1];
			puWeightdown = puweights[2];
		}

                if (isEle){
		for(unsigned lep=0; lep<nele; lep++){
				float *sfvalues = Electron_SF(file_el_sf, Electron_pt[lep], Electron_eta[lep]);
				leptonsf_weight *= sfvalues[0];
				leptonsf_weight_up *= sfvalues[1];
				leptonsf_weight_dn *= sfvalues[2];
				leptonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4]));  // like this for time being
				leptonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6] + sfvalues[7]*sfvalues[7] + sfvalues[8]*sfvalues[8]));  // like this for time being
			} }

		if (isMu){
                for(unsigned lep=0; lep<nmu; lep++){
			        float *sfvalues;
				sfvalues = Muon_SF(file_mu_sf, muon_id_name, Muon_pt[lep], Muon_eta[lep]);
				leptonsf_weight *= *(sfvalues+0);
				leptonsf_weight_up *= *(sfvalues+1);
				leptonsf_weight_dn *= *(sfvalues+2);
				leptonsf_weight_stat *= *(sfvalues+3);
				leptonsf_weight_syst *= *(sfvalues+4);
                        } }

		for(unsigned ph=0; ph<nPhoton; ph++){
                                float *sfvalues = Photon_SF(file_pho_sf, Photon_pt[ph], Photon_eta[ph]);
                                photonsf_weight *= sfvalues[0];
                                photonsf_weight_up *= sfvalues[1];
                                photonsf_weight_dn *= sfvalues[2];
                                photonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4]));  // like this for time being
                                photonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6] + sfvalues[7]*sfvalues[7] + sfvalues[8]*sfvalues[8]));  // like this for time being
				}
      

			if (fabs(leadphovec.Eta()) <= 1.44 && fabs(subleadphovec.Eta()) <= 1.44)
  			{
				haspix_wt = SF_HasPix_EB * SF_HasPix_EB;
				haspix_wt_up = SF_HasPix_EB_up * SF_HasPix_EB_up;
				haspix_wt_dn = SF_HasPix_EB_dn * SF_HasPix_EB_dn;
			}
			if (fabs(leadphovec.Eta()) >= 1.57 && fabs(subleadphovec.Eta()) >= 1.57)
                        {
                                haspix_wt = SF_HasPix_EE * SF_HasPix_EE;
				haspix_wt_up = SF_HasPix_EE_up * SF_HasPix_EE_up;
                                haspix_wt_dn = SF_HasPix_EE_dn * SF_HasPix_EE_dn;
                        }
			if ((fabs(leadphovec.Eta()) <= 1.44 && fabs(subleadphovec.Eta()) >= 1.57) || (fabs(subleadphovec.Eta()) <= 1.44 && fabs(leadphovec.Eta()) >= 1.57))
                        {
                                haspix_wt = SF_HasPix_EB * SF_HasPix_EE;
				haspix_wt_up = SF_HasPix_EB_up * SF_HasPix_EE_up;
                                haspix_wt_dn = SF_HasPix_EB_dn * SF_HasPix_EE_dn;
                        }

                }


	double SF_Trig, SF_Trig_stat, SF_Trig_syst;
	SF_Trig = 1; SF_Trig_stat = 1; SF_Trig_syst = 1;

 	if(!isDATA) {
		      if(isMu) {
		      if(leptonvec.Pt() > 24)
		      {
				int etabin = std::max(1, std::min(h_singlemuon->GetNbinsX(), h_singlemuon->GetXaxis()->FindBin(fabs(leptonvec.Eta()))));
                int ptbin  = std::max(1, std::min(h_singlemuon->GetNbinsY(), h_singlemuon->GetYaxis()->FindBin(leptonvec.Pt())));
                SF_Trig = TMath::Max(float(1.e-3),float(h_singlemuon->GetBinContent(etabin, ptbin))) ;
		        SF_Trig_stat = TMath::Max(float(1.e-3),float(h_singlemuon_stat->GetBinContent(etabin, ptbin))) ;	
		        SF_Trig_syst = TMath::Max(float(1.e-3),float(h_singlemuon_syst->GetBinContent(etabin, ptbin))) ;
	}
	}      

			  if (isEle) {
			  if(leptonvec.Pt() > 32)                              // for 2016 27, for all other yrs 32;
			  {
				int etabin = std::max(1, std::min(h_singleele->GetNbinsX(), h_singleele->GetXaxis()->FindBin(leptonvec.Eta())));
                int ptbin  = std::max(1, std::min(h_singleele->GetNbinsY(), h_singleele->GetYaxis()->FindBin(leptonvec.Pt())));
                SF_Trig = TMath::Max(float(1.e-3),float(h_singleele->GetBinContent(etabin, ptbin))) ;
				SF_Trig_stat = TMath::Max(float(1.e-3),float(h_singleele->GetBinContent(etabin, ptbin))) ;
				SF_Trig_syst = TMath::Max(float(1.e-3),float(h_singleele->GetBinContent(etabin, ptbin))) ;
	}
	}
}

	double SF_Trig_1_up = SF_Trig + sqrt((SF_Trig_stat - SF_Trig)*(SF_Trig_stat - SF_Trig) +  (SF_Trig_syst - SF_Trig)*(SF_Trig_syst - SF_Trig));
        double SF_Trig_1_dn = SF_Trig - sqrt((SF_Trig_stat - SF_Trig)*(SF_Trig_stat - SF_Trig) +  (SF_Trig_syst - SF_Trig)*(SF_Trig_syst - SF_Trig));

        double SF_Trig_2_up = SF_Trig + abs(SF_Trig-1);
        double SF_Trig_2_dn = SF_Trig - abs(SF_Trig-1);



	  //------------------------------------------------------tree fill--------------------------------------------------//


	  if (isEle || isMu) {

	  leppt = leptonvec.Pt();
          lepy = leptonvec.Rapidity();
	  lepeta = leptonvec.Eta();
	  lepphi = leptonvec.Phi();
	  lepe = leptonvec.Energy();

	  Met = miset;
	  MetPhi = misphi;
	  jet_no = nPFJetAK4;

          b1pt = leadbvec_nom.Pt();
          b1y = leadbvec_nom.Rapidity();
          b1eta = leadbvec_nom.Eta();
          b1phi = leadbvec_nom.Phi();
          b1e = leadbvec_nom.Energy();
          b2pt = subleadbvec_nom.Pt();
          b2y = subleadbvec_nom.Rapidity();
          b2eta = subleadbvec_nom.Eta();
          b2phi = subleadbvec_nom.Phi();
          b2e = subleadbvec_nom.Energy();
	  bb_inv_mass = (leadbvec_nom+subleadbvec_nom).M();
	  invmassbbgg = (leadbvec_nom+subleadbvec_nom+leadphovec+subleadphovec).M();
	  
	  b1pt_jecup = leadbvec_jecup.Pt();
	  b1e_jecup = leadbvec_jecup.E();
	  b2pt_jecup = subleadbvec_jecup.Pt();
          b2e_jecup = subleadbvec_jecup.E();

	  b1pt_jecdn = leadbvec_jecdn.Pt();
          b1e_jecdn = leadbvec_jecdn.E();
          b2pt_jecdn = subleadbvec_jecdn.Pt();
          b2e_jecdn = subleadbvec_jecdn.E();

	  b1pt_resup = leadbvec_resup.Pt();
          b1e_resup = leadbvec_resup.E();
          b2pt_resup = subleadbvec_resup.Pt();
          b2e_resup = subleadbvec_resup.E();

	  b1pt_resdn = leadbvec_resdn.Pt();
          b1e_resdn = leadbvec_resdn.E();
          b2pt_resdn = subleadbvec_resdn.Pt();
          b2e_resdn = subleadbvec_resdn.E();

	  b1pt_reg = lbvec_reg.Pt();
	  b1e_reg = lbvec_reg.E();
	  b2pt_reg = slbvec_reg.Pt();
	  b2e_reg = slbvec_reg.E();
	  b1pt_scl = reg_scale_nom_lead.Pt();
          b1e_scl = reg_scale_nom_lead.E();
          b2pt_scl = reg_scale_nom_sublead.Pt();
          b2e_scl = reg_scale_nom_sublead.E();

	  if(isMC){
	  b1pt_final = final_nom_lead.Pt();
          b1y_final = final_nom_lead.Rapidity();
          b1eta_final = final_nom_lead.Eta();
          b1phi_final = final_nom_lead.Phi();
          b1e_final = final_nom_lead.Energy();
          b2pt_final = final_nom_sublead.Pt();
          b2y_final = final_nom_sublead.Rapidity();
          b2eta_final = final_nom_sublead.Eta();
          b2phi_final = final_nom_sublead.Phi();
          b2e_final = final_nom_sublead.Energy();
	  bb_inv_mass_final = (final_nom_lead + final_nom_sublead).M();
	  invmassbbgg_final = (final_nom_lead+final_nom_sublead+leadphovec+subleadphovec).M();
          }

	  b1_DeepFlv = PFJetAK4_btag_DeepFlav[0];
          b2_DeepFlv = PFJetAK4_btag_DeepFlav[1];
	  b1_hadFlv = PFJetAK4_hadronflav[0];
          b2_hadFlv = PFJetAK4_hadronflav[1];
	  b1_partFlv = PFJetAK4_partonflav[0];
	  b2_partFlv = PFJetAK4_partonflav[1];
	  b1_PUid = PFJetAK4_PUID[0];
	  b2_PUid = PFJetAK4_PUID[1];

          if(isMC){

          event_weight = 1.0;
      
          weight_nom = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;

          weight_PU_up = event_weight * puWeightup * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;
          weight_leptonsf_up = event_weight * puWeight * leptonsf_weight_up * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;
          weight_bSF_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF_up[0] * PFJetAK4_btag_DeepFlav_SF_up[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;
	  weight_trig_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig_1_up * SF_Trig_2_up * photonsf_weight * haspix_wt;
	  weight_photonsf_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight_up * haspix_wt;
	  weight_haspixsf_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt_up;
	  weight_pujid_up = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF_Up * SF_Trig * photonsf_weight * haspix_wt_up;

          weight_PU_dn = event_weight * puWeightdown * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;
          weight_leptonsf_dn = event_weight * puWeight * leptonsf_weight_dn * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;
          weight_bSF_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF_dn[0] * PFJetAK4_btag_DeepFlav_SF_dn[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt;
	  weight_trig_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig_1_dn * SF_Trig_2_dn * photonsf_weight * haspix_wt;
	  weight_photonsf_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight_dn * haspix_wt;
	  weight_haspixsf_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF * SF_Trig * photonsf_weight * haspix_wt_dn;
	  weight_pujid_dn = event_weight * puWeight * leptonsf_weight * PFJetAK4_btag_DeepFlav_SF[0] * PFJetAK4_btag_DeepFlav_SF[1] * PuJetIDSF_Down * SF_Trig * photonsf_weight * haspix_wt_up;

	  puwt = puWeight;
	  pujidwt = PuJetIDSF;

          }  

          pho_no = nPhoton;
          pho1pt = leadphovec.Pt();
          pho1y = leadphovec.Rapidity();
          pho1eta = leadphovec.Eta();
          pho1phi = leadphovec.Phi();
          pho1e = leadphovec.Energy();
          pho2pt = subleadphovec.Pt();
          pho2y = subleadphovec.Rapidity();
          pho2eta = subleadphovec.Eta();
          pho2phi = subleadphovec.Phi();
          pho2e = subleadphovec.Energy();
	  pho1MVA = Photon_mvaid_Fall17V2_raw[0];
	  pho2MVA = Photon_mvaid_Fall17V2_raw[1];
	  pho1r9 = Photon_e9by25[0];
	  pho2r9 = Photon_e9by25[1];
          dipho_invmass = (leadphovec+subleadphovec).M();   

	  puvtx = npu_vert_true;
	  evt_rho = Rho;
	  
	  genb1_pt = gb1.Pt();
	  genb1_eta = gb1.Eta();
  	  genb1_phi = gb1.Phi();
	  genb1_mass = gb1.M();
	  genb2_pt = gb2.Pt();
          genb2_eta = gb2.Eta();
          genb2_phi = gb2.Phi();
          genb2_mass = gb2.M();
	  genb_invmass = (gb1+gb2).M();

	  LHEPDF_err = getTheoryEsystematics_PDF(nLHEPDFWeights,LHEPDFWeights,true);

	  Tout->Fill();

}

}

          f->cd();
          delete T1;
          delete f;

}
        infile.close();
        fout->cd();
        fout->Write();
        fout->Close();

        return 0;

}
