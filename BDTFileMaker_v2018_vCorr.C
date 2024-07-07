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


void BDTFileMaker_v2018_vCorr()
{

	float Leppt, Lepeta, Lepphi, Lepe, MetPhi;
	float B1pt, B1eta, B1phi, B1e;
	float B2pt, B2eta, B2phi, B2e;
	float B1pt_jecup, B1e_jecup, B2pt_jecup, B2e_jecup, B1pt_jecdn, B1e_jecdn, B2pt_jecdn, B2e_jecdn;
        float B1pt_resup, B1e_resup, B2pt_resup, B2e_resup, B1pt_resdn, B1e_resdn, B2pt_resdn, B2e_resdn;
	float bb_inv_mass, b1_DeepFlv, b2_DeepFlv;
	float Pho1pt, Pho1eta, Pho1phi, Pho1e, pho1MVA;
	float Pho2pt, Pho2eta, Pho2phi, Pho2e, pho2MVA;
	float dipho_invmass, invmassbbgg, delphi_ggMET, delphi_bb, delphi_gg, delphi_bbgg;
	double weight_nom, weight_PU_up, weight_leptonsf_up, weight_photonsf_up, weight_bSF_up, weight_trig_up, weight_haspixsf_up, weight_pujid_up, weight_PU_dn, weight_leptonsf_dn, weight_photonsf_dn, weight_prefiring_dn, weight_bSF_dn, weight_trig_dn, weight_haspixsf_dn, weight_pujid_dn;
	bool isEle, isMu;
	int ab, bb;


	float leppt, lepeta, lepphi, lepe;
        float b1pt, b1eta, b1phi, b1e;
        float b2pt, b2eta, b2phi, b2e;
	float b1pt_jecup, b1e_jecup, b2pt_jecup, b2e_jecup, b1pt_jecdn, b1e_jecdn, b2pt_jecdn, b2e_jecdn;
	float b1pt_resup, b1e_resup, b2pt_resup, b2e_resup, b1pt_resdn, b1e_resdn, b2pt_resdn, b2e_resdn;
        float leadbdeepjet, subleadbdeepjet;
        float pho1pt, pho1eta, pho1phi, pho1e, leadgMVA, pT1_by_mgg, pT1_by_mbb;
        float pho2pt, pho2eta, pho2phi, pho2e, subleadgMVA, pT2_by_mgg, pT2_by_mbb;
        float m_bb, m_gg, m_bbgg;
	double Weight_nom, Weight_PU_up, Weight_leptonsf_up, Weight_photonsf_up, Weight_bSF_up, Weight_trig_up, Weight_haspixsf_up, Weight_pujid_up, Weight_PU_dn, Weight_leptonsf_dn, Weight_photonsf_dn, Weight_prefiring_dn, Weight_bSF_dn, Weight_trig_dn, Weight_haspixsf_dn, Weight_pujid_dn;
        float evt_wgt;
	double wt;
	int yout;
	int jet_no, Njets;


	double TTGJets = (59730*4.078)/(2.66557e+07);
        double DY = (59730*5343)/(9.14154e+07);
        double TTSL = (59730*365.34)/(1.48113e+11);
        double TT2L2Nu = (59730*88.29)/(1.03533e+10);
	double WH_mA20 = (59730*0.01)/(993835);
        double WH_mA25 = (59730*0.01)/(993842);
        double WH_mA30 = (59730*0.01)/(732869);
        double WH_mA35 = (59730*0.01)/(834837);
        double WH_mA40 = (59730*0.01)/(873839);
        double WH_mA45 = (59730*0.01)/(753850);
        double WH_mA50 = (59730*0.01)/(618893);
        double WH_mA55 = (59730*0.01)/(906844);
        double WH_mA60 = (59730*0.01)/(707896);


	TString str[19] = {"weight_nom", "weight_PU_up", "weight_leptonsf_up", "weight_photonsf_up", "weight_bSF_up", "weight_trig_up", "weight_haspixsf_up", "weight_pujid_up", "weight_jec_up", "weight_jec_dn", "weight_jer_up", "weight_jer_dn", "weight_PU_dn", "weight_leptonsf_dn", "weight_photonsf_dn", "weight_bSF_dn", "weight_trig_dn", "weight_haspixsf_dn", "weight_pujid_dn"};


	TFile *f = new TFile("/home/abala/local_bkp/analysis/root_file/WH_2018/full_reg/corr_v1/WH_mA60_v6.root");
        TTree *tr = (TTree*)f->Get("Tout");

	TFile *fout = new TFile("/home/abala/local_bkp/analysis/BDT/2018/WH60/WH_mA60_vf_weight_nom.root","RECREATE");
	TTree *tree = new TTree("tree","tree");

   	tree->Branch("leppt",&leppt);
	tree->Branch("lepeta",&lepeta);
	
	tree->Branch("b1pt",&b1pt);
        tree->Branch("b1eta",&b1eta);
	tree->Branch("b1phi",&b1phi);
        tree->Branch("b1e",&b1e);
	tree->Branch("b2pt",&b2pt);
        tree->Branch("b2eta",&b2eta);
	tree->Branch("b2phi",&b2phi);
        tree->Branch("b2e",&b2e);

	tree->Branch("b1pt_jecup",&b1pt_jecup);
	tree->Branch("b1e_jecup",&b1e_jecup);
	tree->Branch("b2pt_jecup",&b2pt_jecup);
        tree->Branch("b2e_jecup",&b2e_jecup);
	tree->Branch("b1pt_jecdn",&b1pt_jecdn);
        tree->Branch("b1e_jecdn",&b1e_jecdn);
	tree->Branch("b2pt_jecdn",&b2pt_jecdn);
        tree->Branch("b2e_jecdn",&b2e_jecdn);

	tree->Branch("b1pt_resup",&b1pt_resup);
        tree->Branch("b1e_resup",&b1e_resup);
        tree->Branch("b2pt_resup",&b2pt_resup);
        tree->Branch("b2e_resup",&b2e_resup);
        tree->Branch("b1pt_resdn",&b1pt_resdn);
        tree->Branch("b1e_resdn",&b1e_resdn);
        tree->Branch("b2pt_resdn",&b2pt_resdn);
        tree->Branch("b2e_resdn",&b2e_resdn);

	tree->Branch("pho1pt",&pho1pt);
        tree->Branch("pho1eta",&pho1eta);
	tree->Branch("pho1phi",&pho1phi);
        tree->Branch("pho1e",&pho1e);
        tree->Branch("pho2pt",&pho2pt);
        tree->Branch("pho2eta",&pho2eta);
	tree->Branch("pho2phi",&pho2phi);
        tree->Branch("pho2e",&pho2e);

	tree->Branch("m_bb",&m_bb);
	tree->Branch("m_gg",&m_gg);
	tree->Branch("m_bbgg",&m_bbgg);
	tree->Branch("pT1_by_mgg",&pT1_by_mgg);
	tree->Branch("pT2_by_mgg",&pT2_by_mgg);
	tree->Branch("pT1_by_mbb",&pT1_by_mbb);
        tree->Branch("pT2_by_mbb",&pT2_by_mbb);
	tree->Branch("Njets",&Njets);

	tree->Branch("leadbdeepjet",&leadbdeepjet);
	tree->Branch("subleadbdeepjet",&subleadbdeepjet);
	tree->Branch("leadgMVA",&leadgMVA);
	tree->Branch("subleadgMVA",&subleadgMVA);

	tree->Branch("delphi_ggMET",&delphi_ggMET);
	tree->Branch("delphi_gg",&delphi_gg);
	tree->Branch("delphi_bb",&delphi_bb);
	tree->Branch("delphi_bbgg",&delphi_bbgg);

	tree->Branch("evt_wgt",&evt_wgt);
	tree->Branch("yout",&yout);  

	tree->Branch("weight_nom",&weight_nom);
    	tree->Branch("weight_PU_up",&weight_PU_up);
    	tree->Branch("weight_leptonsf_up",&weight_leptonsf_up);
    	tree->Branch("weight_photonsf_up",&weight_photonsf_up);
    	tree->Branch("weight_bSF_up",&weight_bSF_up);
    	tree->Branch("weight_trig_up",&weight_trig_up);
    	tree->Branch("weight_haspixsf_up",&weight_haspixsf_up);
    	tree->Branch("weight_pujid_up",&weight_pujid_up);
    	tree->Branch("weight_PU_dn",&weight_PU_dn);
    	tree->Branch("weight_leptonsf_dn",&weight_leptonsf_dn);
    	tree->Branch("weight_photonsf_dn",&weight_photonsf_dn);
    	tree->Branch("weight_bSF_dn",&weight_bSF_dn);
    	tree->Branch("weight_trig_dn",&weight_trig_dn);
    	tree->Branch("weight_haspixsf_dn",&weight_haspixsf_dn);
    	tree->Branch("weight_pujid_dn",&weight_pujid_dn);


        tr->SetBranchAddress("leppt",&Leppt);
        tr->SetBranchAddress("lepeta",&Lepeta);

	tr->SetBranchAddress("MetPhi",&MetPhi);

        tr->SetBranchAddress("b1pt",&B1pt);
        tr->SetBranchAddress("b1eta",&B1eta);
	tr->SetBranchAddress("b1phi",&B1phi);
	tr->SetBranchAddress("b1e",&B1e);
        tr->SetBranchAddress("b2pt",&B2pt);
        tr->SetBranchAddress("b2eta",&B2eta);
	tr->SetBranchAddress("b2phi",&B2phi);
        tr->SetBranchAddress("b2e",&B2e);

	tr->SetBranchAddress("b1pt_jecup",&B1pt_jecup);
        tr->SetBranchAddress("b1e_jecup",&B1e_jecup);
        tr->SetBranchAddress("b2pt_jecup",&B2pt_jecup);
        tr->SetBranchAddress("b2e_jecup",&B2e_jecup);
        tr->SetBranchAddress("b1pt_jecdn",&B1pt_jecdn);
        tr->SetBranchAddress("b1e_jecdn",&B1e_jecdn);
        tr->SetBranchAddress("b2pt_jecdn",&B2pt_jecdn);
        tr->SetBranchAddress("b2e_jecdn",&B2e_jecdn);

        tr->SetBranchAddress("b1pt_resup",&B1pt_resup);
        tr->SetBranchAddress("b1e_resup",&B1e_resup);
        tr->SetBranchAddress("b2pt_resup",&B2pt_resup);
        tr->SetBranchAddress("b2e_resup",&B2e_resup);
        tr->SetBranchAddress("b1pt_resdn",&B1pt_resdn);
        tr->SetBranchAddress("b1e_resdn",&B1e_resdn);
        tr->SetBranchAddress("b2pt_resdn",&B2pt_resdn);
        tr->SetBranchAddress("b2e_resdn",&B2e_resdn);

        tr->SetBranchAddress("pho1pt",&Pho1pt);
        tr->SetBranchAddress("pho1eta",&Pho1eta);
        tr->SetBranchAddress("pho1phi",&Pho1phi);
	tr->SetBranchAddress("pho1e",&Pho1e);
        tr->SetBranchAddress("pho2pt",&Pho2pt);
        tr->SetBranchAddress("pho2eta",&Pho2eta);
        tr->SetBranchAddress("pho2phi",&Pho2phi);
	tr->SetBranchAddress("pho2e",&Pho2e);

        tr->SetBranchAddress("dipho_invmass",&dipho_invmass);
        tr->SetBranchAddress("invmassbbgg",&invmassbbgg);
	tr->SetBranchAddress("bb_inv_mass",&bb_inv_mass);

	tr->SetBranchAddress("b1_DeepFlv",&b1_DeepFlv);
	tr->SetBranchAddress("b2_DeepFlv",&b2_DeepFlv);
	tr->SetBranchAddress("pho1MVA",&pho1MVA);
	tr->SetBranchAddress("pho2MVA",&pho2MVA);
	tr->SetBranchAddress("jet_no",&jet_no);

	tr->SetBranchAddress("isMu",&isMu);
	tr->SetBranchAddress("isEle",&isEle);
	tr->SetBranchAddress("ab",&ab);
	tr->SetBranchAddress("bb",&bb);

	tr->SetBranchAddress("weight_nom",&Weight_nom);
	tr->SetBranchAddress("weight_PU_up",&Weight_PU_up);
	tr->SetBranchAddress("weight_leptonsf_up",&Weight_leptonsf_up);
	tr->SetBranchAddress("weight_photonsf_up",&Weight_photonsf_up);
	tr->SetBranchAddress("weight_bSF_up",&Weight_bSF_up);
	tr->SetBranchAddress("weight_trig_up",&Weight_trig_up);
	tr->SetBranchAddress("weight_haspixsf_up",&Weight_haspixsf_up);
	tr->SetBranchAddress("weight_pujid_up",&Weight_pujid_up);
	tr->SetBranchAddress("weight_PU_dn",&Weight_PU_dn);
	tr->SetBranchAddress("weight_leptonsf_dn",&Weight_leptonsf_dn);
	tr->SetBranchAddress("weight_photonsf_dn",&Weight_photonsf_dn);
	tr->SetBranchAddress("weight_bSF_dn",&Weight_bSF_dn);
	tr->SetBranchAddress("weight_trig_dn",&Weight_trig_dn);
	tr->SetBranchAddress("weight_haspixsf_dn",&Weight_haspixsf_dn);
	tr->SetBranchAddress("weight_pujid_dn",&Weight_pujid_dn);


	int nevt;
        nevt=tr->GetEntries();

	for(int ij=0; ij<nevt; ij++){
		tr->GetEntry(ij);

		if((isEle || isMu) && (bb == 1) && (Weight_nom > 0)){

		leppt = Leppt;
		lepeta = Lepeta;

		b1pt = B1pt;
		b1eta = B1eta;
		b1phi = B1phi;
                b1e = B1e;

		b2pt = B2pt;
                b2eta = B2eta;
		b2phi = B2phi;
                b2e = B2e;

		b1pt_jecup = B1pt_jecup;
                b1e_jecup = B1e_jecup;
                b2pt_jecup = B2pt_jecup;
                b2e_jecup = B2e_jecup;

                b1pt_jecdn = B1pt_jecdn;
                b1e_jecdn = B1e_jecdn;
                b2pt_jecdn = B2pt_jecdn;
                b2e_jecdn = B2e_jecdn;

		b1pt_resup = B1pt_resup;
                b1e_resup = B1e_resup;
                b2pt_resup = B2pt_resup;
                b2e_resup = B2e_resup;

                b1pt_resdn = B1pt_resdn;
                b1e_resdn = B1e_resdn;
                b2pt_resdn = B2pt_resdn;
                b2e_resdn = B2e_resdn;

		pho1pt = Pho1pt;
		pho1eta = Pho1eta;
		pho1phi = Pho1phi;
                pho1e = Pho1e;

		pho2pt = Pho2pt;
                pho2eta = Pho2eta;
		pho2phi = Pho2phi;
                pho2e = Pho2e;

		TLorentzVector g1, g2;
                g1.SetPtEtaPhiE(Pho1pt, Pho1eta, Pho1phi, Pho1e);
                g2.SetPtEtaPhiE(Pho2pt, Pho2eta, Pho2phi, Pho2e);
		delphi_gg = fabs(g1.DeltaPhi(g2));
                delphi_ggMET = fabs(PhiInRange((g1+g2).Phi() - MetPhi));

		TLorentzVector b1, b2;
                b1.SetPtEtaPhiE(B1pt, B1eta, B1phi, B1e);
                b2.SetPtEtaPhiE(B2pt, B2eta, B2phi, B2e);
		delphi_bb = fabs(b1.DeltaPhi(b2));
		delphi_bbgg = fabs((b1+b2).DeltaPhi(g1+g2));

		m_bb = bb_inv_mass;
		m_gg = dipho_invmass;
		m_bbgg = invmassbbgg;
		pT1_by_mgg = Pho1pt/dipho_invmass;
		pT2_by_mgg = Pho2pt/dipho_invmass;
		pT1_by_mbb = B1pt/bb_inv_mass;
                pT2_by_mbb = B2pt/bb_inv_mass;
		Njets = jet_no;

		leadbdeepjet = b1_DeepFlv;
		subleadbdeepjet = b2_DeepFlv;
		leadgMVA = pho1MVA;
		subleadgMVA = pho2MVA;

		weight_nom = fabs(Weight_nom);
		weight_PU_up = fabs(Weight_PU_up);
		weight_leptonsf_up = fabs(Weight_leptonsf_up);
		weight_photonsf_up = fabs(Weight_photonsf_up);
		weight_bSF_up = fabs(Weight_bSF_up);
		weight_trig_up = fabs(Weight_trig_up);
		weight_haspixsf_up = fabs(Weight_haspixsf_up);
		weight_pujid_up = fabs(Weight_pujid_up);
		weight_PU_dn = fabs(Weight_PU_dn);
                weight_leptonsf_dn = fabs(Weight_leptonsf_dn);
                weight_photonsf_dn = fabs(Weight_photonsf_dn);
                weight_bSF_dn = fabs(Weight_bSF_dn);
                weight_trig_dn = fabs(Weight_trig_dn);
                weight_haspixsf_dn = fabs(Weight_haspixsf_dn);
                weight_pujid_dn = fabs(Weight_pujid_dn);

		evt_wgt = weight_nom*WH_mA60;
	//	evt_wgt = 1.0;
		yout = 1;


                tree->Fill();}

	}


	f->cd();
        delete tr;
        delete f;

	fout->cd();
	fout->Write();
	fout->Close();
	

}
