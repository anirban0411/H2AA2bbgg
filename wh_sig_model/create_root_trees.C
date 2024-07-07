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


float BDT, CMS_hgg_mass, weight, dZ;
bool isEle, isMu;
int bb;
float bb_inv_mass, dipho_invmass, invmassbbgg;
float b1pt, b1eta, b1phi, b1e, b2pt, b2eta, b2phi, b2e, pho1pt, pho1eta, pho1phi, pho1e, pho2pt, pho2eta, pho2phi, pho2e;
float b1pt_jecup, b1e_jecup, b2pt_jecup, b2e_jecup, b1pt_jecdn , b1e_jecdn, b2pt_jecdn, b2e_jecdn, b1pt_resup, b1e_resup, b2pt_resup, b2e_resup, b1pt_resdn, b1e_resdn, b2pt_resdn, b2e_resdn;
float sigmabb1, sigmah1;
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_photonsf_up, weight_bSF_up, weight_trig_up, weight_haspixsf_up, weight_pujid_up, weight_PU_dn, weight_leptonsf_dn, weight_photonsf_dn, weight_prefiring_dn, weight_bSF_dn, weight_trig_dn, weight_haspixsf_dn, weight_pujid_dn;


void create_root_trees()
{

	double WH_mA20 = (59730*0.01)/(993835);
        double WH_mA25 = (59730*0.01)/(993842);
        double WH_mA30 = (59730*0.01)/(732869);
        double WH_mA35 = (59730*0.01)/(834837);
        double WH_mA40 = (59730*0.01)/(873839);
        double WH_mA45 = (59730*0.01)/(753850);
        double WH_mA50 = (59730*0.01)/(618893);
        double WH_mA55 = (59730*0.01)/(906844);
        double WH_mA60 = (59730*0.01)/(707896);

        double num = 707896;

        sigmabb1 = 12.182;
        sigmah1 = 12.508;

	TFile *f1 = new TFile("/home/abala/local_bkp/analysis/BDT/2018/WH60/WH_mA60_vf_weight_nom.root");
    	TFile *f2 = new TFile("/home/abala/local_bkp/analysis/BDT/2018/tree/final/BDT_UL18_MC_trees_sig_train_60.root");

    	TTree *Tout1 = (TTree*)f1->Get("tree");
    	TTree *Tout2 = (TTree*)f2->Get("Tree");

    	TFile *fout = new TFile("output_WHToAA2B2G_M60_pythia8_ggh_v1.root","RECREATE");

	const int nTrees = 15;

    	TString treeNames[nTrees] = {
        "ggh_60_13TeV_H2AA2bbgg_Tag0", "ggh_60_13TeV_H2AA2bbgg_Tag0_PUUp01sigma",
        "ggh_60_13TeV_H2AA2bbgg_Tag0_leptonIDUp01sigma", "ggh_60_13TeV_H2AA2bbgg_Tag0_photonIDUp01sigma",
   	"ggh_60_13TeV_H2AA2bbgg_Tag0_btagUp01sigma", "ggh_60_13TeV_H2AA2bbgg_Tag0_TrigUp01sigma",
	"ggh_60_13TeV_H2AA2bbgg_Tag0_JECUp01sigma", "ggh_60_13TeV_H2AA2bbgg_Tag0_JERUp01sigma",
	"ggh_60_13TeV_H2AA2bbgg_Tag0_PUDown01sigma", "ggh_60_13TeV_H2AA2bbgg_Tag0_leptonIDDown01sigma",
	"ggh_60_13TeV_H2AA2bbgg_Tag0_photonIDDown01sigma", "ggh_60_13TeV_H2AA2bbgg_Tag0_btagDown01sigma",
	"ggh_60_13TeV_H2AA2bbgg_Tag0_TrigDown01sigma", "ggh_60_13TeV_H2AA2bbgg_Tag0_JECDown01sigma",
	"ggh_60_13TeV_H2AA2bbgg_Tag0_JERDown01sigma"
   	 };

	TTree *trees[nTrees];
    	for (int i = 0; i < nTrees; ++i) {
        trees[i] = new TTree(treeNames[i], treeNames[i]);
        trees[i]->Branch("CMS_hgg_mass", &CMS_hgg_mass, "CMS_hgg_mass/F");
        trees[i]->Branch("weight", &weight, "weight/F");
        trees[i]->Branch("dZ", &dZ, "dZ/F");
   	}

	Tout1->SetBranchAddress("b1pt",&b1pt);
    	Tout1->SetBranchAddress("b1eta",&b1eta);
    	Tout1->SetBranchAddress("b1phi",&b1phi);
    	Tout1->SetBranchAddress("b1e",&b1e);
    	Tout1->SetBranchAddress("b2pt",&b2pt);
    	Tout1->SetBranchAddress("b2eta",&b2eta);
    	Tout1->SetBranchAddress("b2phi",&b2phi);
    	Tout1->SetBranchAddress("b2e",&b2e);
    	Tout1->SetBranchAddress("pho1pt",&pho1pt);
    	Tout1->SetBranchAddress("pho1eta",&pho1eta);
    	Tout1->SetBranchAddress("pho1phi",&pho1phi);
    	Tout1->SetBranchAddress("pho1e",&pho1e);
    	Tout1->SetBranchAddress("pho2pt",&pho2pt);
    	Tout1->SetBranchAddress("pho2eta",&pho2eta);
    	Tout1->SetBranchAddress("pho2phi",&pho2phi);
    	Tout1->SetBranchAddress("pho2e",&pho2e);

    	Tout1->SetBranchAddress("b1pt_jecup",&b1pt_jecup);
    	Tout1->SetBranchAddress("b1e_jecup",&b1e_jecup);
    	Tout1->SetBranchAddress("b2pt_jecup",&b2pt_jecup);
    	Tout1->SetBranchAddress("b2e_jecup",&b2e_jecup);
    	Tout1->SetBranchAddress("b1pt_jecdn",&b1pt_jecdn);
    	Tout1->SetBranchAddress("b1e_jecdn",&b1e_jecdn);
    	Tout1->SetBranchAddress("b2pt_jecdn",&b2pt_jecdn);
    	Tout1->SetBranchAddress("b2e_jecdn",&b2e_jecdn);
    	Tout1->SetBranchAddress("b1pt_resup",&b1pt_resup);
    	Tout1->SetBranchAddress("b1e_resup",&b1e_resup);
    	Tout1->SetBranchAddress("b2pt_resup",&b2pt_resup);
    	Tout1->SetBranchAddress("b2e_resup",&b2e_resup);
    	Tout1->SetBranchAddress("b1pt_resdn",&b1pt_resdn);
    	Tout1->SetBranchAddress("b1e_resdn",&b1e_resdn);
    	Tout1->SetBranchAddress("b2pt_resdn",&b2pt_resdn);
    	Tout1->SetBranchAddress("b2e_resdn",&b2e_resdn);

    	Tout1->SetBranchAddress("weight_nom",&weight_nom);
    	Tout1->SetBranchAddress("weight_PU_up",&weight_PU_up);
    	Tout1->SetBranchAddress("weight_leptonsf_up",&weight_leptonsf_up);
    	Tout1->SetBranchAddress("weight_photonsf_up",&weight_photonsf_up);
    	Tout1->SetBranchAddress("weight_bSF_up",&weight_bSF_up);
    	Tout1->SetBranchAddress("weight_trig_up",&weight_trig_up);
    	Tout1->SetBranchAddress("weight_haspixsf_up",&weight_haspixsf_up);
    	Tout1->SetBranchAddress("weight_pujid_up",&weight_pujid_up);
    	Tout1->SetBranchAddress("weight_PU_dn",&weight_PU_dn);
    	Tout1->SetBranchAddress("weight_leptonsf_dn",&weight_leptonsf_dn);
    	Tout1->SetBranchAddress("weight_photonsf_dn",&weight_photonsf_dn);
    	Tout1->SetBranchAddress("weight_bSF_dn",&weight_bSF_dn);
    	Tout1->SetBranchAddress("weight_trig_dn",&weight_trig_dn);
    	Tout1->SetBranchAddress("weight_haspixsf_dn",&weight_haspixsf_dn);
    	Tout1->SetBranchAddress("weight_pujid_dn",&weight_pujid_dn);

    	Tout2->SetBranchAddress("BDT",&BDT);

	int nevt;
    	nevt=Tout1->GetEntries();

    	double sum_nom = 0.0;
    	double sum_PU_up = 0.0;
    	double sum_leptonsf_up = 0.0;
    	double sum_photonsf_up = 0.0;
    	double sum_bSF_up = 0.0;
    	double sum_trig_up = 0.0;
    	double sum_pujid_up = 0.0;
    	double sum_PU_dn = 0.0;
    	double sum_leptonsf_dn = 0.0;
    	double sum_photonsf_dn = 0.0;
    	double sum_bSF_dn = 0.0;
    	double sum_trig_dn = 0.0;
    	double sum_pujid_dn = 0.0;
    	double sum_jecup = 0.0;
    	double sum_jecdn = 0.0;
    	double sum_resup = 0.0;
    	double sum_resdn = 0.0;


	for (int i = 0; i < nevt; i++)
	{
		Tout1->GetEntry(i);
                Tout2->GetEntry(i);

		TLorentzVector b1, b2, b1JECup, b2JECup, b1JECdn, b2JECdn, b1RESup, b2RESup, b1RESdn, b2RESdn;

                b1.SetPtEtaPhiE(b1pt, b1eta, b1phi, b1e);
                b2.SetPtEtaPhiE(b2pt, b2eta, b2phi, b2e);

                b1JECup.SetPtEtaPhiE(b1pt_jecup, b1eta, b1phi, b1e_jecup);
                b2JECup.SetPtEtaPhiE(b2pt_jecup, b2eta, b2phi, b2e_jecup);

                b1JECdn.SetPtEtaPhiE(b1pt_jecdn, b1eta, b1phi, b1e_jecdn);
                b2JECdn.SetPtEtaPhiE(b2pt_jecdn, b2eta, b2phi, b2e_jecdn);

                b1RESup.SetPtEtaPhiE(b1pt_resup, b1eta, b1phi, b1e_resup);
                b2RESup.SetPtEtaPhiE(b2pt_resup, b2eta, b2phi, b2e_resup);

                b1RESdn.SetPtEtaPhiE(b1pt_resdn, b1eta, b1phi, b1e_resdn);
                b2RESdn.SetPtEtaPhiE(b2pt_resdn, b2eta, b2phi, b2e_resdn);

                float bb_invmass = (b1+b2).M();
                float bb_invmass_jecup = (b1JECup+b2JECup).M();
                float bb_invmass_jecdn = (b1JECdn+b2JECdn).M();
                float bb_invmass_resup = (b1RESup+b2RESup).M();
                float bb_invmass_resdn = (b1RESdn+b2RESdn).M();

                TLorentzVector g1, g2;

                g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                float gg_invmass = (g1+g2).M();

                float bbgg_invmass = (b1+b2+g1+g2).M();
                float bbgg_invmass_jecup = (b1JECup+b2JECup+g1+g2).M();
                float bbgg_invmass_jecdn = (b1JECdn+b2JECdn+g1+g2).M();
                float bbgg_invmass_resup = (b1RESup+b2RESup+g1+g2).M();
                float bbgg_invmass_resdn = (b1RESdn+b2RESdn+g1+g2).M();


                float chi2 = (((bbgg_invmass - 125)/sigmah1)*((bbgg_invmass - 125)/sigmah1) + ((bb_invmass - gg_invmass)/sigmabb1)*((bb_invmass - gg_invmass)/sigmabb1));
                float chi2_jecup = (((bbgg_invmass_jecup - 125)/sigmah1)*((bbgg_invmass_jecup - 125)/sigmah1) + ((bb_invmass_jecup - gg_invmass)/sigmabb1)*((bb_invmass_jecup - gg_invmass)/sigmabb1));
                float chi2_jecdn = (((bbgg_invmass_jecdn - 125)/sigmah1)*((bbgg_invmass_jecdn - 125)/sigmah1) + ((bb_invmass_jecdn - gg_invmass)/sigmabb1)*((bb_invmass_jecdn - gg_invmass)/sigmabb1));
                float chi2_resup = (((bbgg_invmass_resup - 125)/sigmah1)*((bbgg_invmass_resup - 125)/sigmah1) + ((bb_invmass_resup - gg_invmass)/sigmabb1)*((bb_invmass_resup - gg_invmass)/sigmabb1));
                float chi2_resdn = (((bbgg_invmass_resdn - 125)/sigmah1)*((bbgg_invmass_resdn - 125)/sigmah1) + ((bb_invmass_resdn - gg_invmass)/sigmabb1)*((bb_invmass_resdn - gg_invmass)/sigmabb1));

//------------------------------------------------all other syst other than JEC/JER------------------------------------------------//
		
		if (chi2 <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {

		CMS_hgg_mass = gg_invmass;
		dZ = 0.0;

		weight = (weight_nom * 1470 * 0.01)/num;
		trees[0]->Fill();
		sum_nom = sum_nom + weight;

		weight = (weight_PU_up * 1470 * 0.01 * ((1.00929*0.513698)/0.506599))/num;
		trees[1]->Fill();
		sum_PU_up = sum_PU_up + weight;

		weight = (weight_leptonsf_up * 1470 * 0.01 * ((1.01173*0.513698)/0.518231))/num;
		trees[2]->Fill();
		sum_leptonsf_up = sum_leptonsf_up + weight;

		weight = (weight_photonsf_up * 1470 * 0.01 * ((1.01403*0.513698)/0.521166))/num;
		trees[3]->Fill();
		sum_photonsf_up = sum_photonsf_up + weight;

		weight = (weight_bSF_up * 1470 * 0.01 * ((1.05766*0.513698)/0.542902))/num;
		trees[4]->Fill();
		sum_bSF_up = sum_bSF_up + weight;

		weight = (weight_trig_up * 1470 * 0.01 * ((1.00155*0.513698)/0.514517))/num;
		trees[5]->Fill();
		sum_trig_up = sum_trig_up + weight;

		weight = (weight_PU_dn * 1470 * 0.01 * ((0.990383*0.513698)/0.51938))/num;
		trees[8]->Fill();
		sum_PU_dn = sum_PU_dn + weight;

		weight = (weight_leptonsf_dn * 1470 * 0.01 * ((0.988267*0.513698)/0.509164))/num;
		trees[9]->Fill();
		sum_leptonsf_dn = sum_leptonsf_dn + weight;

		weight = (weight_photonsf_dn * 1470 * 0.01 * ((0.986071*0.513698)/0.506284))/num;
		trees[10]->Fill();
		sum_photonsf_dn = sum_photonsf_dn + weight;

		weight = (weight_bSF_dn * 1470 * 0.01 * ((0.943552*0.513698)/0.485115))/num;
		trees[11]->Fill();
		sum_bSF_dn = sum_bSF_dn + weight;

		weight = (weight_trig_dn * 1470 * 0.01 * ((0.99854*0.513698)/0.493199))/num;
		trees[12]->Fill();
		sum_trig_dn = sum_trig_dn + weight;

		}

//-----------------------------------------------------JEC/JER syst--------------------------------------------------//

		if (chi2_jecup <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                CMS_hgg_mass = gg_invmass;
                weight = (weight_nom * 1470 * 0.01 * ((1.10815*0.513698)/0.560017))/num;
                dZ = 0.0;
                trees[6]->Fill();
                sum_jecup = sum_jecup + weight; } 

                if (chi2_jecdn <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                CMS_hgg_mass = gg_invmass;
                weight = (weight_nom * 1470 * 0.01 * ((0.922888*0.513698)/0.47721))/num;
                dZ = 0.0;
                trees[13]->Fill();
                sum_jecdn = sum_jecdn + weight; } 

                if (chi2_resup <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                CMS_hgg_mass = gg_invmass;
                weight = (weight_nom * 1470 * 0.01 * ((1.00681*0.513698)/0.507936))/num;
                dZ = 0.0;
                trees[7]->Fill();
                sum_resup = sum_resup + weight; } 

                if (chi2_resdn <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                CMS_hgg_mass = gg_invmass;
                weight = (weight_nom * 1470 * 0.01 * ((0.992523*0.513698)/0.522956))/num;
                dZ = 0.0;
                trees[14]->Fill();
                sum_resdn = sum_resdn + weight; } 

	}

	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "jecup yield : " << sum_jecup << endl;
    	cout << "jecdn yield : " << sum_jecdn << endl;

    	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "jerup yield : " << sum_resup << endl;
    	cout << "jerdn yield : " << sum_resdn << endl;

    	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "PUup yield : " << sum_PU_up << endl;
   	cout << "PUdn yield : " << sum_PU_dn << endl;

    	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "leptonup yield : " << sum_leptonsf_up << endl;
    	cout << "leptondn yield : " << sum_leptonsf_dn << endl;

    	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "photonup yield : " << sum_photonsf_up << endl;
    	cout << "photondn yield : " << sum_photonsf_dn << endl;

    	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "btagup yield : " << sum_bSF_up << endl;
    	cout << "btagdn yield : " << sum_bSF_dn << endl;

    	cout << "\n";
    	cout << "nom yield : " << sum_nom << endl;
    	cout << "trigup yield : " << sum_trig_up << endl;
    	cout << "trigdn yield : " << sum_trig_dn << endl;


    	cout << "\n";
    	cout << "PU_up:       " << 1 + (sum_PU_up - sum_nom)/sum_nom << endl;
    	cout << "leptonsf_up:       " << 1 + (sum_leptonsf_up - sum_nom)/sum_nom << endl;
    	cout << "photonsf_up:    " << 1 + (sum_photonsf_up - sum_nom)/sum_nom << endl;
    	cout << "bSF_up:    " << 1 + (sum_bSF_up - sum_nom)/sum_nom << endl;
    	cout << "trig_up:   " << 1 + (sum_trig_up - sum_nom)/sum_nom << endl;
    	cout << "pujid_up:  " << 1 + (sum_pujid_up - sum_nom)/sum_nom << endl;
    	cout << "PU_dn:       " << 1 + (sum_PU_dn - sum_nom)/sum_nom << endl;
    	cout << "leptonsf_dn:       " << 1 + (sum_leptonsf_dn - sum_nom)/sum_nom << endl;
    	cout << "photonsf_dn:    " << 1 + (sum_photonsf_dn - sum_nom)/sum_nom << endl;
    	cout << "bSF_dn:    " << 1 + (sum_bSF_dn - sum_nom)/sum_nom << endl;
    	cout << "trig_dn:   " << 1 + (sum_trig_dn - sum_nom)/sum_nom << endl;
    	cout << "pujid_dn:  " << 1 + (sum_pujid_dn - sum_nom)/sum_nom << endl;
    	cout << "jecup:     " << 1 + (sum_jecup - sum_nom)/sum_nom << endl;
    	cout << "jecdn:     " << 1 + (sum_jecdn - sum_nom)/sum_nom << endl;
    	cout << "resup:     " << 1 + (sum_resup - sum_nom)/sum_nom << endl;
    	cout << "resdn:     " << 1 + (sum_resdn - sum_nom)/sum_nom << endl;
	

        f1->Close();
        f2->Close();
        fout->Write();
        fout->Close();

}
