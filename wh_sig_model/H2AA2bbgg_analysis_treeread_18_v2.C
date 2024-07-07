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


float BDT, CMS_hgg_mass, weight, dZ;
bool isEle, isMu;
int bb;
float bb_inv_mass, dipho_invmass, invmassbbgg;
float b1pt, b1eta, b1phi, b1e, b2pt, b2eta, b2phi, b2e, pho1pt, pho1eta, pho1phi, pho1e, pho2pt, pho2eta, pho2phi, pho2e;
float b1pt_jecup, b1e_jecup, b2pt_jecup, b2e_jecup, b1pt_jecdn , b1e_jecdn, b2pt_jecdn, b2e_jecdn, b1pt_resup, b1e_resup, b2pt_resup, b2e_resup, b1pt_resdn, b1e_resdn, b2pt_resdn, b2e_resdn;
float sigmabb1, sigmah1, sigmabb1_jecup, sigmah1_jecup, sigmabb1_jecdn, sigmah1_jecdn, sigmabb1_jerup, sigmah1_jerup, sigmabb1_jerdn, sigmah1_jerdn;
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_photonsf_up, weight_bSF_up, weight_trig_up, weight_haspixsf_up, weight_pujid_up, weight_PU_dn, weight_leptonsf_dn, weight_photonsf_dn, weight_prefiring_dn, weight_bSF_dn, weight_trig_dn, weight_haspixsf_dn, weight_pujid_dn;


void H2AA2bbgg_analysis_treeread_18_v2()
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

	double num = 906844;

	sigmabb1 = 11.193;
        sigmah1 = 12.982;

	sigmabb1_jecup = 12.498;
        sigmah1_jecup = 12.674;

	sigmabb1_jecdn = 12.219;
        sigmah1_jecdn = 12.341;

	sigmabb1_jerup = 12.383;
        sigmah1_jerup = 12.614;

	sigmabb1_jerdn = 12.179;
        sigmah1_jerdn = 12.403;

	TString str[19] = {"weight_nom", "weight_PU_up", "weight_leptonsf_up", "weight_photonsf_up", "weight_bSF_up", "weight_trig_up", "weight_haspixsf_up", "weight_pujid_up", "weight_jec_up", "weight_jec_dn", "weight_jer_up", "weight_jer_dn", "weight_PU_dn", "weight_leptonsf_dn", "weight_photonsf_dn", "weight_bSF_dn", "weight_trig_dn", "weight_haspixsf_dn", "weight_pujid_dn"};


    TFile *f1 = new TFile("/home/abala/local_bkp/analysis/BDT/2018/WH55/WH_mA55_vf_weight_nom.root");
    TFile *f2 = new TFile("/home/abala/local_bkp/analysis/BDT/2018/tree/final/BDT_UL18_MC_trees_sig_train_55.root");

    TTree *Tout1 = (TTree*)f1->Get("tree");
    TTree *Tout2 = (TTree*)f2->Get("Tree");


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
                    		float chi2_jecup = (((bbgg_invmass_jecup - 125)/sigmah1_jecup)*((bbgg_invmass_jecup - 125)/sigmah1_jecup) + ((bb_invmass_jecup - gg_invmass)/sigmabb1_jecup)*((bb_invmass_jecup - gg_invmass)/sigmabb1_jecup));
                    		float chi2_jecdn = (((bbgg_invmass_jecdn - 125)/sigmah1_jecdn)*((bbgg_invmass_jecdn - 125)/sigmah1_jecdn) + ((bb_invmass_jecdn - gg_invmass)/sigmabb1_jecdn)*((bb_invmass_jecdn - gg_invmass)/sigmabb1_jecdn));
                    		float chi2_resup = (((bbgg_invmass_resup - 125)/sigmah1_jerup)*((bbgg_invmass_resup - 125)/sigmah1_jerup) + ((bb_invmass_resup - gg_invmass)/sigmabb1_jerup)*((bb_invmass_resup - gg_invmass)/sigmabb1_jerup));
                    		float chi2_resdn = (((bbgg_invmass_resdn - 125)/sigmah1_jerdn)*((bbgg_invmass_resdn - 125)/sigmah1_jerdn) + ((bb_invmass_resdn - gg_invmass)/sigmabb1_jerdn)*((bb_invmass_resdn - gg_invmass)/sigmabb1_jerdn));


				if (chi2 <= 1 && BDT >= 0.75) {

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_nom * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_nom = sum_nom + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_PU_up * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_PU_up = sum_PU_up + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_leptonsf_up * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_leptonsf_up = sum_leptonsf_up + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_photonsf_up * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_photonsf_up = sum_photonsf_up + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_bSF_up * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_bSF_up = sum_bSF_up + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_trig_up * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_trig_up = sum_trig_up + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_pujid_up * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_pujid_up = sum_pujid_up + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_PU_dn * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_PU_dn = sum_PU_dn + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_leptonsf_dn * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_leptonsf_dn = sum_leptonsf_dn + weight; }

                                if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_photonsf_dn * 1470 * 0.01)/num;
                                dZ = 0.0;
                                sum_photonsf_dn = sum_photonsf_dn + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_bSF_dn * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_bSF_dn = sum_bSF_dn + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_trig_dn * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_trig_dn = sum_trig_dn + weight; }

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_pujid_dn * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_pujid_dn = sum_pujid_dn + weight; } }

				if (chi2_jecup <= 1 && BDT >= 0.75) {

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_nom * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_jecup = sum_jecup + weight; } }

				if (chi2_jecdn <= 1 && BDT >= 0.75) {

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_nom * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_jecdn = sum_jecdn + weight; } }

				if (chi2_resup <= 1 && BDT >= 0.75) {

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_nom * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_resup = sum_resup + weight; } }

				if (chi2_resdn <= 1 && BDT >= 0.75) {

				if (gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass = gg_invmass;
                                weight = (weight_nom * 1470 * 0.01)/num;
                                dZ = 0.0;
				sum_resdn = sum_resdn + weight; } }

	
    }	


    cout << "PU_up:       " << 1 + (sum_PU_up - sum_nom)/sum_nom << endl;
    cout << "leptonsf_up:	" << 1 + (sum_leptonsf_up - sum_nom)/sum_nom << endl;
    cout << "photonsf_up:    " << 1 + (sum_photonsf_up - sum_nom)/sum_nom << endl;
    cout << "bSF_up:    " << 1 + (sum_bSF_up - sum_nom)/sum_nom << endl;
    cout << "trig_up:	" << 1 + (sum_trig_up - sum_nom)/sum_nom << endl;
    cout << "pujid_up:	" << 1 + (sum_pujid_up - sum_nom)/sum_nom << endl;
    cout << "PU_dn:       " << 1 + (sum_PU_dn - sum_nom)/sum_nom << endl;
    cout << "leptonsf_dn:       " << 1 + (sum_leptonsf_dn - sum_nom)/sum_nom << endl;
    cout << "photonsf_dn:    " << 1 + (sum_photonsf_dn - sum_nom)/sum_nom << endl;
    cout << "bSF_dn:	" << 1 + (sum_bSF_dn - sum_nom)/sum_nom << endl;
    cout << "trig_dn:	" << 1 + (sum_trig_dn - sum_nom)/sum_nom << endl;
    cout << "pujid_dn:	" << 1 + (sum_pujid_dn - sum_nom)/sum_nom << endl;
    cout << "jecup:	" << 1 + (sum_jecup - sum_nom)/sum_nom << endl;
    cout << "jecdn:	" << 1 + (sum_jecdn - sum_nom)/sum_nom << endl;
    cout << "resup:	" << 1 + (sum_resup - sum_nom)/sum_nom << endl;
    cout << "resdn:	" << 1 + (sum_resdn - sum_nom)/sum_nom << endl;

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;	
    cout << "jecup yield : " << (sum_jecup*59.73) << endl;
    cout << "jecdn yield : " << (sum_jecdn*59.73) << endl;   

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;
    cout << "jerup yield : " << (1.00821*sum_nom*59.73) << endl;
    cout << "jerdn yield : " << (0.989986*sum_nom*59.73) << endl;

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;
    cout << "PUup yield : " << (1.01116*sum_nom*59.73) << endl;
    cout << "PUdn yield : " << (0.989705*sum_nom*59.73) << endl;

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;
    cout << "leptonup yield : " << (1.00928*sum_nom*59.73) << endl;
    cout << "leptondn yield : " << (0.99072*sum_nom*59.73) << endl;

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;
    cout << "photonup yield : " << (1.01371*sum_nom*59.73) << endl;
    cout << "photondn yield : " << (0.986389*sum_nom*59.73) << endl;

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;
    cout << "btagup yield : " << (1.05983*sum_nom*59.73) << endl;
    cout << "btagdn yield : " << (0.941548*sum_nom*59.73) << endl;

    cout << "\n";
    cout << "nom yield : " << (sum_nom*59.73) << endl;
    cout << "trigup yield : " << (1.00198*sum_nom*59.73) << endl;
    cout << "trigdn yield : " << (0.99675*sum_nom*59.73) << endl;


 
	f1->cd();
        delete Tout1;
        delete f1;

	f2->cd();
        delete Tout2;
        delete f2;

	}
