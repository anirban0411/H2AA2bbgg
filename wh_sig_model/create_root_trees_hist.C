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
float sigmabb1, sigmah1;
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_photonsf_up, weight_bSF_up, weight_trig_up, weight_haspixsf_up, weight_pujid_up, weight_PU_dn, weight_leptonsf_dn, weight_photonsf_dn, weight_prefiring_dn, weight_bSF_dn, weight_trig_dn, weight_haspixsf_dn, weight_pujid_dn;


void create_root_trees_hist()
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


    TH1F* CMS_hgg_mass = new TH1F("CMS_hgg_mass", "CMS_hgg_mass", 100, 0.0, 200.0);
    CMS_hgg_mass->Sumw2();

    TH1F* CMS_hgg_mass_PUUp01sigma = new TH1F("CMS_hgg_mass_PUUp01sigma", "CMS_hgg_mass_PUUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_PUUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_leptonIDUp01sigma = new TH1F("CMS_hgg_mass_leptonIDUp01sigma", "CMS_hgg_mass_leptonIDUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_leptonIDUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_photonIDUp01sigma = new TH1F("CMS_hgg_mass_photonIDUp01sigma", "CMS_hgg_mass_photonIDUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_photonIDUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_btagUp01sigma = new TH1F("CMS_hgg_mass_btagUp01sigma", "CMS_hgg_mass_btagUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_btagUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_TrigUp01sigma = new TH1F("CMS_hgg_mass_TrigUp01sigma", "CMS_hgg_mass_TrigUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_TrigUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_JECUp01sigma = new TH1F("CMS_hgg_mass_JECUp01sigma", "CMS_hgg_mass_JECUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_JECUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_JERUp01sigma = new TH1F("CMS_hgg_mass_JERUp01sigma", "CMS_hgg_mass_JERUp01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_JERUp01sigma->Sumw2();

    TH1F* CMS_hgg_mass_PUDown01sigma = new TH1F("CMS_hgg_mass_PUDown01sigma", "CMS_hgg_mass_PUDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_PUDown01sigma->Sumw2();

    TH1F* CMS_hgg_mass_leptonIDDown01sigma = new TH1F("CMS_hgg_mass_leptonIDDown01sigma", "CMS_hgg_mass_leptonIDDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_leptonIDDown01sigma->Sumw2();

    TH1F* CMS_hgg_mass_photonIDDown01sigma = new TH1F("CMS_hgg_mass_photonIDDown01sigma", "CMS_hgg_mass_photonIDDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_photonIDDown01sigma->Sumw2();

    TH1F* CMS_hgg_mass_btagDown01sigma = new TH1F("CMS_hgg_mass_btagDown01sigma", "CMS_hgg_mass_btagDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_btagDown01sigma->Sumw2();

    TH1F* CMS_hgg_mass_TrigDown01sigma = new TH1F("CMS_hgg_mass_TrigDown01sigma", "CMS_hgg_mass_TrigDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_TrigDown01sigma->Sumw2();

    TH1F* CMS_hgg_mass_JECDown01sigma = new TH1F("CMS_hgg_mass_JECDown01sigma", "CMS_hgg_mass_JECDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_JECDown01sigma->Sumw2();

    TH1F* CMS_hgg_mass_JERDown01sigma = new TH1F("CMS_hgg_mass_JERDown01sigma", "CMS_hgg_mass_JERDown01sigma", 100, 0.0, 200.0);
    CMS_hgg_mass_JERDown01sigma->Sumw2();

    

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


//------------------------------------------------calculation of all weights except JEC/JER--------------------------------------------------------//
                                
				if (chi2 <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {//start of BDT loop//

				CMS_hgg_mass->Fill(gg_invmass,WH_mA60*weight_nom);
				CMS_hgg_mass_PUUp01sigma->Fill(gg_invmass,WH_mA60*weight_PU_up);
				CMS_hgg_mass_leptonIDUp01sigma->Fill(gg_invmass,WH_mA60*weight_leptonsf_up);
				CMS_hgg_mass_photonIDUp01sigma->Fill(gg_invmass,WH_mA60*weight_photonsf_up);
				CMS_hgg_mass_btagUp01sigma->Fill(gg_invmass,WH_mA60*weight_bSF_up);
				CMS_hgg_mass_TrigUp01sigma->Fill(gg_invmass,WH_mA60*weight_trig_up);
				CMS_hgg_mass_PUDown01sigma->Fill(gg_invmass,WH_mA60*weight_PU_dn);
				CMS_hgg_mass_leptonIDDown01sigma->Fill(gg_invmass,WH_mA60*weight_leptonsf_dn);
				CMS_hgg_mass_photonIDDown01sigma->Fill(gg_invmass,WH_mA60*weight_photonsf_dn);
				CMS_hgg_mass_btagDown01sigma->Fill(gg_invmass,WH_mA60*weight_bSF_dn);
				CMS_hgg_mass_TrigDown01sigma->Fill(gg_invmass,WH_mA60*weight_trig_dn);
			
				}// end of BDT loop


//--------------------------------------------------------calculation of JEC/JER----------------------------------------------------//
				
				if (chi2_jecup <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
				CMS_hgg_mass_JECUp01sigma->Fill(gg_invmass,WH_mA60*weight_nom);
				}

				if (chi2_jecdn <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass_JECDown01sigma->Fill(gg_invmass,WH_mA60*weight_nom);
                                }

				if (chi2_resup <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass_JERUp01sigma->Fill(gg_invmass,WH_mA60*weight_nom);
                                }

				if (chi2_resdn <= 1 && BDT >= 0.75 && gg_invmass >= 15 && gg_invmass <= 70) {
                                CMS_hgg_mass_JERDown01sigma->Fill(gg_invmass,WH_mA60*weight_nom);
                                }
	
    }	


    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "jecup yield : " << CMS_hgg_mass_JECUp01sigma->Integral() << endl;
    cout << "jecdn yield : " << CMS_hgg_mass_JECDown01sigma->Integral() << endl;

    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "jerup yield : " << CMS_hgg_mass_JERUp01sigma->Integral() << endl;
    cout << "jerdn yield : " << CMS_hgg_mass_JERDown01sigma->Integral() << endl;

    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "PUup yield : " << CMS_hgg_mass_PUUp01sigma->Integral() << endl;
    cout << "PUdn yield : " << CMS_hgg_mass_PUDown01sigma->Integral() << endl;

    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "leptonup yield : " << CMS_hgg_mass_leptonIDUp01sigma->Integral() << endl;
    cout << "leptondn yield : " << CMS_hgg_mass_leptonIDDown01sigma->Integral() << endl;

    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "photonup yield : " << CMS_hgg_mass_photonIDUp01sigma->Integral() << endl;
    cout << "photondn yield : " << CMS_hgg_mass_photonIDDown01sigma->Integral() << endl;

    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "btagup yield : " << CMS_hgg_mass_btagUp01sigma->Integral() << endl;
    cout << "btagdn yield : " << CMS_hgg_mass_btagDown01sigma->Integral() << endl;

    cout << "\n";
    cout << "nom yield : " << CMS_hgg_mass->Integral() << endl;
    cout << "trigup yield : " << CMS_hgg_mass_TrigUp01sigma->Integral() << endl;
    cout << "trigdn yield : " << CMS_hgg_mass_TrigDown01sigma->Integral() << endl;


    cout << "\n";
    cout << "PU_up:       " << 1 + (CMS_hgg_mass_PUUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "leptonsf_up:	" << 1 + (CMS_hgg_mass_leptonIDUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "photonsf_up:    " << 1 + (CMS_hgg_mass_photonIDUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "bSF_up:    " << 1 + (CMS_hgg_mass_btagUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "trig_up:	" << 1 + (CMS_hgg_mass_TrigUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "PU_dn:       " << 1 + (CMS_hgg_mass_PUDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "leptonsf_dn:       " << 1 + (CMS_hgg_mass_leptonIDDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "photonsf_dn:    " << 1 + (CMS_hgg_mass_photonIDDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "bSF_dn:	" << 1 + (CMS_hgg_mass_btagDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "trig_dn:	" << 1 + (CMS_hgg_mass_TrigDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "jecup:	" << 1 + (CMS_hgg_mass_JECUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "jecdn:	" << 1 + (CMS_hgg_mass_JECDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "resup:	" << 1 + (CMS_hgg_mass_JERUp01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;
    cout << "resdn:	" << 1 + (CMS_hgg_mass_JERDown01sigma->Integral() - CMS_hgg_mass->Integral())/CMS_hgg_mass->Integral() << endl;

	    
	f1->cd();
        delete Tout1;
        delete f1;

	f2->cd();
        delete Tout2;
        delete f2;

	}
