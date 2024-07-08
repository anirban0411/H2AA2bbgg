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


float leppt, lepeta, pho1pt, pho1eta, pho1phi, pho1e, pho2pt, pho2eta, pho2phi, pho2e;
float bb_inv_mass, dipho_invmass, invmassbbgg, genb_invmass;
float b1pt, b1eta, b1phi, b1e, b2pt, b2eta, b2phi, b2e, b1pt_final, b1e_final, b2pt_final, b2e_final, b1pt_scl, b1e_scl, b2pt_scl, b2e_scl;
double weight_nom;
int bb, ab;
bool isEle, isMu;
float sigmabb1, sigmah1;
float chi2;
double evt_rho, puwt, pujidwt;
int puvtx;


void H2AA2bbgg_analysis_treeread_data_mc_3()
{
    TFile *f = new TFile("WH_mA55_v6.root");
    TTree *Tout = (TTree*)f->Get("Tout");
    TFile *fout = new TFile("WH_mA55_v6_isEleMu_testv1.root","RECREATE");
    

    TH1F* h_bb_inv_mass = new TH1F("h_bb_inv_mass", "h_bb_inv_mass", 500, 0.0, 1000.0);
    h_bb_inv_mass->Sumw2();

    TH1F* h_bb_inv_mass_scl = new TH1F("h_bb_inv_mass_scl", "h_bb_inv_mass_scl", 500, 0.0, 1000.0);
    h_bb_inv_mass_scl->Sumw2();

    TH1F* h_bb_inv_mass_sclsmr = new TH1F("h_bb_inv_mass_sclsmr", "h_bb_inv_mass_sclsmr", 500, 0.0, 1000.0);
    h_bb_inv_mass_sclsmr->Sumw2();

    TH1F* h_bb_inv_mass_gen = new TH1F("h_bb_inv_mass_gen", "h_bb_inv_mass_gen", 500, 0.0, 1000.0);
    h_bb_inv_mass_gen->Sumw2();

    TH2F* h_rho_npu = new TH2F("h_rho_npu", "h_rho_npu", 50, 0.0, 100.0, 50, 0.0, 100.0);
    h_rho_npu->Sumw2();

    TH2F* h_rho_puwt = new TH2F("h_rho_puwt", "h_rho_puwt", 50, 0.0, 100.0, 50, 0.0, 2.0);
    h_rho_puwt->Sumw2();

    TH1F* cut_flow = new TH1F("cut_flow", "cut_flow", 10, 0.0, 10.0);
    cut_flow->Sumw2();


    Tout->SetBranchAddress("b1pt",&b1pt);
    Tout->SetBranchAddress("b1eta",&b1eta);
    Tout->SetBranchAddress("b1phi",&b1phi);
    Tout->SetBranchAddress("b1e",&b1e);
    Tout->SetBranchAddress("b2pt",&b2pt);
    Tout->SetBranchAddress("b2eta",&b2eta);
    Tout->SetBranchAddress("b2phi",&b2phi);
    Tout->SetBranchAddress("b2e",&b2e);
    Tout->SetBranchAddress("b1pt_final",&b1pt_final);
    Tout->SetBranchAddress("b1e_final",&b1e_final);
    Tout->SetBranchAddress("b2pt_final",&b2pt_final);
    Tout->SetBranchAddress("b2e_final",&b2e_final);
    Tout->SetBranchAddress("b1pt_scl",&b1pt_scl);
    Tout->SetBranchAddress("b1e_scl",&b1e_scl);
    Tout->SetBranchAddress("b2pt_scl",&b2pt_scl);
    Tout->SetBranchAddress("b2e_scl",&b2e_scl);
    Tout->SetBranchAddress("genb_invmass",&genb_invmass);
    Tout->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout->SetBranchAddress("ab",&ab);
    Tout->SetBranchAddress("bb",&bb);
    Tout->SetBranchAddress("isEle",&isEle);
    Tout->SetBranchAddress("isMu",&isMu);
    Tout->SetBranchAddress("isMu",&isMu);
    Tout->SetBranchAddress("weight_nom",&weight_nom);
    Tout->SetBranchAddress("evt_rho",&evt_rho);
    Tout->SetBranchAddress("puwt",&puwt);
    Tout->SetBranchAddress("puvtx",&puvtx);


    
    int nevt;
    nevt=Tout->GetEntries();
    
    for (int i = 0; i < nevt; i++)
    {
		Tout->GetEntry(i);

		if ((isEle || isMu) && bb == 1)
		{
	
				if (weight_nom <= 0) continue;


				TLorentzVector b1, b2, b1_f, b2_f, b1_sclsmr, b2_sclsmr;

                                b1.SetPtEtaPhiE(b1pt, b1eta, b1phi, b1e);
                                b2.SetPtEtaPhiE(b2pt, b2eta, b2phi, b2e);

				b1_f.SetPtEtaPhiE(b1pt_scl, b1eta, b1phi, b1e_scl);
                                b2_f.SetPtEtaPhiE(b2pt_scl, b2eta, b2phi, b2e_scl);

				b1_sclsmr.SetPtEtaPhiE(b1pt_final, b1eta, b1phi, b1e_final);
                                b2_sclsmr.SetPtEtaPhiE(b2pt_final, b2eta, b2phi, b2e_final);

				
				h_bb_inv_mass->Fill((b1+b2).M());
				h_bb_inv_mass_scl->Fill((b1_f+b2_f).M());
				h_bb_inv_mass_sclsmr->Fill((b1_sclsmr+b2_sclsmr).M());
				h_bb_inv_mass_gen->Fill(genb_invmass);
				h_rho_npu->Fill(evt_rho, puvtx);
				h_rho_puwt->Fill(evt_rho, puwt);

		
		}
	
    }
	

    cout << "a" << endl;

	    
	f->cd();
        delete Tout;
        delete f;

        fout->cd();
        fout->Write();
        fout->Close();
	}
