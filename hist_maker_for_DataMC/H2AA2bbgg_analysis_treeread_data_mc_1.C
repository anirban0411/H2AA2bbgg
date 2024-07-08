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
float bb_inv_mass, dipho_invmass, invmassbbgg;
float b1pt, b1eta, b1phi, b1e, b2pt, b2eta, b2phi, b2e, b1pt_final, b1e_final, b2pt_final, b2e_final, b1pt_scl, b1e_scl, b2pt_scl, b2e_scl;
double weight_nom;
int bb, ab;
bool isEle, isMu;
float sigmabb1, sigmah1;
float chi2;
int bjet_no;
int jet_no;


void H2AA2bbgg_analysis_treeread_data_mc_1()
{
    TFile *f = new TFile("TTbar.root");
    TTree *Tout = (TTree*)f->Get("Tout");
    TFile *fout = new TFile("TTbar_isEleMu.root","RECREATE");
    

    TH1F* h_leppt = new TH1F("h_leppt", "h_leppt", 30, 0.0, 600.0);
    h_leppt->Sumw2();

    TH1F* h_lepeta = new TH1F("h_lepeta", "h_lepeta", 30, -3.0, 3.0);
    h_lepeta->Sumw2();

    TH1F* h_b1pt = new TH1F("h_b1pt", "h_b1pt", 30, 0.0, 600.0);
    h_b1pt->Sumw2();

    TH1F* h_b1eta = new TH1F("h_b1eta", "h_b1eta", 30, -3.0, 3.0);
    h_b1eta->Sumw2();

    TH1F* h_b2pt = new TH1F("h_b2pt", "h_b2pt", 30, 0.0, 600.0);
    h_b2pt->Sumw2();

    TH1F* h_b2eta = new TH1F("h_b2eta", "h_b2eta", 30, -3.0, 3.0);
    h_b2eta->Sumw2();

    TH1F* h_pho1pt = new TH1F("h_pho1pt", "h_pho1pt", 30, 0.0, 600.0);
    h_pho1pt->Sumw2();

    TH1F* h_pho1eta = new TH1F("h_pho1eta", "h_pho1eta", 30, -3.0, 3.0);
    h_pho1eta->Sumw2();

    TH1F* h_pho2pt = new TH1F("h_pho2pt", "h_pho2pt", 30, 0.0, 600.0);
    h_pho2pt->Sumw2();

    TH1F* h_pho2eta = new TH1F("h_pho2eta", "h_pho2eta", 30, -3.0, 3.0);
    h_pho2eta->Sumw2();

    TH1F* h_dipho_invmass = new TH1F("h_dipho_invmass", "h_dipho_invmass", 100, 0.0, 1000.0);
    h_dipho_invmass->Sumw2();

    TH1F* h_invmassbbgg = new TH1F("h_invmassbbgg", "h_invmassbbgg", 100, 0.0, 2000.0);
    h_invmassbbgg->Sumw2();
   
    TH1F* h_bb_inv_mass = new TH1F("h_bb_inv_mass", "h_bb_inv_mass", 50, 0.0, 1000.0);
    h_bb_inv_mass->Sumw2();

    TH1F* h_bb_inv_mass_corr = new TH1F("h_bb_inv_mass_corr", "h_bb_inv_mass_corr", 50, 0.0, 1000.0);
    h_bb_inv_mass_corr->Sumw2();

    TH2F* h_nbjets_ptb1 = new TH2F("h_nbjets_ptb1", "h_nbjets_ptb1", 30, 0.0, 600.0, 10, 0.0, 10.0);
    h_nbjets_ptb1->Sumw2();

    TH2F* h_nbjets_ptb2 = new TH2F("h_nbjets_ptb2", "h_nbjets_ptb2", 30, 0.0, 600.0, 10, 0.0, 10.0);
    h_nbjets_ptb2->Sumw2();

    TH1F* h_bjet_no = new TH1F("h_bjet_no", "h_bjet_no", 10, 0.0, 10.0);
    h_bjet_no->Sumw2();

    TH1F* cut_flow = new TH1F("cut_flow", "cut_flow", 10, 0.0, 10.0);
    cut_flow->Sumw2();


    Tout->SetBranchAddress("leppt",&leppt);
    Tout->SetBranchAddress("lepeta",&lepeta);
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
    Tout->SetBranchAddress("pho1pt",&pho1pt);
    Tout->SetBranchAddress("pho1eta",&pho1eta);
    Tout->SetBranchAddress("pho1phi",&pho1phi);
    Tout->SetBranchAddress("pho1e",&pho1e);
    Tout->SetBranchAddress("pho2pt",&pho2pt);
    Tout->SetBranchAddress("pho2eta",&pho2eta);
    Tout->SetBranchAddress("pho2phi",&pho2phi);
    Tout->SetBranchAddress("pho2e",&pho2e);
    Tout->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout->SetBranchAddress("invmassbbgg",&invmassbbgg);
    Tout->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout->SetBranchAddress("ab",&ab);
    Tout->SetBranchAddress("bb",&bb);
    Tout->SetBranchAddress("isEle",&isEle);
    Tout->SetBranchAddress("isMu",&isMu);
    Tout->SetBranchAddress("isMu",&isMu);
    Tout->SetBranchAddress("weight_nom",&weight_nom);
    Tout->SetBranchAddress("jet_no",&jet_no);


    
    int nevt;
    nevt=Tout->GetEntries();
    
    for (int i = 0; i < nevt; i++)
    {
		Tout->GetEntry(i);

		if ((isEle || isMu) && bb == 1)
		{
	
				if (weight_nom <= 0) continue;

				cut_flow->Fill(0);

				sigmabb1 = 0.20*dipho_invmass + 1.69;
                                sigmah1 = -0.05*dipho_invmass + 16.34;
                                if (sigmabb1 == 0) continue;
                                if (sigmah1 == 0) continue;

                                chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg
                  //              chi2 = (((invmassbbgg - 125)/16.04)*((invmassbbgg - 125)/16.04) + ((bb_inv_mass - dipho_invmass)/5.69)*((bb_inv_mass - dipho_invmass)/5.69));                // for 20 GeV
                  //              chi2 = (((invmassbbgg - 125)/13.98)*((invmassbbgg - 125)/13.98) + ((bb_inv_mass - dipho_invmass)/12.61)*((bb_inv_mass - dipho_invmass)/12.61));                // for 55 GeV

		
				if ((b1pt < 20) || (b2pt < 20)) continue;
				TLorentzVector b1, b2, b1_f, b2_f;

                                b1.SetPtEtaPhiE(b1pt, b1eta, b1phi, b1e);
                                b2.SetPtEtaPhiE(b2pt, b2eta, b2phi, b2e);

				b1_f.SetPtEtaPhiE(b1pt, b1eta, b1phi, b1e);
                                b2_f.SetPtEtaPhiE(b2pt, b2eta, b2phi, b2e);

				TLorentzVector g1, g2;

				g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
				g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

				if (chi2 <= 1) continue;

				cut_flow->Fill(1);
	
				bjet_no = 0;
                                if (jet_no <= 12)
                                {
                                        bjet_no = 2;
                                }
                                if (jet_no > 12 && jet_no <= 18)
                                {
                                        bjet_no = 3;
                                }
                                if (jet_no > 18)
                                {
                                        bjet_no = 4;
                                }

				h_bb_inv_mass->Fill((b1+b2).M(), weight_nom);
				h_bb_inv_mass_corr->Fill((b1_f+b2_f).M(), weight_nom);
				h_leppt->Fill(leppt, weight_nom);
				h_lepeta->Fill(lepeta, weight_nom);
				h_b1pt->Fill(b1_f.Pt(), weight_nom);
				h_b1eta->Fill(b1_f.Eta(), weight_nom);
				h_b2pt->Fill(b2_f.Pt(), weight_nom);
                                h_b2eta->Fill(b2_f.Eta(), weight_nom);
				h_pho1pt->Fill(g1.Pt(), weight_nom);
				h_pho1eta->Fill(g1.Eta(), weight_nom);
				h_pho2pt->Fill(g2.Pt(), weight_nom);
                                h_pho2eta->Fill(g2.Eta(), weight_nom);
				h_dipho_invmass->Fill((g1+g2).M(), weight_nom);
				h_invmassbbgg->Fill((b1_f+b2_f+g1+g2).M(), weight_nom);
				h_bjet_no->Fill(bjet_no);
                                h_nbjets_ptb1->Fill(b1_f.Pt(), bjet_no);
                                h_nbjets_ptb2->Fill(b2_f.Pt(), bjet_no);

		
		}
	
    }
	

	    
	f->cd();
        delete Tout;
        delete f;

        fout->cd();
        fout->Write();
        fout->Close();
	}
