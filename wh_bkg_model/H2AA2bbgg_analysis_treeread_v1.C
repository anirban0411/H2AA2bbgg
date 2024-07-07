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


float BDT, chi2, CMS_hgg_mass, weight, dZ;
bool isEle, isMu;
int bb;
float leppt, lepy, b1pt, b2pt, b1y, b2y, bb_inv_mass, pho1pt, pho2pt, pho1y, pho2y, dipho_invmass, invmassbbgg, invmassbbgg_jesup, invmassbbgg_jesdn, invmassbbgg_resoup, invmassbbgg_resodn, pho1eta, pho1phi, pho1e, pho2eta, pho2phi, pho2e;
float b1_DeepFlv, b2_DeepFlv, pho1MVA, pho2MVA, gg_dR;
float lepeta, b1eta, b2eta, Met, MetPhi, delphi_ggEt;
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_bSF_up, weight_trig_up, weight_PU_dn, weight_leptonsf_dn, weight_bSF_dn, weight_trig_dn;
float sigmabb1, sigmah1;



void H2AA2bbgg_analysis_treeread_v1()
{

    TFile *f1 = new TFile("/home/abala/local_bkp/analysis/BDT/2018/Data/data_bdt_v1_weight_nom.root");
    TFile *f2 = new TFile("/home/abala/local_bkp/analysis/BDT/2018/tree/final/BDT_UL18_MC_trees_bkg_train.root");

    TTree *Tout1 = (TTree*)f1->Get("tree");
    TTree *Tout2 = (TTree*)f2->Get("Tree");

    TFile *fout = new TFile("allData.root","RECREATE");
    TTree *T1 = new TTree("Data_13TeV_AC_Tag0", "Data_13TeV_AC_Tag0");


    T1->Branch("CMS_hgg_mass", &CMS_hgg_mass, "CMS_hgg_mass/F");
    T1->Branch("weight", &weight, "weight/F");
   

    Tout1->SetBranchAddress("pho1pt",&pho1pt);
    Tout1->SetBranchAddress("pho1eta",&pho1eta);
    Tout1->SetBranchAddress("pho1phi",&pho1phi);
    Tout1->SetBranchAddress("pho1e",&pho1e);
    Tout1->SetBranchAddress("pho2pt",&pho2pt);
    Tout1->SetBranchAddress("pho2eta",&pho2eta);
    Tout1->SetBranchAddress("pho2phi",&pho2phi);
    Tout1->SetBranchAddress("pho2e",&pho2e);
    Tout2->SetBranchAddress("BDT",&BDT);

    int count1 = 0;
    int count2 = 0;

    int nevt;
    nevt=Tout1->GetEntries();
   
    for (int i = 0; i < nevt; i++)
    {

		Tout1->GetEntry(i);
		Tout2->GetEntry(i);

				TLorentzVector g1, g2;

                                g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                                g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                                float gg_invmass = (g1+g2).M();

				count1++;

				if (BDT > 0.75) continue;

				count2++;

				if (gg_invmass >= 15 && gg_invmass <= 70) {
				CMS_hgg_mass = gg_invmass;
				weight = 1.0;           

				T1->Fill(); }

	
    }	

        cout << count1 << "	" << count2 << endl;
	    
	f1->cd();
        delete Tout1;
        delete f1;

	f2->cd();
        delete Tout2;
        delete f2;

        fout->cd();
        fout->Write();
        fout->Close();

	}
