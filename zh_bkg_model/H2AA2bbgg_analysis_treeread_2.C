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
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_bSF_up, weight_trig_up, weight_PU_dn, weight_leptonsf_dn, weight_bSF_dn, weight_trig_dn;
float sigmabb1, sigmah1;
float leppt, lepy, b1y, b2y, bb_inv_mass, pho1y, pho2y, dipho_invmass, invmassbbgg, invmassbbgg_jesup, invmassbbgg_jesdn, invmassbbgg_resoup, invmassbbgg_resodn;
float b1pt, b1eta, b1phi, b1e, b2pt, b2eta, b2phi, b2e, pho1pt, pho1eta, pho1phi, pho1e, pho2pt, pho2eta, pho2phi, pho2e;



void H2AA2bbgg_analysis_treeread_2()
{

    TFile *f1 = new TFile("/home/abala/local_bkp/analysis/root_file/ZH_2018/V1/data_v1.root");

    TTree *Tout1 = (TTree*)f1->Get("Tout");

    TFile *fout = new TFile("allData.root","RECREATE");
    TTree *T1 = new TTree("Data_13TeV_H2AA2bbgg_Tag0", "Data_13TeV_H2AA2bbgg_Tag0");


    T1->Branch("CMS_hgg_mass", &CMS_hgg_mass, "CMS_hgg_mass/F");
    T1->Branch("weight", &weight, "weight/F");
   

    Tout1->SetBranchAddress("dipho_invmass",&dipho_invmass);
    Tout1->SetBranchAddress("weight_nom",&weight_nom);
    Tout1->SetBranchAddress("bb",&bb);
    Tout1->SetBranchAddress("isEle",&isEle);
    Tout1->SetBranchAddress("isMu",&isMu);
    Tout1->SetBranchAddress("bb_inv_mass",&bb_inv_mass);
    Tout1->SetBranchAddress("invmassbbgg",&invmassbbgg);
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

    int count1 = 0;
    int count2 = 0;


    int nevt;
    nevt=Tout1->GetEntries();
   
    for (int i = 0; i < nevt; i++)
    {

		Tout1->GetEntry(i);

		if ((isEle || isMu) && bb == 1)
		{

				TLorentzVector b1, b2;

                                b1.SetPtEtaPhiE(b1pt, b1eta, b1phi, b1e);
                                b2.SetPtEtaPhiE(b2pt, b2eta, b2phi, b2e);

                                TLorentzVector g1, g2;

                                g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
                                g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);

                                float delphi_bb = fabs(b1.DeltaPhi(b2));
                                float delphi_gg = fabs(g1.DeltaPhi(g2));                             //delphi variables 
                                float delphi_bbgg = fabs((b1+b2).DeltaPhi(g1+g2));

                                sigmabb1 = 0.191*dipho_invmass + 0.726;
                                sigmah1 = -0.0175*dipho_invmass + 13.7;
                                if (sigmabb1 == 0) continue;
                                if (sigmah1 == 0) continue;

                                chi2 = (((invmassbbgg - 125)/sigmah1)*((invmassbbgg - 125)/sigmah1) + ((bb_inv_mass - dipho_invmass)/sigmabb1)*((bb_inv_mass - dipho_invmass)/sigmabb1));          // for bkg

				count1++;

                                if (chi2 <= 5) continue;
         //                       if (delphi_bb <= M_PI/2 && delphi_gg <= M_PI/2) continue;

                                count2++;

				if (dipho_invmass >= 15 && dipho_invmass <= 70) {
				CMS_hgg_mass = dipho_invmass;
				weight = 1.0;           

				T1->Fill(); }

		
		}
	
    }	


	cout << count1 << "     " << count2 << endl;
	    
	f1->cd();
        delete Tout1;
        delete f1;


        fout->cd();
        fout->Write();
        fout->Close();

	}
