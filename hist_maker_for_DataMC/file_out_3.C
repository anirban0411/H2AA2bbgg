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


float leppt, lepy, b1pt, b2pt, b1y, b2y, bb_inv_mass, pho1pt, pho2pt, pho1y, pho2y, dipho_invmass, invmassbbgg, invmassbbgg_jesup, invmassbbgg_jesdn, invmassbbgg_resoup, invmassbbgg_resodn, pho1eta, pho1phi, pho1e, pho2eta, pho2phi, pho2e;
float b1_DeepFlv, b2_DeepFlv, pho1MVA, pho2MVA, gg_dR;
double weight_nom, weight_PU_up, weight_leptonsf_up, weight_bSF_up, weight_trig_up, weight_PU_dn, weight_leptonsf_dn, weight_bSF_dn, weight_trig_dn;
int bb, ab;
bool isEle, isMu;


void file_out_3()
{

	TFile* f = TFile::Open("combined_2018_wh_v1.root");

	TH1F *h1 = (TH1F*) f->Get("sig_55");
	TH1F *h2 = (TH1F*) f->Get("data_obs");

	TH1F* MA_55 = new TH1F("MA_55", "MA_55", 10, 0.0, 10.0);
        MA_55->Sumw2();

	for(int i = 1; i <= h1->GetNbinsX(); i++)
	{
		float s = h1->GetBinContent(i);
		float b = h2->GetBinContent(i);
	//	cout << s << "   " << b << endl;
		float sig = 0;
		if (s+b <= 0) continue;
		sig = s/sqrt(s+b);
		MA_55->SetBinContent(i,sig);
	}

	MA_55->Draw();

}
