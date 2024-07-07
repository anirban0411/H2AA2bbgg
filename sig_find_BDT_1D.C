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
#include <ctime>


void sig_find_BDT_1D()
{

	TFile* f1 = TFile::Open("/home/abala/local_bkp/analysis/BDT/full_runII/WH202530/BDT_UL18_MC_hists_20.root");
	TFile* f2 = TFile::Open("/home/abala/local_bkp/analysis/BDT/full_runII/bkg_nom/BDT_UL18_MC_hists_bkg.root");

	TH1F *h1 = (TH1F*) f1->Get("BDT_UL18_MC_test_sig_mvaHist");
        TH1F *h2 = (TH1F*) f2->Get("BDT_UL18_MC_test_bkg_mvaHist");

	int dim = h1->GetNbinsX();
	float significance[dim];

	for(int i = 1; i < dim; i++)
	{

		float sum = 0.0;
		float bkg = 0.0;

		for(int j = i+1; j <= dim; j++)
		{
			float s = h1->GetBinContent(j);
                	float b = h2->GetBinContent(j);

			sum = sum + s;
			bkg = bkg + b;
		}

//		if ((sum+bkg) <= 7) continue;

		if ((sum+bkg) <= 0) continue;
		significance[i] = sum/sqrt(sum+bkg);

		cout << (1.0/dim)*i << "	" << significance[i] << endl;

	}

	float max = -999.0;
        int index = -1;

        for (int i = 1; i < dim; i++)
        {
                if(significance[i]>max)
                {
                        max = significance[i];
                        index = i;
                }
        }


        cout << "max sensitivity : " << max << endl;
	cout << "corresponding BDT score : " << (1.0/dim)*index << endl;
        cout << "corresponding bin no : " << index << endl;

}
