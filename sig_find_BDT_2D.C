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


void sig_find_BDT_2D()
{

	TFile* f1 = TFile::Open("full_runII/WH60/BDT_UL18_MC_hists_60.root");
	TFile* f2 = TFile::Open("full_runII/bkg_nom/BDT_UL18_MC_hists_bkg.root");

	TH1F *h1 = (TH1F*) f1->Get("BDT_UL18_MC_test_sig_mvaHist");
        TH1F *h2 = (TH1F*) f2->Get("BDT_UL18_MC_test_bkg_mvaHist");

	TFile *fout = new TFile("twoD_sig_hist_60.root","RECREATE");
	TH2F* h_sig = new TH2F("h_sig", "h_sig", 40, 0.0, 1.0, 40, 0.0, 1.0);
        h_sig->Sumw2();

	int dim = h1->GetNbinsX();
	float significance[dim], significance_1[dim*dim], sig_sum[dim*dim], wt_evt_1[dim], wt_evt_2[dim*dim];


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

                if ((sum+bkg) <= 7) continue;

//		cout << (sum+bkg) << endl;

		wt_evt_1[i] = (sum+bkg);

                significance[i] = sum/sqrt(sum+bkg);

                for (int k = 1; k < i; k++)
                {
                        float sum_1 = 0.0;
                        float bkg_1 = 0.0;

                        for(int l = k+1; l <= i; l++)
                        {
                                float s_1 = h1->GetBinContent(l);
                                float b_1 = h2->GetBinContent(l);

                                sum_1 = sum_1 + s_1;
                                bkg_1 = bkg_1 + b_1;
                        }

                        if ((sum_1+bkg_1) <= 7) continue;

//			cout << (sum_1+bkg_1) << endl;

			wt_evt_2[k] = (sum_1+bkg_1);

                        significance_1[k] = sum_1/sqrt(sum_1+bkg_1);
                }

	}

	int indx = 0;
        for (int i = 1; i < dim; i++)
        {
		for (int k = 1; k < i; k++) 
                {
                        float func1 = sqrt(significance_1[k]*significance_1[k] + significance[i]*significance[i]);
			sig_sum[indx] = func1;
                        indx++;
                }

        }


        float max = -999.0;

        for (int i = 0; i < indx; i++)
        {
                if(sig_sum[i]>max)
                {
                        max = sig_sum[i];
                }

        }


	int a0 = 0;
        int b0 = 0;

        for (int i = 1; i < dim; i++)
        {
                for (int k = 1; k < i; k++)
                {
                        float func1 = sqrt(significance_1[k]*significance_1[k] + significance[i]*significance[i]);

			h_sig->Fill((1.0/dim)*k, (1.0/dim)*i, func1);

//			cout << (1.0/dim)*k << "	" << (1.0/dim)*i << "	" << func1 << endl;

                        if (func1 == max)
                        {
                                a0 = k;
                                b0 = i;
				cout << wt_evt_2[k] << "	" << wt_evt_1[i] << endl;
                                break;
                        }

                }
        }


        cout << "corresponding bin no : " << a0 << "   " << b0 << endl;
        cout << "corresponding BDT score : " << (1.0/dim)*a0 << "        " << (1.0/dim)*b0 << endl;
        cout << "max sensitivity : " << max << endl;

	fout->cd();
        fout->Write();
        fout->Close();

}
