#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TLorentzVector.h"

using namespace std;

// ═══════════════════════════════════════════════════════════════════════════
// Utility: wrap phi into (−π, π]
// ═══════════════════════════════════════════════════════════════════════════
double PhiInRange(const double& phi)
{
    double phiout = phi;
    if (phiout >  2*M_PI || phiout < -2*M_PI) phiout = fmod(phiout, 2*M_PI);
    if (phiout <= -M_PI) phiout += 2*M_PI;
    else if (phiout >   M_PI) phiout -= 2*M_PI;
    return phiout;
}

double delta2R(double eta1, double phi1, double eta2, double phi2)
{
    return sqrt(pow(eta1-eta2, 2) + pow(PhiInRange(phi1-phi2), 2));
}

// ═══════════════════════════════════════════════════════════════════════════
// Global branch variables
// ═══════════════════════════════════════════════════════════════════════════

// Core event variables
float  BDT, CMS_hgg_mass, weight, dZ;
bool   isEle, isMu;
int    bb;
float  bb_inv_mass, dipho_invmass, invmassbbgg;

// Jet & photon 4-vector components (nominal)
float  b1pt,   b1eta,   b1phi,   b1e;
float  b2pt,   b2eta,   b2phi,   b2e;
float  pho1pt, pho1eta, pho1phi, pho1e;
float  pho2pt, pho2eta, pho2phi, pho2e;

// JEC / JER jet pt & energy variations
float  b1pt_jecup, b1e_jecup, b2pt_jecup, b2e_jecup;
float  b1pt_jecdn, b1e_jecdn, b2pt_jecdn, b2e_jecdn;
float  b1pt_resup, b1e_resup, b2pt_resup, b2e_resup;
float  b1pt_resdn, b1e_resdn, b2pt_resdn, b2e_resdn;

// Source-specific JEC uncertainties
float  b1pt_jecup_AbsoluteStat,  b1e_jecup_AbsoluteStat,  b2pt_jecup_AbsoluteStat,  b2e_jecup_AbsoluteStat;
float  b1pt_jecdn_AbsoluteStat,  b1e_jecdn_AbsoluteStat,  b2pt_jecdn_AbsoluteStat,  b2e_jecdn_AbsoluteStat;
float  b1pt_jecup_FlavorQCD,     b1e_jecup_FlavorQCD,     b2pt_jecup_FlavorQCD,     b2e_jecup_FlavorQCD;
float  b1pt_jecdn_FlavorQCD,     b1e_jecdn_FlavorQCD,     b2pt_jecdn_FlavorQCD,     b2e_jecdn_FlavorQCD;
float  b1pt_jecup_PileUpPtBB,    b1e_jecup_PileUpPtBB,    b2pt_jecup_PileUpPtBB,    b2e_jecup_PileUpPtBB;
float  b1pt_jecdn_PileUpPtBB,    b1e_jecdn_PileUpPtBB,    b2pt_jecdn_PileUpPtBB,    b2e_jecdn_PileUpPtBB;
float  b1pt_jecup_PileUpPtEC2,   b1e_jecup_PileUpPtEC2,   b2pt_jecup_PileUpPtEC2,   b2e_jecup_PileUpPtEC2;
float  b1pt_jecdn_PileUpPtEC2,   b1e_jecdn_PileUpPtEC2,   b2pt_jecdn_PileUpPtEC2,   b2e_jecdn_PileUpPtEC2;
float  b1pt_jecup_RelativeBal,   b1e_jecup_RelativeBal,   b2pt_jecup_RelativeBal,   b2e_jecup_RelativeBal;
float  b1pt_jecdn_RelativeBal,   b1e_jecdn_RelativeBal,   b2pt_jecdn_RelativeBal,   b2e_jecdn_RelativeBal;
float  b1pt_jecup_RelativeSample,b1e_jecup_RelativeSample,b2pt_jecup_RelativeSample,b2e_jecup_RelativeSample;
float  b1pt_jecdn_RelativeSample,b1e_jecdn_RelativeSample,b2pt_jecdn_RelativeSample,b2e_jecdn_RelativeSample;

float  sigmabb1, sigmah1;

// Event weights
double weight_nom;
double weight_PU_up,        weight_PU_dn;
double weight_leptonsf_up,  weight_leptonsf_dn;
double weight_photonsf_up,  weight_photonsf_dn;
double weight_prefiring_up, weight_prefiring_dn;
double weight_bSF_up,       weight_bSF_dn;
double weight_trig_up,      weight_trig_dn;
double weight_haspixsf_up,  weight_haspixsf_dn;
double weight_pujid_up,     weight_pujid_dn;
float  weight_PDF_up,       weight_PDF_dn;
float  weight_QCD_up,       weight_QCD_dn;
float  weight_ISR_up,       weight_ISR_dn;
float  weight_FSR_up,       weight_FSR_dn;

// B-tag SF weights
float  btag_sf;
float  btag_jes_up,       btag_jes_dn;
float  btag_lf_up,        btag_lf_dn;
float  btag_hf_up,        btag_hf_dn;
float  btag_hfstats1_up,  btag_hfstats1_dn;
float  btag_hfstats2_up,  btag_hfstats2_dn;
float  btag_lfstats1_up,  btag_lfstats1_dn;
float  btag_lfstats2_up,  btag_lfstats2_dn;
float  btag_cferr1_up,    btag_cferr1_dn;
float  btag_cferr2_up,    btag_cferr2_dn;


// ═══════════════════════════════════════════════════════════════════════════
// Main analysis function
// ═══════════════════════════════════════════════════════════════════════════
void H2AA2bbgg_analysis_treeread_18_for_pas_t2(int mass)
{
    TString mStr = TString::Format("%d", mass);

    // ┌─────────────────────────────────────────────────────────────────────┐
    // │                    A N A L Y S I S   C U T S                       │
    // └─────────────────────────────────────────────────────────────────────┘
    const float CHI2_CUT      = 1.0f;   // χ² selection threshold
    const float BDT_NOM       = 0.90f;   // BDT cut: nominal + MC-scale trees
    const float BDT_SYST      = 0.90f;   // BDT cut: all other systematic trees
    const float GG_MIN        = 10.0f;   // diphoton mass lower bound [GeV]
    const float GG_MAX        = 70.0f;   // diphoton mass upper bound [GeV]
    const float HIGGS_MASS    = 125.0f;  // H mass hypothesis for χ² [GeV]
    const float XSEC          = 1470.0f; // WH production cross section [fb]
    const float BR            = 0.22f;   // branching ratio H→AA→bbγγ
    const float MC_SCALE_UP   = 1.02f;   // photon energy scale-up shift
    const float MC_SCALE_DN   = 0.98f;   // photon energy scale-down shift
    const float NL_SCALE_UP   = 1.0025f; // non-linearity scale-up shift
    const float NL_SCALE_DN   = 0.9975f; // non-linearity scale-down shift
    const float MET_SCALE_UP  = 1.02f;   // MET uncertainty scale-up shift
    const float MET_SCALE_DN  = 0.98f;   // MET uncertainty scale-down shift
    const float THEORY_BR_FAC = 0.125f;  // theory syst. BR factor (PDF/QCD)
    // ──────────────────────────────────────────────────────────────────────

    // ── Generated event counts per mass point ────────────────────────────
    const int    massPoints[] = { 15,      20,      25,      30,      35,
                                  40,      45,      50,      55,      60,     62 };
    const double numEvents[]  = { 978847., 879852., 990842., 855855., 807846.,
                                  993824., 960823., 963832., 978832., 897857., 972851. };
    const int nMass = sizeof(massPoints) / sizeof(massPoints[0]);

    double num = -1.0;
    for (int k = 0; k < nMass; ++k)
        if (massPoints[k] == mass) { num = numEvents[k]; break; }
    if (num < 0) { cerr << "Error: Invalid mass point: " << mass << endl; return; }

    // ── Mass-dependent resolution parameters ─────────────────────────────
    // sigmabb1: bb mass resolution [GeV],  sigmah1: Hγγ mass resolution [GeV]
    struct SigmaEntry { int m; double sbb, sh; };
    const SigmaEntry sigmaTable[] = {
        { 15,  3.580, 13.430 }, { 20,  4.533, 13.649 }, { 25,  5.548, 13.284 },
        { 30,  6.397, 12.859 }, { 35,  7.429, 13.023 }, { 40,  8.358, 13.031 },
        { 45,  9.371, 12.985 }, { 50, 10.330, 12.989 }, { 55, 11.193, 12.982 },
        { 60, 12.182, 12.508 }, { 62, 12.510, 12.490 }
    };
    const int nSigma = sizeof(sigmaTable) / sizeof(sigmaTable[0]);

    sigmabb1 = 6.0f; sigmah1 = 1.0f;  // defaults (overwritten below)
    bool sigmaFound = false;
    for (int k = 0; k < nSigma; ++k) {
        if (sigmaTable[k].m == mass) {
            sigmabb1 = sigmaTable[k].sbb;
            sigmah1  = sigmaTable[k].sh;
            sigmaFound = true; break;
        }
    }
    if (!sigmaFound)
        cerr << "Warning: unknown mass point, using default sigma values." << endl;

    // ── Input / output files ─────────────────────────────────────────────
    TFile *f1    = new TFile("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2018/WH" + mStr + "/WH_mA" + mStr + "_vf_weight_nom_2.root");
    TFile *f2    = new TFile("/home/anirban/data/hh/workspace/analysis/BDT/param_bdt/param_training/2018/tree/param_bdt/v2/BDT_UL18_MC_trees_sig_" + mStr + ".root");
    TTree *Tout1 = (TTree*)f1->Get("tree");
    TTree *Tout2 = (TTree*)f2->Get("Tree");
    TFile *fout  = new TFile("test_1/output_H2AA2bbgg_M" + mStr + "_pythia8_wh.root", "RECREATE");

    // ═══════════════════════════════════════════════════════════════════════
    // Output TTrees  (all share the same three branches)
    // ═══════════════════════════════════════════════════════════════════════
    auto makeTree = [&](const TString& name) -> TTree* {
        TTree* t = new TTree(name, name);
        t->Branch("CMS_hgg_mass", &CMS_hgg_mass, "CMS_hgg_mass/F");
        t->Branch("weight",       &weight,        "weight/F");
        t->Branch("dZ",           &dZ,            "dZ/F");
        return t;
    };

    // Nominal
    TTree *T0  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0");

    // MC photon energy scale systematics
    TTree *T1  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleGain1EBUp01sigma");
    TTree *T2  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleGain1EBDown01sigma");
    TTree *T3  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleGain6EBUp01sigma");
    TTree *T4  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleGain6EBDown01sigma");
    TTree *T5  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleHighR9EBUp01sigma");
    TTree *T6  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleHighR9EBDown01sigma");
    TTree *T7  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleHighR9EEUp01sigma");
    TTree *T8  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleHighR9EEDown01sigma");
    TTree *T9  = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleLowR9EBUp01sigma");
    TTree *T10 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleLowR9EBDown01sigma");
    TTree *T11 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleLowR9EEUp01sigma");
    TTree *T12 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MCScaleLowR9EEDown01sigma");

    // JEC systematics
    TTree *T13 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECAbsolutecombinedUp01sigma");
    TTree *T14 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECAbsolutecombinedDown01sigma");
    TTree *T15 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECFlavorQCDUp01sigma");
    TTree *T16 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECFlavorQCDDown01sigma");
    TTree *T17 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECBBEC1combinedUp01sigma");
    TTree *T18 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECBBEC1combinedDown01sigma");

    // Non-linearity nuisance
    TTree *T19 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_NonLinearityUp01sigma");
    TTree *T20 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_NonLinearityDown01sigma");

    // JEC EC2 + RelativeBal + RelativeSample
    TTree *T21 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECEC2combinedUp01sigma");
    TTree *T22 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECEC2combinedDown01sigma");
    TTree *T23 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECRelativeBalUp01sigma");
    TTree *T24 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECRelativeBalDown01sigma");
    TTree *T25 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECRelativeSamplecombinedUp01sigma");
    TTree *T26 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JECRelativeSamplecombinedDown01sigma");

    // JER
    TTree *T27 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JERcombinedUp01sigma");
    TTree *T28 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JERcombinedDown01sigma");

    // Trigger weight
    TTree *T29 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_TriggerWeightcombinedUp01sigma");
    TTree *T30 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_TriggerWeightcombinedDown01sigma");

    // B-tag reshape (JES)
    TTree *T31 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JetBTagReshapeWeightUp01sigma");
    TTree *T32 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_JetBTagReshapeWeightDown01sigma");

    // Photon MVA shift
    TTree *T33 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MvaShiftcombinedUp01sigma");
    TTree *T34 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_MvaShiftcombinedDown01sigma");

    // PUJID shift
    TTree *T35 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_PUJIDShiftcombinedUp01sigma");
    TTree *T36 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_PUJIDShiftcombinedDown01sigma");

    // MET uncertainty
    TTree *T37 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_metUncUncertaintycombinedUp01sigma");
    TTree *T38 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_metUncUncertaintycombinedDown01sigma");

    // Theory: PDF / QCD scale / alphaS
    TTree *T39 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_pdf_Higgs_VHUp01sigma");
    TTree *T40 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_pdf_Higgs_VHDown01sigma");
    TTree *T41 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_alphaS_VHUp01sigma");
    TTree *T42 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_alphaS_VHDown01sigma");

    // PU
    TTree *T43 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_PUcombinedUp01sigma");
    TTree *T44 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_PUcombinedDown01sigma");

    // B-tag nuisances
    TTree *T45 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_JESUp01sigma");
    TTree *T46 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_JESDown01sigma");
    TTree *T47 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_HFUp01sigma");
    TTree *T48 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_HFDown01sigma");
    TTree *T49 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_LFUp01sigma");
    TTree *T50 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_LFDown01sigma");
    TTree *T51 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_CFERR1Up01sigma");
    TTree *T52 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_CFERR1Down01sigma");
    TTree *T53 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_CFERR2Up01sigma");
    TTree *T54 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_CFERR2Down01sigma");
    TTree *T55 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_HFSTAT1combinedUp01sigma");
    TTree *T56 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_HFSTAT1combinedDown01sigma");
    TTree *T57 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_HFSTAT2combinedUp01sigma");
    TTree *T58 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_HFSTAT2combinedDown01sigma");
    TTree *T59 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_LFSTAT1combinedUp01sigma");
    TTree *T60 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_LFSTAT1combinedDown01sigma");
    TTree *T61 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_LFSTAT2combinedUp01sigma");
    TTree *T62 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_btag_LFSTAT2combinedDown01sigma");

    // Lepton ID / prefiring / ISR / FSR
    TTree *T63 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_LeptonIDWeightcombinedUp01sigma");
    TTree *T64 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_LeptonIDWeightcombinedDown01sigma");
    TTree *T65 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_prefireWeightcombinedUp01sigma");
    TTree *T66 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_prefireWeightcombinedDown01sigma");
    TTree *T67 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_ISRUp01sigma");
    TTree *T68 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_ISRDown01sigma");
    TTree *T69 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_FSRUp01sigma");
    TTree *T70 = makeTree("wh_" + mStr + "_13TeV_H2AA2bbgg_Tag0_FSRDown01sigma");

    // ═══════════════════════════════════════════════════════════════════════
    // Connect input branches
    // ═══════════════════════════════════════════════════════════════════════

    // Jet 4-vectors (nominal)
    Tout1->SetBranchAddress("b1pt",  &b1pt);   Tout1->SetBranchAddress("b1eta", &b1eta);
    Tout1->SetBranchAddress("b1phi", &b1phi);  Tout1->SetBranchAddress("b1e",   &b1e);
    Tout1->SetBranchAddress("b2pt",  &b2pt);   Tout1->SetBranchAddress("b2eta", &b2eta);
    Tout1->SetBranchAddress("b2phi", &b2phi);  Tout1->SetBranchAddress("b2e",   &b2e);

    // Photon 4-vectors
    Tout1->SetBranchAddress("pho1pt",  &pho1pt);  Tout1->SetBranchAddress("pho1eta", &pho1eta);
    Tout1->SetBranchAddress("pho1phi", &pho1phi); Tout1->SetBranchAddress("pho1e",   &pho1e);
    Tout1->SetBranchAddress("pho2pt",  &pho2pt);  Tout1->SetBranchAddress("pho2eta", &pho2eta);
    Tout1->SetBranchAddress("pho2phi", &pho2phi); Tout1->SetBranchAddress("pho2e",   &pho2e);

    // JEC / JER variations
    Tout1->SetBranchAddress("b1pt_jecup", &b1pt_jecup); Tout1->SetBranchAddress("b1e_jecup", &b1e_jecup);
    Tout1->SetBranchAddress("b2pt_jecup", &b2pt_jecup); Tout1->SetBranchAddress("b2e_jecup", &b2e_jecup);
    Tout1->SetBranchAddress("b1pt_jecdn", &b1pt_jecdn); Tout1->SetBranchAddress("b1e_jecdn", &b1e_jecdn);
    Tout1->SetBranchAddress("b2pt_jecdn", &b2pt_jecdn); Tout1->SetBranchAddress("b2e_jecdn", &b2e_jecdn);
    Tout1->SetBranchAddress("b1pt_resup", &b1pt_resup); Tout1->SetBranchAddress("b1e_resup", &b1e_resup);
    Tout1->SetBranchAddress("b2pt_resup", &b2pt_resup); Tout1->SetBranchAddress("b2e_resup", &b2e_resup);
    Tout1->SetBranchAddress("b1pt_resdn", &b1pt_resdn); Tout1->SetBranchAddress("b1e_resdn", &b1e_resdn);
    Tout1->SetBranchAddress("b2pt_resdn", &b2pt_resdn); Tout1->SetBranchAddress("b2e_resdn", &b2e_resdn);

    // Source-specific JEC variations
    Tout1->SetBranchAddress("b1pt_jecup_AbsoluteStat",   &b1pt_jecup_AbsoluteStat);  Tout1->SetBranchAddress("b1e_jecup_AbsoluteStat",  &b1e_jecup_AbsoluteStat);
    Tout1->SetBranchAddress("b2pt_jecup_AbsoluteStat",   &b2pt_jecup_AbsoluteStat);  Tout1->SetBranchAddress("b2e_jecup_AbsoluteStat",  &b2e_jecup_AbsoluteStat);
    Tout1->SetBranchAddress("b1pt_jecdn_AbsoluteStat",   &b1pt_jecdn_AbsoluteStat);  Tout1->SetBranchAddress("b1e_jecdn_AbsoluteStat",  &b1e_jecdn_AbsoluteStat);
    Tout1->SetBranchAddress("b2pt_jecdn_AbsoluteStat",   &b2pt_jecdn_AbsoluteStat);  Tout1->SetBranchAddress("b2e_jecdn_AbsoluteStat",  &b2e_jecdn_AbsoluteStat);
    Tout1->SetBranchAddress("b1pt_jecup_FlavorQCD",      &b1pt_jecup_FlavorQCD);     Tout1->SetBranchAddress("b1e_jecup_FlavorQCD",     &b1e_jecup_FlavorQCD);
    Tout1->SetBranchAddress("b2pt_jecup_FlavorQCD",      &b2pt_jecup_FlavorQCD);     Tout1->SetBranchAddress("b2e_jecup_FlavorQCD",     &b2e_jecup_FlavorQCD);
    Tout1->SetBranchAddress("b1pt_jecdn_FlavorQCD",      &b1pt_jecdn_FlavorQCD);     Tout1->SetBranchAddress("b1e_jecdn_FlavorQCD",     &b1e_jecdn_FlavorQCD);
    Tout1->SetBranchAddress("b2pt_jecdn_FlavorQCD",      &b2pt_jecdn_FlavorQCD);     Tout1->SetBranchAddress("b2e_jecdn_FlavorQCD",     &b2e_jecdn_FlavorQCD);
    Tout1->SetBranchAddress("b1pt_jecup_PileUpPtBB",     &b1pt_jecup_PileUpPtBB);    Tout1->SetBranchAddress("b1e_jecup_PileUpPtBB",    &b1e_jecup_PileUpPtBB);
    Tout1->SetBranchAddress("b2pt_jecup_PileUpPtBB",     &b2pt_jecup_PileUpPtBB);    Tout1->SetBranchAddress("b2e_jecup_PileUpPtBB",    &b2e_jecup_PileUpPtBB);
    Tout1->SetBranchAddress("b1pt_jecdn_PileUpPtBB",     &b1pt_jecdn_PileUpPtBB);    Tout1->SetBranchAddress("b1e_jecdn_PileUpPtBB",    &b1e_jecdn_PileUpPtBB);
    Tout1->SetBranchAddress("b2pt_jecdn_PileUpPtBB",     &b2pt_jecdn_PileUpPtBB);    Tout1->SetBranchAddress("b2e_jecdn_PileUpPtBB",    &b2e_jecdn_PileUpPtBB);
    Tout1->SetBranchAddress("b1pt_jecup_PileUpPtEC2",    &b1pt_jecup_PileUpPtEC2);   Tout1->SetBranchAddress("b1e_jecup_PileUpPtEC2",   &b1e_jecup_PileUpPtEC2);
    Tout1->SetBranchAddress("b2pt_jecup_PileUpPtEC2",    &b2pt_jecup_PileUpPtEC2);   Tout1->SetBranchAddress("b2e_jecup_PileUpPtEC2",   &b2e_jecup_PileUpPtEC2);
    Tout1->SetBranchAddress("b1pt_jecdn_PileUpPtEC2",    &b1pt_jecdn_PileUpPtEC2);   Tout1->SetBranchAddress("b1e_jecdn_PileUpPtEC2",   &b1e_jecdn_PileUpPtEC2);
    Tout1->SetBranchAddress("b2pt_jecdn_PileUpPtEC2",    &b2pt_jecdn_PileUpPtEC2);   Tout1->SetBranchAddress("b2e_jecdn_PileUpPtEC2",   &b2e_jecdn_PileUpPtEC2);
    Tout1->SetBranchAddress("b1pt_jecup_RelativeBal",    &b1pt_jecup_RelativeBal);   Tout1->SetBranchAddress("b1e_jecup_RelativeBal",   &b1e_jecup_RelativeBal);
    Tout1->SetBranchAddress("b2pt_jecup_RelativeBal",    &b2pt_jecup_RelativeBal);   Tout1->SetBranchAddress("b2e_jecup_RelativeBal",   &b2e_jecup_RelativeBal);
    Tout1->SetBranchAddress("b1pt_jecdn_RelativeBal",    &b1pt_jecdn_RelativeBal);   Tout1->SetBranchAddress("b1e_jecdn_RelativeBal",   &b1e_jecdn_RelativeBal);
    Tout1->SetBranchAddress("b2pt_jecdn_RelativeBal",    &b2pt_jecdn_RelativeBal);   Tout1->SetBranchAddress("b2e_jecdn_RelativeBal",   &b2e_jecdn_RelativeBal);
    Tout1->SetBranchAddress("b1pt_jecup_RelativeSample", &b1pt_jecup_RelativeSample);Tout1->SetBranchAddress("b1e_jecup_RelativeSample",&b1e_jecup_RelativeSample);
    Tout1->SetBranchAddress("b2pt_jecup_RelativeSample", &b2pt_jecup_RelativeSample);Tout1->SetBranchAddress("b2e_jecup_RelativeSample",&b2e_jecup_RelativeSample);
    Tout1->SetBranchAddress("b1pt_jecdn_RelativeSample", &b1pt_jecdn_RelativeSample);Tout1->SetBranchAddress("b1e_jecdn_RelativeSample",&b1e_jecdn_RelativeSample);
    Tout1->SetBranchAddress("b2pt_jecdn_RelativeSample", &b2pt_jecdn_RelativeSample);Tout1->SetBranchAddress("b2e_jecdn_RelativeSample",&b2e_jecdn_RelativeSample);

    // Event weights
    Tout1->SetBranchAddress("weight_nom",         &weight_nom);
    Tout1->SetBranchAddress("weight_PU_up",       &weight_PU_up);       Tout1->SetBranchAddress("weight_PU_dn",       &weight_PU_dn);
    Tout1->SetBranchAddress("weight_leptonsf_up", &weight_leptonsf_up); Tout1->SetBranchAddress("weight_leptonsf_dn", &weight_leptonsf_dn);
    Tout1->SetBranchAddress("weight_photonsf_up", &weight_photonsf_up); Tout1->SetBranchAddress("weight_photonsf_dn", &weight_photonsf_dn);
    Tout1->SetBranchAddress("weight_prefiring_up",&weight_prefiring_up);Tout1->SetBranchAddress("weight_prefiring_dn",&weight_prefiring_dn);
    Tout1->SetBranchAddress("weight_bSF_up",      &weight_bSF_up);      Tout1->SetBranchAddress("weight_bSF_dn",      &weight_bSF_dn);
    Tout1->SetBranchAddress("weight_trig_up",     &weight_trig_up);     Tout1->SetBranchAddress("weight_trig_dn",     &weight_trig_dn);
    Tout1->SetBranchAddress("weight_haspixsf_up", &weight_haspixsf_up); Tout1->SetBranchAddress("weight_haspixsf_dn", &weight_haspixsf_dn);
    Tout1->SetBranchAddress("weight_pujid_up",    &weight_pujid_up);    Tout1->SetBranchAddress("weight_pujid_dn",    &weight_pujid_dn);
    Tout1->SetBranchAddress("weight_PDF_up",      &weight_PDF_up);      Tout1->SetBranchAddress("weight_PDF_dn",      &weight_PDF_dn);
    Tout1->SetBranchAddress("weight_QCD_up",      &weight_QCD_up);      Tout1->SetBranchAddress("weight_QCD_dn",      &weight_QCD_dn);
    Tout1->SetBranchAddress("weight_ISR_up",      &weight_ISR_up);      Tout1->SetBranchAddress("weight_ISR_dn",      &weight_ISR_dn);
    Tout1->SetBranchAddress("weight_FSR_up",      &weight_FSR_up);      Tout1->SetBranchAddress("weight_FSR_dn",      &weight_FSR_dn);

    // B-tag SF weights
    Tout1->SetBranchAddress("btag_sf",           &btag_sf);
    Tout1->SetBranchAddress("btag_jes_up",       &btag_jes_up);       Tout1->SetBranchAddress("btag_jes_dn",       &btag_jes_dn);
    Tout1->SetBranchAddress("btag_lf_up",        &btag_lf_up);        Tout1->SetBranchAddress("btag_lf_dn",        &btag_lf_dn);
    Tout1->SetBranchAddress("btag_hf_up",        &btag_hf_up);        Tout1->SetBranchAddress("btag_hf_dn",        &btag_hf_dn);
    Tout1->SetBranchAddress("btag_hfstats1_up",  &btag_hfstats1_up);  Tout1->SetBranchAddress("btag_hfstats1_dn",  &btag_hfstats1_dn);
    Tout1->SetBranchAddress("btag_hfstats2_up",  &btag_hfstats2_up);  Tout1->SetBranchAddress("btag_hfstats2_dn",  &btag_hfstats2_dn);
    Tout1->SetBranchAddress("btag_lfstats1_up",  &btag_lfstats1_up);  Tout1->SetBranchAddress("btag_lfstats1_dn",  &btag_lfstats1_dn);
    Tout1->SetBranchAddress("btag_lfstats2_up",  &btag_lfstats2_up);  Tout1->SetBranchAddress("btag_lfstats2_dn",  &btag_lfstats2_dn);
    Tout1->SetBranchAddress("btag_cferr1_up",    &btag_cferr1_up);    Tout1->SetBranchAddress("btag_cferr1_dn",    &btag_cferr1_dn);
    Tout1->SetBranchAddress("btag_cferr2_up",    &btag_cferr2_up);    Tout1->SetBranchAddress("btag_cferr2_dn",    &btag_cferr2_dn);

    // BDT score
    Tout2->SetBranchAddress("BDT", &BDT);

    int nevt = Tout1->GetEntries();

    // ═══════════════════════════════════════════════════════════════════════
    // Yield accumulators for systematic ratio printout
    // ═══════════════════════════════════════════════════════════════════════
    double sum_nom         = 0.0;
    double sum_PU_up       = 0.0, sum_PU_dn       = 0.0;
    double sum_leptonsf_up = 0.0, sum_leptonsf_dn = 0.0;
    double sum_photonsf_up = 0.0, sum_photonsf_dn = 0.0;
    double sum_prefiring_up= 0.0, sum_prefiring_dn= 0.0;
    double sum_bSF_up      = 0.0, sum_bSF_dn      = 0.0;
    double sum_trig_up     = 0.0, sum_trig_dn     = 0.0;
    double sum_pujid_up    = 0.0, sum_pujid_dn    = 0.0;
    double sum_PDF_up      = 0.0, sum_PDF_dn      = 0.0;
    double sum_QCD_up      = 0.0, sum_QCD_dn      = 0.0;
    double sum_ISR_up      = 0.0, sum_ISR_dn      = 0.0;
    double sum_FSR_up      = 0.0, sum_FSR_dn      = 0.0;
    double sum_jecup       = 0.0, sum_jecdn       = 0.0;
    double sum_resup       = 0.0, sum_resdn       = 0.0;
    double sum_jes_up      = 0.0, sum_jes_dn      = 0.0;
    double sum_lf_up       = 0.0, sum_lf_dn       = 0.0;
    double sum_hf_up       = 0.0, sum_hf_dn       = 0.0;
    double sum_hfstats1_up = 0.0, sum_hfstats1_dn = 0.0;
    double sum_hfstats2_up = 0.0, sum_hfstats2_dn = 0.0;
    double sum_lfstats1_up = 0.0, sum_lfstats1_dn = 0.0;
    double sum_lfstats2_up = 0.0, sum_lfstats2_dn = 0.0;
    double sum_cferr1_up   = 0.0, sum_cferr1_dn   = 0.0;
    double sum_cferr2_up   = 0.0, sum_cferr2_dn   = 0.0;

    double c_0 = 0.0, c_1 = 0.0, c_2 = 0.0;
    double c_int_0 = 0.0, c_int_1 = 0.0, c_int_2 = 0.0;

    // ═══════════════════════════════════════════════════════════════════════
    // PASS 1: Compute b-tag SF normalization ratios
    //
    // Prescription:  r_X = sum_before / sum_after_X
    //   sum_before  = Σ (weight_nom / btag_sf)          over all events
    //   sum_after_X = Σ (weight_nom / btag_sf * btag_X) over all events
    // Multiplying by r_X in Pass 2 preserves the total event yield.
    // ═══════════════════════════════════════════════════════════════════════
    double sum_before           = 0.0;
    double sum_after_jes_up     = 0.0, sum_after_jes_dn     = 0.0;
    double sum_after_hf_up      = 0.0, sum_after_hf_dn      = 0.0;
    double sum_after_hfstats1_up= 0.0, sum_after_hfstats1_dn= 0.0;
    double sum_after_hfstats2_up= 0.0, sum_after_hfstats2_dn= 0.0;
    double sum_after_lf_up      = 0.0, sum_after_lf_dn      = 0.0;
    double sum_after_lfstats1_up= 0.0, sum_after_lfstats1_dn= 0.0;
    double sum_after_lfstats2_up= 0.0, sum_after_lfstats2_dn= 0.0;
    double sum_after_cferr1_up  = 0.0, sum_after_cferr1_dn  = 0.0;
    double sum_after_cferr2_up  = 0.0, sum_after_cferr2_dn  = 0.0;

    for (int i = 0; i < nevt; ++i)
    {
        Tout1->GetEntry(i);
        if (btag_sf <= 0.0) continue;  // guard against zero / bad SF

        double w_base = weight_nom / btag_sf;
        sum_before            += w_base;
        sum_after_jes_up      += w_base * btag_jes_up;
        sum_after_jes_dn      += w_base * btag_jes_dn;
        sum_after_hf_up       += w_base * btag_hf_up;
        sum_after_hf_dn       += w_base * btag_hf_dn;
        sum_after_hfstats1_up += w_base * btag_hfstats1_up;
        sum_after_hfstats1_dn += w_base * btag_hfstats1_dn;
        sum_after_hfstats2_up += w_base * btag_hfstats2_up;
        sum_after_hfstats2_dn += w_base * btag_hfstats2_dn;
        sum_after_lf_up       += w_base * btag_lf_up;
        sum_after_lf_dn       += w_base * btag_lf_dn;
        sum_after_lfstats1_up += w_base * btag_lfstats1_up;
        sum_after_lfstats1_dn += w_base * btag_lfstats1_dn;
        sum_after_lfstats2_up += w_base * btag_lfstats2_up;
        sum_after_lfstats2_dn += w_base * btag_lfstats2_dn;
        sum_after_cferr1_up   += w_base * btag_cferr1_up;
        sum_after_cferr1_dn   += w_base * btag_cferr1_dn;
        sum_after_cferr2_up   += w_base * btag_cferr2_up;
        sum_after_cferr2_dn   += w_base * btag_cferr2_dn;
    }

    // r_X = sum_before / sum_after_X  (fall back to 1.0 if denominator is zero)
    auto safeRatio = [&](double after) -> double {
        return (after > 0.0) ? sum_before / after : 1.0;
    };
    double r_btag_jes_up      = safeRatio(sum_after_jes_up);
    double r_btag_jes_dn      = safeRatio(sum_after_jes_dn);
    double r_btag_hf_up       = safeRatio(sum_after_hf_up);
    double r_btag_hf_dn       = safeRatio(sum_after_hf_dn);
    double r_btag_hfstats1_up = safeRatio(sum_after_hfstats1_up);
    double r_btag_hfstats1_dn = safeRatio(sum_after_hfstats1_dn);
    double r_btag_hfstats2_up = safeRatio(sum_after_hfstats2_up);
    double r_btag_hfstats2_dn = safeRatio(sum_after_hfstats2_dn);
    double r_btag_lf_up       = safeRatio(sum_after_lf_up);
    double r_btag_lf_dn       = safeRatio(sum_after_lf_dn);
    double r_btag_lfstats1_up = safeRatio(sum_after_lfstats1_up);
    double r_btag_lfstats1_dn = safeRatio(sum_after_lfstats1_dn);
    double r_btag_lfstats2_up = safeRatio(sum_after_lfstats2_up);
    double r_btag_lfstats2_dn = safeRatio(sum_after_lfstats2_dn);
    double r_btag_cferr1_up   = safeRatio(sum_after_cferr1_up);
    double r_btag_cferr1_dn   = safeRatio(sum_after_cferr1_dn);
    double r_btag_cferr2_up   = safeRatio(sum_after_cferr2_up);
    double r_btag_cferr2_dn   = safeRatio(sum_after_cferr2_dn);

    // ═══════════════════════════════════════════════════════════════════════
    // PASS 2: Main event loop
    // ═══════════════════════════════════════════════════════════════════════
    for (int i = 0; i < nevt; ++i)
    {
        Tout1->GetEntry(i);
        Tout2->GetEntry(i);

        // ── Jet 4-vectors for all JEC/JER variations ──────────────────────
        TLorentzVector b1,  b2;    // [0] nominal
        TLorentzVector b3,  b4;    // [1] JEC up
        TLorentzVector b5,  b6;    // [2] JER up
        TLorentzVector b7,  b8;    // [3] JEC down
        TLorentzVector b9,  b10;   // [4] JER down
        TLorentzVector b11, b12;   // [5] AbsoluteStat up
        TLorentzVector b13, b14;   // [6] AbsoluteStat dn
        TLorentzVector b15, b16;   // [7] FlavorQCD up
        TLorentzVector b17, b18;   // [8] FlavorQCD dn
        TLorentzVector b19, b20;   // [9] PileUpPtBB up
        TLorentzVector b21, b22;   // [10] PileUpPtBB dn
        TLorentzVector b23, b24;   // [11] PileUpPtEC2 up
        TLorentzVector b25, b26;   // [12] PileUpPtEC2 dn
        TLorentzVector b27, b28;   // [13] RelativeBal up
        TLorentzVector b29, b30;   // [14] RelativeBal dn
        TLorentzVector b31, b32;   // [15] RelativeSample up
        TLorentzVector b33, b34;   // [16] RelativeSample dn

        b1.SetPtEtaPhiE(b1pt, b1eta, b1phi, b1e);
        b2.SetPtEtaPhiE(b2pt, b2eta, b2phi, b2e);
        b3.SetPtEtaPhiE(b1pt_jecup, b1eta, b1phi, b1e_jecup);
        b4.SetPtEtaPhiE(b2pt_jecup, b2eta, b2phi, b2e_jecup);
        b5.SetPtEtaPhiE(b1pt_resup, b1eta, b1phi, b1e_resup);
        b6.SetPtEtaPhiE(b2pt_resup, b2eta, b2phi, b2e_resup);
        b7.SetPtEtaPhiE(b1pt_jecdn, b1eta, b1phi, b1e_jecdn);
        b8.SetPtEtaPhiE(b2pt_jecdn, b2eta, b2phi, b2e_jecdn);
        b9.SetPtEtaPhiE(b1pt_resdn, b1eta, b1phi, b1e_resdn);
        b10.SetPtEtaPhiE(b2pt_resdn, b2eta, b2phi, b2e_resdn);
        b11.SetPtEtaPhiE(b1pt_jecup_AbsoluteStat,  b1eta, b1phi, b1e_jecup_AbsoluteStat);
        b12.SetPtEtaPhiE(b2pt_jecup_AbsoluteStat,  b2eta, b2phi, b2e_jecup_AbsoluteStat);
        b13.SetPtEtaPhiE(b1pt_jecdn_AbsoluteStat,  b1eta, b1phi, b1e_jecdn_AbsoluteStat);
        b14.SetPtEtaPhiE(b2pt_jecdn_AbsoluteStat,  b2eta, b2phi, b2e_jecdn_AbsoluteStat);
        b15.SetPtEtaPhiE(b1pt_jecup_FlavorQCD,     b1eta, b1phi, b1e_jecup_FlavorQCD);
        b16.SetPtEtaPhiE(b2pt_jecup_FlavorQCD,     b2eta, b2phi, b2e_jecup_FlavorQCD);
        b17.SetPtEtaPhiE(b1pt_jecdn_FlavorQCD,     b1eta, b1phi, b1e_jecdn_FlavorQCD);
        b18.SetPtEtaPhiE(b2pt_jecdn_FlavorQCD,     b2eta, b2phi, b2e_jecdn_FlavorQCD);
        b19.SetPtEtaPhiE(b1pt_jecup_PileUpPtBB,    b1eta, b1phi, b1e_jecup_PileUpPtBB);
        b20.SetPtEtaPhiE(b2pt_jecup_PileUpPtBB,    b2eta, b2phi, b2e_jecup_PileUpPtBB);
        b21.SetPtEtaPhiE(b1pt_jecdn_PileUpPtBB,    b1eta, b1phi, b1e_jecdn_PileUpPtBB);
        b22.SetPtEtaPhiE(b2pt_jecdn_PileUpPtBB,    b2eta, b2phi, b2e_jecdn_PileUpPtBB);
        b23.SetPtEtaPhiE(b1pt_jecup_PileUpPtEC2,   b1eta, b1phi, b1e_jecup_PileUpPtEC2);
        b24.SetPtEtaPhiE(b2pt_jecup_PileUpPtEC2,   b2eta, b2phi, b2e_jecup_PileUpPtEC2);
        b25.SetPtEtaPhiE(b1pt_jecdn_PileUpPtEC2,   b1eta, b1phi, b1e_jecdn_PileUpPtEC2);
        b26.SetPtEtaPhiE(b2pt_jecdn_PileUpPtEC2,   b2eta, b2phi, b2e_jecdn_PileUpPtEC2);
        b27.SetPtEtaPhiE(b1pt_jecup_RelativeBal,   b1eta, b1phi, b1e_jecup_RelativeBal);
        b28.SetPtEtaPhiE(b2pt_jecup_RelativeBal,   b2eta, b2phi, b2e_jecup_RelativeBal);
        b29.SetPtEtaPhiE(b1pt_jecdn_RelativeBal,   b1eta, b1phi, b1e_jecdn_RelativeBal);
        b30.SetPtEtaPhiE(b2pt_jecdn_RelativeBal,   b2eta, b2phi, b2e_jecdn_RelativeBal);
        b31.SetPtEtaPhiE(b1pt_jecup_RelativeSample,b1eta, b1phi, b1e_jecup_RelativeSample);
        b32.SetPtEtaPhiE(b2pt_jecup_RelativeSample,b2eta, b2phi, b2e_jecup_RelativeSample);
        b33.SetPtEtaPhiE(b1pt_jecdn_RelativeSample,b1eta, b1phi, b1e_jecdn_RelativeSample);
        b34.SetPtEtaPhiE(b2pt_jecdn_RelativeSample,b2eta, b2phi, b2e_jecdn_RelativeSample);

        // ── Photon 4-vectors & diphoton mass ──────────────────────────────
        TLorentzVector g1, g2;
        g1.SetPtEtaPhiE(pho1pt, pho1eta, pho1phi, pho1e);
        g2.SetPtEtaPhiE(pho2pt, pho2eta, pho2phi, pho2e);
        float gg = (g1 + g2).M();

        // ── bb invariant masses for each jet variation ─────────────────────
        float bb1  = (b1+b2).M();   float bb2  = (b3+b4).M();   float bb3  = (b5+b6).M();
        float bb4  = (b7+b8).M();   float bb5  = (b9+b10).M();  float bb6  = (b11+b12).M();
        float bb7  = (b13+b14).M(); float bb8  = (b15+b16).M(); float bb9  = (b17+b18).M();
        float bb10 = (b19+b20).M(); float bb11 = (b21+b22).M(); float bb12 = (b23+b24).M();
        float bb13 = (b25+b26).M(); float bb14 = (b27+b28).M(); float bb15 = (b29+b30).M();
        float bb16 = (b31+b32).M(); float bb17 = (b33+b34).M();

        // ── bbγγ invariant masses ─────────────────────────────────────────
        float bbgg1  = (b1+b2+g1+g2).M();   float bbgg2  = (b3+b4+g1+g2).M();
        float bbgg3  = (b5+b6+g1+g2).M();   float bbgg4  = (b7+b8+g1+g2).M();
        float bbgg5  = (b9+b10+g1+g2).M();  float bbgg6  = (b11+b12+g1+g2).M();
        float bbgg7  = (b13+b14+g1+g2).M(); float bbgg8  = (b15+b16+g1+g2).M();
        float bbgg9  = (b17+b18+g1+g2).M(); float bbgg10 = (b19+b20+g1+g2).M();
        float bbgg11 = (b21+b22+g1+g2).M(); float bbgg12 = (b23+b24+g1+g2).M();
        float bbgg13 = (b25+b26+g1+g2).M(); float bbgg14 = (b27+b28+g1+g2).M();
        float bbgg15 = (b29+b30+g1+g2).M(); float bbgg16 = (b31+b32+g1+g2).M();
        float bbgg17 = (b33+b34+g1+g2).M();

        // ── χ² = ((M_bbgg − M_H)/σ_H)² + ((M_bb − M_gg)/σ_bb)² ─────────
        float chi2[17];
        chi2[0]  = pow((bbgg1  - HIGGS_MASS)/sigmah1, 2) + pow((bb1  - gg)/sigmabb1, 2);
        chi2[1]  = pow((bbgg2  - HIGGS_MASS)/sigmah1, 2) + pow((bb2  - gg)/sigmabb1, 2);
        chi2[2]  = pow((bbgg3  - HIGGS_MASS)/sigmah1, 2) + pow((bb3  - gg)/sigmabb1, 2);
        chi2[3]  = pow((bbgg4  - HIGGS_MASS)/sigmah1, 2) + pow((bb4  - gg)/sigmabb1, 2);
        chi2[4]  = pow((bbgg5  - HIGGS_MASS)/sigmah1, 2) + pow((bb5  - gg)/sigmabb1, 2);
        chi2[5]  = pow((bbgg6  - HIGGS_MASS)/sigmah1, 2) + pow((bb6  - gg)/sigmabb1, 2);
        chi2[6]  = pow((bbgg7  - HIGGS_MASS)/sigmah1, 2) + pow((bb7  - gg)/sigmabb1, 2);
        chi2[7]  = pow((bbgg8  - HIGGS_MASS)/sigmah1, 2) + pow((bb8  - gg)/sigmabb1, 2);
        chi2[8]  = pow((bbgg9  - HIGGS_MASS)/sigmah1, 2) + pow((bb9  - gg)/sigmabb1, 2);
        chi2[9]  = pow((bbgg10 - HIGGS_MASS)/sigmah1, 2) + pow((bb10 - gg)/sigmabb1, 2);
        chi2[10] = pow((bbgg11 - HIGGS_MASS)/sigmah1, 2) + pow((bb11 - gg)/sigmabb1, 2);
        chi2[11] = pow((bbgg12 - HIGGS_MASS)/sigmah1, 2) + pow((bb12 - gg)/sigmabb1, 2);
        chi2[12] = pow((bbgg13 - HIGGS_MASS)/sigmah1, 2) + pow((bb13 - gg)/sigmabb1, 2);
        chi2[13] = pow((bbgg14 - HIGGS_MASS)/sigmah1, 2) + pow((bb14 - gg)/sigmabb1, 2);
        chi2[14] = pow((bbgg15 - HIGGS_MASS)/sigmah1, 2) + pow((bb15 - gg)/sigmabb1, 2);
        chi2[15] = pow((bbgg16 - HIGGS_MASS)/sigmah1, 2) + pow((bb16 - gg)/sigmabb1, 2);
        chi2[16] = pow((bbgg17 - HIGGS_MASS)/sigmah1, 2) + pow((bb17 - gg)/sigmabb1, 2);

        // ── Pre-compute common selection flags ────────────────────────────
        bool inGGWindow  = (gg >= GG_MIN && gg <= GG_MAX);
        bool passChi2    = (chi2[0] <= CHI2_CUT);
        bool sel_nom     = (passChi2 && BDT >= BDT_NOM  && inGGWindow);
        bool sel_syst    = (passChi2 && BDT >= BDT_SYST && inGGWindow);

        // Event counter
        if (inGGWindow) c_0++;

        // Helper: set branches and fill tree
        auto fillTree = [&](TTree* t, float massVal, double wgtFactor) {
            CMS_hgg_mass = massVal;
            weight       = static_cast<float>((wgtFactor * XSEC * BR) / num);
            dZ           = 0.0f;
            t->Fill();
        };

        // ── Nominal ───────────────────────────────────────────────────────
        if (sel_nom) { fillTree(T0, gg, weight_nom); c_1++; }

        // ── MC photon energy scale (BDT_NOM) ─────────────────────────────
        if (sel_nom) { fillTree(T1,  gg * MC_SCALE_UP, weight_nom); }
        if (sel_nom) { fillTree(T2,  gg * MC_SCALE_DN, weight_nom); }
        if (sel_nom) { fillTree(T3,  gg * MC_SCALE_UP, weight_nom); }
        if (sel_nom) { fillTree(T4,  gg * MC_SCALE_DN, weight_nom); }
        if (sel_nom) { fillTree(T5,  gg * MC_SCALE_UP, weight_nom); }
        if (sel_nom) { fillTree(T6,  gg * MC_SCALE_DN, weight_nom); }
        if (sel_nom) { fillTree(T7,  gg * MC_SCALE_UP, weight_nom); }
        if (sel_nom) { fillTree(T8,  gg * MC_SCALE_DN, weight_nom); }
        if (sel_nom) { fillTree(T9,  gg * MC_SCALE_UP, weight_nom); }
        if (sel_nom) { fillTree(T10, gg * MC_SCALE_DN, weight_nom); }

        // ── MC photon energy scale (BDT_SYST) ────────────────────────────
        if (sel_syst) { fillTree(T11, gg * MC_SCALE_UP, weight_nom); }
        if (sel_syst) { fillTree(T12, gg * MC_SCALE_DN, weight_nom); }

        // ── Non-linearity nuisance ────────────────────────────────────────
        if (sel_syst) { fillTree(T19, gg * NL_SCALE_UP, weight_nom); }
        if (sel_syst) { fillTree(T20, gg * NL_SCALE_DN, weight_nom); }

        // ── Pileup reweighting ────────────────────────────────────────────
        if (sel_syst) { fillTree(T43, gg, weight_PU_up); }
        if (sel_syst) { fillTree(T44, gg, weight_PU_dn); }

        // ── Lepton ID SF ──────────────────────────────────────────────────
        if (sel_syst) { fillTree(T63, gg, weight_leptonsf_up); }
        if (sel_syst) { fillTree(T64, gg, weight_leptonsf_dn); }

        // ── Prefiring weight ──────────────────────────────────────────────
        if (sel_syst) { fillTree(T65, gg, weight_prefiring_up); }
        if (sel_syst) { fillTree(T66, gg, weight_prefiring_dn); }

        // ── MET uncertainty ───────────────────────────────────────────────
        if (sel_syst) { fillTree(T37, gg, weight_nom * MET_SCALE_UP); }
        if (sel_syst) { fillTree(T38, gg, weight_nom * MET_SCALE_DN); }

        // ── Theory systematics (PDF / QCD / ISR / FSR) ───────────────────
        if (sel_syst) { fillTree(T39, gg, weight_PDF_up * THEORY_BR_FAC); }
        if (sel_syst) { fillTree(T40, gg, weight_PDF_dn * THEORY_BR_FAC); }
        if (sel_syst) { fillTree(T41, gg, weight_QCD_up * THEORY_BR_FAC); }
        if (sel_syst) { fillTree(T42, gg, weight_QCD_dn * THEORY_BR_FAC); }
        if (sel_syst) { fillTree(T67, gg, weight_ISR_up); }
        if (sel_syst) { fillTree(T68, gg, weight_ISR_dn); }
        if (sel_syst) { fillTree(T69, gg, weight_FSR_up); }
        if (sel_syst) { fillTree(T70, gg, weight_FSR_dn); }

        // ── Photon MVA shift ──────────────────────────────────────────────
        if (sel_syst) { fillTree(T33, gg, weight_photonsf_up); }
        if (sel_syst) { fillTree(T34, gg, weight_photonsf_dn); }

        // ── PUJID shift ───────────────────────────────────────────────────
        if (sel_syst) { fillTree(T35, gg, weight_pujid_up); }
        if (sel_syst) { fillTree(T36, gg, weight_pujid_dn); }

        // ── Trigger weight ────────────────────────────────────────────────
        if (sel_syst) { fillTree(T29, gg, weight_trig_up); }
        if (sel_syst) { fillTree(T30, gg, weight_trig_dn); }

        // ── B-tag reshape (JES) ───────────────────────────────────────────
        if (sel_syst) { fillTree(T31, gg, btag_jes_up); }
        if (sel_syst) { fillTree(T32, gg, btag_jes_dn); }

        // ── B-tag nuisances (normalized by r_X) ──────────────────────────
        if (sel_syst) { fillTree(T45, gg, weight_nom * (btag_jes_up     /btag_sf) * r_btag_jes_up);      }
        if (sel_syst) { fillTree(T46, gg, weight_nom * (btag_jes_dn     /btag_sf) * r_btag_jes_dn);      }
        if (sel_syst) { fillTree(T47, gg, weight_nom * (btag_hf_up      /btag_sf) * r_btag_hf_up);       }
        if (sel_syst) { fillTree(T48, gg, weight_nom * (btag_hf_dn      /btag_sf) * r_btag_hf_dn);       }
        if (sel_syst) { fillTree(T49, gg, weight_nom * (btag_lf_up      /btag_sf) * r_btag_lf_up);       }
        if (sel_syst) { fillTree(T50, gg, weight_nom * (btag_lf_dn      /btag_sf) * r_btag_lf_dn);       }
        if (sel_syst) { fillTree(T51, gg, weight_nom * (btag_cferr1_up  /btag_sf) * r_btag_cferr1_up);   }
        if (sel_syst) { fillTree(T52, gg, weight_nom * (btag_cferr1_dn  /btag_sf) * r_btag_cferr1_dn);   }
        if (sel_syst) { fillTree(T53, gg, weight_nom * (btag_cferr2_up  /btag_sf) * r_btag_cferr2_up);   }
        if (sel_syst) { fillTree(T54, gg, weight_nom * (btag_cferr2_dn  /btag_sf) * r_btag_cferr2_dn);   }
        if (sel_syst) { fillTree(T55, gg, weight_nom * (btag_hfstats1_up/btag_sf) * r_btag_hfstats1_up); }
        if (sel_syst) { fillTree(T56, gg, weight_nom * (btag_hfstats1_dn/btag_sf) * r_btag_hfstats1_dn); }
        if (sel_syst) { fillTree(T57, gg, weight_nom * (btag_hfstats2_up/btag_sf) * r_btag_hfstats2_up); }
        if (sel_syst) { fillTree(T58, gg, weight_nom * (btag_hfstats2_dn/btag_sf) * r_btag_hfstats2_dn); }
        if (sel_syst) { fillTree(T59, gg, weight_nom * (btag_lfstats1_up/btag_sf) * r_btag_lfstats1_up); }
        if (sel_syst) { fillTree(T60, gg, weight_nom * (btag_lfstats1_dn/btag_sf) * r_btag_lfstats1_dn); }
        if (sel_syst) { fillTree(T61, gg, weight_nom * (btag_lfstats2_up/btag_sf) * r_btag_lfstats2_up); }
        if (sel_syst) { fillTree(T62, gg, weight_nom * (btag_lfstats2_dn/btag_sf) * r_btag_lfstats2_dn); }

        // ── JEC / JER varied chi2 trees ───────────────────────────────────
        if (chi2[5]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T13, gg, weight_nom); }
        if (chi2[6]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T14, gg, weight_nom); }
        if (chi2[7]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T15, gg, weight_nom); }
        if (chi2[8]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T16, gg, weight_nom); }
        if (chi2[9]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T17, gg, weight_nom); }
        if (chi2[10] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T18, gg, weight_nom); }
        if (chi2[11] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T21, gg, weight_nom); }
        if (chi2[12] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T22, gg, weight_nom); }
        if (chi2[13] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T23, gg, weight_nom); }
        if (chi2[14] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T24, gg, weight_nom); }
        if (chi2[15] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T25, gg, weight_nom); }
        if (chi2[16] <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T26, gg, weight_nom); }
        if (chi2[2]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T27, gg, weight_nom); }
        if (chi2[4]  <= CHI2_CUT && BDT >= BDT_SYST && inGGWindow) { fillTree(T28, gg, weight_nom); }

        // ── Yield sums for systematic ratio printout ──────────────────────
        const double weightFactors[41] = {
            weight_nom,       weight_PU_up,        weight_leptonsf_up,  weight_photonsf_up,  weight_prefiring_up,
            weight_bSF_up,    weight_trig_up,       weight_pujid_up,     weight_PU_dn,        weight_leptonsf_dn,
            weight_photonsf_dn,weight_prefiring_dn, weight_bSF_dn,       weight_trig_dn,      weight_pujid_dn,
            weight_PDF_up,    weight_PDF_dn,        weight_QCD_up,       weight_QCD_dn,       weight_ISR_up,
            weight_ISR_dn,    weight_FSR_up,        weight_FSR_dn,
            btag_jes_up,      btag_jes_dn,          btag_lf_up,          btag_lf_dn,
            btag_hf_up,       btag_hf_dn,           btag_hfstats1_up,    btag_hfstats1_dn,
            btag_hfstats2_up, btag_hfstats2_dn,     btag_lfstats1_up,    btag_lfstats1_dn,
            btag_lfstats2_up, btag_lfstats2_dn,     btag_cferr1_up,      btag_cferr1_dn,
            btag_cferr2_up,   btag_cferr2_dn
        };
        double* sums[41] = {
            &sum_nom,          &sum_PU_up,         &sum_leptonsf_up,  &sum_photonsf_up,  &sum_prefiring_up,
            &sum_bSF_up,       &sum_trig_up,        &sum_pujid_up,     &sum_PU_dn,        &sum_leptonsf_dn,
            &sum_photonsf_dn,  &sum_prefiring_dn,   &sum_bSF_dn,       &sum_trig_dn,      &sum_pujid_dn,
            &sum_PDF_up,       &sum_PDF_dn,          &sum_QCD_up,       &sum_QCD_dn,       &sum_ISR_up,
            &sum_ISR_dn,       &sum_FSR_up,          &sum_FSR_dn,
            &sum_jes_up,       &sum_jes_dn,          &sum_lf_up,        &sum_lf_dn,
            &sum_hf_up,        &sum_hf_dn,           &sum_hfstats1_up,  &sum_hfstats1_dn,
            &sum_hfstats2_up,  &sum_hfstats2_dn,     &sum_lfstats1_up,  &sum_lfstats1_dn,
            &sum_lfstats2_up,  &sum_lfstats2_dn,     &sum_cferr1_up,    &sum_cferr1_dn,
            &sum_cferr2_up,    &sum_cferr2_dn
        };


    } // end event loop

    // ═══════════════════════════════════════════════════════════════════════
    // Print systematic ratios  1 + (syst − nom) / nom
    // ═══════════════════════════════════════════════════════════════════════
    const char* sys_names_up[] = {
        "PU_up",       "leptonsf_up", "photonsf_up",  "prefiring_up", "bSF_up",
        "trig_up",     "pujid_up",    "PDF_up",        "QCD_up",       "ISR_up",
        "FSR_up",      "jes_up",      "lf_up",         "hf_up",        "hfstats1_up",
        "hfstats2_up", "lfstats1_up", "lfstats2_up",   "cferr1_up",    "cferr2_up"
    };
    const char* sys_names_dn[] = {
        "PU_dn",       "leptonsf_dn", "photonsf_dn",  "prefiring_dn", "bSF_dn",
        "trig_dn",     "pujid_dn",    "PDF_dn",        "QCD_dn",       "ISR_dn",
        "FSR_dn",      "jes_dn",      "lf_dn",         "hf_dn",        "hfstats1_dn",
        "hfstats2_dn", "lfstats1_dn", "lfstats2_dn",   "cferr1_dn",    "cferr2_dn"
    };
    double* sys_sums_up[] = {
        &sum_PU_up,       &sum_leptonsf_up, &sum_photonsf_up,  &sum_prefiring_up, &sum_bSF_up,
        &sum_trig_up,     &sum_pujid_up,    &sum_PDF_up,        &sum_QCD_up,       &sum_ISR_up,
        &sum_FSR_up,      &sum_jes_up,      &sum_lf_up,         &sum_hf_up,        &sum_hfstats1_up,
        &sum_hfstats2_up, &sum_lfstats1_up, &sum_lfstats2_up,   &sum_cferr1_up,    &sum_cferr2_up
    };
    double* sys_sums_dn[] = {
        &sum_PU_dn,       &sum_leptonsf_dn, &sum_photonsf_dn,  &sum_prefiring_dn, &sum_bSF_dn,
        &sum_trig_dn,     &sum_pujid_dn,    &sum_PDF_dn,        &sum_QCD_dn,       &sum_ISR_dn,
        &sum_FSR_dn,      &sum_jes_dn,      &sum_lf_dn,         &sum_hf_dn,        &sum_hfstats1_dn,
        &sum_hfstats2_dn, &sum_lfstats1_dn, &sum_lfstats2_dn,   &sum_cferr1_dn,    &sum_cferr2_dn
    };

    cout << c_0 << " " << c_1 << endl;

    // ── Cleanup & write ───────────────────────────────────────────────────
    f1->cd();  delete Tout1; delete f1;
    f2->cd();  delete Tout2; delete f2;
    fout->cd(); fout->Write(); fout->Close();
}
