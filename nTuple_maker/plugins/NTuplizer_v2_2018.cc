//Twikis used:
/*
Prefiring weights: https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe#2018_UL 
Electron MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2
Photon MVA ID: https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2 
Main EGamma: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaIDRecipesRun2#MVA_based_electron_Identificatio
JEC: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources
DeepAKX: https://twiki.cern.ch/twiki/bin/viewauth/CMS/DeepAKXTagging
Btag SF (recipe): https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagCalibration
Btag SF (2018UL): https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
Rochester correction: https://gitlab.cern.ch/akhukhun/roccor
*/

// system include files
#include <memory>
// user include files
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TAxis.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TRandom.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "GeneratorInterface/Pythia8Interface/plugins/ReweightUserHooks.h"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include <string>
#include <iostream>
#include <fstream>
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "FWCore/Utilities/interface/typelookup.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include  "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "RecoEgamma/ElectronIdentification/interface/ElectronMVAEstimatorRun2.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReportEntry.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtTrigReport.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include "FWCore/Utilities/interface/typelookup.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondFormats/BTauObjects/interface/BTagEntry.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <fastjet/GhostedAreaSpec.hh>
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"

// Rochester correction for muons //
#include "RoccoR.h"

//for storing vectors in tree//

# include <vector>


#ifdef __MAKECINT__
    
    #pragma link C++ class std::vector+;
    #pragma link C++ class std::vector<float>+;
    #pragma link C++ class std::vector<std::vector<float> >+;
    
#endif


using namespace std;
using namespace edm;
using namespace reco;  
using namespace CLHEP;
using namespace trigger;
using namespace math;
using namespace fastjet;
using namespace fastjet::contrib;

const float mu_mass = 0.105658;
const float el_mass = 0.000511;
const float pival = acos(-1.);

static const int nsrc = 24;
const char* jecsrcnames[nsrc] = {
	 "AbsoluteStat", "AbsoluteScale","AbsoluteMPFBias", 
	 "FlavorQCD", "Fragmentation", 
	 "PileUpDataMC",  "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", //"PileUpPtHF",
	 "PileUpPtRef",
	 "RelativeFSR", "RelativeJEREC1", "RelativeJEREC2", //"RelativeJERHF",
	 "RelativePtBB", "RelativePtEC1", "RelativePtEC2", //"RelativePtHF", 
	 "RelativeBal", "RelativeSample", "RelativeStatEC", "RelativeStatFSR", //"RelativeStatHF", 
	 "SinglePionECAL", "SinglePionHCAL","TimePtEta",
	 "Total"
	};
const int njecmcmx = 2*nsrc + 1 ;

struct triggervar{
  TLorentzVector  trg4v;
  bool		  both;
  bool            level1;
  bool            highl;
  int             ihlt;
  int             prescl;
  int             pdgId;
  int			  type;
};

int getbinid(double val, int nbmx, double* array) {
  if (val<array[0]) return -2;
  for (int ix=0; ix<=nbmx; ix++) {
    if (val < array[ix]) return ix-1;
  }
  return -3;
}

double theta_to_eta(double theta) { return -log(tan(theta/2.)); }

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

double diff_func(double f1, double f2){
  double ff = pow(f1-f2,2)*1./pow(f1+f2,2);
  return ff;
}


TLorentzVector productX(TLorentzVector X, TLorentzVector Y, float pro1, float pro2)
{
  float b1, b2, b3;
  float c1, c2, c3;
  
  b1 = X.Px();
  b2 = X.Py();
  b3 = X.Pz();
  
  c1 = Y.Px();
  c2 = Y.Py();
  c3 = Y.Pz();
  
  float d1, d2, e1, e2, X1, X2;
  
  X1 = pro1;
  X2 = pro2;
  
  d1 = (c2*X1 - b2*X2)*1./(b1*c2 - b2*c1);
  d2 = (c1*X1 - b1*X2)*1./(b2*c1 - b1*c2);
  e1 = (b2*c3 - b3*c2)*1./(b1*c2 - b2*c1);
  e2 = (b1*c3 - b3*c1)*1./(b2*c1 - b1*c2);
  
  float A, B, C;
  A = (e1*e1 + e2*e2+ 1);
  B = 2*(d1*e1 + d2*e2);
  C = d1*d1 + d2*d2 - 1;
  
  float sol;
  
  if((pow(B,2) - (4*A*C)) < 0){
    sol = -1*B/(2*A);
    
    float A1, A2, A3;
    A3 = sol;
    A1 = d1 + e1*A3;
    A2 = d2 + e2*A3;
    
    TLorentzVector vec4;
    vec4.SetPxPyPzE(A1,A2,A3,0);
    return vec4;
  }
  else{
    float sol1 = (-1*B+sqrt((pow(B,2) - (4*A*C))))*1./(2*A);
    float sol2 =  (-1*B-sqrt((pow(B,2) - (4*A*C))))*1./(2*A);
    (sol1>sol2)?sol=sol1:sol=sol2;
    
    float A1, A2, A3;
    A3 = sol;
    A1 = d1 + e1*A3;
    A2 = d2 + e2*A3;
    
    TLorentzVector vec4;
    vec4.SetPxPyPzE(A1,A2,A3,0);
    return vec4;;
  }
}

struct JetIDVars
{
  float NHF, NEMF, MUF, CHF, CEMF;
  int NumConst, NumNeutralParticle, CHM;
};

bool getJetID(JetIDVars vars, string jettype="CHS", int year=2018, double eta=0, bool tightLepVeto=true, bool UltraLegacy=false){
  
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  
  if (jettype!="CHS" && jettype!="PUPPI"){
    cout<<"Don't know your jet type! I know only CHS & PUPPI :D"<<endl;
    return false;
  }
  
  float NHF, NEMF, MUF, CHF, CEMF;
  int NumConst, NumNeutralParticle, CHM;
  
  NHF = vars.NHF; 
  NEMF = vars.NEMF;
  MUF = vars.MUF;
  CHF = vars.CHF;
  CEMF = vars.CEMF;
  NumConst = vars.NumConst;
  NumNeutralParticle = vars.NumNeutralParticle;
  CHM = vars.CHM;
  
  bool JetID = false;
  
  if(!UltraLegacy){
    
    if(year==2018 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10));
    }
    
    if(year==2018 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.6 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)<=3.0 && NHF<0.99) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>2 && NumNeutralParticle<15));
    }
    
    if(year==2017 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.02 && NEMF<0.99 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>10));
    }
    
    if(year==2017 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) ||
 (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 &&  NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NHF<0.99) || (fabs(eta)>3.0 && NEMF<0.90 && NHF>0.02 && NumNeutralParticle>2 && NumNeutralParticle<15));
    }

    if(year==2016 && jettype=="CHS"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.90 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CEMF<0.99 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9  && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.01 && NHF<0.98 && NumNeutralParticle>2) || (fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>10));
	}
    
    if(year==2016 && jettype=="PUPPI"){
      
      JetID = ( (fabs(eta)<=2.4 && CEMF<0.9 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)<=2.4 && CEMF<0.99 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || (fabs(eta)>2.4 && fabs(eta)<=2.7 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ));
      if(fabs(eta)>2.7) { JetID = false; }
	}
  }
  
  else {
    
    if(year==2017||year==2018){
      
      if(jettype=="CHS"){
	
	JetID = ( fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && CHM>0 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CHM>0 && NEMF<0.99 && NHF < 0.9 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.01 && NEMF<0.99 && NumNeutralParticle>1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10) ;
      }
      
      if(jettype=="PUPPI"){
	
	JetID =  ( fabs(eta)<=2.6 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.6 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && CEMF<0.8 && NEMF<0.99 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)>2.6 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NHF<0.9999 ) ||( fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>2 ) ;
      }
      // there is a inconsistency between table & lines in https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
      // table is chosen as it is consistent with the slides https://indico.cern.ch/event/937597/contributions/3940302/attachments/2073315/3481068/ULJetID_UL17_UL18_AK4PUPPI.pdf 
    }
    
    if(year==2016){
      
      if(jettype=="CHS"){
	
	JetID =  ( fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.4 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.9 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10) ;
	
      }
      
      if(jettype=="PUPPI"){
	
	JetID = ( fabs(eta)<=2.4 && CEMF<0.8 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && MUF <0.8 && NHF < 0.9 && tightLepVeto ) || ( fabs(eta)<=2.4 && CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 && !tightLepVeto ) || ( fabs(eta)>2.4 && fabs(eta)<=2.7 && NEMF<0.99 && NHF < 0.98 ) || ( fabs(eta)>2.7 && fabs(eta)<=3.0 && NumNeutralParticle>=1 ) || ( fabs(eta)>3.0 && NEMF<0.90 && NumNeutralParticle>2  ) ;
      }
    }	
  }
  
  return JetID;
  
}

bool Muon_Tight_ID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst, float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Tight_Muon
	bool tightid = false;
	if(muonisGL && muonisPF){
		if(muonchi<10 && muonhit>0 && muonmst>1){
			if(fabs(muontrkvtx)<0.2 && fabs(muondz)<0.5){
				if(muonpixhit>0 && muontrklay>5){
					tightid = true;
				}
			}
		}
	}
	return tightid;
}

bool StoreMuon(pat::Muon muon1, float ptcut, float etacut){
	
	if (((muon1.isTrackerMuon() || muon1.isGlobalMuon()) && (muon1.isPFMuon())) && (muon1.pt()>=ptcut) && (fabs(muon1.eta())<=etacut)) {                                                                
			return true;
	}
	else{
			return false;
		}
}

bool StoreElectron(pat::Electron electron1, float ptcut, float etacut){
	
	GsfTrackRef gsftrk1 = electron1.gsfTrack();                                                                                                      
    if ((!gsftrk1.isNull()) && (electron1.pt()>=ptcut) && (fabs(electron1.eta())<=etacut) && (gsftrk1->ndof()>=9)) {
			return true;
		}
    else{
			return false;
		}
}

TLorentzVector LeptonJet_subtraction(vector<auto> leps, pat::Jet jet, TLorentzVector jet4v){
	
	TLorentzVector newjet4v;
	newjet4v = jet4v;
	
	if (leps.size()>0) {                                                                                           
		for (unsigned int ilep = 0; ilep<leps.size(); ilep++) {
	  
			bool lepmember = false;
			
			for(unsigned int jd = 0 ; jd < leps[ilep].numberOfSourceCandidatePtrs() ; ++jd) {
				
				if(leps[ilep].sourceCandidatePtr(jd).isNonnull() && leps[ilep].sourceCandidatePtr(jd).isAvailable()){
					const reco::Candidate* jcand = leps[ilep].sourceCandidatePtr(jd).get();
				
					for(unsigned int ic = 0 ; ic < jet.numberOfSourceCandidatePtrs() ; ++ic) {  
					
						if(jet.sourceCandidatePtr(ic).isNonnull() && jet.sourceCandidatePtr(ic).isAvailable()){
							const reco::Candidate* icand = jet.sourceCandidatePtr(ic).get();
							if (delta2R(jcand->eta(),jcand->phi(),icand->eta(),icand->phi()) < 0.00001)    
							{
								TLorentzVector tmpvec(jcand->px(),jcand->py(),jcand->pz(),jcand->energy());
								newjet4v = jet4v - tmpvec;
								lepmember = true; 
								break;
							}
						}
					}		
				
				}
								
			}//jd
			
			if(lepmember) break;
			
		}//ilep
	}
    
    return newjet4v;
}

void Read_JEC(double &total_JEC,  double &tmprecpt, 
			  double jeteta, double Rho, bool isData,
			  pat::Jet jet,
			  FactorizedJetCorrector *jecL1Fast, FactorizedJetCorrector *jecL2Relative, FactorizedJetCorrector *jecL3Absolute, FactorizedJetCorrector*jecL2L3Residual)
{
	
    double total_cor =1;
      
    jecL1Fast->setJetPt(tmprecpt); jecL1Fast->setJetA(jet.jetArea()); jecL1Fast->setRho(Rho);jecL1Fast->setJetEta(jeteta);
    double corFactorL1Fast = jecL1Fast->getCorrection();
    total_cor *= corFactorL1Fast;
    tmprecpt = tmprecpt * corFactorL1Fast;
      
    jecL2Relative->setJetPt(tmprecpt); jecL2Relative->setJetEta(jeteta);
    double corFactorL2Relative = jecL2Relative->getCorrection();
    total_cor *= corFactorL2Relative ;
    tmprecpt = tmprecpt * corFactorL2Relative;
      
    jecL3Absolute->setJetPt(tmprecpt); jecL3Absolute->setJetEta(jeteta);
    double corFactorL3Absolute = jecL3Absolute->getCorrection();
    total_cor *= corFactorL3Absolute ;
    tmprecpt = tmprecpt * corFactorL3Absolute;
      
    double corFactorL2L3Residual=1.;
      
    if(isData){
		jecL2L3Residual->setJetPt(tmprecpt); jecL2L3Residual->setJetEta(jeteta);
		corFactorL2L3Residual = jecL2L3Residual->getCorrection();
		total_cor*= corFactorL2L3Residual;
		tmprecpt *=corFactorL2L3Residual;
	}
	
	total_JEC = total_cor;
	
	return;     
}

void Read_JER(std::string mPtResoFile, std::string mPtSFFile, double tmprecpt, TLorentzVector pfjet4v, double Rho, edm::Handle<reco::GenJetCollection>  genjets, double dRcut, vector<double> &SFs)
{
 
	JME::JetResolution resolution;
	resolution = JME::JetResolution(mPtResoFile.c_str());
	JME::JetResolutionScaleFactor res_sf;
	res_sf = JME::JetResolutionScaleFactor(mPtSFFile.c_str());
	
	JME::JetParameters parameters_5 = {{JME::Binning::JetPt, tmprecpt}, {JME::Binning::JetEta, pfjet4v.Eta()}, {JME::Binning::Rho, Rho}};
	double rp = resolution.getResolution(parameters_5);
	double gaus_rp = gRandom->Gaus(0.,rp);
	double sf = res_sf.getScaleFactor(parameters_5, Variation::NOMINAL);
	double sf_up = res_sf.getScaleFactor(parameters_5, Variation::UP);
	double sf_dn = res_sf.getScaleFactor(parameters_5, Variation::DOWN);
	
	bool match = false;
	int match_gen = -1;
		
	for (unsigned get = 0; get<(genjets->size()); get++) {
		TLorentzVector genjet4v((*genjets)[get].px(),(*genjets)[get].py(),(*genjets)[get].pz(), (*genjets)[get].energy());
		if((delta2R(pfjet4v.Rapidity(),pfjet4v.Phi(),genjet4v.Rapidity(),genjet4v.Phi()) < (dRcut)) &&(fabs(tmprecpt-genjet4v.Pt())<(3*fabs(rp)*tmprecpt))){
			match = true;
			match_gen = get;
			break;
		}
	}
		
	if(match && (match_gen>=0)){
	  
		SFs.push_back((sf-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_up-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
		SFs.push_back((sf_dn-1.)*(tmprecpt-(*genjets)[match_gen].pt())*1./tmprecpt);
	  
	}else{
	  
		SFs.push_back(sqrt(max(0.,(sf*sf-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_up*sf_up-1))) * gaus_rp);
		SFs.push_back(sqrt(max(0.,(sf_dn*sf_dn-1))) * gaus_rp);
	}
      	
}

float getEtaForEA(auto obj){
	float eta;
	if(abs(obj->pdgId())==11||abs(obj->pdgId())==22) { eta = obj->superCluster()->eta(); }     
	else { eta = obj->eta(); }
	return eta;    
}

std::unique_ptr<EffectiveAreas> ea_mu_miniiso_, ea_el_miniiso_;

void Read_MiniIsolation(auto obj, double Rho, vector<float> &isovalues)
{
	pat::PFIsolation iso = obj->miniPFIsolation();                                                                                                                                                                                                   
	float chg = iso.chargedHadronIso();                                                                                                                     
	float neu = iso.neutralHadronIso();                                                                                                                     
	float pho = iso.photonIso();                                                                                       
	                                                                                   
	float ea;
	if(abs(obj->pdgId())==13) { ea = ea_mu_miniiso_->getEffectiveArea(fabs(getEtaForEA(obj))); }
	else { ea = ea_el_miniiso_->getEffectiveArea(fabs(getEtaForEA(obj))); }  
	                                                                                    
	float R = 10.0/std::min(std::max(obj->pt(), 50.0),200.0);                                                                      
	ea *= std::pow(R / 0.3, 2);                                                                                                                  	
	float tot = (chg+std::max(0.0,neu+pho-(Rho)*ea));
	
	isovalues.push_back(tot);
	isovalues.push_back(chg);
	isovalues.push_back(neu);
	isovalues.push_back(pho);	
	
	for(unsigned ij=0; ij<isovalues.size(); ij++){
		isovalues[ij] *= 1./obj->pt();
	}
}

std::unique_ptr<EffectiveAreas> ea_el_pfiso_;

void Read_ElePFIsolation(auto obj, double Rho, vector<float> &isovalues)
{
	auto iso = obj->pfIsolationVariables();   
	auto  ea = ea_el_pfiso_->getEffectiveArea(fabs(getEtaForEA(obj)));                                                    
    float val = iso.sumChargedHadronPt + max(0., iso.sumNeutralHadronEt + iso.sumPhotonEt - (Rho)*ea); 
    float val04 = (obj->chargedHadronIso()+std::max(0.0,obj->neutralHadronIso()+obj->photonIso()-(Rho)*ea*16./9.));
    isovalues.push_back(val);
    isovalues.push_back(val04);
    
    for(unsigned ij=0; ij<isovalues.size(); ij++){
		isovalues[ij] *= 1./obj->pt();
	}    
}

//class declaration
//
class Leptop : public edm::EDAnalyzer {
public:
  explicit Leptop(const edm::ParameterSet&);
  ~Leptop();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  void fillmetarray();
  // ----------member data ---------------------------
  int Nevt;
  bool isData;
  bool isMC;
  bool isFastSIM;
  int year;
  bool isUltraLegacy;
  bool isSoftDrop;
  bool add_prefireweights;
  bool store_electron_scalnsmear, store_electron_addvariabs;
  bool store_fatjet_constituents;
  bool read_btagSF;
  bool subtractLepton_fromAK4, subtractLepton_fromAK8;
  
  uint nPDFsets;
  
  std::string theRootFileName;
  std::string theHLTTag;
  std::string softdropmass;
  std::string Nsubjettiness_tau1;
  std::string Nsubjettiness_tau2;
  std::string Nsubjettiness_tau3;
  std::string subjets;
  std::string toptagger_DAK8;
  std::string Wtagger_DAK8;
  std::string Ztagger_DAK8;
  std::string Htagger_DAK8;
  std::string bbtagger_DAK8;
  std::string toptagger_PNet;
  std::string Wtagger_PNet;
  std::string Ztagger_PNet;
  std::string Xbbtagger_PNet;
  std::string Xcctagger_PNet;
  std::string Xqqtagger_PNet;
  
  edm::EDGetTokenT<double> tok_Rho_;
  edm::EDGetTokenT<reco::BeamSpot> tok_beamspot_;
  edm::EDGetTokenT<reco::VertexCollection> tok_primaryVertices_;
  edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> tok_sv;
  edm::EDGetTokenT<pat::METCollection>tok_mets_, tok_mets_PUPPI_;
  edm::EDGetTokenT<edm::View<pat::Jet>>tok_pfjetAK8s_;
  edm::EDGetTokenT<edm::View<pat::Jet>>tok_pfjetAK4s_;
  edm::EDGetTokenT<edm::View<pat::Jet>>tok_pfjetAK4sB_;
  edm::EDGetTokenT<edm::View<pat::Muon>> tok_muons_;
  edm::EDGetTokenT<edm::View<pat::Electron>> tok_electrons_;
  edm::EDGetTokenT<edm::View<pat::Photon>>tok_photons_;
  edm::EDGetTokenT<edm::View<pat::Tau>>tok_taus_;
  
  //std::unique_ptr<EffectiveAreas> ea_miniiso_;
  
  edm::EDGetTokenT<reco::GenMETCollection>tok_genmets_;
  edm::EDGetTokenT<reco::GenJetCollection>tok_genjetAK8s_;
  edm::EDGetTokenT<reco::GenJetCollection>tok_genjetAK4s_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>>tok_genparticles_;
  edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
  
  edm::EDGetTokenT<HepMCProduct> tok_HepMC ;
  edm::EDGetTokenT<GenEventInfoProduct> tok_wt_;
  edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileup_;
  
  //edm::EDGetTokenT <edm::ValueMap <bool> > tok_mvaPhoID_FallV2_WP90;
  //edm::EDGetTokenT <edm::ValueMap <bool> > tok_mvaPhoID_FallV2_WP80;
  edm::EDGetTokenT <edm::ValueMap <float> > tok_mvaPhoID_FallV2_raw;

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
  
  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
  
  // object cuts //
  int iTag;
  int iTagMET;
  
  double minjPt;
  double minGenPt;
  double maxEta;
  double maxgenEta;
  double AK8PtCut;
  double AK8GenPtCut;
  double minmuPt;
  double minePt;
  double mingmPt;
  double mintauPt;
  
  double beta ;
  double z_cut;     
  
  // Root file & tree //
  
  TFile* theFile;
  
  TTree* T1;
  TTree* T2;
  
  // HLTConfigProvider hltConfig_;
  
  unsigned ievt;
  
  static const int njetmx = 100; 
  static const int njetmxAK8 =100;
  static const int npartmx = 100; 
  static const int nconsmax = 1000; 
  static const int njetconsmax = 3; 
    
  int irunold;
  int irun, ilumi, ifltr, ibrnch;
  
  int nprim, npvert, PV_npvsGood, PV_ndof;
  float PV_x, PV_y, PV_z, PV_chi2;
  
  double Generator_weight;
  double weights[njetmx];
  
  double Rho ;
  
  double prefiringweight, prefiringweightup, prefiringweightdown;
  
  int nPFJetAK8;
  float PFJetAK8_pt[njetmxAK8], PFJetAK8_y[njetmxAK8], PFJetAK8_eta[njetmxAK8], PFJetAK8_phi[njetmxAK8], PFJetAK8_mass[njetmxAK8];
  
  bool PFJetAK8_jetID_tightlepveto[njetmxAK8], PFJetAK8_jetID[njetmxAK8];
  
  float PFJetAK8_btag_DeepCSV[njetmxAK8];
  float PFJetAK8_DeepTag_DAK8_TvsQCD[njetmxAK8], PFJetAK8_DeepTag_DAK8_WvsQCD[njetmxAK8], PFJetAK8_DeepTag_DAK8_ZvsQCD[njetmxAK8], PFJetAK8_DeepTag_DAK8_HvsQCD[njetmxAK8], PFJetAK8_DeepTag_DAK8_bbvsQCD[njetmxAK8]; 
  float PFJetAK8_DeepTag_PNet_TvsQCD[njetmxAK8], PFJetAK8_DeepTag_PNet_WvsQCD[njetmxAK8], PFJetAK8_DeepTag_PNet_ZvsQCD[njetmxAK8], PFJetAK8_DeepTag_PNet_XbbvsQCD[njetmxAK8], PFJetAK8_DeepTag_PNet_XccvsQCD[njetmxAK8], PFJetAK8_DeepTag_PNet_XqqvsQCD[njetmxAK8], PFJetAK8_DeepTag_PNet_QCD[njetmxAK8]; 
  
  float PFJetAK8_CHF[njetmxAK8], PFJetAK8_NHF[njetmxAK8], PFJetAK8_MUF[njetmxAK8], PFJetAK8_PHF[njetmxAK8], PFJetAK8_CEMF[njetmxAK8], PFJetAK8_NEMF[njetmxAK8], PFJetAK8_EEF[njetmxAK8], PFJetAK8_HFHF[njetmxAK8], /*PFJetAK8_HFEMF[njetmxAK8],*/ PFJetAK8_HOF[njetmxAK8];
  int PFJetAK8_CHM[njetmxAK8], PFJetAK8_NHM[njetmxAK8], PFJetAK8_MUM[njetmxAK8], PFJetAK8_PHM[njetmxAK8], PFJetAK8_Neucons[njetmxAK8], PFJetAK8_Chcons[njetmxAK8], PFJetAK8_EEM[njetmxAK8], PFJetAK8_HFHM[njetmxAK8];// PFJetAK8_HFEMM[njetmxAK8];
  
  float PFJetAK8_chrad[njetmxAK8], PFJetAK8_pTD[njetmxAK8]; 
  float PFJetAK8_sdmass[njetmxAK8], PFJetAK8_tau1[njetmxAK8], PFJetAK8_tau2[njetmxAK8], PFJetAK8_tau3[njetmxAK8];
  
  float PFJetAK8_sub1pt[njetmxAK8], PFJetAK8_sub1eta[njetmxAK8], PFJetAK8_sub1phi[njetmxAK8], PFJetAK8_sub1mass[njetmxAK8], PFJetAK8_sub1JEC[njetmxAK8], PFJetAK8_sub1btag[njetmxAK8]; 
  float PFJetAK8_sub2pt[njetmxAK8], PFJetAK8_sub2eta[njetmxAK8], PFJetAK8_sub2phi[njetmxAK8], PFJetAK8_sub2mass[njetmxAK8], PFJetAK8_sub2JEC[njetmxAK8], PFJetAK8_sub2btag[njetmxAK8];
  
  float PFJetAK8_JEC[njetmxAK8];
  float PFJetAK8_reso[njetmxAK8], PFJetAK8_resoup[njetmxAK8], PFJetAK8_resodn[njetmxAK8];
  float PFJetAK8_jesup_pu[njetmx], PFJetAK8_jesup_rel[njetmx], PFJetAK8_jesup_scale[njetmx], PFJetAK8_jesup_total[njetmx], PFJetAK8_jesdn_pu[njetmx], PFJetAK8_jesdn_rel[njetmx], PFJetAK8_jesdn_scale[njetmx], PFJetAK8_jesdn_total[njetmx];

  float PFJetAK8_jesup_AbsoluteStat[njetmxAK8], PFJetAK8_jesdn_AbsoluteStat[njetmxAK8];
  float PFJetAK8_jesup_AbsoluteScale[njetmxAK8], PFJetAK8_jesdn_AbsoluteScale[njetmxAK8];
  float PFJetAK8_jesup_AbsoluteMPFBias[njetmxAK8], PFJetAK8_jesdn_AbsoluteMPFBias[njetmxAK8];
  float PFJetAK8_jesup_FlavorQCD[njetmxAK8], PFJetAK8_jesdn_FlavorQCD[njetmxAK8];
  float PFJetAK8_jesup_Fragmentation[njetmxAK8], PFJetAK8_jesdn_Fragmentation[njetmxAK8];
  float PFJetAK8_jesup_PileUpDataMC[njetmxAK8], PFJetAK8_jesdn_PileUpDataMC[njetmxAK8];
  float PFJetAK8_jesup_PileUpPtBB[njetmxAK8], PFJetAK8_jesdn_PileUpPtBB[njetmxAK8];
  float PFJetAK8_jesup_PileUpPtEC1[njetmxAK8], PFJetAK8_jesdn_PileUpPtEC1[njetmxAK8];
  float PFJetAK8_jesup_PileUpPtEC2[njetmxAK8], PFJetAK8_jesdn_PileUpPtEC2[njetmxAK8];
  float PFJetAK8_jesup_PileUpPtRef[njetmxAK8], PFJetAK8_jesdn_PileUpPtRef[njetmxAK8];
  float PFJetAK8_jesup_RelativeFSR[njetmxAK8], PFJetAK8_jesdn_RelativeFSR[njetmxAK8];
  float PFJetAK8_jesup_RelativeJEREC1[njetmxAK8], PFJetAK8_jesdn_RelativeJEREC1[njetmxAK8];
  float PFJetAK8_jesup_RelativeJEREC2[njetmxAK8], PFJetAK8_jesdn_RelativeJEREC2[njetmxAK8];
  float PFJetAK8_jesup_RelativePtBB[njetmxAK8], PFJetAK8_jesdn_RelativePtBB[njetmxAK8];
  float PFJetAK8_jesup_RelativePtEC1[njetmxAK8], PFJetAK8_jesdn_RelativePtEC1[njetmxAK8];
  float PFJetAK8_jesup_RelativePtEC2[njetmxAK8], PFJetAK8_jesdn_RelativePtEC2[njetmxAK8];
  float PFJetAK8_jesup_RelativeBal[njetmxAK8], PFJetAK8_jesdn_RelativeBal[njetmxAK8];
  float PFJetAK8_jesup_RelativeSample[njetmxAK8], PFJetAK8_jesdn_RelativeSample[njetmxAK8];
  float PFJetAK8_jesup_RelativeStatEC[njetmxAK8], PFJetAK8_jesdn_RelativeStatEC[njetmxAK8];
  float PFJetAK8_jesup_RelativeStatFSR[njetmxAK8], PFJetAK8_jesdn_RelativeStatFSR[njetmxAK8];
  float PFJetAK8_jesup_SinglePionECAL[njetmxAK8], PFJetAK8_jesdn_SinglePionECAL[njetmxAK8];
  float PFJetAK8_jesup_SinglePionHCAL[njetmxAK8], PFJetAK8_jesdn_SinglePionHCAL[njetmxAK8];
  float PFJetAK8_jesup_TimePtEta[njetmxAK8], PFJetAK8_jesdn_TimePtEta[njetmxAK8];
  float PFJetAK8_jesup_Total[njetmxAK8], PFJetAK8_jesdn_Total[njetmxAK8];
  
  int nPFJetAK8_cons;
  float PFJetAK8_cons_pt[nconsmax], PFJetAK8_cons_eta[nconsmax], PFJetAK8_cons_phi[nconsmax], PFJetAK8_cons_mass[nconsmax];
  int PFJetAK8_cons_jetIndex[nconsmax], PFJetAK8_cons_pdgId[nconsmax];
  
  int nPFJetAK4;
  float PFJetAK4_pt[njetmx], PFJetAK4_pt_reg[njetmx], PFJetAK4_eta[njetmx], PFJetAK4_y[njetmx], PFJetAK4_phi[njetmx], PFJetAK4_mass[njetmx], PFJetAK4_energy[njetmx], PFJetAK4_eta_reg[njetmx], PFJetAK4_y_reg[njetmx], PFJetAK4_phi_reg[njetmx], PFJetAK4_mass_reg[njetmx], PFJetAK4_energy_reg[njetmx];
  float PFJetAK4_btag_DeepCSV[njetmx], PFJetAK4_btag_DeepFlav[njetmx]; 
  bool PFJetAK4_jetID[njetmx], PFJetAK4_jetID_tightlepveto[njetmx];
  
  float PFJetAK4_btag_DeepCSV_SF[njetmx], PFJetAK4_btag_DeepCSV_SF_up[njetmx], PFJetAK4_btag_DeepCSV_SF_dn[njetmx];
  float PFJetAK4_btag_DeepFlav_SF[njetmx], PFJetAK4_btag_DeepFlav_SF_up[njetmx], PFJetAK4_btag_DeepFlav_SF_dn[njetmx];
  
  float PFJetAK4_reso[njetmx], PFJetAK4_resoup[njetmx], PFJetAK4_resodn[njetmx];
  float PFJetAK4_bcorr[njetmx], PFJetAK4_breso[njetmx];
  float PFJetAK4_JEC[njetmx];
  
  int PFJetAK4_hadronflav[njetmx], PFJetAK4_partonflav[njetmx];
  int PFJetAK4_Ncons[njetmx];
  float PFJetAK4_qgl[njetmx], PFJetAK4_PUID[njetmx];
  
  float PFJetAK4_jesup_AbsoluteStat[njetmx], PFJetAK4_jesdn_AbsoluteStat[njetmx];
  float PFJetAK4_jesup_AbsoluteScale[njetmx], PFJetAK4_jesdn_AbsoluteScale[njetmx];
  float PFJetAK4_jesup_AbsoluteMPFBias[njetmx], PFJetAK4_jesdn_AbsoluteMPFBias[njetmx];
  float PFJetAK4_jesup_FlavorQCD[njetmx], PFJetAK4_jesdn_FlavorQCD[njetmx];
  float PFJetAK4_jesup_Fragmentation[njetmx], PFJetAK4_jesdn_Fragmentation[njetmx];
  float PFJetAK4_jesup_PileUpDataMC[njetmx], PFJetAK4_jesdn_PileUpDataMC[njetmx];
  float PFJetAK4_jesup_PileUpPtBB[njetmx], PFJetAK4_jesdn_PileUpPtBB[njetmx];
  float PFJetAK4_jesup_PileUpPtEC1[njetmx], PFJetAK4_jesdn_PileUpPtEC1[njetmx];
  float PFJetAK4_jesup_PileUpPtEC2[njetmx], PFJetAK4_jesdn_PileUpPtEC2[njetmx];
  float PFJetAK4_jesup_PileUpPtRef[njetmx], PFJetAK4_jesdn_PileUpPtRef[njetmx];
  float PFJetAK4_jesup_RelativeFSR[njetmx], PFJetAK4_jesdn_RelativeFSR[njetmx];
  float PFJetAK4_jesup_RelativeJEREC1[njetmx], PFJetAK4_jesdn_RelativeJEREC1[njetmx];
  float PFJetAK4_jesup_RelativeJEREC2[njetmx], PFJetAK4_jesdn_RelativeJEREC2[njetmx];
  float PFJetAK4_jesup_RelativePtBB[njetmx], PFJetAK4_jesdn_RelativePtBB[njetmx];
  float PFJetAK4_jesup_RelativePtEC1[njetmx], PFJetAK4_jesdn_RelativePtEC1[njetmx];
  float PFJetAK4_jesup_RelativePtEC2[njetmx], PFJetAK4_jesdn_RelativePtEC2[njetmx];
  float PFJetAK4_jesup_RelativeBal[njetmx], PFJetAK4_jesdn_RelativeBal[njetmx];
  float PFJetAK4_jesup_RelativeSample[njetmx], PFJetAK4_jesdn_RelativeSample[njetmx];
  float PFJetAK4_jesup_RelativeStatEC[njetmx], PFJetAK4_jesdn_RelativeStatEC[njetmx];
  float PFJetAK4_jesup_RelativeStatFSR[njetmx], PFJetAK4_jesdn_RelativeStatFSR[njetmx];
  float PFJetAK4_jesup_SinglePionECAL[njetmx], PFJetAK4_jesdn_SinglePionECAL[njetmx];
  float PFJetAK4_jesup_SinglePionHCAL[njetmx], PFJetAK4_jesdn_SinglePionHCAL[njetmx];
  float PFJetAK4_jesup_TimePtEta[njetmx], PFJetAK4_jesdn_TimePtEta[njetmx];
  float PFJetAK4_jesup_Total[njetmx], PFJetAK4_jesdn_Total[njetmx];
  
  static const int ngenjetAK8mx =100;
  
  int nGenJetAK8;
  float GenJetAK8_pt[njetmxAK8], GenJetAK8_eta[njetmxAK8], GenJetAK8_phi[njetmxAK8], GenJetAK8_mass[njetmxAK8], GenJetAK8_sdmass[njetmxAK8]; 
  int GenJetAK8_hadronflav[njetmxAK8], GenJetAK8_partonflav[njetmxAK8];
  
  int nGenJetAK8_cons;
  float GenJetAK8_cons_pt[nconsmax], GenJetAK8_cons_eta[nconsmax], GenJetAK8_cons_phi[nconsmax], GenJetAK8_cons_mass[nconsmax];
  int GenJetAK8_cons_jetIndex[nconsmax], GenJetAK8_cons_pdgId[nconsmax];
  
  int nGenJetAK4;
  float GenJetAK4_pt[njetmx], GenJetAK4_eta[njetmx], GenJetAK4_phi[njetmx], GenJetAK4_mass[njetmx];
  int GenJetAK4_hadronflav[njetmx], GenJetAK4_partonflav[njetmx];
  
  int nGenPart;
  int GenPart_status[npartmx], GenPart_pdg[npartmx], GenPart_mompdg[npartmx], GenPart_momstatus[npartmx], GenPart_grmompdg[npartmx], GenPart_momid[npartmx], GenPart_daugno[npartmx];
  float GenPart_pt[npartmx], GenPart_eta[npartmx], GenPart_phi[npartmx], GenPart_mass[npartmx]; //GenPart_q[npartmx];
  bool GenPart_fromhard[npartmx], GenPart_fromhardbFSR[npartmx], GenPart_isPromptFinalState[npartmx], GenPart_isLastCopyBeforeFSR[npartmx], GenPart_isDirectPromptTauDecayProductFinalState[npartmx];
  
  static const int nlhemax = 10;
  int nLHEPart;
  float LHEPart_pt[nlhemax], LHEPart_eta[nlhemax], LHEPart_phi[nlhemax], LHEPart_m[nlhemax];
  int LHEPart_pdg[nlhemax];
  
  static const int nlhescalemax = 50;
  int nLHEScaleWeights;
  float LHEScaleWeights[nlhescalemax];
  
  static const int nlhepdfmax = 103; // be consistent with nPDFsets (nlhepdfmax should be >= nPDFsets)
  int nLHEPDFWeights;
  float LHEPDFWeights[nlhepdfmax];
  
  static const int nalpsmax = 3;
  int nLHEAlpsWeights;
  float LHEAlpsWeights[nalpsmax];
  
  static const int nlhepsmax = 8;
  int nLHEPSWeights;
  float LHEPSWeights[nlhepsmax];
  float wgt_isr_up, wgt_isr_dn, wgt_fsr_up, wgt_fsr_dn;
  double LHE_weight;
 

  float BTAG_SF, BTAG_jes_up, BTAG_jes_dn, BTAG_lf_up, BTAG_lf_dn, BTAG_hf_up, BTAG_hf_dn, BTAG_hfstats1_up, BTAG_hfstats1_dn, BTAG_hfstats2_up, BTAG_hfstats2_dn, BTAG_lfstats1_up, BTAG_lfstats1_dn, BTAG_lfstats2_up, BTAG_lfstats2_dn, BTAG_cferr1_up, BTAG_cferr1_dn, BTAG_cferr2_up, BTAG_cferr2_dn; 
  float miset , misphi , sumEt, misetsig;
  float miset_UnclusEup, miset_UnclusEdn;
  float misphi_UnclusEup, misphi_UnclusEdn;
  
  float miset_PUPPI , misphi_PUPPI , sumEt_PUPPI, misetsig_PUPPI;
  float miset_PUPPI_JESup, miset_PUPPI_JESdn, miset_PUPPI_JERup, miset_PUPPI_JERdn, miset_PUPPI_UnclusEup, miset_PUPPI_UnclusEdn;
  float misphi_PUPPI_JESup, misphi_PUPPI_JESdn, misphi_PUPPI_JERup, misphi_PUPPI_JERdn, misphi_PUPPI_UnclusEup, misphi_PUPPI_UnclusEdn;
  float genmiset, genmisphi, genmisetsig;
  
  int nMuon;
  
  float Muon_minchiso[njetmx], Muon_minnhiso[njetmx], Muon_minphiso[njetmx], Muon_minisoall[njetmx]; 
  float Muon_charge[njetmx], Muon_p[njetmx], Muon_pt[njetmx], Muon_eta[njetmx], Muon_phi[njetmx], Muon_e[njetmx], Muon_dz[njetmx], Muon_ip3d[njetmx], Muon_ptErr[njetmx], Muon_chi[njetmx], Muon_ecal[njetmx], Muon_hcal[njetmx]; //Muon_emiso[njetmx], Muon_hadiso[njetmx], Muon_tkpt03[njetmx], Muon_tkpt05[njetmx];
  
  float Muon_posmatch[njetmx], Muon_trkink[njetmx], Muon_segcom[njetmx], Muon_pfiso[njetmx], Muon_dxy[njetmx], Muon_dxyErr[njetmx], Muon_hit[njetmx], Muon_pixhit[njetmx], Muon_mst[njetmx], Muon_trklay[njetmx], Muon_valfrac[njetmx],Muon_dxy_sv[njetmx];
  int Muon_ndf[njetmx];
  
  bool Muon_isPF[njetmx], Muon_isGL[njetmx], Muon_isTRK[njetmx];
  bool Muon_isGoodGL[njetmx], Muon_isTight[njetmx], Muon_isHighPt[njetmx], Muon_isHighPttrk[njetmx], Muon_isMed[njetmx], Muon_isMedPr[njetmx], Muon_isLoose[njetmx], Muon_TightID[njetmx];
  
  float Muon_corrected_pt[njetmx], Muon_correctedUp_pt[njetmx], Muon_correctedDown_pt[njetmx];
  
  int nElectron;
  bool Electron_mvaid_Fallv2WP90[njetmx], Electron_mvaid_Fallv2WP90_noIso[njetmx];
  bool Electron_mvaid_Fallv2WP80[njetmx], Electron_mvaid_Fallv2WP80_noIso[njetmx];

  float Electron_charge[njetmx], Electron_pt[njetmx], Electron_eta[njetmx], Electron_phi[njetmx], Electron_e[njetmx], Electron_e_ECAL[njetmx], Electron_p[njetmx];
  float Electron_dxy[njetmx],  Electron_dxyErr[njetmx], Electron_dxy_sv[njetmx], Electron_dz[njetmx], Electron_dzErr[njetmx], Electron_ip3d[njetmx];
  float Electron_hovere[njetmx], Electron_qovrper[njetmx], Electron_chi[njetmx]; //Electron_emiso03[njetmx], Electron_hadiso03[njetmx], Electron_emiso04[njetmx], Electron_hadiso04[njetmx];
  float Electron_eoverp[njetmx], Electron_ietaieta[njetmx], Electron_etain[njetmx], Electron_phiin[njetmx], Electron_fbrem[njetmx]; 
  float Electron_nohits[njetmx], Electron_misshits[njetmx];
  float Electron_pfiso_drcor[njetmx];
  float Electron_pfiso_eacor[njetmx];
  float Electron_pfiso04_eacor[njetmx];
  int Electron_ndf[njetmx];
  
  float Electron_eccalTrkEnergyPostCorr[njetmx];
  float Electron_energyScaleValue[njetmx];
  float Electron_energyScaleUp[njetmx];
  float Electron_energyScaleDown[njetmx];
  float Electron_energySigmaValue[njetmx];
  float Electron_energySigmaUp[njetmx];
  float Electron_energySigmaDown[njetmx];
  float Electron_energyScaleStatUp[njetmx];
  float Electron_energyScaleStatDown[njetmx];
  float Electron_energyScaleSystUp[njetmx];
  float Electron_energyScaleSystDown[njetmx];
  float Electron_energyScaleGainUp[njetmx];
  float Electron_energyScaleGainDown[njetmx];
  float Electron_energySigmaRhoUp[njetmx];
  float Electron_energySigmaRhoDown[njetmx];
  float Electron_energySigmaPhiUp[njetmx];
  float Electron_energySigmaPhiDown[njetmx];
  float Electron_energySmearNrSigma[njetmx];
  float Electron_ecalEnergyPreCorr[njetmx];
  float Electron_ecalEnergyErrPreCorr[njetmx];
  float Electron_ecalEnergyPostCorr[njetmx];
  float Electron_ecalEnergyErrPostCorr[njetmx];
  
  float Electron_supcl_eta[njetmx]; 
  float Electron_supcl_phi[njetmx]; 
  float Electron_supcl_e[njetmx]; 
  float Electron_supcl_rawE[njetmx]; 
  float Electron_sigmaieta[njetmx];
  float Electron_sigmaiphi[njetmx];
  float Electron_r9full[njetmx];
  float Electron_supcl_etaw[njetmx];
  float Electron_supcl_phiw[njetmx];
  float Electron_hcaloverecal[njetmx];
  float Electron_cloctftrkn[njetmx];
  float Electron_cloctftrkchi2[njetmx];
  float Electron_e1x5bye5x5[njetmx];
  float Electron_normchi2[njetmx];
  float Electron_hitsmiss[njetmx];
  float Electron_trkmeasure[njetmx];
  float Electron_convtxprob[njetmx];
  float Electron_ecloverpout[njetmx];
  float Electron_ecaletrkmomentum[njetmx];
  float Electron_deltaetacltrkcalo[njetmx];
  float Electron_supcl_preshvsrawe[njetmx];
  bool Electron_convVeto[njetmx]; 
  
  float Electron_pfisolsumphet[njetmx];
  float Electron_pfisolsumchhadpt[njetmx];
  float Electron_pfsiolsumneuhadet[njetmx];
  float Electron_minchiso[njetmx];
  float Electron_minnhiso[njetmx];
  float Electron_minphiso[njetmx];
  float Electron_minisoall[njetmx]; 
  
  
  int nPhoton;
  bool Photon_passEveto[njetmx];
  bool Photon_PixelSeed[njetmx];
  bool Photon_mvaid_Fall17V2_WP90[njetmx];
  bool Photon_mvaid_Fall17V2_WP80[njetmx];
  float Photon_mvaid_Fall17V2_raw[njetmx];
  bool Photon_mvaid_Spring16V1_WP90[njetmx];
  bool Photon_mvaid_Spring16V1_WP80[njetmx];
  float Photon_pt[njetmx];
  float Photon_e[njetmx];
  float Photon_eta[njetmx];
  float Photon_phi[njetmx];
  float Photon_e1by9[njetmx];
  float Photon_e9by25[njetmx];
  float Photon_hadbyem[njetmx];
  float Photon_trkiso[njetmx];
  float Photon_emiso[njetmx];
  float Photon_hadiso[njetmx];
  float Photon_chhadiso[njetmx];
  float Photon_neuhadiso[njetmx];
  float Photon_PUiso[njetmx];
  float Photon_phoiso[njetmx];
  float Photon_ietaieta[njetmx];
  float Photon_energyScaleValue[njetmx];
  float Photon_energyScaleUp[njetmx];
  float Photon_energyScaleDown[njetmx];
  float Photon_energySigmaValue[njetmx];
  float Photon_energySigmaUp[njetmx];
  float Photon_energySigmaDown[njetmx];
  float Photon_energyScaleStatUp[njetmx];
  float Photon_energyScaleStatDown[njetmx];
  float Photon_energyScaleSystUp[njetmx];
  float Photon_energyScaleSystDown[njetmx];
  float Photon_energyScaleGainUp[njetmx];
  float Photon_energyScaleGainDown[njetmx];
  float Photon_energySigmaRhoUp[njetmx];
  float Photon_energySigmaRhoDown[njetmx];
  float Photon_energySigmaPhiUp[njetmx];
  float Photon_energySigmaPhiDown[njetmx];
  float Photon_energySmearNrSigma[njetmx];
  float Photon_ecalEnergyPreCorr[njetmx];
  float Photon_ecalEnergyErrPreCorr[njetmx];
  float Photon_ecalEnergyPostCorr[njetmx];
  float Photon_ecalEnergyErrPostCorr[njetmx];
  
  int nTau;
  float Tau_pt[njetmx];
  float Tau_eta[njetmx];
  float Tau_phi[njetmx];
  float Tau_e[njetmx];
  bool Tau_isPF[njetmx];
  float Tau_dxy[njetmx];
  float Tau_dz[njetmx];
  int Tau_charge[njetmx];
  int Tau_decayMode[njetmx];
  bool Tau_decayModeinding[njetmx];
  bool Tau_decayModeindingNewDMs[njetmx];

  int Tau_muiso[njetmx];
  float Tau_eiso_raw[njetmx];
  int Tau_eiso[njetmx];
  float Tau_eiso2018_raw[njetmx];
  int Tau_eiso2018[njetmx];

  float Tau_jetiso_deeptau2017v2p1_raw[njetmx];
  int Tau_jetiso_deeptau2017v2p1[njetmx];
  float Tau_eiso_deeptau2017v2p1_raw[njetmx];
  int Tau_eiso_deeptau2017v2p1[njetmx];
  float Tau_muiso_deeptau2017v2p1_raw[njetmx];
  int Tau_muiso_deeptau2017v2p1[njetmx];

  float Tau_rawiso[njetmx];
  float Tau_rawisodR03[njetmx];
  float Tau_puCorr[njetmx];
  float Tau_leadtrkpt[njetmx];
  float Tau_leadtrketa[njetmx];
  float Tau_leadtrkphi[njetmx];
  float Tau_leadtrkdxy[njetmx];
  float Tau_leadtrkdz[njetmx];

  // Trigger Info //
  
  int nTrigObj;
  float TrigObj_pt[njetmx], TrigObj_eta[njetmx],TrigObj_phi[njetmx], TrigObj_mass[njetmx];
  bool TrigObj_HLT[njetmx], TrigObj_L1[njetmx],  TrigObj_Both[njetmx];
  int  TrigObj_Ihlt[njetmx], TrigObj_pdgId[njetmx], TrigObj_type[njetmx];
  
  // Collision Info //
  
  float Generator_qscale, Generator_x1, Generator_x2, Generator_xpdf1, Generator_xpdf2, Generator_scalePDF;
  int Generator_id1, Generator_id2;
  
  int npu_vert;
  int npu_vert_true;
  
  // HL triggers //
  
  //----------------------------------------2018 triggers--------------------------------------------//
 
  static const int nHLTmx = 21;
  const char *hlt_name[nHLTmx] = {"HLT_IsoMu24_v","HLT_Mu50_v",  // single-muon triggers
								  "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v", "HLT_Ele40_WPTight_Gsf_v",  "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v", // single-electron triggers
								  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_DoubleEle33_CaloIdL_MW_v", "HLT_DoubleEle25_CaloIdL_MW_v", // double-lepton triggers
								  "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v", "HLT_PFJet500_v", "HLT_PFHT1050_v", "HLT_AK8PFJet400_TrimMass30_v", "HLT_AK8PFHT800_TrimMass50_v", // jet triggers & diphoton triggers
								  "HLT_Photon200_v", // photon trigger
								  "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v", "HLT_PFMETTypeOne140_PFMHT140_IDTight_v" // MET trigger
								  }; 
  
  bool hlt_IsoMu24;
  bool hlt_Mu50; 
  bool hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165; 
  bool hlt_Ele115_CaloIdVT_GsfTrkIdT;
  bool hlt_Ele40_WPTight_Gsf;
  bool hlt_Ele32_WPTight_Gsf;
  bool hlt_Ele28_eta2p1_WPTight_Gsf_HT150;
  bool hlt_Mu17_Mu8;
  bool hlt_Ele23_Ele12;
  bool hlt_DoubleEle33;
  bool hlt_DoubleEle25_CaloIdL_MW;
  bool hlt_Diphoton_30_18;
  bool hlt_PFJet500;
  bool hlt_HT1050;
  bool hlt_AK8PFJet400_TrimMass30;
  bool hlt_AK8PFHT800_TrimMass50;
  bool hlt_Photon200;
  bool hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  bool hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
  bool hlt_PFMETNoMu140_PFMHTNoMu140_IDTight;
  bool hlt_PFMETTypeOne140_PFMHT140_IDTight;

//---------------------------------------------2017 triggers-------------------------------------------//
/*
  static const int nHLTmx = 21;
  const char *hlt_name[nHLTmx] = {"HLT_IsoMu24_v","HLT_IsoMu27_v",  // single-muon triggers
                                                                  "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v", "HLT_Ele40_WPTight_Gsf_v",  "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v", // single-electron triggers
                                                                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_DoubleEle33_CaloIdL_MW_v", "HLT_DoubleEle25_CaloIdL_MW_v", // double-lepton triggers
                                                                  "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v", "HLT_PFJet500_v", "HLT_PFHT1050_v", "HLT_AK8PFJet400_TrimMass30_v", "HLT_AK8PFHT800_TrimMass50_v", // jet triggers & diphoton triggers
                                                                  "HLT_Photon200_v", // photon trigger
                                                                  "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v", "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v", "HLT_PFMETTypeOne140_PFMHT140_IDTight_v" // MET trigger
                                                                  };

  bool hlt_IsoMu24;
  bool hlt_IsoMu27;
  bool hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165;
  bool hlt_Ele115_CaloIdVT_GsfTrkIdT;
  bool hlt_Ele40_WPTight_Gsf;
  bool hlt_Ele32_WPTight_Gsf;
  bool hlt_Ele28_eta2p1_WPTight_Gsf_HT150;
  bool hlt_Mu17_Mu8;
  bool hlt_Ele23_Ele12;
  bool hlt_DoubleEle33;
  bool hlt_DoubleEle25_CaloIdL_MW;
  bool hlt_Diphoton_30_18;
  bool hlt_PFJet500;
  bool hlt_HT1050;
  bool hlt_AK8PFJet400_TrimMass30;
  bool hlt_AK8PFHT800_TrimMass50;
  bool hlt_Photon200;
  bool hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  bool hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60;
  bool hlt_PFMETNoMu140_PFMHTNoMu140_IDTight;
  bool hlt_PFMETTypeOne140_PFMHT140_IDTight;
*/
//---------------------------------------------2016 triggers--------------------------------------------// 
  
/*  static const int nHLTmx = 21;
  const char *hlt_name[nHLTmx] = {"HLT_IsoMu24_v","HLT_IsoTkMu24_v",  // single-muon triggers
                                                                  "HLT_Ele27_WPTight_Gsf_v","HLT_Ele115_CaloIdVT_GsfTrkIdT_v", "HLT_Ele40_WPTight_Gsf_v",  "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v", // single-electron triggers
                                                                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", "HLT_DoubleEle25_CaloIdL_MW_v", // double-lepton triggers
                                                                  "HLT_Diphoton30_18_R9IdL_AND_HE_AND_IsoCaloId_NoPixelVeto_v", "HLT_PFJet500_v", "HLT_PFHT1050_v", "HLT_AK8PFJet400_TrimMass30_v", "HLT_AK8PFHT800_TrimMass50_v", // jet triggers & diphoton triggers
                                                                  "HLT_Photon200_v", // photon trigger
                                                                  "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v" // double-lepton triggers
                                                                  };


  bool hlt_IsoMu24;
  bool hlt_IsoTkMu24;
  bool hlt_Ele27_WPTight_Gsf;
  bool hlt_Ele115_CaloIdVT_GsfTrkIdT;
  bool hlt_Ele40_WPTight_Gsf;
  bool hlt_Ele32_WPTight_Gsf;
  bool hlt_Ele28_eta2p1_WPTight_Gsf_HT150;
  bool hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8;
  bool hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
  bool hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
  bool hlt_DoubleEle25_CaloIdL_MW;
  bool hlt_Diphoton_30_18;
  bool hlt_PFJet500;
  bool hlt_HT1050;
  bool hlt_AK8PFJet400_TrimMass30;
  bool hlt_AK8PFHT800_TrimMass50;
  bool hlt_Photon200;
  bool hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
  bool hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
  bool hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
  bool hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
  
*/

  int trig_value;
  
  HLTPrescaleProvider hltPrescaleProvider_;
    
  // ---- Jet Corrector Parameter ---- //
  
  JetCorrectorParameters *L1FastAK4, *L2RelativeAK4, *L3AbsoluteAK4, *L2L3ResidualAK4;
  vector<JetCorrectorParameters> vecL1FastAK4, vecL2RelativeAK4, vecL3AbsoluteAK4, vecL2L3ResidualAK4;
  FactorizedJetCorrector *jecL1FastAK4, *jecL2RelativeAK4, *jecL3AbsoluteAK4, *jecL2L3ResidualAK4;
  
  JetCorrectorParameters *L1FastAK8, *L2RelativeAK8, *L3AbsoluteAK8, *L2L3ResidualAK8;
  vector<JetCorrectorParameters> vecL1FastAK8, vecL2RelativeAK8, vecL3AbsoluteAK8, vecL2L3ResidualAK8;
  FactorizedJetCorrector *jecL1FastAK8, *jecL2RelativeAK8, *jecL3AbsoluteAK8, *jecL2L3ResidualAK8;
  
  // ---- Jet Corrector Parameter End---- //

  // BTagCalibration Begin //

  BTagCalibration calib_deepcsv, calib_deepflav;
  BTagCalibrationReader reader_deepcsv, reader_deepflav;
  BTagCalibrationReader reader; 
  // BTagCalibration End //

  // ---- Jet Resolution Parameter ---- //
  
  std::string mJECL1FastFileAK4, mJECL2RelativeFileAK4, mJECL3AbsoluteFileAK4, mJECL2L3ResidualFileAK4, mJECL1FastFileAK8, mJECL2RelativeFileAK8, mJECL3AbsoluteFileAK8, mJECL2L3ResidualFileAK8;
  std::string mPtResoFileAK4, mPtResoFileAK8, mPtSFFileAK4, mPtSFFileAK8;
  
  std::string mJECUncFileAK4;
  std::vector<JetCorrectionUncertainty*> vsrc ;
  
  std::string mJECUncFileAK8;
  std::vector<JetCorrectionUncertainty*> vsrcAK8 ;
  
  // ---- Jet Resolution Parameter End---- //
  
  // ---- B tagging scale factor files --- //
  
  std::string mBtagSF_DeepCSV;
  std::string mBtagSF_DeepFlav;
  std::string mBtagSF_DeepFlav_itr; 
  // ---- B tagging scale factor files End --- //
  
  // ---- Rochester correction files --- //
  
  std::string mRochcorFolder;
  
  // Electron MVA ID //
  
  std::string melectronID_isowp90, melectronID_noisowp90;
  std::string melectronID_isowp80, melectronID_noisowp80;
  
  // Photon MVA ID //
  
  std::string mPhoID_FallV2_WP90, mPhoID_FallV2_WP80;
  std::string mPhoID_SpringV1_WP90, mPhoID_SpringV1_WP80;
  
  // Rochester correction for muons//
  
  RoccoR roch_cor; 
  
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

Leptop::Leptop(const edm::ParameterSet& pset):
  hltPrescaleProvider_(pset, consumesCollector(), *this)  
{
  //now do what ever initialization is needed
  
  edm::Service<TFileService> fs;
  
  isData    = pset.getUntrackedParameter<bool>("Data",false);
  isMC      = pset.getUntrackedParameter<bool>("MonteCarlo", false);
  isFastSIM      = pset.getUntrackedParameter<bool>("FastSIM", false);
  year		= pset.getUntrackedParameter<int>("YEAR", 2018);
  isUltraLegacy = pset.getUntrackedParameter<bool>("UltraLegacy", false);
  isSoftDrop      = pset.getUntrackedParameter<bool>("SoftDrop_ON",false);
  theRootFileName = pset.getUntrackedParameter<string>("RootFileName");
  theHLTTag = pset.getUntrackedParameter<string>("HLTTag", "HLT");
  add_prefireweights = pset.getUntrackedParameter<bool>("add_prefireweights", false);
  store_electron_scalnsmear = pset.getUntrackedParameter<bool>("store_electron_scalnsmear", false);
  store_electron_addvariabs = pset.getUntrackedParameter<bool>("store_electron_addvariabs", false);
  store_fatjet_constituents = pset.getUntrackedParameter<bool>("store_fatjet_constituents", false);
  read_btagSF = pset.getUntrackedParameter<bool>("Read_btagging_SF", false);
  subtractLepton_fromAK8 = pset.getUntrackedParameter<bool>("Subtract_Lepton_fromAK8", false);
  subtractLepton_fromAK4 = pset.getUntrackedParameter<bool>("Subtract_Lepton_fromAK4", false);
  
  minjPt = pset.getUntrackedParameter<double>("minjPt",20.);
  minGenPt = pset.getUntrackedParameter<double>("minGenPt",15.);
  AK8PtCut = pset.getUntrackedParameter<double>("AK8PtCut",0.);
  AK8GenPtCut = pset.getUntrackedParameter<double>("AK8GenPtCut",0.);
  minmuPt = pset.getUntrackedParameter<double>("minmuPt",10.);
  minePt = pset.getUntrackedParameter<double>("minePt",10.);
  mingmPt = pset.getUntrackedParameter<double>("mingmPt",20.);
  mintauPt = pset.getUntrackedParameter<double>("mintauPt",25.);
  
  maxEta = pset.getUntrackedParameter<double>("maxEta",3.);
  maxgenEta = pset.getUntrackedParameter<double>("maxgenEta",5.);
 
  softdropmass = pset.getUntrackedParameter<string>("softdropmass");
  beta = pset.getUntrackedParameter<double>("beta",0);
  z_cut = pset.getUntrackedParameter<double>("z_cut",0.1); 
  Nsubjettiness_tau1 = pset.getUntrackedParameter<string>("tau1");
  Nsubjettiness_tau2 = pset.getUntrackedParameter<string>("tau2");
  Nsubjettiness_tau3 = pset.getUntrackedParameter<string>("tau3");
  subjets = pset.getUntrackedParameter<string>("subjets");
  toptagger_DAK8 = pset.getUntrackedParameter<string>("toptagger_DAK8");
  Wtagger_DAK8 = pset.getUntrackedParameter<string>("Wtagger_DAK8");
  Ztagger_DAK8 = pset.getUntrackedParameter<string>("Ztagger_DAK8");
  Htagger_DAK8 = pset.getUntrackedParameter<string>("Htagger_DAK8");
  bbtagger_DAK8 = pset.getUntrackedParameter<string>("bbtagger_DAK8");  
  toptagger_PNet = pset.getUntrackedParameter<string>("toptagger_PNet");
  Wtagger_PNet = pset.getUntrackedParameter<string>("Wtagger_PNet");
  Ztagger_PNet = pset.getUntrackedParameter<string>("Ztagger_PNet");
  Xbbtagger_PNet = pset.getUntrackedParameter<string>("Xbbtagger_PNet");
  Xcctagger_PNet = pset.getUntrackedParameter<string>("Xcctagger_PNet");  
  Xqqtagger_PNet = pset.getUntrackedParameter<string>("Xqqtagger_PNet");  
  
  ea_mu_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_MuonMiniIso")).fullPath()));
  ea_el_miniiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_EleMiniIso")).fullPath()));
  ea_el_pfiso_.reset(new EffectiveAreas((pset.getParameter<edm::FileInPath>("EAFile_ElePFIso")).fullPath()));
 
  tok_beamspot_ = consumes<reco::BeamSpot> (pset.getParameter<edm::InputTag>("Beamspot"));
  tok_primaryVertices_ =consumes<reco::VertexCollection>( pset.getParameter<edm::InputTag>("PrimaryVertices"));
  tok_sv =consumes<reco::VertexCompositePtrCandidateCollection>( pset.getParameter<edm::InputTag>("SecondaryVertices"));
  tok_Rho_ = consumes<double>(pset.getParameter<edm::InputTag>("PFRho"));
     
  tok_mets_= consumes<pat::METCollection> ( pset.getParameter<edm::InputTag>("PFMet"));
  tok_mets_PUPPI_ = consumes<pat::METCollection> ( pset.getParameter<edm::InputTag>("PuppiMet"));
  tok_genmets_= consumes<reco::GenMETCollection> ( pset.getParameter<edm::InputTag>("GENMet"));
    
  tok_muons_ = consumes<edm::View<pat::Muon>> ( pset.getParameter<edm::InputTag>("Muons"));
  tok_electrons_ = consumes<edm::View<pat::Electron>> ( pset.getParameter<edm::InputTag>("Electrons"));
  tok_photons_ = consumes<edm::View<pat::Photon>>  ( pset.getParameter<edm::InputTag>("Photons"));
  tok_taus_ = consumes<edm::View<pat::Tau>>  ( pset.getParameter<edm::InputTag>("Taus"));
  
  tok_pfjetAK8s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK8"));
  tok_pfjetAK4s_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFJetsAK4"));
  tok_pfjetAK4sB_= consumes<edm::View<pat::Jet>>( pset.getParameter<edm::InputTag>("PFAK4JetsB"));
  
  if(isMC){
	  
    tok_genjetAK8s_= consumes<reco::GenJetCollection>( pset.getParameter<edm::InputTag>("GENJetAK8"));
    tok_genjetAK4s_= consumes<reco::GenJetCollection>( pset.getParameter<edm::InputTag>("GENJetAK4"));
    tok_genparticles_ = consumes<std::vector<reco::GenParticle>>( pset.getParameter<edm::InputTag>("GenParticles"));
    jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>(pset.getParameter<edm::InputTag>("jetFlavourInfos"));
    
    tok_HepMC = consumes<HepMCProduct>(pset.getParameter<edm::InputTag>("Generator"));
    tok_wt_ = consumes<GenEventInfoProduct>(pset.getParameter<edm::InputTag>("Generator")) ;
    lheEventProductToken_ = consumes<LHEEventProduct>(pset.getParameter<edm::InputTag>("LHEEventProductInputTag")) ;
    pileup_ = consumes<std::vector<PileupSummaryInfo> >(pset.getParameter<edm::InputTag>("slimmedAddPileupInfo"));
    nPDFsets      = pset.getUntrackedParameter<uint>("nPDFsets", 103);
    
  }
  
  melectronID_isowp90       = pset.getParameter<std::string>("electronID_isowp90");
  melectronID_noisowp90     = pset.getParameter<std::string>("electronID_noisowp90");
  melectronID_isowp80       = pset.getParameter<std::string>("electronID_isowp80");
  melectronID_noisowp80     = pset.getParameter<std::string>("electronID_noisowp80");
  
  mPhoID_FallV2_WP90       = pset.getParameter<std::string>("PhoID_FallV2_WP90");
  mPhoID_FallV2_WP80       = pset.getParameter<std::string>("PhoID_FallV2_WP80");
  mPhoID_SpringV1_WP90       = pset.getParameter<std::string>("PhoID_SpringV1_WP90");
  mPhoID_SpringV1_WP80       = pset.getParameter<std::string>("PhoID_SpringV1_WP80");
  
  //tok_mvaPhoID_FallV2_WP90 = consumes<edm::ValueMap <bool> >(pset.getParameter<edm::InputTag>("label_mvaPhoID_FallV2_WP90"));
  //tok_mvaPhoID_FallV2_WP80 = consumes<edm::ValueMap <bool> >(pset.getParameter<edm::InputTag>("label_mvaPhoID_FallV2_WP80"));
  tok_mvaPhoID_FallV2_raw = consumes<edm::ValueMap <float> >(pset.getParameter<edm::InputTag>("label_mvaPhoID_FallV2_Value"));
  
  mJECL1FastFileAK4         = pset.getParameter<std::string>("jecL1FastFileAK4");
  mJECL1FastFileAK8         = pset.getParameter<std::string>("jecL1FastFileAK8");
  mJECL2RelativeFileAK4     = pset.getParameter<std::string>("jecL2RelativeFileAK4");
  mJECL2RelativeFileAK8     = pset.getParameter<std::string>("jecL2RelativeFileAK8");
  mJECL3AbsoluteFileAK4     = pset.getParameter<std::string>("jecL3AbsoluteFileAK4");
  mJECL3AbsoluteFileAK8     = pset.getParameter<std::string> ("jecL3AbsoluteFileAK8");
  mJECL2L3ResidualFileAK4   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK4");
  mJECL2L3ResidualFileAK8   = pset.getParameter<std::string> ("jecL2L3ResidualFileAK8");
  
  mPtResoFileAK4  = pset.getParameter<std::string>("PtResoFileAK4");
  mPtResoFileAK8  = pset.getParameter<std::string>("PtResoFileAK8");
  mPtSFFileAK4  = pset.getParameter<std::string>("PtSFFileAK4");
  mPtSFFileAK8  = pset.getParameter<std::string>("PtSFFileAK8");
  
  mJECUncFileAK4 = pset.getParameter<std::string>("JECUncFileAK4");
  mJECUncFileAK8 = pset.getParameter<std::string>("JECUncFileAK8");               
  
  mBtagSF_DeepCSV = pset.getParameter<std::string>("BtagSFFile_DeepCSV");
  mBtagSF_DeepFlav = pset.getParameter<std::string>("BtagSFFile_DeepFlav");
  mBtagSF_DeepFlav_itr = pset.getParameter<std::string>("BtagSFFile_DeepFlav_Iter"),

  mRochcorFolder = pset.getParameter<std::string>("RochcorFolder");
  
  if(!isFastSIM){
	triggerBits_ = consumes<edm::TriggerResults> ( pset.getParameter<edm::InputTag>("bits"));
	triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(pset.getParameter<edm::InputTag>("objects"));
	triggerPrescales_ = consumes<pat::PackedTriggerPrescales>(pset.getParameter<edm::InputTag>("prescales"));
  }
  
  if(add_prefireweights){
	prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
	prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
	prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));
	}
  
  theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  theFile->cd();
 
  T1 = new TTree("Events", "XtoYH");
 
  T1->Branch("irun", &irun, "irun/I");  
  T1->Branch("ilumi", &ilumi, "ilumi/I");  
  T1->Branch("ievt", &ievt, "ievt/i");
  
  // primary vertices //
  
  T1->Branch("nprim", &nprim, "nprim/I");
  T1->Branch("npvert", &npvert, "npvert/I");
  T1->Branch("PV_npvsGood", &PV_npvsGood, "PV_npvsGood/I");
  T1->Branch("PV_ndof", &PV_ndof, "PV_ndof/I");
  T1->Branch("PV_chi2", &PV_chi2, "PV_chi2/F");
  T1->Branch("PV_x", &PV_x, "PV_x/F");
  T1->Branch("PV_y", &PV_y, "PV_y/F");
  T1->Branch("PV_z", &PV_z, "PV_z/F");
  
  // energy density //
  
  T1->Branch("Rho", &Rho, "Rho/D") ;
  
  // trigger info //
  
  T1->Branch("trig_value",&trig_value,"trig_value/I");  
  
//-----------------------------------------------2018 triggers-------------------------------------------------//
  T1->Branch("hlt_IsoMu24",&hlt_IsoMu24,"hlt_IsoMu24/O");
  T1->Branch("hlt_Mu50",&hlt_Mu50,"hlt_Mu50/O");
  T1->Branch("hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",&hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165,"hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165/O");
  T1->Branch("hlt_Ele115_CaloIdVT_GsfTrkIdT",&hlt_Ele115_CaloIdVT_GsfTrkIdT,"hlt_Ele115_CaloIdVT_GsfTrkIdT/O");
  T1->Branch("hlt_Ele40_WPTight_Gsf",&hlt_Ele40_WPTight_Gsf,"hlt_Ele40_WPTight_Gsf/O");
  T1->Branch("hlt_Ele32_WPTight_Gsf",&hlt_Ele32_WPTight_Gsf,"hlt_Ele32_WPTight_Gsf/O");
  T1->Branch("hlt_Ele28_eta2p1_WPTight_Gsf_HT150",&hlt_Ele28_eta2p1_WPTight_Gsf_HT150,"hlt_Ele28_eta2p1_WPTight_Gsf_HT150/O");
  T1->Branch("hlt_Mu17_Mu8",&hlt_Mu17_Mu8,"hlt_Mu17_Mu8/O");
  T1->Branch("hlt_Ele23_Ele12",&hlt_Ele23_Ele12,"hlt_Ele23_Ele12/O");
  T1->Branch("hlt_DoubleEle33",&hlt_DoubleEle33,"hlt_DoubleEle33/O");
  T1->Branch("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW,"hlt_DoubleEle25_CaloIdL_MW/O");
  T1->Branch("hlt_Diphoton_30_18",&hlt_Diphoton_30_18,"hlt_Diphoton_30_18/O");
  T1->Branch("hlt_PFJet500",&hlt_PFJet500,"hlt_PFJet500/O");
  T1->Branch("hlt_HT1050",&hlt_HT1050,"hlt_HT1050/O");
  T1->Branch("hlt_AK8PFJet400_TrimMass30",&hlt_AK8PFJet400_TrimMass30,"hlt_AK8PFJet400_TrimMass30/O");
  T1->Branch("hlt_AK8PFHT800_TrimMass50",&hlt_AK8PFHT800_TrimMass50,"hlt_AK8PFHT800_TrimMass50/O");
  T1->Branch("hlt_Photon200",&hlt_Photon200,"hlt_Photon200/O");
  T1->Branch("hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,"hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
  T1->Branch("hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",&hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60,"hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60/O");
  T1->Branch("hlt_PFMETNoMu140_PFMHTNoMu140_IDTight",&hlt_PFMETNoMu140_PFMHTNoMu140_IDTight,"hlt_PFMETNoMu140_PFMHTNoMu140_IDTight/O");
  T1->Branch("hlt_PFMETTypeOne140_PFMHT140_IDTight",&hlt_PFMETTypeOne140_PFMHT140_IDTight,"hlt_PFMETTypeOne140_PFMHT140_IDTight/O");

//--------------------------------------------------2017 triggers---------------------------------------------//
/*  T1->Branch("hlt_IsoMu24",&hlt_IsoMu24,"hlt_IsoMu24/O");
  T1->Branch("hlt_IsoMu27",&hlt_IsoMu27,"hlt_IsoMu27/O");
  T1->Branch("hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",&hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165,"hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165/O");
  T1->Branch("hlt_Ele115_CaloIdVT_GsfTrkIdT",&hlt_Ele115_CaloIdVT_GsfTrkIdT,"hlt_Ele115_CaloIdVT_GsfTrkIdT/O");
  T1->Branch("hlt_Ele40_WPTight_Gsf",&hlt_Ele40_WPTight_Gsf,"hlt_Ele40_WPTight_Gsf/O");
  T1->Branch("hlt_Ele32_WPTight_Gsf",&hlt_Ele32_WPTight_Gsf,"hlt_Ele32_WPTight_Gsf/O");
  T1->Branch("hlt_Ele28_eta2p1_WPTight_Gsf_HT150",&hlt_Ele28_eta2p1_WPTight_Gsf_HT150,"hlt_Ele28_eta2p1_WPTight_Gsf_HT150/O");
  T1->Branch("hlt_Mu17_Mu8",&hlt_Mu17_Mu8,"hlt_Mu17_Mu8/O");
  T1->Branch("hlt_Ele23_Ele12",&hlt_Ele23_Ele12,"hlt_Ele23_Ele12/O");
  T1->Branch("hlt_DoubleEle33",&hlt_DoubleEle33,"hlt_DoubleEle33/O");
  T1->Branch("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW,"hlt_DoubleEle25_CaloIdL_MW/O");
  T1->Branch("hlt_Diphoton_30_18",&hlt_Diphoton_30_18,"hlt_Diphoton_30_18/O");
  T1->Branch("hlt_PFJet500",&hlt_PFJet500,"hlt_PFJet500/O");
  T1->Branch("hlt_HT1050",&hlt_HT1050,"hlt_HT1050/O");
  T1->Branch("hlt_AK8PFJet400_TrimMass30",&hlt_AK8PFJet400_TrimMass30,"hlt_AK8PFJet400_TrimMass30/O");
  T1->Branch("hlt_AK8PFHT800_TrimMass50",&hlt_AK8PFHT800_TrimMass50,"hlt_AK8PFHT800_TrimMass50/O");
  T1->Branch("hlt_Photon200",&hlt_Photon200,"hlt_Photon200/O");
  T1->Branch("hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",&hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,"hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
  T1->Branch("hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60",&hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60,"hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60/O");
  T1->Branch("hlt_PFMETNoMu140_PFMHTNoMu140_IDTight",&hlt_PFMETNoMu140_PFMHTNoMu140_IDTight,"hlt_PFMETNoMu140_PFMHTNoMu140_IDTight/O");
  T1->Branch("hlt_PFMETTypeOne140_PFMHT140_IDTight",&hlt_PFMETTypeOne140_PFMHT140_IDTight,"hlt_PFMETTypeOne140_PFMHT140_IDTight/O");
*/
//-----------------------------------------2016 triggers--------------------------------------------//
/*  T1->Branch("hlt_IsoMu24",&hlt_IsoMu24,"hlt_IsoMu24/O");
  T1->Branch("hlt_IsoTkMu24",&hlt_IsoTkMu24,"hlt_IsoTkMu24/O");
  T1->Branch("hlt_Ele27_WPTight_Gsf",&hlt_Ele27_WPTight_Gsf,"hlt_Ele27_WPTight_Gsf/O");
  T1->Branch("hlt_Ele115_CaloIdVT_GsfTrkIdT",&hlt_Ele115_CaloIdVT_GsfTrkIdT,"hlt_Ele115_CaloIdVT_GsfTrkIdT/O");
  T1->Branch("hlt_Ele40_WPTight_Gsf",&hlt_Ele40_WPTight_Gsf,"hlt_Ele40_WPTight_Gsf/O");
  T1->Branch("hlt_Ele32_WPTight_Gsf",&hlt_Ele32_WPTight_Gsf,"hlt_Ele32_WPTight_Gsf/O");
  T1->Branch("hlt_Ele28_eta2p1_WPTight_Gsf_HT150",&hlt_Ele28_eta2p1_WPTight_Gsf_HT150,"hlt_Ele28_eta2p1_WPTight_Gsf_HT150/O");
  T1->Branch("hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",&hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,"hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8/O");
  T1->Branch("hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,"hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL/O");
  T1->Branch("hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,"hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
  T1->Branch("hlt_DoubleEle25_CaloIdL_MW",&hlt_DoubleEle25_CaloIdL_MW,"hlt_DoubleEle25_CaloIdL_MW/O");
  T1->Branch("hlt_Diphoton_30_18",&hlt_Diphoton_30_18,"hlt_Diphoton_30_18/O");
  T1->Branch("hlt_PFJet500",&hlt_PFJet500,"hlt_PFJet500/O");
  T1->Branch("hlt_HT1050",&hlt_HT1050,"hlt_HT1050/O");
  T1->Branch("hlt_AK8PFJet400_TrimMass30",&hlt_AK8PFJet400_TrimMass30,"hlt_AK8PFJet400_TrimMass30/O");
  T1->Branch("hlt_AK8PFHT800_TrimMass50",&hlt_AK8PFHT800_TrimMass50,"hlt_AK8PFHT800_TrimMass50/O");
  T1->Branch("hlt_Photon200",&hlt_Photon200,"hlt_Photon200/O");
  T1->Branch("hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,"hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/O");
  T1->Branch("hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,"hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ/O");
  T1->Branch("hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",&hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL,"hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL/O");
  T1->Branch("hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,"hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ/O");
*/

  T1->Branch("nTrigObj",&nTrigObj,"nTrigObj/I");
  T1->Branch("TrigObj_pt",TrigObj_pt,"TrigObj_pt[nTrigObj]/F");
  T1->Branch("TrigObj_eta",TrigObj_eta,"TrigObj_eta[nTrigObj]/F");
  T1->Branch("TrigObj_phi",TrigObj_phi,"TrigObj_phi[nTrigObj]/F");
  T1->Branch("TrigObj_mass",TrigObj_mass,"TrigObj_mass[nTrigObj]/F");
  T1->Branch("TrigObj_HLT",TrigObj_HLT,"TrigObj_HLT[nTrigObj]/O");
  T1->Branch("TrigObj_L1",TrigObj_L1,"TrigObj_L1[nTrigObj]/O");
  //T1->Branch("TrigObj_Both",TrigObj_Both,"TrigObj_Both[nTrigObj]/O");
  T1->Branch("TrigObj_Ihlt",TrigObj_Ihlt,"TrigObj_Ihlt[nTrigObj]/I");
  T1->Branch("TrigObj_pdgId",TrigObj_pdgId,"TrigObj_pdgId[nTrigObj]/I");
  T1->Branch("TrigObj_type",TrigObj_type,"TrigObj_type[nTrigObj]/I");
  
  // Prefire weights //
  
  T1->Branch("prefiringweight",&prefiringweight,"prefiringweight/D");
  T1->Branch("prefiringweightup",&prefiringweightup,"prefiringweightup/D");
  T1->Branch("prefiringweightdown",&prefiringweightdown,"prefiringweightdown/D");
  
  // MET info //
  
  T1->Branch("CHSMET_pt",&miset,"miset/F") ;
  T1->Branch("CHSMET_phi",&misphi,"misphi/F") ;
  T1->Branch("CHSMET_sig",&misetsig,"misetsig/F");
  T1->Branch("CHSMET_sumEt",&sumEt,"sumEt/F");
  
  T1->Branch("CHSMET_pt_UnclusEup",&miset_UnclusEup,"miset_CHS_UnclusEup/F") ;
  T1->Branch("CHSMET_pt_UnclusEdn",&miset_UnclusEdn,"miset_CHS_UnclusEdn/F") ;
  T1->Branch("CHSMET_phi_UnclusEup",&misphi_UnclusEup,"CHSMET_phi_UnclusEup/F") ;
  T1->Branch("CHSMET_phi_UnclusEdn",&misphi_UnclusEdn,"CHSMET_phi_UnclusEdn/F") ;
  
  
  T1->Branch("PuppiMET_pt",&miset_PUPPI,"miset_PUPPI/F") ;
  T1->Branch("PuppiMET_phi",&misphi_PUPPI,"misphi_PUPPI/F") ;
  T1->Branch("PuppiMET_sig",&misetsig_PUPPI,"misetsig_PUPPI/F");
  T1->Branch("PuppiMET_sumEt",&sumEt_PUPPI,"sumEt_PUPPI/F");
  
  T1->Branch("PuppiMET_pt_JESup",&miset_PUPPI_JESup,"miset_PUPPI_JESup/F") ;
  T1->Branch("PuppiMET_pt_JESdn",&miset_PUPPI_JESdn,"miset_PUPPI_JESdn/F") ;
  T1->Branch("PuppiMET_pt_JERup",&miset_PUPPI_JERup,"miset_PUPPI_JERup/F") ;
  T1->Branch("PuppiMET_pt_JERdn",&miset_PUPPI_JERdn,"miset_PUPPI_JERdn/F") ;
  T1->Branch("PuppiMET_pt_UnclusEup",&miset_PUPPI_UnclusEup,"miset_PUPPI_UnclusEup/F") ;
  T1->Branch("PuppiMET_pt_UnclusEdn",&miset_PUPPI_UnclusEdn,"miset_PUPPI_UnclusEdn/F") ;
  T1->Branch("PuppiMET_phi_JESup",&misphi_PUPPI_JESup,"misphi_PUPPI_JESup/F") ;
  T1->Branch("PuppiMET_phi_JESdn",&misphi_PUPPI_JESdn,"misphi_PUPPI_JESdn/F") ;
  T1->Branch("PuppiMET_phi_JERup",&misphi_PUPPI_JERup,"misphi_PUPPI_JERup/F") ;
  T1->Branch("PuppiMET_phi_JERdn",&misphi_PUPPI_JERdn,"misphi_PUPPI_JERdn/F") ;
  T1->Branch("PuppiMET_phi_UnclusEup",&misphi_PUPPI_UnclusEup,"misphi_PUPPI_UnclusEup/F") ;
  T1->Branch("PuppiMET_phi_UnclusEdn",&misphi_PUPPI_UnclusEdn,"misphi_PUPPI_UnclusEdn/F") ;
  
  // AK8 jet info //
  
  T1->Branch("nPFJetAK8",&nPFJetAK8, "nPFJetAK8/I"); 
  T1->Branch("PFJetAK8_pt",PFJetAK8_pt,"PFJetAK8_pt[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_y",PFJetAK8_y,"PFJetAK8_y[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_eta",PFJetAK8_eta,"PFJetAK8_eta[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_phi",PFJetAK8_phi,"PFJetAK8_phi[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_mass",PFJetAK8_mass,"PFJetAK8_mass[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jetID",PFJetAK8_jetID,"PFJetAK8_jetID[nPFJetAK8]/O");
  T1->Branch("PFJetAK8_jetID_tightlepveto",PFJetAK8_jetID_tightlepveto,"PFJetAK8_jetID_tightlepveto[nPFJetAK8]/O");
  T1->Branch("PFJetAK8_JEC",PFJetAK8_JEC,"PFJetAK8_JEC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_CHF",PFJetAK8_CHF,"PFJetAK8_CHF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_NHF",PFJetAK8_NHF,"PFJetAK8_NHF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_CEMF",PFJetAK8_CEMF,"PFJetAK8_CEMF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_NEMF",PFJetAK8_NEMF,"PFJetAK8_NEMF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_MUF",PFJetAK8_MUF,"PFJetAK8_MUF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_PHF",PFJetAK8_PHF,"PFJetAK8_PHF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_EEF",PFJetAK8_EEF,"PFJetAK8_EEF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_HFHF",PFJetAK8_HFHF,"PFJetAK8_HFHF[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_CHM",PFJetAK8_CHM,"PFJetAK8_CHM[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_NHM",PFJetAK8_NHM,"PFJetAK8_NHM[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_MUM",PFJetAK8_MUM,"PFJetAK8_MUM[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_PHM",PFJetAK8_PHM,"PFJetAK8_PHM[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_EEM",PFJetAK8_EEM,"PFJetAK8_EEM[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_HFHM",PFJetAK8_HFHM,"PFJetAK8_HFHM[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_Neucons",PFJetAK8_Neucons,"PFJetAK8_Neucons[nPFJetAK8]/I");
  T1->Branch("PFJetAK8_Chcons",PFJetAK8_Chcons,"PFJetAK8_Chcons[nPFJetAK8]/I");
  
  T1->Branch("PFJetAK8_JER",PFJetAK8_reso,"PFJetAK8_reso[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_JERup",PFJetAK8_resoup,"PFJetAK8_resoup[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_JERdn",PFJetAK8_resodn,"PFJetAK8_resodn[nPFJetAK8]/F");
  
  T1->Branch("PFJetAK8_msoftdrop",PFJetAK8_sdmass,"PFJetAK8_sdmass[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_tau1",PFJetAK8_tau1,"PFJetAK8_tau1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_tau2",PFJetAK8_tau2,"PFJetAK8_tau2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_tau3",PFJetAK8_tau3,"PFJetAK8_tau3[nPFJetAK8]/F");
  
  T1->Branch("PFJetAK8_btag_DeepCSV",PFJetAK8_btag_DeepCSV,"PFJetAK8_btag_DeepCSV[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_DAK8MD_TvsQCD",PFJetAK8_DeepTag_DAK8_TvsQCD,"PFJetAK8_DeepTag_DAK8_TvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_DAK8MD_WvsQCD",PFJetAK8_DeepTag_DAK8_WvsQCD,"PFJetAK8_DeepTag_DAK8_WvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_DAK8MD_ZvsQCD",PFJetAK8_DeepTag_DAK8_ZvsQCD,"PFJetAK8_DeepTag_DAK8_ZvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_DAK8MD_HvsQCD",PFJetAK8_DeepTag_DAK8_HvsQCD,"PFJetAK8_DeepTag_DAK8_HvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_DAK8MD_bbvsQCD",PFJetAK8_DeepTag_DAK8_bbvsQCD,"PFJetAK8_DeepTag_DAK8_bbvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNet_TvsQCD",PFJetAK8_DeepTag_PNet_TvsQCD,"PFJetAK8_DeepTag_PNet_TvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNet_WvsQCD",PFJetAK8_DeepTag_PNet_WvsQCD,"PFJetAK8_DeepTag_PNet_WvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNet_ZvsQCD",PFJetAK8_DeepTag_PNet_ZvsQCD,"PFJetAK8_DeepTag_PNet_ZvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNetMD_XbbvsQCD",PFJetAK8_DeepTag_PNet_XbbvsQCD,"PFJetAK8_DeepTag_PNet_XbbvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNetMD_XccvsQCD",PFJetAK8_DeepTag_PNet_XccvsQCD,"PFJetAK8_DeepTag_PNet_XccvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNetMD_XqqvsQCD",PFJetAK8_DeepTag_PNet_XqqvsQCD,"PFJetAK8_DeepTag_PNet_XqqvsQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_DeepTag_PNetMD_QCD",PFJetAK8_DeepTag_PNet_QCD,"PFJetAK8_DeepTag_PNet_QCD[nPFJetAK8]/F");
  
  T1->Branch("PFJetAK8_sub1pt",PFJetAK8_sub1pt,"PFJetAK8_sub1pt[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub1eta",PFJetAK8_sub1eta,"PFJetAK8_sub1eta[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub1phi",PFJetAK8_sub1phi,"PFJetAK8_sub1phi[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub1mass",PFJetAK8_sub1mass,"PFJetAK8_sub1mass[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub1JEC",PFJetAK8_sub1JEC,"PFJetAK8_sub1JEC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub1btag",PFJetAK8_sub1btag,"PFJetAK8_sub1btag[nPFJetAK8]/F");
  
  T1->Branch("PFJetAK8_sub2pt",PFJetAK8_sub2pt,"PFJetAK8_sub2pt[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub2eta",PFJetAK8_sub2eta,"PFJetAK8_sub2eta[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub2phi",PFJetAK8_sub2phi,"PFJetAK8_sub2phi[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub2mass",PFJetAK8_sub2mass,"PFJetAK8_sub2mass[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub2JEC",PFJetAK8_sub2JEC,"PFJetAK8_sub2JEC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_sub2btag",PFJetAK8_sub2btag,"PFJetAK8_sub2btag[nPFJetAK8]/F");
  
  T1->Branch("PFJetAK8_jesup_AbsoluteStat",PFJetAK8_jesup_AbsoluteStat,"PFJetAK8_jesup_AbsoluteStat[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_AbsoluteScale",PFJetAK8_jesup_AbsoluteScale,"PFJetAK8_jesup_AbsoluteScale[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_AbsoluteMPFBias",PFJetAK8_jesup_AbsoluteMPFBias,"PFJetAK8_jesup_AbsoluteMPFBias[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_FlavorQCD",PFJetAK8_jesup_FlavorQCD,"PFJetAK8_jesup_FlavorQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_Fragmentation",PFJetAK8_jesup_Fragmentation,"PFJetAK8_jesup_Fragmentation[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_PileUpDataMC",PFJetAK8_jesup_PileUpDataMC,"PFJetAK8_jesup_PileUpDataMC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_PileUpPtBB",PFJetAK8_jesup_PileUpPtBB,"PFJetAK8_jesup_PileUpPtBB[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_PileUpPtEC1",PFJetAK8_jesup_PileUpPtEC1,"PFJetAK8_jesup_PileUpPtEC1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_PileUpPtEC2",PFJetAK8_jesup_PileUpPtEC2,"PFJetAK8_jesup_PileUpPtEC2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_PileUpPtRef",PFJetAK8_jesup_PileUpPtRef,"PFJetAK8_jesup_PileUpPtRef[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeFSR",PFJetAK8_jesup_RelativeFSR,"PFJetAK8_jesup_RelativeFSR[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeJEREC1",PFJetAK8_jesup_RelativeJEREC1,"PFJetAK8_jesup_RelativeJEREC1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeJEREC2",PFJetAK8_jesup_RelativeJEREC2,"PFJetAK8_jesup_RelativeJEREC2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativePtBB",PFJetAK8_jesup_RelativePtBB,"PFJetAK8_jesup_RelativePtBB[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativePtEC1",PFJetAK8_jesup_RelativePtEC1,"PFJetAK8_jesup_RelativePtEC1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativePtEC2",PFJetAK8_jesup_RelativePtEC2,"PFJetAK8_jesup_RelativePtEC2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeBal",PFJetAK8_jesup_RelativeBal,"PFJetAK8_jesup_RelativeBal[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeSample",PFJetAK8_jesup_RelativeSample,"PFJetAK8_jesup_RelativeSample[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeStatEC",PFJetAK8_jesup_RelativeStatEC,"PFJetAK8_jesup_RelativeStatEC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_RelativeStatFSR",PFJetAK8_jesup_RelativeStatFSR,"PFJetAK8_jesup_RelativeStatFSR[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_SinglePionECAL",PFJetAK8_jesup_SinglePionECAL,"PFJetAK8_jesup_SinglePionECAL[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_SinglePionHCAL",PFJetAK8_jesup_SinglePionHCAL,"PFJetAK8_jesup_SinglePionHCAL[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_TimePtEta",PFJetAK8_jesup_TimePtEta,"PFJetAK8_jesup_TimePtEta[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesup_Total",PFJetAK8_jesup_Total,"PFJetAK8_jesup_Total[nPFJetAK8]/F");
  
  T1->Branch("PFJetAK8_jesdn_AbsoluteStat",PFJetAK8_jesdn_AbsoluteStat,"PFJetAK8_jesdn_AbsoluteStat[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_AbsoluteScale",PFJetAK8_jesdn_AbsoluteScale,"PFJetAK8_jesdn_AbsoluteScale[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_AbsoluteMPFBias",PFJetAK8_jesdn_AbsoluteMPFBias,"PFJetAK8_jesdn_AbsoluteMPFBias[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_FlavorQCD",PFJetAK8_jesdn_FlavorQCD,"PFJetAK8_jesdn_FlavorQCD[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_Fragmentation",PFJetAK8_jesdn_Fragmentation,"PFJetAK8_jesdn_Fragmentation[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_PileUpDataMC",PFJetAK8_jesdn_PileUpDataMC,"PFJetAK8_jesdn_PileUpDataMC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_PileUpPtBB",PFJetAK8_jesdn_PileUpPtBB,"PFJetAK8_jesdn_PileUpPtBB[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_PileUpPtEC1",PFJetAK8_jesdn_PileUpPtEC1,"PFJetAK8_jesdn_PileUpPtEC1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_PileUpPtEC2",PFJetAK8_jesdn_PileUpPtEC2,"PFJetAK8_jesdn_PileUpPtEC2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_PileUpPtRef",PFJetAK8_jesdn_PileUpPtRef,"PFJetAK8_jesdn_PileUpPtRef[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeFSR",PFJetAK8_jesdn_RelativeFSR,"PFJetAK8_jesdn_RelativeFSR[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeJEREC1",PFJetAK8_jesdn_RelativeJEREC1,"PFJetAK8_jesdn_RelativeJEREC1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeJEREC2",PFJetAK8_jesdn_RelativeJEREC2,"PFJetAK8_jesdn_RelativeJEREC2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativePtBB",PFJetAK8_jesdn_RelativePtBB,"PFJetAK8_jesdn_RelativePtBB[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativePtEC1",PFJetAK8_jesdn_RelativePtEC1,"PFJetAK8_jesdn_RelativePtEC1[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativePtEC2",PFJetAK8_jesdn_RelativePtEC2,"PFJetAK8_jesdn_RelativePtEC2[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeBal",PFJetAK8_jesdn_RelativeBal,"PFJetAK8_jesdn_RelativeBal[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeSample",PFJetAK8_jesdn_RelativeSample,"PFJetAK8_jesdn_RelativeSample[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeStatEC",PFJetAK8_jesdn_RelativeStatEC,"PFJetAK8_jesdn_RelativeStatEC[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_RelativeStatFSR",PFJetAK8_jesdn_RelativeStatFSR,"PFJetAK8_jesdn_RelativeStatFSR[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_SinglePionECAL",PFJetAK8_jesdn_SinglePionECAL,"PFJetAK8_jesdn_SinglePionECAL[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_SinglePionHCAL",PFJetAK8_jesdn_SinglePionHCAL,"PFJetAK8_jesdn_SinglePionHCAL[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_TimePtEta",PFJetAK8_jesdn_TimePtEta,"PFJetAK8_jesdn_TimePtEta[nPFJetAK8]/F");
  T1->Branch("PFJetAK8_jesdn_Total",PFJetAK8_jesdn_Total,"PFJetAK8_jesdn_Total[nPFJetAK8]/F");
  
  //gROOT->ProcessLine(".L CustomRootDict.cc+");
  
  if(store_fatjet_constituents){
    T1->Branch("nPFJetAK8_cons",&nPFJetAK8_cons,"nPFJetAK8_cons/I");
    T1->Branch("PFJetAK8_cons_pt",PFJetAK8_cons_pt, "PFJetAK8_cons_pt[nPFJetAK8_cons]/F");
    T1->Branch("PFJetAK8_cons_eta",PFJetAK8_cons_eta, "PFJetAK8_cons_eta[nPFJetAK8_cons]/F");
    T1->Branch("PFJetAK8_cons_phi",PFJetAK8_cons_phi, "PFJetAK8_cons_phi[nPFJetAK8_cons]/F");
    T1->Branch("PFJetAK8_cons_mass",PFJetAK8_cons_mass, "PFJetAK8_cons_mass[nPFJetAK8_cons]/F");
    T1->Branch("PFJetAK8_cons_pdgId",PFJetAK8_cons_pdgId, "PFJetAK8_cons_pdgId[nPFJetAK8_cons]/I");
    T1->Branch("PFJetAK8_cons_jetIndex",PFJetAK8_cons_jetIndex, "PFJetAK8_cons_jetIndex[nPFJetAK8_cons]/I");
  }
  

  // AK4 jet info //
 
  T1->Branch("nPFJetAK4",&nPFJetAK4,"nPFJetAK4/I"); 

  T1->Branch("PFJetAK4_pt",PFJetAK4_pt,"PFJetAK4_pt[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_eta",PFJetAK4_eta,"PFJetAK4_eta[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_y",PFJetAK4_y,"PFJetAK4_y[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_phi",PFJetAK4_phi,"PFJetAK4_phi[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_mass",PFJetAK4_mass,"PFJetAK4_mass[nPFJetAK4]/F");  
  T1->Branch("PFJetAK4_energy",PFJetAK4_energy,"PFJetAK4_energy[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_pt_reg",PFJetAK4_pt_reg,"PFJetAK4_pt_reg[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_eta_reg",PFJetAK4_eta_reg,"PFJetAK4_eta_reg[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_y_reg",PFJetAK4_y_reg,"PFJetAK4_y_reg[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_phi_reg",PFJetAK4_phi_reg,"PFJetAK4_phi_reg[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_mass_reg",PFJetAK4_mass_reg,"PFJetAK4_mass_reg[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_energy_reg",PFJetAK4_energy_reg,"PFJetAK4_energy_reg[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_bcorr",PFJetAK4_bcorr,"PFJetAK4_bcorr[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_breso",PFJetAK4_breso,"PFJetAK4_breso[nPFJetAK4]/F");

  T1->Branch("PFJetAK4_jetID",PFJetAK4_jetID,"PFJetAK4_jetID[nPFJetAK4]/O");
  T1->Branch("PFJetAK4_jetID_tightlepveto",PFJetAK4_jetID_tightlepveto,"PFJetAK4_jetID_tightlepveto[nPFJetAK4]/O");
  
  T1->Branch("PFJetAK4_JEC",PFJetAK4_JEC,"PFJetAK4_JEC[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepCSV",PFJetAK4_btag_DeepCSV,"PFJetAK4_btag_DeepCSV[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepFlav",PFJetAK4_btag_DeepFlav,"PFJetAK4_btag_DeepFlav[nPFJetAK4]/F");
 
  T1->Branch("PFJetAK4_JER",PFJetAK4_reso,"PFJetAK4_reso[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_JERup",PFJetAK4_resoup,"PFJetAK4_resoup[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_JERdn",PFJetAK4_resodn,"PFJetAK4_resodn[nPFJetAK4]/F"); 
  
  T1->Branch("PFJetAK4_hadronflav",PFJetAK4_hadronflav,"PFJetAK4_hadronflav[nPFJetAK4]/I");
  T1->Branch("PFJetAK4_partonflav",PFJetAK4_partonflav,"PFJetAK4_partonflav[nPFJetAK4]/I");
  T1->Branch("PFJetAK4_qgl",PFJetAK4_qgl,"PFJetAK4_qgl[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_PUID",PFJetAK4_PUID,"PFJetAK4_PUID[nPFJetAK4]/F");
  
  T1->Branch("PFJetAK4_jesup_AbsoluteStat",PFJetAK4_jesup_AbsoluteStat,"PFJetAK4_jesup_AbsoluteStat[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_AbsoluteScale",PFJetAK4_jesup_AbsoluteScale,"PFJetAK4_jesup_AbsoluteScale[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_AbsoluteMPFBias",PFJetAK4_jesup_AbsoluteMPFBias,"PFJetAK4_jesup_AbsoluteMPFBias[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_FlavorQCD",PFJetAK4_jesup_FlavorQCD,"PFJetAK4_jesup_FlavorQCD[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_Fragmentation",PFJetAK4_jesup_Fragmentation,"PFJetAK4_jesup_Fragmentation[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_PileUpDataMC",PFJetAK4_jesup_PileUpDataMC,"PFJetAK4_jesup_PileUpDataMC[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_PileUpPtBB",PFJetAK4_jesup_PileUpPtBB,"PFJetAK4_jesup_PileUpPtBB[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_PileUpPtEC1",PFJetAK4_jesup_PileUpPtEC1,"PFJetAK4_jesup_PileUpPtEC1[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_PileUpPtEC2",PFJetAK4_jesup_PileUpPtEC2,"PFJetAK4_jesup_PileUpPtEC2[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_PileUpPtRef",PFJetAK4_jesup_PileUpPtRef,"PFJetAK4_jesup_PileUpPtRef[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeFSR",PFJetAK4_jesup_RelativeFSR,"PFJetAK4_jesup_RelativeFSR[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeJEREC1",PFJetAK4_jesup_RelativeJEREC1,"PFJetAK4_jesup_RelativeJEREC1[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeJEREC2",PFJetAK4_jesup_RelativeJEREC2,"PFJetAK4_jesup_RelativeJEREC2[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativePtBB",PFJetAK4_jesup_RelativePtBB,"PFJetAK4_jesup_RelativePtBB[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativePtEC1",PFJetAK4_jesup_RelativePtEC1,"PFJetAK4_jesup_RelativePtEC1[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativePtEC2",PFJetAK4_jesup_RelativePtEC2,"PFJetAK4_jesup_RelativePtEC2[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeBal",PFJetAK4_jesup_RelativeBal,"PFJetAK4_jesup_RelativeBal[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeSample",PFJetAK4_jesup_RelativeSample,"PFJetAK4_jesup_RelativeSample[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeStatEC",PFJetAK4_jesup_RelativeStatEC,"PFJetAK4_jesup_RelativeStatEC[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_RelativeStatFSR",PFJetAK4_jesup_RelativeStatFSR,"PFJetAK4_jesup_RelativeStatFSR[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_SinglePionECAL",PFJetAK4_jesup_SinglePionECAL,"PFJetAK4_jesup_SinglePionECAL[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_SinglePionHCAL",PFJetAK4_jesup_SinglePionHCAL,"PFJetAK4_jesup_SinglePionHCAL[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_TimePtEta",PFJetAK4_jesup_TimePtEta,"PFJetAK4_jesup_TimePtEta[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesup_Total",PFJetAK4_jesup_Total,"PFJetAK4_jesup_Total[nPFJetAK4]/F");
  
  T1->Branch("PFJetAK4_jesdn_AbsoluteStat",PFJetAK4_jesdn_AbsoluteStat,"PFJetAK4_jesdn_AbsoluteStat[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_AbsoluteScale",PFJetAK4_jesdn_AbsoluteScale,"PFJetAK4_jesdn_AbsoluteScale[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_AbsoluteMPFBias",PFJetAK4_jesdn_AbsoluteMPFBias,"PFJetAK4_jesdn_AbsoluteMPFBias[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_FlavorQCD",PFJetAK4_jesdn_FlavorQCD,"PFJetAK4_jesdn_FlavorQCD[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_Fragmentation",PFJetAK4_jesdn_Fragmentation,"PFJetAK4_jesdn_Fragmentation[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_PileUpDataMC",PFJetAK4_jesdn_PileUpDataMC,"PFJetAK4_jesdn_PileUpDataMC[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_PileUpPtBB",PFJetAK4_jesdn_PileUpPtBB,"PFJetAK4_jesdn_PileUpPtBB[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_PileUpPtEC1",PFJetAK4_jesdn_PileUpPtEC1,"PFJetAK4_jesdn_PileUpPtEC1[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_PileUpPtEC2",PFJetAK4_jesdn_PileUpPtEC2,"PFJetAK4_jesdn_PileUpPtEC2[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_PileUpPtRef",PFJetAK4_jesdn_PileUpPtRef,"PFJetAK4_jesdn_PileUpPtRef[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeFSR",PFJetAK4_jesdn_RelativeFSR,"PFJetAK4_jesdn_RelativeFSR[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeJEREC1",PFJetAK4_jesdn_RelativeJEREC1,"PFJetAK4_jesdn_RelativeJEREC1[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeJEREC2",PFJetAK4_jesdn_RelativeJEREC2,"PFJetAK4_jesdn_RelativeJEREC2[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativePtBB",PFJetAK4_jesdn_RelativePtBB,"PFJetAK4_jesdn_RelativePtBB[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativePtEC1",PFJetAK4_jesdn_RelativePtEC1,"PFJetAK4_jesdn_RelativePtEC1[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativePtEC2",PFJetAK4_jesdn_RelativePtEC2,"PFJetAK4_jesdn_RelativePtEC2[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeBal",PFJetAK4_jesdn_RelativeBal,"PFJetAK4_jesdn_RelativeBal[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeSample",PFJetAK4_jesdn_RelativeSample,"PFJetAK4_jesdn_RelativeSample[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeStatEC",PFJetAK4_jesdn_RelativeStatEC,"PFJetAK4_jesdn_RelativeStatEC[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_RelativeStatFSR",PFJetAK4_jesdn_RelativeStatFSR,"PFJetAK4_jesdn_RelativeStatFSR[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_SinglePionECAL",PFJetAK4_jesdn_SinglePionECAL,"PFJetAK4_jesdn_SinglePionECAL[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_SinglePionHCAL",PFJetAK4_jesdn_SinglePionHCAL,"PFJetAK4_jesdn_SinglePionHCAL[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_TimePtEta",PFJetAK4_jesdn_TimePtEta,"PFJetAK4_jesdn_TimePtEta[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_jesdn_Total",PFJetAK4_jesdn_Total,"PFJetAK4_jesdn_Total[nPFJetAK4]/F");
  
  T1->Branch("PFJetAK4_btag_DeepCSV_SF",PFJetAK4_btag_DeepCSV_SF,"PFJetAK4_btag_DeepCSV_SF[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepCSV_SF_up",PFJetAK4_btag_DeepCSV_SF_up,"PFJetAK4_btag_DeepCSV_SF_up[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepCSV_SF_dn",PFJetAK4_btag_DeepCSV_SF_dn,"PFJetAK4_btag_DeepCSV_SF_dn[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepFlav_SF",PFJetAK4_btag_DeepFlav_SF,"PFJetAK4_btag_DeepFlav_SF[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepFlav_SF_up",PFJetAK4_btag_DeepFlav_SF_up,"PFJetAK4_btag_DeepFlav_SF_up[nPFJetAK4]/F");
  T1->Branch("PFJetAK4_btag_DeepFlav_SF_dn",PFJetAK4_btag_DeepFlav_SF_dn,"PFJetAK4_btag_DeepFlav_SF_dn[nPFJetAK4]/F");
  
  // Muon info //
  
  T1->Branch("nMuon",&nMuon,"nMuon/I");
  T1->Branch("Muon_isPF",Muon_isPF,"Muon_isPF[nMuon]/O");
  T1->Branch("Muon_isGL",Muon_isGL,"Muon_isGL[nMuon]/O");
  T1->Branch("Muon_isTRK",Muon_isTRK,"Muon_isTRK[nMuon]/O");
  T1->Branch("Muon_isLoose",Muon_isLoose,"Muon_isLoose[nMuon]/O");
  T1->Branch("Muon_isGoodGL",Muon_isGoodGL,"Muon_isGoodGL[nMuon]/O");
  T1->Branch("Muon_isMed",Muon_isMed,"Muon_isMed[nMuon]/O");
  T1->Branch("Muon_isMedPr",Muon_isMedPr,"Muon_isMedPr[nMuon]/O");
  T1->Branch("Muon_isTight",Muon_isTight,"Muon_isTight[nMuon]/O");
  T1->Branch("Muon_isHighPt",Muon_isHighPt,"Muon_isHighPt[nMuon]/O"); 
  T1->Branch("Muon_isHighPttrk",Muon_isHighPttrk,"Muon_isHighPttrk[nMuon]/O");
  T1->Branch("Muon_TightID",Muon_TightID,"Muon_TightID[nMuon]/O");
  
  T1->Branch("Muon_pt",Muon_pt,"Muon_pt[nMuon]/F");
  T1->Branch("Muon_p",Muon_p,"Muon_p[nMuon]/F");
  T1->Branch("Muon_eta",Muon_eta,"Muon_eta[nMuon]/F");
  T1->Branch("Muon_phi",Muon_phi,"Muon_phi[nMuon]/F");
  T1->Branch("Muon_e",Muon_e,"Muon_e[nMuon]/F");
  
  T1->Branch("Muon_minisoch", Muon_minchiso, "Muon_minchiso[nMuon]/F");
  T1->Branch("Muon_minisonh", Muon_minnhiso, "Muon_minnhiso[nMuon]/F");
  T1->Branch("Muon_minisoph", Muon_minphiso, "Muon_minphiso[nMuon]/F");
  T1->Branch("Muon_minisoall", Muon_minisoall, "Muon_minisoall[nMuon]/F");
  T1->Branch("Muon_dxy",Muon_dxy,"Muon_dxy[nMuon]/F");
  T1->Branch("Muon_dz",Muon_dz,"Muon_dz[nMuon]/F");
  T1->Branch("Muon_dxyErr",Muon_dxyErr,"Muon_dxyErr[nMuon]/F");
  T1->Branch("Muon_ip3d",Muon_ip3d,"Muon_ip3d[nMuon]/F");
  T1->Branch("Muon_ptErr",Muon_ptErr,"Muon_ptErr[nMuon]/F");
  T1->Branch("Muon_chi",Muon_chi,"Muon_chi[nMuon]/F");
  T1->Branch("Muon_ndf",Muon_ndf,"Muon_ndf[nMuon]/I");
  T1->Branch("Muon_ecal",Muon_ecal,"Muon_ecal[nMuon]/F");
  T1->Branch("Muon_hcal",Muon_hcal,"Muon_hcal[nMuon]/F");
  T1->Branch("Muon_pfiso",Muon_pfiso,"Muon_pfiso[nMuon]/F");
  T1->Branch("Muon_posmatch",Muon_posmatch,"Muon_posmatch[nMuon]/F");
  T1->Branch("Muon_trkink",Muon_trkink,"Muon_trkink[nMuon]/F");
  T1->Branch("Muon_segcom",Muon_segcom,"Muon_segcom[nMuon]/F");
  T1->Branch("Muon_hit",Muon_hit,"Muon_hit[nMuon]/F");
  T1->Branch("Muon_pixhit",Muon_pixhit,"Muon_pixhit[nMuon]/F");
  T1->Branch("Muon_mst",Muon_mst,"Muon_mst[nMuon]/F");
  T1->Branch("Muon_trklay",Muon_trklay,"Muon_trklay[nMuon]/F"); 
  T1->Branch("Muon_valfrac",Muon_valfrac,"Muon_valfrac[nMuon]/F"); 
  T1->Branch("Muon_dxy_sv",Muon_dxy_sv,"Muon_dxy_sv[nMuon]/F");
  
  T1->Branch("Muon_corrected_pt",Muon_corrected_pt,"Muon_corrected_pt[nMuon]/F");
  T1->Branch("Muon_correctedUp_pt",Muon_correctedUp_pt,"Muon_correctedUp_pt[nMuon]/F");
  T1->Branch("Muon_correctedDown_pt",Muon_correctedDown_pt,"Muon_correctedDown_pt[nMuon]/F");
  
  // Electron info //
  
  T1->Branch("nElectron",&nElectron,"nElectron/I");
  T1->Branch("Electron_pt",Electron_pt,"Electron_pt[nElectron]/F");
  T1->Branch("Electron_eta",Electron_eta,"Electron_eta[nElectron]/F");
  T1->Branch("Electron_phi",Electron_phi,"Electron_phi[nElectron]/F");
  T1->Branch("Electron_p",Electron_p,"Electron_p[nElectron]/F");
  T1->Branch("Electron_e",Electron_e,"Electron_e[nElectron]/F");
  T1->Branch("Electron_e_ECAL",Electron_e_ECAL,"Electron_e_ECAL[nElectron]/F");
  
  T1->Branch("Electron_mvaid_Fallv2WP90",Electron_mvaid_Fallv2WP90,"Electron_mvaid_Fallv2WP90[nElectron]/O");
  T1->Branch("Electron_mvaid_Fallv2WP90_noIso",Electron_mvaid_Fallv2WP90_noIso,"Electron_mvaid_Fallv2WP90_noIso[nElectron]/O");
  T1->Branch("Electron_mvaid_Fallv2WP80",Electron_mvaid_Fallv2WP80,"Electron_mvaid_Fallv2WP80[nElectron]/O");
  T1->Branch("Electron_mvaid_Fallv2WP80_noIso",Electron_mvaid_Fallv2WP80_noIso,"Electron_mvaid_Fallv2WP80_noIso[nElectron]/O");
  
  T1->Branch("Electron_dxy",Electron_dxy,"Electron_dxy[nElectron]/F");
  T1->Branch("Electron_dxyErr",Electron_dxyErr,"Electron_dxyErr[nElectron]/F");
  T1->Branch("Electron_dz",Electron_dz,"Electron_dz[nElectron]/F");
  T1->Branch("Electron_dzErr",Electron_dzErr,"Electron_dzErr[nElectron]/F");
  T1->Branch("Electron_ip3d",Electron_ip3d,"Electron_ip3d[nElectron]/F");
  T1->Branch("Electron_dxy_sv",Electron_dxy_sv,"Electron_dxy_sv[nElectron]/F");
  
  T1->Branch("Electron_hovere",Electron_hovere,"Electron_hovere[nElectron]/F");
  T1->Branch("Electron_chi",Electron_chi,"Electron_chi[nElectron]/F");
  T1->Branch("Electron_ndf",Electron_ndf,"Electron_ndf[nElectron]/I");
  T1->Branch("Electron_eoverp",Electron_eoverp,"Electron_eoverp[nElectron]/F");
  T1->Branch("Electron_ietaieta",Electron_ietaieta,"Electron_ietaieta[nElectron]/F");
  T1->Branch("Electron_misshits",Electron_misshits,"Electron_misshits[nElectron]/F");
  T1->Branch("Electron_pfiso_drcor",Electron_pfiso_drcor,"Electron_pfiso_drcor[nElectron]/F");
  T1->Branch("Electron_pfiso_eacor",Electron_pfiso_eacor,"Electron_pfiso_eacor[nElectron]/F");
  T1->Branch("Electron_pfiso04_eacor",Electron_pfiso04_eacor,"Electron_pfiso04_eacor[nElectron]/F");
  
  if(store_electron_scalnsmear){
	T1->Branch("Electron_eccalTrkEnergyPostCorr",Electron_eccalTrkEnergyPostCorr,"Electron_eccalTrkEnergyPostCorr[nElectron]/F");
	T1->Branch("Electron_energyScaleValue",Electron_energyScaleValue,"Electron_energyScaleValue[nElectron]/F");
	T1->Branch("Electron_energyScaleUp",Electron_energyScaleUp,"Electron_energyScaleUp[nElectron]/F");
	T1->Branch("Electron_energyScaleDown",Electron_energyScaleDown,"Electron_energyScaleDown[nElectron]/F");
	T1->Branch("Electron_energySigmaValue",Electron_energySigmaValue,"Electron_energySigmaValue[nElectron]/F");
	T1->Branch("Electron_energySigmaUp",Electron_energySigmaUp,"Electron_energySigmaUp[nElectron]/F");
	T1->Branch("Electron_energySigmaDown",Electron_energySigmaDown,"Electron_energySigmaDown[nElectron]/F");
	T1->Branch("Electron_energyScaleStatUp",Electron_energyScaleStatUp,"Electron_energyScaleStatUp[nElectron]/F");
        T1->Branch("Electron_energyScaleStatDown",Electron_energyScaleStatDown,"Electron_energyScaleStatDown[nElectron]/F");
        T1->Branch("Electron_energyScaleSystUp",Electron_energyScaleSystUp,"Electron_energyScaleSystUp[nElectron]/F");
        T1->Branch("Electron_energyScaleSystDown",Electron_energyScaleSystDown,"Electron_energyScaleSystDown[nElectron]/F");
        T1->Branch("Electron_energyScaleGainUp",Electron_energyScaleGainUp,"Electron_energyScaleGainUp[nElectron]/F");
        T1->Branch("Electron_energyScaleGainDown",Electron_energyScaleGainDown,"Electron_energyScaleGainDown[nElectron]/F");
	T1->Branch("Electron_energySigmaRhoUp",Electron_energySigmaRhoUp,"Electron_energySigmaRhoUp[nElectron]/F");
        T1->Branch("Electron_energySigmaRhoDown",Electron_energySigmaRhoDown,"Electron_energySigmaRhoDown[nElectron]/F");
        T1->Branch("Electron_energySigmaPhiUp",Electron_energySigmaPhiUp,"Electron_energySigmaPhiUp[nElectron]/F");
        T1->Branch("Electron_energySigmaPhiDown",Electron_energySigmaPhiDown,"Electron_energySigmaPhiDown[nElectron]/F");
        T1->Branch("Electron_energySmearNrSigma",Electron_energySmearNrSigma,"Electron_energySmearNrSigma[nElectron]/F");
	T1->Branch("Electron_ecalEnergyPreCorr",Electron_ecalEnergyPreCorr,"Electron_ecalEnergyPreCorr[nElectron]/F");
  	T1->Branch("Electron_ecalEnergyErrPreCorr",Electron_ecalEnergyErrPreCorr,"Electron_ecalEnergyErrPreCorr[nElectron]/F");
  	T1->Branch("Electron_ecalEnergyPostCorr",Electron_ecalEnergyPostCorr,"Electron_ecalEnergyPostCorr[nElectron]/F");
  	T1->Branch("Electron_ecalEnergyErrPostCorr",Electron_ecalEnergyErrPostCorr,"Electron_ecalEnergyErrPostCorr[nElectron]/F");
  }
  
  T1->Branch("Electron_supcl_eta",Electron_supcl_eta,"Electron_supcl_eta[nElectron]/F");
  T1->Branch("Electron_supcl_phi",Electron_supcl_phi,"Electron_supcl_phi[nElectron]/F");
  T1->Branch("Electron_supcl_e",Electron_supcl_e,"Electron_supcl_e[nElectron]/F");
  T1->Branch("Electron_supcl_rawE",Electron_supcl_rawE,"Electron_supcl_rawE[nElectron]/F");
  T1->Branch("Electron_sigmaieta", Electron_sigmaieta, "Electron_sigmaieta[nElectron]/F");
  T1->Branch("Electron_sigmaiphi", Electron_sigmaiphi, "Electron_sigmaiphi[nElectron]/F");
 
  T1->Branch("Electron_r9full", Electron_r9full, "Electron_r9full[nElectron]/F");
  T1->Branch("Electron_hcaloverecal", Electron_hcaloverecal, "Electron_hcaloverecal[nElectron]/F");
  T1->Branch("Electron_hitsmiss", Electron_hitsmiss, "Electron_hitsmiss[nElectron]/F");
  T1->Branch("Electron_ecloverpout", Electron_ecloverpout, "Electron_ecloverpout[nElectron]/F"); 
  T1->Branch("Electron_convVeto", Electron_convVeto, "Electron_convVeto[nElectron]/O");

  T1->Branch("Electron_pfisolsumphet", Electron_pfisolsumphet, "Electron_pfisolsumphet[nElectron]/F");
  T1->Branch("Electron_pfisolsumchhadpt", Electron_pfisolsumchhadpt, "Electron_pfisolsumchhadpt[nElectron]/F");
  T1->Branch("Electron_pfsiolsumneuhadet", Electron_pfsiolsumneuhadet, "Electron_pfsiolsumneuhadet[nElectron]/F");
  
  T1->Branch("Electron_minisoch", Electron_minchiso, "Electron_minchiso[nElectron]/F");
  T1->Branch("Electron_minisonh", Electron_minnhiso, "Electron_minnhiso[nElectron]/F");
  T1->Branch("Electron_minisoph", Electron_minphiso, "Electron_minphiso[nElectron]/F");
  T1->Branch("Electron_minisoall", Electron_minisoall, "Electron_minisoall[nElectron]/F");
  
  if(store_electron_addvariabs){
	
	T1->Branch("Electron_etain",Electron_etain,"Electron_etain[nElectron]/F");
	T1->Branch("Electron_phiin",Electron_phiin,"Electron_phiin[nElectron]/F");
	T1->Branch("Electron_fbrem",Electron_fbrem,"Electron_fbrem[nElectron]/F");
	T1->Branch("Electron_supcl_etaw", Electron_supcl_etaw, "Electron_supcl_etaw[nElectron]/F");
	T1->Branch("Electron_supcl_phiw", Electron_supcl_phiw, "Electron_supcl_phiw[nElectron]/F");
	T1->Branch("Electron_cloctftrkn", Electron_cloctftrkn, "Electron_cloctftrkn[nElectron]/F");
	T1->Branch("Electron_cloctftrkchi2", Electron_cloctftrkchi2, "Electron_cloctftrkchi2[nElectron]/F");
	T1->Branch("Electron_e1x5bye5x5", Electron_e1x5bye5x5, "Electron_e1x5bye5x5[nElectron]/F");
	T1->Branch("Electron_normchi2", Electron_normchi2, "Electron_normchi2[nElectron]/F");
	T1->Branch("Electron_trkmeasure", Electron_trkmeasure, "Electron_trkmeasure[nElectron]/F");
	T1->Branch("Electron_convtxprob", Electron_convtxprob, "Electron_convtxprob[nElectron]/F");
	T1->Branch("Electron_deltaetacltrkcalo", Electron_deltaetacltrkcalo, "Electron_deltaetacltrkcalo[nElectron]/F");
	T1->Branch("Electron_supcl_preshvsrawe", Electron_supcl_preshvsrawe, "Electron_supcl_preshvsrawe[nElectron]/F");
	T1->Branch("Electron_ecaletrkmomentum", Electron_ecaletrkmomentum, "Electron_ecaletrkmomentum[nElectron]/F");
  
  }
  
  // Photon Info //
  
  T1->Branch("nPhoton",&nPhoton,"nPhoton/I");
  T1->Branch("Photon_pt",Photon_pt,"Photon_pt[nPhoton]/F");
  T1->Branch("Photon_e",Photon_e,"Photon_e[nPhoton]/F");
  T1->Branch("Photon_eta",Photon_eta,"Photon_eta[nPhoton]/F");
  T1->Branch("Photon_phi",Photon_phi,"Photon_phi[nPhoton]/F");
  T1->Branch("Photon_mvaid_Fall17V2_raw",Photon_mvaid_Fall17V2_raw,"Photon_mvaid_Fall17V2_raw[nPhoton]/F");
  T1->Branch("Photon_mvaid_Fall17V2_WP90",Photon_mvaid_Fall17V2_WP90,"Photon_mvaid_Fall17V2_WP90[nPhoton]/O");
  T1->Branch("Photon_mvaid_Fall17V2_WP80",Photon_mvaid_Fall17V2_WP80,"Photon_mvaid_Fall17V2_WP80[nPhoton]/O");
  T1->Branch("Photon_mvaid_Spring16V1_WP90",Photon_mvaid_Spring16V1_WP90,"Photon_mvaid_Spring16V1_WP90[nPhoton]/O");
  T1->Branch("Photon_mvaid_Spring16V1_WP80",Photon_mvaid_Spring16V1_WP80,"Photon_mvaid_Spring16V1_WP80[nPhoton]/O");
  T1->Branch("Photon_e1by9",Photon_e1by9,"Photon_e1by9[nPhoton]/F");
  T1->Branch("Photon_e9by25",Photon_e9by25,"Photon_e9by25[nPhoton]/F");
  T1->Branch("Photon_trkiso",Photon_trkiso,"Photon_trkiso[nPhoton]/F");
  T1->Branch("Photon_emiso",Photon_emiso,"Photon_emiso[nPhoton]/F");
  T1->Branch("Photon_hadiso",Photon_hadiso,"Photon_hadiso[nPhoton]/F");
  T1->Branch("Photon_chhadiso",Photon_chhadiso,"Photon_chhadiso[nPhoton]/F");
  T1->Branch("Photon_neuhadiso",Photon_neuhadiso,"Photon_neuhadiso[nPhoton]/F");
  T1->Branch("Photon_phoiso",Photon_phoiso,"Photon_phoiso[nPhoton]/F");
  T1->Branch("Photon_PUiso",Photon_PUiso,"Photon_PUiso[nPhoton]/F");
  T1->Branch("Photon_hadbyem",Photon_hadbyem,"Photon_hadbyem[nPhoton]/F");
  T1->Branch("Photon_ietaieta",Photon_ietaieta,"Photon_ietaieta[nPhoton]/F");
  T1->Branch("Photon_passEveto",Photon_passEveto,"Photon_passEveto[nPhoton]/O");
  T1->Branch("Photon_PixelSeed",Photon_PixelSeed,"Photon_PixelSeed[nPhoton]/O");

  T1->Branch("Photon_energyScaleValue",Photon_energyScaleValue,"Photon_energyScaleValue[nPhoton]/F");
  T1->Branch("Photon_energyScaleUp",Photon_energyScaleUp,"Photon_energyScaleUp[nPhoton]/F");
  T1->Branch("Photon_energyScaleDown",Photon_energyScaleDown,"Photon_energyScaleDown[nPhoton]/F");
  T1->Branch("Photon_energySigmaValue",Photon_energySigmaValue,"Photon_energySigmaValue[nPhoton]/F");
  T1->Branch("Photon_energySigmaUp",Photon_energySigmaUp,"Photon_energySigmaUp[nPhoton]/F");
  T1->Branch("Photon_energySigmaDown",Photon_energySigmaDown,"Photon_energySigmaDown[nPhoton]/F");
  T1->Branch("Photon_energyScaleStatUp",Photon_energyScaleStatUp,"Photon_energyScaleStatUp[nPhoton]/F");
  T1->Branch("Photon_energyScaleStatDown",Photon_energyScaleStatDown,"Photon_energyScaleStatDown[nPhoton]/F");
  T1->Branch("Photon_energyScaleSystUp",Photon_energyScaleSystUp,"Photon_energyScaleSystUp[nPhoton]/F");
  T1->Branch("Photon_energyScaleSystDown",Photon_energyScaleSystDown,"Photon_energyScaleSystDown[nPhoton]/F");
  T1->Branch("Photon_energyScaleGainUp",Photon_energyScaleGainUp,"Photon_energyScaleGainUp[nPhoton]/F");
  T1->Branch("Photon_energyScaleGainDown",Photon_energyScaleGainDown,"Photon_energyScaleGainDown[nPhoton]/F");
  T1->Branch("Photon_energySigmaRhoUp",Photon_energySigmaRhoUp,"Photon_energySigmaRhoUp[nPhoton]/F");
  T1->Branch("Photon_energySigmaRhoDown",Photon_energySigmaRhoDown,"Photon_energySigmaRhoDown[nPhoton]/F");
  T1->Branch("Photon_energySigmaPhiUp",Photon_energySigmaPhiUp,"Photon_energySigmaPhiUp[nPhoton]/F");
  T1->Branch("Photon_energySigmaPhiDown",Photon_energySigmaPhiDown,"Photon_energySigmaPhiDown[nPhoton]/F");
  T1->Branch("Photon_energySmearNrSigma",Photon_energySmearNrSigma,"Photon_energySmearNrSigma[nPhoton]/F");
  T1->Branch("Photon_ecalEnergyPreCorr",Photon_ecalEnergyPreCorr,"Photon_ecalEnergyPreCorr[nPhoton]/F");
  T1->Branch("Photon_ecalEnergyErrPreCorr",Photon_ecalEnergyErrPreCorr,"Photon_ecalEnergyErrPreCorr[nPhoton]/F");
  T1->Branch("Photon_ecalEnergyPostCorr",Photon_ecalEnergyPostCorr,"Photon_ecalEnergyPostCorr[nPhoton]/F");
  T1->Branch("Photon_ecalEnergyErrPostCorr",Photon_ecalEnergyErrPostCorr,"Photon_ecalEnergyErrPostCorr[nPhoton]/F");

  // Tau Info //
  
  T1->Branch("nTau",&nTau,"nTau/I");
  T1->Branch("Tau_isPF",Tau_isPF,"Tau_isPF[nTau]/O");
  T1->Branch("Tau_pt",Tau_pt,"Tau_pt[nTau]/F");
  T1->Branch("Tau_eta",Tau_eta,"Tau_eta[nTau]/F");
  T1->Branch("Tau_phi",Tau_phi,"Tau_phi[nTau]/F");
  T1->Branch("Tau_e",Tau_e,"Tau_e[nTau]/F");
  T1->Branch("Tau_charge",Tau_charge,"Tau_charge[nTau]/I");
  T1->Branch("Tau_dxy",Tau_dxy,"Tau_dxy[nTau]/F");
  //T1->Branch("Tau_dz",Tau_dz,"Tau_dz[nTau]/F");
  T1->Branch("Tau_leadtrkdxy",Tau_leadtrkdxy,"Tau_leadtrkdxy[nTau]/F");
  T1->Branch("Tau_leadtrkdz",Tau_leadtrkdz,"Tau_leadtrkdz[nTau]/F");
  T1->Branch("Tau_leadtrkpt",Tau_leadtrkpt,"Tau_leadtrkpt[nTau]/F");
  T1->Branch("Tau_leadtrketa",Tau_leadtrketa,"Tau_leadtrketa[nTau]/F");
  T1->Branch("Tau_leadtrkphi",Tau_leadtrkphi,"Tau_leadtrkphi[nTau]/F");
  T1->Branch("Tau_decayMode",Tau_decayMode,"Tau_decayMode[nTau]/I");
  T1->Branch("Tau_decayModeinding",Tau_decayModeinding,"Tau_decayModeinding[nTau]/O");
  T1->Branch("Tau_decayModeindingNewDMs",Tau_decayModeindingNewDMs,"Tau_decayModeindingNewDMs[nTau]/O");
  T1->Branch("Tau_eiso2018_raw",Tau_eiso2018_raw,"Tau_eiso2018_raw[nTau]/F");
  T1->Branch("Tau_eiso2018",Tau_eiso2018,"Tau_eiso2018[nTau]/I");
  T1->Branch("Tau_jetiso_deeptau2017v2p1_raw",Tau_jetiso_deeptau2017v2p1_raw,"Tau_jetiso_deeptau2017v2p1_raw[nTau]/F");
  T1->Branch("Tau_jetiso_deeptau2017v2p1",Tau_jetiso_deeptau2017v2p1,"Tau_jetiso_deeptau2017v2p1[nTau]/I");
  T1->Branch("Tau_eiso_deeptau2017v2p1_raw",Tau_eiso_deeptau2017v2p1_raw,"Tau_eiso_deeptau2017v2p1_raw[nTau]/F");
  T1->Branch("Tau_eiso_deeptau2017v2p1",Tau_eiso_deeptau2017v2p1,"Tau_eiso_deeptau2017v2p1[nTau]/I");
  T1->Branch("Tau_muiso_deeptau2017v2p1_raw",Tau_muiso_deeptau2017v2p1_raw,"Tau_muiso_deeptau2017v2p1_raw[nTau]/F");
  T1->Branch("Tau_muiso_deeptau2017v2p1",Tau_muiso_deeptau2017v2p1,"Tau_muiso_deeptau2017v2p1[nTau]/I");
  T1->Branch("Tau_rawiso",Tau_rawiso,"Tau_rawiso[nTau]/F");
  T1->Branch("Tau_rawisodR03",Tau_rawisodR03,"Tau_rawisodR03[nTau]/F");
  T1->Branch("Tau_puCorr",Tau_puCorr,"Tau_puCorr[nTau]/F");
  
  // MC Info //
  
  if(isMC){
	  
  // generator-related info //
  
  T1->Branch("Generator_weight", &Generator_weight, "Generator_weight/D") ;
  T1->Branch("Generator_qscale",&Generator_qscale,"Generator_qscale/F");
  T1->Branch("Generator_x1",&Generator_x1,"Generator_x1/F");
  T1->Branch("Generator_x2",&Generator_x2,"Generator_x2/F");
  T1->Branch("Generator_xpdf1",&Generator_xpdf1,"Generator_xpdf1/F");
  T1->Branch("Generator_xpdf2",&Generator_xpdf2,"Generator_xpdf2/F");
  T1->Branch("Generator_id1",&Generator_id1,"Generator_id1/I");
  T1->Branch("Generator_id2",&Generator_id2,"Generator_id2/I");
  T1->Branch("Generator_scalePDF",&Generator_scalePDF,"Generator_scalePDF/F");
  
  T1->Branch("npu_vert",&npu_vert,"npu_vert/I");
  T1->Branch("npu_vert_true",&npu_vert_true,"npu_vert_true/I");
	  
  // GEN MET info //    
  
  T1->Branch("GENMET_pt",&genmiset,"genmiset/F") ;
  T1->Branch("GENMET_phi",&genmisphi,"genmisphi/F") ;
  
  // GEN AK8 jet info //  
  
  T1->Branch("nGenJetAK8",&nGenJetAK8, "nGenJetAK8/I");
  T1->Branch("GenJetAK8_pt",GenJetAK8_pt,"GenJetAK8_pt[nGenJetAK8]/F");
  T1->Branch("GenJetAK8_eta",GenJetAK8_eta,"GenJetAK8_eta[nGenJetAK8]/F");
  T1->Branch("GenJetAK8_phi",GenJetAK8_phi,"GenJetAK8_phi[nGenJetAK8]/F");
  T1->Branch("GenJetAK8_mass",GenJetAK8_mass,"GenJetAK8_mass[nGenJetAK8]/F"); 
  T1->Branch("GenJetAK8_msoftdrop",GenJetAK8_sdmass,"GenJetAK8_sdmass[nGenJetAK8]/F");
  T1->Branch("GenJetAK8_hadronflav",GenJetAK8_hadronflav,"GenJetAK8_hadronflav[nGenJetAK8]/I");
  T1->Branch("GenJetAK8_partonflav",GenJetAK8_partonflav,"GenJetAK8_partonflav[nGenJetAK8]/I");

  if(store_fatjet_constituents){
    T1->Branch("nGenJetAK8_cons",&nGenJetAK8_cons,"nGenJetAK8_cons/I");
    T1->Branch("GenJetAK8_cons_pt",GenJetAK8_cons_pt, "GenJetAK8_cons_pt[nGenJetAK8_cons]/F");
    T1->Branch("GenJetAK8_cons_eta",GenJetAK8_cons_eta, "GenJetAK8_cons_eta[nGenJetAK8_cons]/F");
    T1->Branch("GenJetAK8_cons_phi",GenJetAK8_cons_phi, "GenJetAK8_cons_phi[nGenJetAK8_cons]/F");
    T1->Branch("GenJetAK8_cons_mass",GenJetAK8_cons_mass, "GenJetAK8_cons_mass[nGenJetAK8_cons]/F");
    T1->Branch("GenJetAK8_cons_pdgId",GenJetAK8_cons_pdgId, "GenJetAK8_cons_pdgId[nGenJetAK8_cons]/I");
    T1->Branch("GenJetAK8_cons_jetIndex",GenJetAK8_cons_jetIndex, "GenJetAK8_cons_jetIndex[nGenJetAK8_cons]/I");
  }
  
  // GEN AK4 jet info //  
 
  T1->Branch("nGenJetAK4",&nGenJetAK4, "nGenJetAK4/I");
  T1->Branch("GenJetAK4_pt",GenJetAK4_pt,"GenJetAK4_pt[nGenJetAK4]/F");
  T1->Branch("GenJetAK4_eta",GenJetAK4_eta,"GenJetAK4_eta[nGenJetAK4]/F");
  T1->Branch("GenJetAK4_phi",GenJetAK4_phi,"GenJetAK4_phi[nGenJetAK4]/F");
  T1->Branch("GenJetAK4_mass",GenJetAK4_mass,"GenJetAK4_mass[nGenJetAK4]/F");
  T1->Branch("GenJetAK4_hadronflav",GenJetAK4_hadronflav,"GenJetAK4_hadronflav[nGenJetAK4]/I");
  T1->Branch("GenJetAK4_partonflav",GenJetAK4_partonflav,"GenJetAK4_partonflav[nGenJetAK4]/I");
  
  // GEN particles info //  
  
  T1->Branch("nGenPart",&nGenPart, "nGenPart/I");
  T1->Branch("GenPart_pt",GenPart_pt,"GenPart_pt[nGenPart]/F");
  T1->Branch("GenPart_eta",GenPart_eta,"GenPart_eta[nGenPart]/F");
  T1->Branch("GenPart_phi",GenPart_phi,"GenPart_phi[nGenPart]/F");
  T1->Branch("GenPart_mass",GenPart_mass,"GenPart_mass[nGenPart]/F");
  T1->Branch("GenPart_status",GenPart_status,"GenPart_status[nGenPart]/I");
  T1->Branch("GenPart_pdgId",GenPart_pdg,"GenPart_pdg[nGenPart]/I");
  T1->Branch("GenPart_mompdgId",GenPart_mompdg,"GenPart_mompdg[nGenPart]/I");
  T1->Branch("GenPart_momstatus",GenPart_momstatus,"GenPart_momstatus[nGenPart]/I");
  T1->Branch("GenPart_grmompdgId",GenPart_grmompdg,"GenPart_grmompdg[nGenPart]/I");
  T1->Branch("GenPart_daugno",GenPart_daugno,"GenPart_daugno[nGenPart]/I");
  T1->Branch("GenPart_fromhard",GenPart_fromhard,"GenPart_fromhard[nGenPart]/O");
  T1->Branch("GenPart_fromhardbFSR",GenPart_fromhardbFSR,"GenPart_fromhardbFSR[nGenPart]/O");
  T1->Branch("GenPart_isPromptFinalState",GenPart_isPromptFinalState,"GenPart_isPromptFinalState[nGenPart]/O");
  T1->Branch("GenPart_isLastCopyBeforeFSR",GenPart_isLastCopyBeforeFSR,"GenPart_isLastCopyBeforeFSR[nGenPart]/O");
  T1->Branch("GenPart_isDirectPromptTauDecayProductFinalState",GenPart_isDirectPromptTauDecayProductFinalState,"GenPart_isDirectPromptTauDecayProductFinalState[nGenPart]/O");
  
  
  // LHE Info //
  
  T1->Branch("nLHEPart",&nLHEPart, "nLHEPart/I");
  T1->Branch("LHEPart_pdg",LHEPart_pdg,"LHEPart_pdg[nLHEPart]/I");
  T1->Branch("LHEPart_pt",LHEPart_pt,"LHEPart_pt[nLHEPart]/F");
  T1->Branch("LHEPart_eta",LHEPart_eta,"LHEPart_eta[nLHEPart]/F");
  T1->Branch("LHEPart_phi",LHEPart_phi,"LHEPart_phi[nLHEPart]/F");
  T1->Branch("LHEPart_m",LHEPart_m,"LHEPart_m[nLHEPart]/F");
  
  T1->Branch("LHE_weight",&LHE_weight, "LHE_weight/D");
  T1->Branch("nLHEScaleWeights",&nLHEScaleWeights, "nLHEScaleWeights/I");
  T1->Branch("LHEScaleWeights",LHEScaleWeights,"LHEScaleWeights[nLHEScaleWeights]/F");
  T1->Branch("nLHEPDFWeights",&nLHEPDFWeights, "nLHEPDFWeights/I");
  T1->Branch("LHEPDFWeights",LHEPDFWeights,"LHEPDFWeights[nLHEPDFWeights]/F");
  T1->Branch("nLHEAlpsWeights",&nLHEAlpsWeights, "nLHEAlpsWeights/I");
  T1->Branch("LHEAlpsWeights",LHEAlpsWeights,"LHEAlpsWeights[nLHEAlpsWeights]/F");
  T1->Branch("nLHEPSWeights",&nLHEPSWeights, "nLHEPSWeights/I");
  T1->Branch("LHEPSWeights",LHEPSWeights,"LHEPSWeights[nLHEPSWeights]/F");
  T1->Branch("wgt_isr_up", &wgt_isr_up, "wgt_isr_up/F");
  T1->Branch("wgt_isr_dn", &wgt_isr_dn, "wgt_isr_dn/F");
  T1->Branch("wgt_fsr_up", &wgt_fsr_up, "wgt_fsr_up/F");
  T1->Branch("wgt_fsr_dn", &wgt_fsr_dn, "wgt_fsr_dn/F");

  //BTag SF
  T1->Branch("BTAG_SF" , &BTAG_SF, "BTAG_SF/F");
  T1->Branch("BTAG_jes_up" , &BTAG_jes_up, "BTAG_jes_up/F");
  T1->Branch("BTAG_jes_dn" , &BTAG_jes_dn, "BTAG_jes_dn/F");
  T1->Branch("BTAG_lf_up", &BTAG_lf_up , "BTAG_lf_up/F");
  T1->Branch("BTAG_lf_dn", &BTAG_lf_dn , "BTAG_lf_dn/F");
  T1->Branch("BTAG_hf_up", &BTAG_hf_up , "BTAG_hf_up/F");
  T1->Branch("BTAG_hf_dn", &BTAG_hf_dn , "BTAG_hf_dn/F");
  T1->Branch("BTAG_hfstats1_up", &BTAG_hfstats1_up , "BTAG_hfstats1_up/F");
  T1->Branch("BTAG_hfstats1_dn", &BTAG_hfstats1_dn , "BTAG_hfstats1_dn/F");
  T1->Branch("BTAG_hfstats2_up", &BTAG_hfstats2_up , "BTAG_hfstats2_up/F");
  T1->Branch("BTAG_hfstats2_dn", &BTAG_hfstats2_dn , "BTAG_hfstats2_dn/F");
  T1->Branch("BTAG_lfstats1_up", &BTAG_lfstats1_up , "BTAG_lfstats1_up/F");
  T1->Branch("BTAG_lfstats1_dn", &BTAG_lfstats1_dn , "BTAG_lfstats1_dn/F");
  T1->Branch("BTAG_lfstats2_up", &BTAG_lfstats2_up , "BTAG_lfstats2_up/F");
  T1->Branch("BTAG_lfstats2_dn", &BTAG_lfstats2_dn , "BTAG_lfstats2_dn/F");
  T1->Branch("BTAG_cferr1_up", &BTAG_cferr1_up , "BTAG_cferr1_up/F");
  T1->Branch("BTAG_cferr1_dn", &BTAG_cferr1_dn , "BTAG_cferr1_dn/F");
  T1->Branch("BTAG_cferr2_up", &BTAG_cferr2_up , "BTAG_cferr2_up/F");
  T1->Branch("BTAG_cferr2_dn", &BTAG_cferr2_dn , "BTAG_cferr2_dn/F");
  } //isMC
  
  T2 = new TTree("Events_All", "XtoYH");
  
  T2->Branch("ievt", &ievt, "ievt/i");
  T2->Branch("PV_npvsGood", &PV_npvsGood, "PV_npvsGood/I");
  
  if(isMC){
	  
  T2->Branch("Generator_weight", &Generator_weight, "Generator_weight/D") ;
  T2->Branch("npu_vert",&npu_vert,"npu_vert/I");
  T2->Branch("npu_vert_true",&npu_vert_true,"npu_vert_true/I");
  T2->Branch("LHE_weight",&LHE_weight, "LHE_weight/D");
  
  }
  
  Nevt=0;
}


Leptop::~Leptop()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Leptop::analyze(const edm::Event& iEvent, const edm::EventSetup& pset) {
  
  using namespace edm;
  Nevt++;
  cout << Nevt << endl;
  
  irun = iEvent.id().run();
  ilumi = iEvent.luminosityBlock();
  ievt = iEvent.id().event();
  
  if (Nevt%100==1)cout <<"Leptop::analyze "<<Nevt<<" "<<iEvent.id().run()<<" "<<iEvent.id().event()<<endl;
    
  // First store all MC information //
  
  edm::Handle<reco::GenJetCollection> genjetAK8s;
  edm::Handle<reco::GenJetCollection> genjetAK4s;
 
  wgt_isr_up = 1.0, wgt_isr_dn = 1.0, wgt_fsr_up = 1.0, wgt_fsr_dn = 1.0;

  if(isMC){
	
	// MC weights //
	  
    edm::Handle<GenEventInfoProduct>eventinfo ;  
    iEvent.getByToken(tok_wt_,eventinfo) ;

 // implementation of PS weight //    


/*    nLHEPSWeights = 0;
    
    if (eventinfo.isValid()){
		
		Generator_weight = eventinfo->weight();
		
		// Generator information //
		
		Generator_qscale = eventinfo->qScale();
		Generator_x1 = (*eventinfo->pdf()).x.first;
        Generator_x2 = (*eventinfo->pdf()).x.second;
        Generator_id1 = (*eventinfo->pdf()).id.first;
        Generator_id2 = (*eventinfo->pdf()).id.second;
        Generator_xpdf1 = (*eventinfo->pdf()).xPDF.first;
        Generator_xpdf2 = (*eventinfo->pdf()).xPDF.second;
        Generator_scalePDF = (*eventinfo->pdf()).scalePDF;
        
        //cout<<"eventinfo->weights().size() "<<eventinfo->weights().size()<<" GEN weight "<<Generator_weight<<endl;
        
        // Parton shower weights //
        if(eventinfo->weights().size()>2){
			for(unsigned int i=2; i<eventinfo->weights().size(); ++i){
				LHEPSWeights[nLHEPSWeights] = eventinfo->weights()[i]/eventinfo->weights()[1];
				nLHEPSWeights++;
				if(nLHEPSWeights >= nlhepsmax) break;
			}
		}
        
    }
    else
      {
		Generator_weight = Generator_qscale = Generator_x1 = Generator_x2 = Generator_id1 = Generator_id2 = Generator_xpdf1 = Generator_xpdf2 = Generator_scalePDF = -10000;
      }
*/

    nLHEPSWeights = 8;
      if (eventinfo.isValid()){
        Generator_weight = eventinfo->weight();
        Generator_qscale = eventinfo->qScale();
        Generator_x1 = (*eventinfo->pdf()).x.first;
        Generator_x2 = (*eventinfo->pdf()).x.second;
        Generator_id1 = (*eventinfo->pdf()).id.first;
        Generator_id2 = (*eventinfo->pdf()).id.second;
        Generator_xpdf1 = (*eventinfo->pdf()).xPDF.first;
        Generator_xpdf2 = (*eventinfo->pdf()).xPDF.second;
        Generator_scalePDF = (*eventinfo->pdf()).scalePDF;
        if(eventinfo->weights().size()>2){
                LHEPSWeights[0] = eventinfo->weights()[2]/eventinfo->weights()[1];
                LHEPSWeights[1] = eventinfo->weights()[3]/eventinfo->weights()[1];
                LHEPSWeights[2] = eventinfo->weights()[4]/eventinfo->weights()[1];
                LHEPSWeights[3] = eventinfo->weights()[5]/eventinfo->weights()[1];
                LHEPSWeights[4] = eventinfo->weights()[24]/eventinfo->weights()[1];
                LHEPSWeights[5] = eventinfo->weights()[25]/eventinfo->weights()[1];
                LHEPSWeights[6] = eventinfo->weights()[26]/eventinfo->weights()[1];
                LHEPSWeights[7] = eventinfo->weights()[27]/eventinfo->weights()[1];
                wgt_isr_up = eventinfo->weights()[25]/eventinfo->weights()[1]; wgt_isr_dn = eventinfo->weights()[24]/eventinfo->weights()[1];
                wgt_fsr_up = eventinfo->weights()[3] /eventinfo->weights()[1]; wgt_fsr_dn = eventinfo->weights()[2] /eventinfo->weights()[1];
                }
        }
       else
        {
                Generator_weight = Generator_qscale = Generator_x1 = Generator_x2 = Generator_id1 = Generator_id2 = Generator_xpdf1 = Generator_xpdf2 = Generator_scalePDF = -10000;
        }


    
    // LHE-level particles //
      
    edm::Handle<LHEEventProduct>lheeventinfo ;
    iEvent.getByToken(lheEventProductToken_,lheeventinfo) ;
    
    nLHEScaleWeights = 0;
	nLHEPDFWeights = 0;
	nLHEAlpsWeights = 0;
    
    nLHEPart = 0;
    
    if(lheeventinfo.isValid()){
      
      // LHE-level particles //
      
      const auto & hepeup = lheeventinfo->hepeup();
      const auto & pup = hepeup.PUP;
      
      for (unsigned int ij = 0; ij  < pup.size(); ++ij) {
		if(hepeup.ISTUP[ij]==1){// status==1 --> particles which stay up to final state                          
			TLorentzVector p4(pup[ij][0], pup[ij][1], pup[ij][2], pup[ij][3]);
			LHEPart_pt[nLHEPart] = p4.Pt();
			LHEPart_eta[nLHEPart] = p4.Eta();
			LHEPart_phi[nLHEPart] = p4.Phi();
			LHEPart_m[nLHEPart] = p4.M();
			LHEPart_pdg[nLHEPart] = (hepeup.IDUP[ij]);
			nLHEPart++;
			if(nLHEPart>=nlhemax) break;
			}
		}
	
	 // LHE-level weights //
	
	  LHE_weight = lheeventinfo->originalXWGTUP();
	  	  
	  for ( unsigned int index = 0; index < lheeventinfo->weights().size(); ++index ) {	
		//cout<<"Index "<<index+1<<" Id "<<lheeventinfo->weights()[index].id<<" weight "<<lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP()<<endl;//" muR "<<lheeventinfo->weights()[index].MUR<<" muF "<<lheeventinfo->weights()[index].MUF<<" DYN Scale "<<lheeventinfo->weights()[index].DYN_SCALE<<endl;
		if(index<nlhescalemax && nLHEScaleWeights<nlhescalemax){
			LHEScaleWeights[nLHEScaleWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
			nLHEScaleWeights++;
		}
		if(index>=nlhescalemax && index<(nlhescalemax+nPDFsets)  && nLHEPDFWeights<nlhepdfmax){
			LHEPDFWeights[nLHEPDFWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
			nLHEPDFWeights++;
		}
		if(index>=(nlhescalemax+nPDFsets) && index<(nlhescalemax+nPDFsets+nalpsmax) && nLHEAlpsWeights<nalpsmax){
			LHEAlpsWeights[nLHEAlpsWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
			nLHEAlpsWeights++;
			}
	  }
	  		
	}// lheeventinfo.isValid()

    // Flavor tagging of GEN jets using ghost-matching //                                                             

	edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos);

	// AK8 GEN jet //

    nGenJetAK8 = 0;
    nGenJetAK8_cons = 0;
    
    iEvent.getByToken(tok_genjetAK8s_, genjetAK8s);
    
    if(genjetAK8s.isValid()){
		
		std::vector<int> partonFlavour_AK8;
		std::vector<int> hadronFlavour_AK8;

		for (const reco::GenJet & jet : *genjetAK8s) {
			
			bool matched = false;
			for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
				if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < 0.1) {
					partonFlavour_AK8.push_back(jetFlavourInfoMatching.second.getPartonFlavour());
					hadronFlavour_AK8.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
					matched = true;
					break;
				}
			}
		
			if (!matched) {
				partonFlavour_AK8.push_back(-100);
				hadronFlavour_AK8.push_back(-100);
			}
		}
      
        JetDefinition pfjetAK8_Def(antikt_algorithm,0.8,E_scheme);
		SoftDrop sd(beta,z_cut,0.8);
      
		for(unsigned gjet = 0; gjet<genjetAK8s->size(); gjet++)	{
	
			TLorentzVector genjetAK8_4v((*genjetAK8s)[gjet].px(),(*genjetAK8s)[gjet].py(),(*genjetAK8s)[gjet].pz(), (*genjetAK8s)[gjet].energy());
			if(genjetAK8_4v.Pt()<AK8GenPtCut) continue;
			if(abs(genjetAK8_4v.Eta())>maxgenEta) continue;
	
			GenJetAK8_pt[nGenJetAK8] = genjetAK8_4v.Pt();
			GenJetAK8_eta[nGenJetAK8] = genjetAK8_4v.Eta();
			GenJetAK8_phi[nGenJetAK8] = genjetAK8_4v.Phi();
			GenJetAK8_mass[nGenJetAK8] = (*genjetAK8s)[gjet].mass();
			GenJetAK8_hadronflav[nGenJetAK8] = (int)hadronFlavour_AK8[gjet];
			GenJetAK8_partonflav[nGenJetAK8] = partonFlavour_AK8[gjet];
			
			std::vector<reco::CandidatePtr> daught((*genjetAK8s)[gjet].daughterPtrVector());
	
			vector <fastjet::PseudoJet> fjInputs;
			fjInputs.resize(0);
			for (unsigned int i2 = 0; i2< daught.size(); ++i2) {
				PseudoJet psjet = PseudoJet( (*daught[i2]).px(),(*daught[i2]).py(),(*daught[i2]).pz(),(*daught[i2]).energy() );
				psjet.set_user_index(i2);
				fjInputs.push_back(psjet); 
        
        // Storing 4-momenta of jet constituents//
        if(store_fatjet_constituents && nGenJetAK8<njetconsmax){
          GenJetAK8_cons_pt[nGenJetAK8_cons] = daught[i2]->pt();
          GenJetAK8_cons_eta[nGenJetAK8_cons] = daught[i2]->eta();
          GenJetAK8_cons_phi[nGenJetAK8_cons] = daught[i2]->phi();
          GenJetAK8_cons_mass[nGenJetAK8_cons] = daught[i2]->mass();
          GenJetAK8_cons_pdgId[nGenJetAK8_cons] = daught[i2]->pdgId();
          GenJetAK8_cons_jetIndex[nGenJetAK8_cons] = nGenJetAK8;   
          nGenJetAK8_cons++;
        }
        // end of candidate storage //
        
			} //i2
	
			vector <fastjet::PseudoJet> sortedJets;
			fastjet::ClusterSequence clustSeq(fjInputs, pfjetAK8_Def);
			fjInputs.clear();
			sortedJets    = sorted_by_pt(clustSeq.inclusive_jets());
	
			if(sortedJets.size()>0){
				GenJetAK8_sdmass[nGenJetAK8] = (sd(sortedJets[0])).m();
			}
			sortedJets.clear();
	
			if (++nGenJetAK8>=njetmxAK8) break;
	
		}
	}
      
    // AK4 GEN jet //
      
	nGenJetAK4 = 0;
	iEvent.getByToken(tok_genjetAK4s_, genjetAK4s);
	
	if(genjetAK4s.isValid()){
		
		std::vector<int> partonFlavour_AK4;
		std::vector<int> hadronFlavour_AK4;
      
		for (const reco::GenJet & jet : *genjetAK4s) {
		  
			bool matched = false;
			for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
				if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) < 0.1) {
					partonFlavour_AK4.push_back(jetFlavourInfoMatching.second.getPartonFlavour());
					hadronFlavour_AK4.push_back(jetFlavourInfoMatching.second.getHadronFlavour());
					matched = true;
					break;
				}
			}
		
			if (!matched) {
				partonFlavour_AK4.push_back(-100);
				hadronFlavour_AK4.push_back(-100);
			}	
		}
	
		for(unsigned gjet = 0; gjet<genjetAK4s->size(); gjet++)	{
	
			TLorentzVector genjetAK44v((*genjetAK4s)[gjet].px(),(*genjetAK4s)[gjet].py(),(*genjetAK4s)[gjet].pz(), (*genjetAK4s)[gjet].energy());
			if(genjetAK44v.Pt()<minGenPt) continue;
			if(abs(genjetAK44v.Eta())>maxgenEta) continue;
	
			GenJetAK4_pt[nGenJetAK4] = genjetAK44v.Pt();
			GenJetAK4_eta[nGenJetAK4] = genjetAK44v.Eta();
			GenJetAK4_phi[nGenJetAK4] = genjetAK44v.Phi();
			GenJetAK4_mass[nGenJetAK4] = (*genjetAK4s)[gjet].mass();
			GenJetAK4_hadronflav[nGenJetAK4] = (int)hadronFlavour_AK4[gjet];
			GenJetAK4_partonflav[nGenJetAK4] = partonFlavour_AK4[gjet];

			if (++nGenJetAK4>=njetmx) break;
      
		}
		
    }
    
    // Gen particles //
    
	nGenPart = 0;
    edm::Handle<std::vector<reco::GenParticle>> genparticles;  
	iEvent.getByToken(tok_genparticles_,genparticles);
    
	if(genparticles.isValid()){
	
		for(unsigned ig=0; ig<(genparticles->size()); ig++){
			
			if(!(((*genparticles)[ig].status()==1)||(abs((*genparticles)[ig].status())==22)||((*genparticles)[ig].status()==23)|((*genparticles)[ig].status()==2))) continue;
			if(!((*genparticles)[ig].isHardProcess()||(*genparticles)[ig].fromHardProcessBeforeFSR()||(*genparticles)[ig].isLastCopyBeforeFSR()||(*genparticles)[ig].isDirectPromptTauDecayProductFinalState())) continue;
	  
			if(!((abs((*genparticles)[ig].pdgId())>=1 && abs((*genparticles)[ig].pdgId())<=6) || (abs((*genparticles)[ig].pdgId())>=11 && abs((*genparticles)[ig].pdgId())<=16) || (abs((*genparticles)[ig].pdgId())>=22 && abs((*genparticles)[ig].pdgId())<=24))) continue;
			// important condition on pdg id -> May be changed in other analyses //
	  
			GenPart_status[nGenPart] = (*genparticles)[ig].status();
			GenPart_pdg[nGenPart] = (*genparticles)[ig].pdgId();
			GenPart_daugno[nGenPart] = (*genparticles)[ig].numberOfDaughters();
			GenPart_fromhard[nGenPart] = (*genparticles)[ig].isHardProcess();
			GenPart_fromhardbFSR[nGenPart] = (*genparticles)[ig].fromHardProcessBeforeFSR();
			GenPart_isLastCopyBeforeFSR[nGenPart] = (*genparticles)[ig].isLastCopyBeforeFSR();
			GenPart_isPromptFinalState[nGenPart] = (*genparticles)[ig].isPromptFinalState();
			GenPart_isDirectPromptTauDecayProductFinalState[nGenPart] = (*genparticles)[ig].isDirectPromptTauDecayProductFinalState();
			GenPart_pt[nGenPart] = (*genparticles)[ig].pt();
			GenPart_eta[nGenPart] = (*genparticles)[ig].eta();
			GenPart_phi[nGenPart] = (*genparticles)[ig].phi();
			GenPart_mass[nGenPart] = (*genparticles)[ig].mass();
			
			int mompdg, momstatus, grmompdg;
			mompdg = momstatus = grmompdg = 0;
			
			if((*genparticles)[ig].numberOfMothers()>0){
				
				// mother pdg id & status //
			
				const Candidate * mom = (*genparticles)[ig].mother();
				
				while(mom->pdgId() == (*genparticles)[ig].pdgId())
				{
					mom = mom->mother();
				}
				
				if(!(*genparticles)[ig].isPromptFinalState() && !(*genparticles)[ig].isDirectPromptTauDecayProductFinalState()){
					while(mom->status()==2){
						mom = mom->mother();	
					}
				}
				
				mompdg = mom->pdgId();
				momstatus = mom->status();
	  
				// grand-mother pdg id //
						
				if(mom->numberOfMothers()>0){
	
					const Candidate * grmom = mom->mother();
					
					while(grmom->pdgId() == mompdg)
					{
						if(grmom->numberOfMothers()>0){
							grmom = grmom->mother();
						}
						else{ break; }
					}
					
					grmompdg  = grmom->pdgId();	
				} 
				
			}
			
			GenPart_mompdg[nGenPart] = mompdg;
			GenPart_momstatus[nGenPart] = momstatus;
			GenPart_grmompdg[nGenPart] = grmompdg; 
			
			nGenPart++;
			if(nGenPart>=npartmx) break;
		}
    }
    
    // Pileup information //
    
    npu_vert = 0;
	npu_vert_true = 0;
    
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(pileup_, PupInfo);
    if (PupInfo.isValid()) {
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
		if (PVI->getBunchCrossing()==0) {
			npu_vert = PVI->getPU_NumInteractions();
			npu_vert_true = PVI->getTrueNumInteractions();
			break;
			}
		}
    }

  }//isMC
  
  // Primary vertex info //
  
  Handle<VertexCollection> primaryVertices;
  iEvent.getByToken(tok_primaryVertices_, primaryVertices);
  reco::Vertex vertex;
  
  if (primaryVertices.isValid()) {
	  
	if(primaryVertices->size() > 0){  
		vertex = primaryVertices->at(0); 
		PV_ndof = vertex.ndof();
		PV_chi2 = vertex.normalizedChi2();
		PV_x = vertex.position().x();
		PV_y = vertex.position().y();
		PV_z = vertex.position().z();
	} 
	 
    int ndofct_org=0;
    int nchict_org=0;
    int nvert_org = 0;
    int nprimi_org = 0;
    
    for (reco::VertexCollection::const_iterator vert=primaryVertices->begin(); vert<primaryVertices->end(); vert++) {
      nvert_org++;
      if (vert->isValid() && !vert->isFake()) {
		if (vert->ndof() > 4 && fabs(vert->position().z()) <= 24 && fabs(vert->position().Rho()) <= 2) {
			nprimi_org++;
			}
		if (vert->ndof()>7) {
			ndofct_org++;
			if (vert->normalizedChi2()<5) nchict_org++;
			}
		}
    }
    
    nprim = min(999,nvert_org) + 1000*min(999,ndofct_org) + 1000000*min(999,nchict_org);
    npvert = nchict_org;
    PV_npvsGood = nprimi_org;
    
  } else { 
    nprim = -100;
    npvert = -100;
    PV_npvsGood = -100;
    PV_ndof = -100;
    PV_chi2 = PV_x = PV_y = PV_z = -1000;
  }
  
  reco::TrackBase::Point beamPoint(0,0, 0);
  edm::Handle<reco::BeamSpot> beamSpotH;
  
  iEvent.getByToken(tok_beamspot_, beamSpotH);  //Label("offlineBeamSpot",beamSpotH);
  if (beamSpotH.isValid()){
    beamPoint = beamSpotH->position();
  }

  edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVertices;
  iEvent.getByToken(tok_sv,secondaryVertices);
  
  // Energy density info //
  
  edm::Handle<double> Rho_PF;
  iEvent.getByToken(tok_Rho_,Rho_PF);
  Rho = *Rho_PF;
  
  // Store trigger information //
  
  bool booltrg[nHLTmx]= {false};
  
  if(!isFastSIM){
  
	const char* variab1;
  
	edm::Handle<edm::TriggerResults> trigRes;
	iEvent.getByToken(triggerBits_, trigRes);
  
	const edm::TriggerNames &names = iEvent.triggerNames(*trigRes);
  
	edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
	iEvent.getByToken(triggerObjects_, triggerObjects);
  
	edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
	iEvent.getByToken(triggerPrescales_, triggerPrescales);
  
	int ihlttrg[nHLTmx+1]= {0};
  
	for (int jk=0; jk<nHLTmx; jk++) {
		for(unsigned ij = 0; ij<trigRes->size(); ++ij) {
			std::string name = names.triggerName(ij);
			variab1 = name.c_str(); 
			if (strstr(variab1,hlt_name[jk]) && ((strlen(variab1)-strlen(hlt_name[jk]))<5))
			{
				if ((trigRes->accept(ij))){   //||(isMC)) {
					ihlttrg[jk] = ihlttrg[nHLTmx] = 1;
					booltrg[jk] = true;
					break;
				}
			}
		}//ij     
	}//jk
  
	trig_value = 1; 
  
	for (int jk=1; jk<(nHLTmx+1); jk++) {
		if(ihlttrg[nHLTmx-jk]>0) {
			trig_value+=(1<<jk);
		}
	}
  
  // Trigger objects //
    
	vector<triggervar> alltrgobj;
	if (trigRes.isValid()) { 
    
		const char* variab2 ;
    
		alltrgobj.clear(); 
    
		for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      
			obj.unpackPathNames(names);
			std::vector<std::string> pathNamesAll  = obj.pathNames(false);
      
			for (unsigned ih = 0, n = pathNamesAll.size(); ih < n; ih++) {
	
				variab2 = pathNamesAll[ih].c_str(); 
	
				for (int jk=0; jk<nHLTmx; jk++) {
					if (strstr(variab2,hlt_name[jk]) && (strlen(variab2)-strlen(hlt_name[jk])<5)) {
	    	    
						if(obj.pt()>20 && fabs(obj.eta())<3.0) {
	      
							triggervar tmpvec1;
	      
							tmpvec1.both = obj.hasPathName( pathNamesAll[ih], true, true );
							tmpvec1.highl  = obj.hasPathName( pathNamesAll[ih], false, true );
							tmpvec1.level1 = obj.hasPathName( pathNamesAll[ih], true, false );
							tmpvec1.trg4v = TLorentzVector(obj.px(), obj.py(), obj.pz(), obj.energy());
							tmpvec1.pdgId = obj.pdgId();
							tmpvec1.prescl = 1;    //triggerPrescales->getPrescaleForIndex(ih);
							tmpvec1.ihlt = jk;
							tmpvec1.type = (obj.type(92) + 2*(obj.coll("hltEgammaCandidates")) + 4*obj.type(83) + 8*(obj.coll("hltIterL3MuonCandidates")) + 16*obj.type(84) + 32*(obj.coll("*Tau*"))
								  + 64*obj.type(87) + 128*(obj.coll("L1ETM")) + 256*obj.type(85));
								  //order: e/gamma + muon + tau + met + jet
							alltrgobj.push_back(tmpvec1);
							break;
						}
					}
				}//jk 
			}//ih
		}
	}
    
	int xht=0;
	nTrigObj = alltrgobj.size();
	if(nTrigObj>njetmx) { nTrigObj = njetmx; }
	if(nTrigObj > 0){
		for(unsigned int iht=0; iht<(unsigned int)nTrigObj; iht++){
			if(alltrgobj[iht].trg4v.Pt()>20 && fabs(alltrgobj[iht].trg4v.Eta())<3.0) {
				TrigObj_pt[xht] = alltrgobj[iht].trg4v.Pt();
				TrigObj_eta[xht] = alltrgobj[iht].trg4v.Eta();
				TrigObj_phi[xht] = alltrgobj[iht].trg4v.Phi();
				TrigObj_mass[xht] = alltrgobj[iht].trg4v.M();
				TrigObj_HLT[xht] = alltrgobj[iht].highl;
				TrigObj_L1[xht] = alltrgobj[iht].level1;
				TrigObj_Both[xht] = alltrgobj[iht].both;
				TrigObj_Ihlt[xht] = alltrgobj[iht].ihlt;
				TrigObj_pdgId[xht] = alltrgobj[iht].pdgId;
				TrigObj_type[xht] = alltrgobj[iht].type;
				xht++;
				if(xht>=njetmx) break;
			}
		//if(iht == (njetmx-1)) break;
		}
	}

  } // !isFastSIM
  
   // fill trigger booleans //
 
//-------------------------------------------2018 triggers---------------------------------------//
  for(int jk=0; jk<nHLTmx; jk++) {
	  if(jk==0) 	  {  hlt_IsoMu24 = booltrg[jk]; }
	  else if(jk==1)  {  hlt_Mu50 = booltrg[jk]; }
	  else if(jk==2)  {  hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = booltrg[jk]; }
	  else if(jk==3)  {  hlt_Ele115_CaloIdVT_GsfTrkIdT = booltrg[jk]; }
	  else if(jk==4)  {  hlt_Ele40_WPTight_Gsf = booltrg[jk]; }
	  else if(jk==5)  {  hlt_Ele32_WPTight_Gsf = booltrg[jk]; }
	  else if(jk==6)  {  hlt_Ele28_eta2p1_WPTight_Gsf_HT150 = booltrg[jk]; }
	  else if(jk==7)  {  hlt_Mu17_Mu8 = booltrg[jk]; }
	  else if(jk==8)  {  hlt_Ele23_Ele12 = booltrg[jk]; }
	  else if(jk==9)  {  hlt_DoubleEle33 = booltrg[jk]; }
	  else if(jk==10) {  hlt_DoubleEle25_CaloIdL_MW = booltrg[jk]; }
	  else if(jk==11) {  hlt_Diphoton_30_18 = booltrg[jk]; }
	  else if(jk==12) {  hlt_PFJet500 = booltrg[jk]; }
	  else if(jk==13) {  hlt_HT1050 = booltrg[jk]; }
	  else if(jk==14) {  hlt_AK8PFJet400_TrimMass30 = booltrg[jk]; }
	  else if(jk==15) {  hlt_AK8PFHT800_TrimMass50 = booltrg[jk]; }
	  else if(jk==16) {  hlt_Photon200 = booltrg[jk]; }
	  else if(jk==17) {  hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = booltrg[jk]; }
	  else if(jk==18) {  hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 = booltrg[jk]; }
	  else if(jk==19) {  hlt_PFMETNoMu140_PFMHTNoMu140_IDTight = booltrg[jk]; }
	  else if(jk==20) {  hlt_PFMETTypeOne140_PFMHT140_IDTight = booltrg[jk]; }
  } 

//-------------------------------------------2017 triggers-----------------------------------------//
/* for(int jk=0; jk<nHLTmx; jk++) {
          if(jk==0)       {  hlt_IsoMu24 = booltrg[jk]; }
          else if(jk==1)  {  hlt_IsoMu27 = booltrg[jk]; }
          else if(jk==2)  {  hlt_Ele50_CaloIdVT_GsfTrkIdT_PFJet165 = booltrg[jk]; }
          else if(jk==3)  {  hlt_Ele115_CaloIdVT_GsfTrkIdT = booltrg[jk]; }
          else if(jk==4)  {  hlt_Ele40_WPTight_Gsf = booltrg[jk]; }
          else if(jk==5)  {  hlt_Ele32_WPTight_Gsf = booltrg[jk]; }
          else if(jk==6)  {  hlt_Ele28_eta2p1_WPTight_Gsf_HT150 = booltrg[jk]; }
          else if(jk==7)  {  hlt_Mu17_Mu8 = booltrg[jk]; }
          else if(jk==8)  {  hlt_Ele23_Ele12 = booltrg[jk]; }
          else if(jk==9)  {  hlt_DoubleEle33 = booltrg[jk]; }
          else if(jk==10) {  hlt_DoubleEle25_CaloIdL_MW = booltrg[jk]; }
          else if(jk==11) {  hlt_Diphoton_30_18 = booltrg[jk]; }
          else if(jk==12) {  hlt_PFJet500 = booltrg[jk]; }
          else if(jk==13) {  hlt_HT1050 = booltrg[jk]; }
          else if(jk==14) {  hlt_AK8PFJet400_TrimMass30 = booltrg[jk]; }
          else if(jk==15) {  hlt_AK8PFHT800_TrimMass50 = booltrg[jk]; }
          else if(jk==16) {  hlt_Photon200 = booltrg[jk]; }
          else if(jk==17) {  hlt_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = booltrg[jk]; }
          else if(jk==18) {  hlt_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60 = booltrg[jk]; }
          else if(jk==19) {  hlt_PFMETNoMu140_PFMHTNoMu140_IDTight = booltrg[jk]; }
          else if(jk==20) {  hlt_PFMETTypeOne140_PFMHT140_IDTight = booltrg[jk]; }
  }
*/
//-----------------------------------------------2016 triggers--------------------------------------------//
/*  for(int jk=0; jk<nHLTmx; jk++) {
          if(jk==0)       {  hlt_IsoMu24 = booltrg[jk]; }
          else if(jk==1)  {  hlt_IsoTkMu24 = booltrg[jk]; }
          else if(jk==2)  {  hlt_Ele27_WPTight_Gsf = booltrg[jk]; }
          else if(jk==3)  {  hlt_Ele115_CaloIdVT_GsfTrkIdT = booltrg[jk]; }
          else if(jk==4)  {  hlt_Ele40_WPTight_Gsf = booltrg[jk]; }
          else if(jk==5)  {  hlt_Ele32_WPTight_Gsf = booltrg[jk]; }
          else if(jk==6)  {  hlt_Ele28_eta2p1_WPTight_Gsf_HT150 = booltrg[jk]; }
          else if(jk==7)  {  hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 = booltrg[jk]; }
          else if(jk==8)  {  hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL = booltrg[jk]; }
          else if(jk==9)  {  hlt_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ = booltrg[jk]; }
          else if(jk==10) {  hlt_DoubleEle25_CaloIdL_MW = booltrg[jk]; }
          else if(jk==11) {  hlt_Diphoton_30_18 = booltrg[jk]; }
          else if(jk==12) {  hlt_PFJet500 = booltrg[jk]; }
          else if(jk==13) {  hlt_HT1050 = booltrg[jk]; }
          else if(jk==14) {  hlt_AK8PFJet400_TrimMass30 = booltrg[jk]; }
          else if(jk==15) {  hlt_AK8PFHT800_TrimMass50 = booltrg[jk]; }
          else if(jk==16) {  hlt_Photon200 = booltrg[jk]; }
          else if(jk==17) {  hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL = booltrg[jk]; }
          else if(jk==18) {  hlt_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ = booltrg[jk]; }
          else if(jk==19) {  hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL = booltrg[jk]; }
          else if(jk==20) {  hlt_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ = booltrg[jk]; }
 }
*/
  
  // trigger filling end //
  
  // End of trigger info //
  
  // Prefire weights //
  
  if(add_prefireweights){
  
	edm::Handle< double > theprefweight;
	iEvent.getByToken(prefweight_token, theprefweight ) ;	
	prefiringweight =(*theprefweight);

	edm::Handle< double > theprefweightup;
	iEvent.getByToken(prefweightup_token, theprefweightup ) ;
	prefiringweightup =(*theprefweightup);
    
	edm::Handle< double > theprefweightdown;
	iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;   
	prefiringweightdown =(*theprefweightdown);
 
  }
  
  // End of prefire weights //
  
  // ====== RECO-objects now  ==========//
  
  // MET //
    
  miset = misphi = misetsig = sumEt = genmiset = genmisphi = genmisetsig = -1000 ;
  miset_UnclusEup = miset_UnclusEdn = -100;
  misphi_UnclusEup = misphi_UnclusEdn = -100;
  
  // MET uncertainty ids are taken from: https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/DataFormats/PatCandidates/interface/MET.h#L152-L158 //
  
  // CHS MET //
  
  edm::Handle<pat::METCollection> pfmet_ ;
  iEvent.getByToken(tok_mets_,pfmet_) ;
  
  if(pfmet_.isValid()){
	  
    const pat::MET &met = pfmet_->front();
    
    miset = met.corPt(); //met.pt();
    misphi = met.corPhi();//met.phi();
    misetsig = met.significance();
    sumEt = met.corSumEt();//sumEt();
    
    miset_UnclusEup = met.shiftedPt(pat::MET::METUncertainty(10));
	miset_UnclusEdn = met.shiftedPt(pat::MET::METUncertainty(11));
	
	misphi_UnclusEup = met.shiftedPhi(pat::MET::METUncertainty(10));
	misphi_UnclusEdn = met.shiftedPhi(pat::MET::METUncertainty(11));
	    
    if(isMC){
      genmiset = met.genMET()->pt();
      genmisphi = met.genMET()->phi();
      genmisetsig = met.genMET()->significance();
    }
  }
  
  // PUPPI MET //
  
  edm::Handle<pat::METCollection> pfmet_PUPPI_ ;
  iEvent.getByToken(tok_mets_PUPPI_,pfmet_PUPPI_) ;
  
  miset_PUPPI = misphi_PUPPI = misetsig_PUPPI = sumEt_PUPPI = -100;
  miset_PUPPI_JESup = miset_PUPPI_JESdn = miset_PUPPI_JERup = miset_PUPPI_JERdn = miset_PUPPI_UnclusEup = miset_PUPPI_UnclusEdn = -100;
  misphi_PUPPI_JESup = misphi_PUPPI_JESdn = misphi_PUPPI_JERup = misphi_PUPPI_JERdn = misphi_PUPPI_UnclusEup = misphi_PUPPI_UnclusEdn = -100;
  
  if(pfmet_PUPPI_.isValid()){
	
	const pat::MET &met = pfmet_PUPPI_->front();
	  
	miset_PUPPI = met.corPt(); 
    misphi_PUPPI = met.corPhi();
    misetsig_PUPPI = met.significance();
    sumEt_PUPPI = met.corSumEt();
    
    miset_PUPPI_JESup = met.shiftedPt(pat::MET::METUncertainty(2));
	miset_PUPPI_JESdn = met.shiftedPt(pat::MET::METUncertainty(3));
	miset_PUPPI_JERup = met.shiftedPt(pat::MET::METUncertainty(0));
	miset_PUPPI_JERdn = met.shiftedPt(pat::MET::METUncertainty(1));
	miset_PUPPI_UnclusEup = met.shiftedPt(pat::MET::METUncertainty(10));
	miset_PUPPI_UnclusEdn = met.shiftedPt(pat::MET::METUncertainty(11));
	
	misphi_PUPPI_JESup = met.shiftedPhi(pat::MET::METUncertainty(2));
	misphi_PUPPI_JESdn = met.shiftedPhi(pat::MET::METUncertainty(3));
	misphi_PUPPI_JERup = met.shiftedPhi(pat::MET::METUncertainty(0));
	misphi_PUPPI_JERdn = met.shiftedPhi(pat::MET::METUncertainty(1));
	misphi_PUPPI_UnclusEup = met.shiftedPhi(pat::MET::METUncertainty(10));
	misphi_PUPPI_UnclusEdn = met.shiftedPhi(pat::MET::METUncertainty(11));
	
	//See DataFormats/PatCandidates/interface/MET.h for the names of uncertainty sources //
	
  }
  
  // Muons //
    
  nMuon = 0;                                                                                                                                        
  std::vector<pat::Muon> tlvmu;
  edm::Handle<edm::View<pat::Muon>> muons;                                                                                                          
  iEvent.getByToken(tok_muons_, muons);                                                                                                             
    
  if(muons.isValid() && muons->size()>0) {                                                                                                           
    
	edm::View<pat::Muon>::const_iterator muon1;                                                                                                      

    for( muon1 = muons->begin(); muon1 < muons->end(); muon1++ ) {                                                                                   

		if (StoreMuon(*muon1,minmuPt,maxEta)) {                                                                
			
			Muon_pt[nMuon] = muon1->pt();                                                                         
			TrackRef trktrk = muon1->innerTrack();                                                                                                       
			Muon_p[nMuon] = trktrk->p()*muon1->charge();                                                                                                                 
			Muon_eta[nMuon] = muon1->eta();                                                                                                              
			Muon_phi[nMuon] = muon1->phi();    
			Muon_e[nMuon] = muon1->energy();
			Muon_ptErr[nMuon] = trktrk->ptError();                                                                                                                                                                                                                                
                                                                                                                                                       
			//MiniIsolation: begin//                                                                                      
			vector<float> isovalues;
			Read_MiniIsolation(muon1,Rho,isovalues);
			Muon_minisoall[nMuon] = isovalues[0];
			Muon_minchiso[nMuon] = isovalues[1];
			Muon_minnhiso[nMuon] = isovalues[2];
			Muon_minphiso[nMuon] = isovalues[3];
			//MiniIsolation: end//  
			                                         
			// Basic id variables //    
			                                  
			Muon_isPF[nMuon] = muon1->isPFMuon();                                                                                                        
			Muon_isGL[nMuon] = muon1->isGlobalMuon();                                                                                                    
			Muon_isTRK[nMuon] = muon1->isTrackerMuon();                                                                                                  
			Muon_isLoose[nMuon] = (muon::isLooseMuon(*muon1));                                                                                           
			Muon_isMed[nMuon] = (muon::isMediumMuon(*muon1));                                                                                            
			Muon_isMedPr[nMuon] = false;                                                                          
			if(muon::isMediumMuon(*muon1)) {                                                                                                             
				if ((std::abs(muon1->muonBestTrack()->dz(vertex.position())) < 0.1) && (std::abs(muon1->muonBestTrack()->dxy(vertex.position())) < 0.02)){                                                                                                                  
					Muon_isMedPr[nMuon] = true;                                                                                                              
				}                                                                                                                                          
			}                                                                                                                                      
			Muon_isGoodGL[nMuon] = (muon1->isGlobalMuon() && muon1->globalTrack()->normalizedChi2() < 3 && muon1->combinedQuality().chi2LocalPosition < 12 && muon1->combinedQuality().trkKink < 20 && (muon::segmentCompatibility(*muon1)) > 0.303);                     
			Muon_isTight[nMuon] = (muon::isTightMuon(*muon1,vertex));                                                                                    
			Muon_isHighPt[nMuon] = (muon::isHighPtMuon(*muon1,vertex));                                                                                  
			Muon_isHighPttrk[nMuon] = (muon::isTrackerHighPtMuon(*muon1,vertex));   
			
			// Displacement //
			
			Muon_dxy[nMuon] = muon1->muonBestTrack()->dxy(vertex.position());                                                                         
			Muon_dz[nMuon] = muon1->muonBestTrack()->dz(vertex.position());  
			Muon_dxyErr[nMuon] = muon1->edB(pat::Muon::PV2D);   
			Muon_ip3d[nMuon] =  muon1->dB(pat::Muon::PV3D)/muon1->edB(pat::Muon::PV3D);     
			 
			// energy & track info //
			                                                                     
			Muon_ecal[nMuon] = (muon1->calEnergy()).em;                                                                                                  
			Muon_hcal[nMuon] = (muon1->calEnergy()).had;                                                                         
			Muon_posmatch[nMuon] = muon1->combinedQuality().chi2LocalPosition;                                                                           
			Muon_trkink[nMuon] = muon1->combinedQuality().trkKink;                                                                                       
			Muon_segcom[nMuon] = muon::segmentCompatibility(*muon1);                                                                                     
			  
			// Displacement w.r.t secondary vertex //
			                                                                    
			float dzmumin = 1000;                                                                                                                        
			float dxymumin = 1000;                                                                                                                       
			if(secondaryVertices.isValid()){                                                                                                                          
				for(unsigned int isv=0; isv<(secondaryVertices->size()); isv++){                                                                                        
					const auto &sv = (*secondaryVertices)[isv];                                                                                                           
					reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
					float dztmp = fabs(muon1->muonBestTrack()->dz(svpoint));
					if(dztmp < dzmumin){
						dzmumin = dztmp;                                                                                   
						dxymumin = muon1->muonBestTrack()->dxy(svpoint);                                                                                       
					}                                                                                                                                        
				}                                                                                                                                          
			}                                                                                                                                            
			Muon_dxy_sv[nMuon] = dxymumin; 
			
			// Track info //
			                                                                                                                                                                                                                    
			TrackRef trkglb =muon1->globalTrack();                                                                                                       
			if ((!muon1->isGlobalMuon())) {                                                                                                              
				if (muon1->isTrackerMuon()) {                                                                                                              
					trkglb =muon1->innerTrack();                                                                                                             
				} else {                                                                                                                                   
					trkglb =muon1->outerTrack();                                                                                                             
				}                                                                                                                                          
			}
			                                                                                                                                            
			Muon_chi[nMuon] = trkglb->normalizedChi2();                                                                                                  
			Muon_ndf[nMuon] = (int)trkglb->ndof();                                                                                                       
			Muon_hit[nMuon] = trkglb->hitPattern().numberOfValidMuonHits();                                                                              
			Muon_mst[nMuon] = muon1->numberOfMatchedStations();                                                                                          
			Muon_pixhit[nMuon] = trktrk->hitPattern().numberOfValidPixelHits();                                                                          
			Muon_trklay[nMuon] = trktrk->hitPattern().trackerLayersWithMeasurement();                                                                    
			Muon_valfrac[nMuon] = trktrk->validFraction();                                                        
			Muon_pfiso[nMuon] = (muon1->pfIsolationR04().sumChargedHadronPt + max(0., muon1->pfIsolationR04().sumNeutralHadronEt + muon1->pfIsolationR04().sumPhotonEt - 0.5*muon1->pfIsolationR04().sumPUPt))/muon1->pt();                                               
			
			bool mu_id = Muon_Tight_ID(Muon_isGL[nMuon],Muon_isPF[nMuon],
										   Muon_chi[nMuon],Muon_hit[nMuon],Muon_mst[nMuon],
										   Muon_dxy[nMuon],Muon_dz[nMuon],Muon_pixhit[nMuon],Muon_trklay[nMuon]);
		    Muon_TightID[nMuon] = mu_id;
		    
			if (Muon_pt[nMuon]>15 && fabs(Muon_eta[nMuon])<2.5 && mu_id && Muon_dxy[nMuon]<0.2 && Muon_dz[nMuon]<0.5) {
				tlvmu.push_back(*muon1);
			}
			
			// Application of Rochester correction //
			
			float rcSF, rcSF_error;
			
			if(!isMC){
				// Data
				rcSF = roch_cor.kScaleDT(muon1->charge(), Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon]); 
				rcSF_error = roch_cor.kScaleDTerror(muon1->charge(), Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon]); 
			}
			else{
				// MC
				bool gen_match = false;
				float match_genpt = -100;
				float dR_cut = 0.1;
				for(int ipart=0; ipart<nGenPart; ipart++)
				{
					if((GenPart_status[ipart]==1) && (GenPart_pdg[ipart]==(-1*muon1->charge()*13)) && (delta2R(GenPart_eta[ipart],GenPart_phi[ipart],Muon_eta[nMuon], Muon_phi[nMuon])<dR_cut))
					{
						dR_cut = delta2R(GenPart_eta[ipart],GenPart_phi[ipart],Muon_eta[nMuon], Muon_phi[nMuon]);
						gen_match = true;
						match_genpt = GenPart_pt[ipart];
					}
				}
				if(gen_match){
					// when matched gen muon is available
					rcSF = roch_cor.kSpreadMC(muon1->charge(), Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon], match_genpt); 
					rcSF_error = roch_cor.kSpreadMCerror(muon1->charge(), Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon], match_genpt);
				} 
				else{
					// when matched gen muon is not available
					rcSF = roch_cor.kSmearMC(muon1->charge(), Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon], Muon_trklay[nMuon], gRandom->Rndm()); 
					rcSF_error = roch_cor.kSmearMCerror(muon1->charge(), Muon_pt[nMuon], Muon_eta[nMuon], Muon_phi[nMuon], Muon_trklay[nMuon], gRandom->Rndm());
				}
			}
						
			Muon_corrected_pt[nMuon] = Muon_pt[nMuon]*rcSF;
			Muon_correctedUp_pt[nMuon] = Muon_pt[nMuon]*max(rcSF+rcSF_error,float(0.));
			Muon_correctedDown_pt[nMuon] = Muon_pt[nMuon]*max(rcSF-rcSF_error,float(0.));
			
			// End of Rochester correction //
			
			if (++nMuon>=njetmx) break;                                                                                                                 
		
		}                                                                                                                                              
      }                                                                                                                                               
  }// muon loop 
  
  // Electrons //
      
  nElectron = 0;             
  std::vector<pat::Electron> tlvel;
    
  for(const auto& electron1 : iEvent.get(tok_electrons_) ) {                                                                                          
 
    if (!StoreElectron(electron1,minePt,maxEta)) continue;
                                                   
	GsfTrackRef gsftrk1 = electron1.gsfTrack();   																														
    TrackRef ctftrk = electron1.closestCtfTrackRef();    
    
    Electron_pt[nElectron] = electron1.pt();   	                                                                                 
    Electron_eta[nElectron] = electron1.eta();                                                                                                                 
    Electron_phi[nElectron] = electron1.phi();                                                                                                                 
    Electron_e[nElectron] = electron1.energy(); 
    Electron_e_ECAL[nElectron] = electron1.ecalEnergy();                                                                                                 
    Electron_p[nElectron] = electron1.trackMomentumAtVtx().R()*electron1.charge();     
    
    // MVA id //
    
    Electron_mvaid_Fallv2WP90[nElectron] = electron1.electronID(melectronID_isowp90);                                                                                 
    Electron_mvaid_Fallv2WP90_noIso[nElectron] = electron1.electronID(melectronID_noisowp90);                                                                             
    Electron_mvaid_Fallv2WP80[nElectron] = electron1.electronID(melectronID_isowp80);                                                                                 
    Electron_mvaid_Fallv2WP80_noIso[nElectron] = electron1.electronID(melectronID_noisowp80);   
    
    // displacement //
                                                                                 
    Electron_dxy[nElectron] = gsftrk1->dxy(vertex.position());  
    Electron_dxyErr[nElectron] = electron1.edB(pat::Electron::PV2D);                                                                                           
    Electron_dz[nElectron] = gsftrk1->dz(vertex.position()); 
    Electron_dzErr[nElectron] = electron1.edB(pat::Electron::PVDZ);     
    Electron_ip3d[nElectron] =  electron1.dB(pat::Electron::PV3D)/electron1.edB(pat::Electron::PV3D);                                                                                           
    
    // supercluste info //
    
    Electron_supcl_e[nElectron] = electron1.superCluster()->energy();  
    Electron_supcl_rawE[nElectron] = electron1.superCluster()->rawEnergy();      
    Electron_supcl_eta[nElectron] = electron1.superCluster()->eta();                                                                                           
    Electron_supcl_phi[nElectron] = electron1.superCluster()->phi();                                                                                           
    
    // scaling & smearing factors //
    
    if(store_electron_scalnsmear){
		Electron_eccalTrkEnergyPostCorr[nElectron] = electron1.userFloat("ecalTrkEnergyPostCorr");
		Electron_energyScaleValue[nElectron] = electron1.userFloat("energyScaleValue");
		Electron_energySigmaValue[nElectron] = electron1.userFloat("energySigmaValue");
		Electron_energyScaleUp[nElectron] = electron1.userFloat("energyScaleUp");
		Electron_energyScaleDown[nElectron] = electron1.userFloat("energyScaleDown");
		Electron_energySigmaUp[nElectron] = electron1.userFloat("energySigmaUp");
		Electron_energySigmaDown[nElectron] = electron1.userFloat("energySigmaDown");
		Electron_energyScaleStatUp[nElectron] = electron1.userFloat("energyScaleStatUp");
		Electron_energyScaleStatDown[nElectron] = electron1.userFloat("energyScaleStatDown");
		Electron_energyScaleSystUp[nElectron] = electron1.userFloat("energyScaleSystUp");
		Electron_energyScaleSystDown[nElectron] = electron1.userFloat("energyScaleSystDown");
		Electron_energyScaleGainUp[nElectron] = electron1.userFloat("energyScaleGainUp");
		Electron_energyScaleGainDown[nElectron] = electron1.userFloat("energyScaleGainDown");
		Electron_energySigmaRhoUp[nElectron] = electron1.userFloat("energySigmaRhoUp");
		Electron_energySigmaRhoDown[nElectron] = electron1.userFloat("energySigmaRhoDown");
		Electron_energySigmaPhiUp[nElectron] = electron1.userFloat("energySigmaPhiUp");
		Electron_energySigmaPhiDown[nElectron] = electron1.userFloat("energySigmaPhiDown");
		Electron_energySmearNrSigma[nElectron] = electron1.userFloat("energySmearNrSigma");
		Electron_ecalEnergyPreCorr[nElectron] = electron1.userFloat("ecalEnergyPreCorr");
		Electron_ecalEnergyErrPreCorr[nElectron] = electron1.userFloat("ecalEnergyErrPreCorr");
		Electron_ecalEnergyPostCorr[nElectron] = electron1.userFloat("ecalEnergyPostCorr");
		Electron_ecalEnergyErrPostCorr[nElectron] = electron1.userFloat("ecalEnergyErrPostCorr");

	}
    // end of scaling & smearing factors //                                                                                                
                 
    // shape of energy deposition //             
                                                                                                                                                                                                                                              
	Electron_sigmaieta[nElectron] = electron1.full5x5_sigmaIetaIeta();                                                                                         
    Electron_sigmaiphi[nElectron] = electron1.full5x5_sigmaIphiIphi();                                                                                         
    Electron_r9full[nElectron] = electron1.full5x5_r9();                                                                                                                                                                                        
    Electron_hcaloverecal[nElectron] = electron1.full5x5_hcalOverEcal();                                                                                                                                                                                                                                            
    Electron_hitsmiss[nElectron] =  electron1.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);   
    Electron_eoverp[nElectron] = electron1.eSuperClusterOverP();                                                                                               
    Electron_hovere[nElectron] = electron1.hadronicOverEm();                                                                                                   
    Electron_ietaieta[nElectron] = electron1.sigmaIetaIeta();                                                                                                  
    Electron_ecloverpout[nElectron] = electron1.eEleClusterOverPout();  
    Electron_convVeto[nElectron] = electron1.passConversionVeto();
    
    if(store_electron_addvariabs){
		
		Electron_etain[nElectron] = electron1.deltaEtaSuperClusterTrackAtVtx();                                                                                    
		Electron_phiin[nElectron] = electron1.deltaPhiSuperClusterTrackAtVtx();                                                                                    
		Electron_fbrem[nElectron] = electron1.fbrem();   
		Electron_supcl_preshvsrawe[nElectron] = electron1.superCluster()->preshowerEnergy()/electron1.superCluster()->rawEnergy();
		Electron_cloctftrkn[nElectron] = electron1.closestCtfTrackNLayers();                                                                                       
		Electron_cloctftrkchi2[nElectron] = electron1.closestCtfTrackNormChi2();   
		Electron_e1x5bye5x5[nElectron] = 1.-electron1.full5x5_e1x5()/electron1.full5x5_e5x5();                                                                     
		Electron_normchi2[nElectron] =  electron1.gsfTrack()->normalizedChi2(); 
		Electron_supcl_etaw[nElectron] = electron1.superCluster()->etaWidth();                                                                                     
		Electron_supcl_phiw[nElectron] = electron1.superCluster()->phiWidth(); 
		Electron_trkmeasure[nElectron] = electron1.gsfTrack()->hitPattern().trackerLayersWithMeasurement();  
		Electron_convtxprob[nElectron] = electron1.convVtxFitProb();                                                                                                                                                   
		Electron_ecaletrkmomentum[nElectron] = 1.0/(electron1.ecalEnergy())-1.0/(electron1.trackMomentumAtVtx().R());                                              
		Electron_deltaetacltrkcalo[nElectron] = electron1.deltaEtaSeedClusterTrackAtCalo();                                                                        
                                           	
	}
    
    // isolation variables //                              
                                                          
    Electron_pfisolsumphet[nElectron] = electron1.pfIsolationVariables().sumPhotonEt;                                                                          
    Electron_pfisolsumchhadpt[nElectron] = electron1.pfIsolationVariables().sumChargedHadronPt;                                                                
    Electron_pfsiolsumneuhadet[nElectron] = electron1.pfIsolationVariables().sumNeutralHadronEt;                                                               
	                                                                                      
    const reco::GsfElectron::PflowIsolationVariables& pfIso = electron1.pfIsolationVariables();                                                      
    Electron_pfiso_drcor[nElectron] = (pfIso.sumChargedHadronPt + max(0., pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5*pfIso.sumPUPt))*1./electron1.pt();      
    vector<float> pfisovalues;                                                                                     
    Read_ElePFIsolation(&electron1,Rho,pfisovalues);
    Electron_pfiso_eacor[nElectron] = pfisovalues[0];
    Electron_pfiso04_eacor[nElectron] = pfisovalues[1];
    
    // Displacement w.r.t secondary vertex //
    
    float dzmin = 1000;                                                                                                                              
    float dxymin = 1000;
    if(secondaryVertices.isValid()){                                                                                                                              
		for(unsigned int isv=0; isv<(secondaryVertices->size()); isv++){                                                                                            
		const auto &sv = (*secondaryVertices)[isv];                                                                                                               
     	  reco::TrackBase::Point svpoint(sv.vx(),sv.vy(),sv.vz());
		  float dztmp =fabs(gsftrk1->dz(svpoint));
		  if(dztmp < dzmin){                                                                                                      
			dzmin = dztmp;                                                                                                        
			dxymin = gsftrk1->dxy(svpoint);    
			}                                                                                                                                            
		}                                                                                                                                              
    }     
                                                                                                                                                 
    Electron_dxy_sv[nElectron] = dxymin;  
     
    // track info //
    Electron_chi[nElectron] = gsftrk1->chi2();                                                                                                                 
    Electron_ndf[nElectron] = (int)gsftrk1->ndof();                                                                                                            
    Electron_misshits[nElectron] = (int)gsftrk1->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
    
    //MiniIsolation: begin//                                                                                      
	vector<float> isovalues;
	Read_MiniIsolation(&electron1,Rho,isovalues);
	Electron_minisoall[nElectron] = isovalues[0];
	Electron_minchiso[nElectron] = isovalues[1];
	Electron_minnhiso[nElectron] = isovalues[2];
	Electron_minphiso[nElectron] = isovalues[3];
	//MiniIsolation: end//  
	
	bool impact_pass = 	((fabs(Electron_supcl_eta[nElectron])<1.4442 && fabs(Electron_dxy[nElectron])<0.05 && fabs(Electron_dz[nElectron])<0.1)
					   ||(fabs(Electron_supcl_eta[nElectron])>1.5660 && fabs(Electron_dxy[nElectron])<(2*0.05) && fabs(Electron_dz[nElectron])<(2*0.1)));

	
	if(Electron_pt[nElectron]>15 && fabs(Electron_eta[nElectron])<2.5 && Electron_mvaid_Fallv2WP90_noIso[nElectron] && impact_pass){
		tlvel.push_back(electron1);
	}
  
    if(++nElectron>=njetmx) break;                                                                                                                      
  }
  
  // AK8 jets //
  
  nPFJetAK8 = 0;
  nPFJetAK8_cons = 0;
  
  edm::Handle<edm::View<pat::Jet>> pfjetAK8s;
  iEvent.getByToken(tok_pfjetAK8s_, pfjetAK8s);	
  
  if(pfjetAK8s.isValid()){
    
    for (unsigned jet = 0; jet< pfjetAK8s->size(); jet++) {
      
      const auto &ak8jet = (*pfjetAK8s)[jet];

      TLorentzVector pfjetAK8_4v(ak8jet.correctedP4("Uncorrected").px(),ak8jet.correctedP4("Uncorrected").py(),ak8jet.correctedP4("Uncorrected").pz(), ak8jet.correctedP4("Uncorrected").energy());
     
	  if(subtractLepton_fromAK8){
		pfjetAK8_4v = LeptonJet_subtraction(tlvmu,ak8jet,pfjetAK8_4v);
		pfjetAK8_4v = LeptonJet_subtraction(tlvel,ak8jet,pfjetAK8_4v);
	  }
	  
      double tmprecpt = pfjetAK8_4v.Pt();
      
      double total_cor =1;
      Read_JEC(total_cor,tmprecpt,pfjetAK8_4v.Eta(),Rho,isData,ak8jet,jecL1FastAK8,jecL2RelativeAK8,jecL3AbsoluteAK8,jecL2L3ResidualAK8);
      PFJetAK8_JEC[nPFJetAK8] = total_cor;
            
      if(tmprecpt<AK8PtCut) continue;
      if(abs(pfjetAK8_4v.Eta())>maxEta) continue;
      
      PFJetAK8_pt[nPFJetAK8] = 	pfjetAK8_4v.Pt();
      PFJetAK8_y[nPFJetAK8] = pfjetAK8_4v.Rapidity();
      PFJetAK8_eta[nPFJetAK8] = pfjetAK8_4v.Eta();
      PFJetAK8_phi[nPFJetAK8] = pfjetAK8_4v.Phi();
      PFJetAK8_mass[nPFJetAK8] = ak8jet.correctedP4("Uncorrected").mass();
      PFJetAK8_btag_DeepCSV[nPFJetAK8] = ak8jet.bDiscriminator("pfDeepCSVJetTags:probb")+ak8jet.bDiscriminator("pfDeepCSVJetTags:probbb");
      
      // DNN-based tagger info //
      
      PFJetAK8_DeepTag_DAK8_TvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(toptagger_DAK8);
      PFJetAK8_DeepTag_DAK8_WvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Wtagger_DAK8);
      PFJetAK8_DeepTag_DAK8_ZvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Ztagger_DAK8);
      PFJetAK8_DeepTag_DAK8_HvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Htagger_DAK8);
      PFJetAK8_DeepTag_DAK8_bbvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(bbtagger_DAK8);
      
      PFJetAK8_DeepTag_PNet_TvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(toptagger_PNet);
      PFJetAK8_DeepTag_PNet_WvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Wtagger_PNet);
      PFJetAK8_DeepTag_PNet_ZvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Ztagger_PNet);
      PFJetAK8_DeepTag_PNet_XbbvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Xbbtagger_PNet);
      PFJetAK8_DeepTag_PNet_XccvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Xcctagger_PNet);
      PFJetAK8_DeepTag_PNet_XqqvsQCD[nPFJetAK8] = ak8jet.bDiscriminator(Xqqtagger_PNet);
      PFJetAK8_DeepTag_PNet_QCD[nPFJetAK8] = (ak8jet.bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDbb")+ak8jet.bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDcc")+
										    ak8jet.bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDb")+ak8jet.bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDc")+
										    ak8jet.bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDothers"));
      
      //  JER //
       
      if(isMC){
	
		TLorentzVector tmp4v;
		tmp4v.SetPtEtaPhiM(PFJetAK8_JEC[nPFJetAK8]*pfjetAK8_4v.Pt(),pfjetAK8_4v.Eta(),pfjetAK8_4v.Phi(),PFJetAK8_JEC[nPFJetAK8]*pfjetAK8_4v.M());
		
      	vector<double> SFs;
      	Read_JER(mPtResoFileAK8, mPtSFFileAK8, tmprecpt, tmp4v, Rho, genjetAK8s, 0.5*0.8, SFs);
      	
      	PFJetAK8_reso[nPFJetAK8] = SFs[0];
      	PFJetAK8_resoup[nPFJetAK8] = SFs[1];
      	PFJetAK8_resodn[nPFJetAK8] = SFs[2];
      	      	
      }//isMC
      
      // JES uncertainty //
      
      for(int isrc =0 ; isrc<njecmcmx; isrc++){
	
        double sup = 1.0 ;
	
		if((isrc>0)&&(isrc<=nsrc)){
		  
		  JetCorrectionUncertainty *jecUnc = vsrcAK8[isrc-1];
		  jecUnc->setJetEta(ak8jet.eta());
		  jecUnc->setJetPt(tmprecpt);
		  
		  sup += jecUnc->getUncertainty(true);         
		  if(isrc==1){ PFJetAK8_jesup_AbsoluteStat[nPFJetAK8] = sup; }
		  if(isrc==2){ PFJetAK8_jesup_AbsoluteScale[nPFJetAK8] = sup; }
		  if(isrc==3){ PFJetAK8_jesup_AbsoluteMPFBias[nPFJetAK8] = sup; }
		  if(isrc==4){ PFJetAK8_jesup_FlavorQCD[nPFJetAK8] = sup; }
		  if(isrc==5){ PFJetAK8_jesup_Fragmentation[nPFJetAK8] = sup; }
		  if(isrc==6){ PFJetAK8_jesup_PileUpDataMC[nPFJetAK8] = sup; }
		  if(isrc==7){ PFJetAK8_jesup_PileUpPtBB[nPFJetAK8] = sup; }
		  if(isrc==8){ PFJetAK8_jesup_PileUpPtEC1[nPFJetAK8] = sup; }
		  if(isrc==9){ PFJetAK8_jesup_PileUpPtEC2[nPFJetAK8] = sup; }
		  if(isrc==10){ PFJetAK8_jesup_PileUpPtRef[nPFJetAK8] = sup; }
		  if(isrc==11){ PFJetAK8_jesup_RelativeFSR[nPFJetAK8] = sup; }
		  if(isrc==12){ PFJetAK8_jesup_RelativeJEREC1[nPFJetAK8] = sup; }
		  if(isrc==13){ PFJetAK8_jesup_RelativeJEREC2[nPFJetAK8] = sup; }
		  if(isrc==14){ PFJetAK8_jesup_RelativePtBB[nPFJetAK8] = sup; }
		  if(isrc==15){ PFJetAK8_jesup_RelativePtEC1[nPFJetAK8] = sup; }
		  if(isrc==16){ PFJetAK8_jesup_RelativePtEC2[nPFJetAK8] = sup; }
		  if(isrc==17){ PFJetAK8_jesup_RelativeBal[nPFJetAK8] = sup; }
		  if(isrc==18){ PFJetAK8_jesup_RelativeSample[nPFJetAK8] = sup; }
		  if(isrc==19){ PFJetAK8_jesup_RelativeStatEC[nPFJetAK8] = sup; }
		  if(isrc==20){ PFJetAK8_jesup_RelativeStatFSR[nPFJetAK8] = sup; }
		  if(isrc==21){ PFJetAK8_jesup_SinglePionECAL[nPFJetAK8] = sup; }
		  if(isrc==22){ PFJetAK8_jesup_SinglePionHCAL[nPFJetAK8] = sup; }
		  if(isrc==23){ PFJetAK8_jesup_TimePtEta[nPFJetAK8] = sup; }
		  if(isrc==24){ PFJetAK8_jesup_Total[nPFJetAK8] = sup; }
		
		}
		else if(isrc>nsrc){
		  
		  JetCorrectionUncertainty *jecUnc = vsrcAK8[isrc-1-nsrc];
		  jecUnc->setJetEta(ak8jet.eta());
		  jecUnc->setJetPt(tmprecpt);
		  
		  sup -= jecUnc->getUncertainty(false);
		  if(isrc==(nsrc+1)){ PFJetAK8_jesdn_AbsoluteStat[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+2)){ PFJetAK8_jesdn_AbsoluteScale[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+3)){ PFJetAK8_jesdn_AbsoluteMPFBias[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+4)){ PFJetAK8_jesdn_FlavorQCD[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+5)){ PFJetAK8_jesdn_Fragmentation[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+6)){ PFJetAK8_jesdn_PileUpDataMC[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+7)){ PFJetAK8_jesdn_PileUpPtBB[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+8)){ PFJetAK8_jesdn_PileUpPtEC1[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+9)){ PFJetAK8_jesdn_PileUpPtEC2[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+10)){ PFJetAK8_jesdn_PileUpPtRef[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+11)){ PFJetAK8_jesdn_RelativeFSR[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+12)){ PFJetAK8_jesdn_RelativeJEREC1[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+13)){ PFJetAK8_jesdn_RelativeJEREC2[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+14)){ PFJetAK8_jesdn_RelativePtBB[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+15)){ PFJetAK8_jesdn_RelativePtEC1[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+16)){ PFJetAK8_jesdn_RelativePtEC2[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+17)){ PFJetAK8_jesdn_RelativeBal[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+18)){ PFJetAK8_jesdn_RelativeSample[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+19)){ PFJetAK8_jesdn_RelativeStatEC[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+20)){ PFJetAK8_jesdn_RelativeStatFSR[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+21)){ PFJetAK8_jesdn_SinglePionECAL[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+22)){ PFJetAK8_jesdn_SinglePionHCAL[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+23)){ PFJetAK8_jesdn_TimePtEta[nPFJetAK8] = sup; }
		  if(isrc==(nsrc+24)){ PFJetAK8_jesdn_Total[nPFJetAK8] = sup; }
		}
		
      }
      
      // Jet id //
       
      PFJetAK8_CHF[nPFJetAK8] = ak8jet.chargedHadronEnergyFraction();
      PFJetAK8_NHF[nPFJetAK8] = ak8jet.neutralHadronEnergyFraction();
      PFJetAK8_CEMF[nPFJetAK8] = ak8jet.chargedEmEnergyFraction();
      PFJetAK8_NEMF[nPFJetAK8] = ak8jet.neutralEmEnergyFraction();
      PFJetAK8_MUF[nPFJetAK8] = ak8jet.muonEnergyFraction();
      PFJetAK8_PHF[nPFJetAK8] = ak8jet.photonEnergyFraction();
      PFJetAK8_EEF[nPFJetAK8] = ak8jet.electronEnergyFraction();
      PFJetAK8_HFHF[nPFJetAK8] = ak8jet.HFHadronEnergyFraction();
      
      PFJetAK8_CHM[nPFJetAK8] = ak8jet.chargedHadronMultiplicity();
      PFJetAK8_NHM[nPFJetAK8] = ak8jet.neutralHadronMultiplicity();
      PFJetAK8_MUM[nPFJetAK8] = ak8jet.muonMultiplicity();
      PFJetAK8_PHM[nPFJetAK8] = ak8jet.photonMultiplicity();
      PFJetAK8_EEM[nPFJetAK8] = ak8jet.electronMultiplicity();
      PFJetAK8_HFHM[nPFJetAK8] = ak8jet.HFHadronMultiplicity();
      
      PFJetAK8_Chcons[nPFJetAK8] = ak8jet.chargedMultiplicity();
      PFJetAK8_Neucons[nPFJetAK8] = ak8jet.neutralMultiplicity();
      
      JetIDVars idvars; 
      idvars.NHF = PFJetAK8_NHF[nPFJetAK8];
      idvars.NEMF = PFJetAK8_NEMF[nPFJetAK8];
      idvars.MUF = PFJetAK8_MUF[nPFJetAK8];
      idvars.CHF = PFJetAK8_CHF[nPFJetAK8];
      idvars.CEMF = PFJetAK8_CEMF[nPFJetAK8];
      idvars.NumConst = (PFJetAK8_Chcons[nPFJetAK8]+PFJetAK8_Neucons[nPFJetAK8]);
      idvars.NumNeutralParticle = PFJetAK8_Neucons[nPFJetAK8];
      idvars.CHM = PFJetAK8_CHM[nPFJetAK8];
      
      PFJetAK8_jetID[nPFJetAK8] = getJetID(idvars,"PUPPI",year,PFJetAK8_eta[nPFJetAK8],false,isUltraLegacy);
      PFJetAK8_jetID_tightlepveto[nPFJetAK8] = getJetID(idvars,"PUPPI",year,PFJetAK8_eta[nPFJetAK8],true,isUltraLegacy);  
     
      // subjet info & classical observables //
     
      PFJetAK8_sub1pt[nPFJetAK8] = PFJetAK8_sub1eta[nPFJetAK8] = PFJetAK8_sub1phi[nPFJetAK8] = PFJetAK8_sub1mass[nPFJetAK8] = PFJetAK8_sub1btag[nPFJetAK8] = -100;              
      PFJetAK8_sub2pt[nPFJetAK8] = PFJetAK8_sub2eta[nPFJetAK8] = PFJetAK8_sub2phi[nPFJetAK8] = PFJetAK8_sub2mass[nPFJetAK8] = PFJetAK8_sub2btag[nPFJetAK8] = -100;                                                        
      PFJetAK8_sdmass[nPFJetAK8] = PFJetAK8_tau1[nPFJetAK8] = PFJetAK8_tau2[nPFJetAK8] = PFJetAK8_tau3[nPFJetAK8] = -100;                                                                      
      
      if(isSoftDrop){
	
		PFJetAK8_tau1[nPFJetAK8] = ak8jet.userFloat(Nsubjettiness_tau1);
		PFJetAK8_tau2[nPFJetAK8] = ak8jet.userFloat(Nsubjettiness_tau2);
		PFJetAK8_tau3[nPFJetAK8] = ak8jet.userFloat(Nsubjettiness_tau3);
		
		PFJetAK8_sdmass[nPFJetAK8] = (ak8jet.groomedMass(subjets) > 0)? ak8jet.groomedMass(subjets) : 0;
		
		for(unsigned int isub=0; isub<((ak8jet.subjets(subjets)).size()); isub++){
		
			const auto ak8subjet = (ak8jet.subjets(subjets))[isub];
	    
			if(isub==0){
				PFJetAK8_sub1pt[nPFJetAK8] = ak8subjet->correctedP4("Uncorrected").pt();
				PFJetAK8_sub1eta[nPFJetAK8] = ak8subjet->eta();
				PFJetAK8_sub1phi[nPFJetAK8] = ak8subjet->phi();
				PFJetAK8_sub1mass[nPFJetAK8] = ak8subjet->correctedP4("Uncorrected").mass();	 
				PFJetAK8_sub1JEC[nPFJetAK8] = ak8subjet->pt()*1./ak8subjet->correctedP4("Uncorrected").pt();
				PFJetAK8_sub1btag[nPFJetAK8] = ak8subjet->bDiscriminator("pfDeepCSVJetTags:probb")+ak8subjet->bDiscriminator("pfDeepCSVJetTags:probbb");
			}
			else if(isub==1){
				PFJetAK8_sub2pt[nPFJetAK8] = ak8subjet->correctedP4("Uncorrected").pt();
				PFJetAK8_sub2eta[nPFJetAK8] = ak8subjet->eta();
				PFJetAK8_sub2phi[nPFJetAK8] = ak8subjet->phi();
				PFJetAK8_sub2mass[nPFJetAK8] = ak8subjet->correctedP4("Uncorrected").mass();	
				PFJetAK8_sub2JEC[nPFJetAK8] = ak8subjet->pt()*1./ak8subjet->correctedP4("Uncorrected").pt(); 
				PFJetAK8_sub2btag[nPFJetAK8] = ak8subjet->bDiscriminator("pfDeepCSVJetTags:probb")+ak8subjet->bDiscriminator("pfDeepCSVJetTags:probbb");
			}
		}
		
	
      }//isSoftDrop
      
      // Storing 4-momenta of jet constituents//
      if(store_fatjet_constituents && nGenJetAK8<njetconsmax)     {
        for(unsigned int ic = 0 ; ic < ak8jet.numberOfSourceCandidatePtrs() ; ++ic) {  
          if(ak8jet.sourceCandidatePtr(ic).isNonnull() && ak8jet.sourceCandidatePtr(ic).isAvailable()){
            
            if(nPFJetAK8_cons>= nconsmax) break;	
            const reco::Candidate* jcand = ak8jet.sourceCandidatePtr(ic).get();
            PFJetAK8_cons_pt[nPFJetAK8_cons] = jcand->pt();
            PFJetAK8_cons_eta[nPFJetAK8_cons] = jcand->eta();
            PFJetAK8_cons_phi[nPFJetAK8_cons] = jcand->phi();
            PFJetAK8_cons_mass[nPFJetAK8_cons] = jcand->mass();
            PFJetAK8_cons_pdgId[nPFJetAK8_cons] = jcand->pdgId();
            PFJetAK8_cons_jetIndex[nPFJetAK8_cons] = nPFJetAK8;   
            nPFJetAK8_cons++;
          }
        }
      }
            
      // end of candidate storage //
            
      nPFJetAK8++;	
      if(nPFJetAK8 >= njetmxAK8) { break;}
      
    }
  }
  
  nPFJetAK4 = 0;
  edm::Handle<edm::View<pat::Jet>> pfjetAK4s;
  edm::Handle<edm::View<pat::Jet>> pfjetAK4sB;
  iEvent.getByToken(tok_pfjetAK4s_, pfjetAK4s);
  iEvent.getByToken(tok_pfjetAK4sB_, pfjetAK4sB);


  for (unsigned jet = 0; jet< pfjetAK4s->size(); jet++) {
      
	const auto &ak4jet = (*pfjetAK4s)[jet];
    TLorentzVector pfjetAK4_4v(ak4jet.correctedP4("Uncorrected").px(),ak4jet.correctedP4("Uncorrected").py(),ak4jet.correctedP4("Uncorrected").pz(), ak4jet.correctedP4("Uncorrected").energy());
    
    if(subtractLepton_fromAK4){
		pfjetAK4_4v = LeptonJet_subtraction(tlvmu,ak4jet,pfjetAK4_4v);
		pfjetAK4_4v = LeptonJet_subtraction(tlvel,ak4jet,pfjetAK4_4v);
	}
    
    double tmprecpt = pfjetAK4_4v.Pt();
   
    double total_cor =1;
    Read_JEC(total_cor,tmprecpt,pfjetAK4_4v.Eta(),Rho,isData,ak4jet,jecL1FastAK4,jecL2RelativeAK4,jecL3AbsoluteAK4,jecL2L3ResidualAK4);  
    PFJetAK4_JEC[nPFJetAK4] = total_cor;
    
    tmprecpt = pfjetAK4_4v.Pt();
    if(tmprecpt<minjPt) continue;
    if(abs(pfjetAK4_4v.Rapidity())>maxEta) continue;
      
    PFJetAK4_pt[nPFJetAK4] = 	tmprecpt;
    PFJetAK4_eta[nPFJetAK4] = 	pfjetAK4_4v.Eta();
    PFJetAK4_y[nPFJetAK4] = pfjetAK4_4v.Rapidity();
    PFJetAK4_phi[nPFJetAK4] = pfjetAK4_4v.Phi();
    PFJetAK4_mass[nPFJetAK4] = pfjetAK4_4v.M(); 
    PFJetAK4_energy[nPFJetAK4] = pfjetAK4_4v.E();
   
    for (unsigned jetB = 0; jetB< pfjetAK4sB->size(); jetB++) { 

	 const auto &ak4jetB = (*pfjetAK4sB)[jetB];
         TLorentzVector pfjetAK4_4vB(ak4jetB.correctedP4("Uncorrected").px(),ak4jetB.correctedP4("Uncorrected").py(),ak4jetB.correctedP4("Uncorrected").pz(), ak4jetB.correctedP4("Uncorrected").energy());
         if (pfjetAK4_4v.Eta() == pfjetAK4_4vB.Eta()) {
         PFJetAK4_bcorr[nPFJetAK4] =  ak4jetB.userFloat("bJetRegCorr");
         PFJetAK4_breso[nPFJetAK4] =  ak4jetB.userFloat("bJetRegRes");
         break;
        }
    }

   TLorentzVector AK4j_4v;
   AK4j_4v.SetPtEtaPhiE(PFJetAK4_pt[nPFJetAK4], PFJetAK4_eta[nPFJetAK4], PFJetAK4_phi[nPFJetAK4], PFJetAK4_energy[nPFJetAK4]);
   TLorentzVector AK4j_4v_reg;
   AK4j_4v_reg.SetPtEtaPhiE(PFJetAK4_pt[nPFJetAK4]*PFJetAK4_bcorr[nPFJetAK4], PFJetAK4_eta[nPFJetAK4], PFJetAK4_phi[nPFJetAK4], PFJetAK4_energy[nPFJetAK4]*PFJetAK4_bcorr[nPFJetAK4]);

   PFJetAK4_pt_reg[nPFJetAK4] = AK4j_4v_reg.Pt();
   PFJetAK4_y_reg[nPFJetAK4] = AK4j_4v_reg.Rapidity();
   PFJetAK4_eta_reg[nPFJetAK4] = AK4j_4v_reg.Eta();
   PFJetAK4_phi_reg[nPFJetAK4] = AK4j_4v_reg.Phi();
   PFJetAK4_mass_reg[nPFJetAK4] = AK4j_4v_reg.M();
   PFJetAK4_energy_reg[nPFJetAK4] = AK4j_4v_reg.E();


    // JER //
     
    if(isMC){
		
		TLorentzVector tmp4v;
		tmp4v.SetPtEtaPhiM(PFJetAK4_JEC[nPFJetAK4]*pfjetAK4_4v.Pt(),pfjetAK4_4v.Eta(),pfjetAK4_4v.Phi(),PFJetAK4_JEC[nPFJetAK4]*pfjetAK4_4v.M());
		
		vector<double> SFs;
      	Read_JER(mPtResoFileAK4, mPtSFFileAK4, tmprecpt, tmp4v, Rho, genjetAK4s, 0.5*0.4, SFs);
      	
      	PFJetAK4_reso[nPFJetAK4] = SFs[0];
      	PFJetAK4_resoup[nPFJetAK4] = SFs[1];
      	PFJetAK4_resodn[nPFJetAK4] = SFs[2];
		
    }//isMC
      
     // JES uncertainty //
      
    for(int isrc =0 ; isrc<njecmcmx; isrc++){
	
		double sup = 1.0 ;
	
		if((isrc>0)&&(isrc<=nsrc)){
	  
			JetCorrectionUncertainty *jecUnc = vsrc[isrc-1];
			jecUnc->setJetEta(ak4jet.eta());
			jecUnc->setJetPt(tmprecpt);
	  
			sup += jecUnc->getUncertainty(true);         
			if(isrc==1){ PFJetAK4_jesup_AbsoluteStat[nPFJetAK4] = sup; }
			if(isrc==2){ PFJetAK4_jesup_AbsoluteScale[nPFJetAK4] = sup; }
			if(isrc==3){ PFJetAK4_jesup_AbsoluteMPFBias[nPFJetAK4] = sup; }
			if(isrc==4){ PFJetAK4_jesup_FlavorQCD[nPFJetAK4] = sup; }
			if(isrc==5){ PFJetAK4_jesup_Fragmentation[nPFJetAK4] = sup; }
			if(isrc==6){ PFJetAK4_jesup_PileUpDataMC[nPFJetAK4] = sup; }
			if(isrc==7){ PFJetAK4_jesup_PileUpPtBB[nPFJetAK4] = sup; }
			if(isrc==8){ PFJetAK4_jesup_PileUpPtEC1[nPFJetAK4] = sup; }
			if(isrc==9){ PFJetAK4_jesup_PileUpPtEC2[nPFJetAK4] = sup; }
			if(isrc==10){ PFJetAK4_jesup_PileUpPtRef[nPFJetAK4] = sup; }
			if(isrc==11){ PFJetAK4_jesup_RelativeFSR[nPFJetAK4] = sup; }
			if(isrc==12){ PFJetAK4_jesup_RelativeJEREC1[nPFJetAK4] = sup; }
			if(isrc==13){ PFJetAK4_jesup_RelativeJEREC2[nPFJetAK4] = sup; }
			if(isrc==14){ PFJetAK4_jesup_RelativePtBB[nPFJetAK4] = sup; }
			if(isrc==15){ PFJetAK4_jesup_RelativePtEC1[nPFJetAK4] = sup; }
			if(isrc==16){ PFJetAK4_jesup_RelativePtEC2[nPFJetAK4] = sup; }
			if(isrc==17){ PFJetAK4_jesup_RelativeBal[nPFJetAK4] = sup; }
			if(isrc==18){ PFJetAK4_jesup_RelativeSample[nPFJetAK4] = sup; }
			if(isrc==19){ PFJetAK4_jesup_RelativeStatEC[nPFJetAK4] = sup; }
			if(isrc==20){ PFJetAK4_jesup_RelativeStatFSR[nPFJetAK4] = sup; }
			if(isrc==21){ PFJetAK4_jesup_SinglePionECAL[nPFJetAK4] = sup; }
			if(isrc==22){ PFJetAK4_jesup_SinglePionHCAL[nPFJetAK4] = sup; }
			if(isrc==23){ PFJetAK4_jesup_TimePtEta[nPFJetAK4] = sup; }
			if(isrc==24){ PFJetAK4_jesup_Total[nPFJetAK4] = sup; }
		}
	
		else if(isrc>nsrc){
	  
			JetCorrectionUncertainty *jecUnc = vsrc[isrc-1-nsrc];
		    jecUnc->setJetEta(ak4jet.eta());
		    jecUnc->setJetPt(tmprecpt);
	  
			sup -= jecUnc->getUncertainty(false);
			if(isrc==(nsrc+1)){ PFJetAK4_jesdn_AbsoluteStat[nPFJetAK4] = sup; }
			if(isrc==(nsrc+2)){ PFJetAK4_jesdn_AbsoluteScale[nPFJetAK4] = sup; }
			if(isrc==(nsrc+3)){ PFJetAK4_jesdn_AbsoluteMPFBias[nPFJetAK4] = sup; }
			if(isrc==(nsrc+4)){ PFJetAK4_jesdn_FlavorQCD[nPFJetAK4] = sup; }
			if(isrc==(nsrc+5)){ PFJetAK4_jesdn_Fragmentation[nPFJetAK4] = sup; }
			if(isrc==(nsrc+6)){ PFJetAK4_jesdn_PileUpDataMC[nPFJetAK4] = sup; }
			if(isrc==(nsrc+7)){ PFJetAK4_jesdn_PileUpPtBB[nPFJetAK4] = sup; }
			if(isrc==(nsrc+8)){ PFJetAK4_jesdn_PileUpPtEC1[nPFJetAK4] = sup; }
			if(isrc==(nsrc+9)){ PFJetAK4_jesdn_PileUpPtEC2[nPFJetAK4] = sup; }
			if(isrc==(nsrc+10)){ PFJetAK4_jesdn_PileUpPtRef[nPFJetAK4] = sup; }
			if(isrc==(nsrc+11)){ PFJetAK4_jesdn_RelativeFSR[nPFJetAK4] = sup; }
			if(isrc==(nsrc+12)){ PFJetAK4_jesdn_RelativeJEREC1[nPFJetAK4] = sup; }
			if(isrc==(nsrc+13)){ PFJetAK4_jesdn_RelativeJEREC2[nPFJetAK4] = sup; }
			if(isrc==(nsrc+14)){ PFJetAK4_jesdn_RelativePtBB[nPFJetAK4] = sup; }
			if(isrc==(nsrc+15)){ PFJetAK4_jesdn_RelativePtEC1[nPFJetAK4] = sup; }
			if(isrc==(nsrc+16)){ PFJetAK4_jesdn_RelativePtEC2[nPFJetAK4] = sup; }
			if(isrc==(nsrc+17)){ PFJetAK4_jesdn_RelativeBal[nPFJetAK4] = sup; }
			if(isrc==(nsrc+18)){ PFJetAK4_jesdn_RelativeSample[nPFJetAK4] = sup; }
			if(isrc==(nsrc+19)){ PFJetAK4_jesdn_RelativeStatEC[nPFJetAK4] = sup; }
			if(isrc==(nsrc+20)){ PFJetAK4_jesdn_RelativeStatFSR[nPFJetAK4] = sup; }
			if(isrc==(nsrc+21)){ PFJetAK4_jesdn_SinglePionECAL[nPFJetAK4] = sup; }
			if(isrc==(nsrc+22)){ PFJetAK4_jesdn_SinglePionHCAL[nPFJetAK4] = sup; }
			if(isrc==(nsrc+23)){ PFJetAK4_jesdn_TimePtEta[nPFJetAK4] = sup; }
			if(isrc==(nsrc+24)){ PFJetAK4_jesdn_Total[nPFJetAK4] = sup; }
			
		}
	
    }
      
    // JES uncertainty Ends //
    
    // Jet id //
      
    JetIDVars AK4idvars;
      
    AK4idvars.NHF = ak4jet.neutralHadronEnergyFraction();
    AK4idvars.NEMF = ak4jet.neutralEmEnergyFraction();
    AK4idvars.MUF = ak4jet.muonEnergyFraction();
    AK4idvars.CHF = ak4jet.chargedHadronEnergyFraction();
    AK4idvars.CEMF = ak4jet.chargedEmEnergyFraction();
    AK4idvars.NumConst = (ak4jet.chargedMultiplicity()+ak4jet.neutralMultiplicity());
    AK4idvars.NumNeutralParticle = ak4jet.neutralMultiplicity();
    AK4idvars.CHM = ak4jet.chargedHadronMultiplicity();
     
    PFJetAK4_jetID[nPFJetAK4] = getJetID(AK4idvars,"CHS",year,PFJetAK4_eta[nPFJetAK4],false,isUltraLegacy);
    PFJetAK4_jetID_tightlepveto[nPFJetAK4] = getJetID(AK4idvars,"CHS",year,PFJetAK4_eta[nPFJetAK4],true,isUltraLegacy);
      
    PFJetAK4_hadronflav[nPFJetAK4] = ak4jet.hadronFlavour();
    PFJetAK4_partonflav[nPFJetAK4] = ak4jet.partonFlavour();
      
    PFJetAK4_qgl[nPFJetAK4] = ak4jet.userFloat("QGTagger:qgLikelihood");
    PFJetAK4_PUID[nPFJetAK4] = ak4jet.userFloat("pileupJetId:fullDiscriminant");
    PFJetAK4_PUID[nPFJetAK4] = ak4jet.userFloat("pileupJetId:fullDiscriminant");
    PFJetAK4_PUID[nPFJetAK4] = ak4jet.userFloat("pileupJetId:fullDiscriminant");
    
    // B tagging stuffs //
    
    PFJetAK4_btag_DeepCSV[nPFJetAK4] = ak4jet.bDiscriminator("pfDeepCSVJetTags:probb")+ak4jet.bDiscriminator("pfDeepCSVJetTags:probbb");
    PFJetAK4_btag_DeepFlav[nPFJetAK4] = ak4jet.bDiscriminator("pfDeepFlavourJetTags:probb") + ak4jet.bDiscriminator("pfDeepFlavourJetTags:probbb")+ak4jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
    
    BTagEntry::JetFlavor btv_flav;
    if(abs(PFJetAK4_hadronflav[nPFJetAK4])==5){ btv_flav = BTagEntry::FLAV_B; }
    else if (abs(PFJetAK4_hadronflav[nPFJetAK4])==4){ btv_flav = BTagEntry::FLAV_C; }
    else { btv_flav = BTagEntry::FLAV_UDSG; }
    
    if(read_btagSF){
		
		PFJetAK4_btag_DeepCSV_SF[nPFJetAK4] = reader_deepcsv.eval_auto_bounds("central",btv_flav,fabs(pfjetAK4_4v.Eta()),tmprecpt); 
		PFJetAK4_btag_DeepCSV_SF_up[nPFJetAK4] = reader_deepcsv.eval_auto_bounds("up",btv_flav,fabs(pfjetAK4_4v.Eta()),tmprecpt);
		PFJetAK4_btag_DeepCSV_SF_dn[nPFJetAK4] = reader_deepcsv.eval_auto_bounds("down",btv_flav,fabs(pfjetAK4_4v.Eta()),tmprecpt);
	
		PFJetAK4_btag_DeepFlav_SF[nPFJetAK4] = reader_deepflav.eval_auto_bounds("central",btv_flav,fabs(pfjetAK4_4v.Eta()),tmprecpt); 
		PFJetAK4_btag_DeepFlav_SF_up[nPFJetAK4] = reader_deepflav.eval_auto_bounds("up",btv_flav,fabs(pfjetAK4_4v.Eta()),tmprecpt);
		PFJetAK4_btag_DeepFlav_SF_dn[nPFJetAK4] = reader_deepflav.eval_auto_bounds("down",btv_flav,fabs(pfjetAK4_4v.Eta()),tmprecpt);
	}
	// Note that btag SF is derived after applying JEC //

    nPFJetAK4++;	
    if(nPFJetAK4 >= njetmx) { break;}

    
  }
   //BTAG iterattion method
   double btgSF[9];
   double btgSF_up[9];
   double btgSF_dn[9];
   for(int qq =0 ; qq < 9; qq++)
   { 
      btgSF[qq] = 1.0;
      btgSF_up[qq] = 1.0;
      btgSF_dn[qq] = 1.0;
   }
  
   if(isMC){
    std::string btag_unc_up[9]= {"up_jes"  ,"up_lf"  ,"up_hf"  ,"up_hfstats1"  ,"up_hfstats2"  ,"up_lfstats1"  ,"up_lfstats2"  ,"up_cferr1"  ,"up_cferr2"};
    std::string btag_unc_dn[9]= {"down_jes","down_lf","down_hf","down_hfstats1","down_hfstats2","down_lfstats1","down_lfstats2","down_cferr1","down_cferr2"};
    for(int qq =0 ; qq < 9; qq++)
    {

      double jetSF = 1;
      double jetSF_Up   = 1;
      double jetSF_Down = 1;

      double btgSF_temp = 1;
      double btgSF_up_temp   = 1;
      double btgSF_dn_temp = 1;

      for(int jj = 0; jj < nPFJetAK4; jj++)
      {
        if(fabs(PFJetAK4_eta[jj]) > 2.5 ) continue;
        if(PFJetAK4_pt[jj] < 20.0 ) continue;
        jetSF = 1;
        jetSF_Up   = 1;
        jetSF_Down = 1;
        if( fabs(PFJetAK4_hadronflav[jj])==5 )
         {
           jetSF      = reader.eval(BTagEntry::FLAV_B, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
           jetSF_Up   = reader.eval_auto_bounds(btag_unc_up[qq],BTagEntry::FLAV_B, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
           jetSF_Down = reader.eval_auto_bounds(btag_unc_dn[qq],BTagEntry::FLAV_B, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
         }
        else if( fabs(PFJetAK4_hadronflav[jj])==4 )
         {
           jetSF      = reader.eval(BTagEntry::FLAV_C, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
           jetSF_Up   = reader.eval_auto_bounds(btag_unc_up[qq],BTagEntry::FLAV_C, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
           jetSF_Down = reader.eval_auto_bounds(btag_unc_dn[qq],BTagEntry::FLAV_C, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
         }
        else 
         {
           jetSF      = reader.eval(BTagEntry::FLAV_UDSG, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
           jetSF_Up   = reader.eval_auto_bounds(btag_unc_up[qq],BTagEntry::FLAV_UDSG, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
           jetSF_Down = reader.eval_auto_bounds(btag_unc_dn[qq],BTagEntry::FLAV_UDSG, fabs(PFJetAK4_eta[jj]),PFJetAK4_pt[jj],PFJetAK4_btag_DeepFlav[jj]);
         }
        if (jetSF!=0) btgSF_temp*=jetSF;
        if (jetSF_Up!=0)   btgSF_up_temp*=jetSF_Up;
        if (jetSF_Down!=0) btgSF_dn_temp*=jetSF_Down;  
      }
      btgSF[qq] = btgSF_temp;
      if(btgSF_up_temp != 1)  
             { btgSF_up[qq] = btgSF_up_temp;}
      else {btgSF_up[qq] = btgSF_temp;}

      if(btgSF_dn_temp != 1)  
             { btgSF_dn[qq] = btgSF_dn_temp;}
      else {btgSF_dn[qq] = btgSF_temp;}
    }
    BTAG_SF = btgSF[0];
    BTAG_jes_up      = btgSF_up[0]; BTAG_jes_dn      = btgSF_dn[0]; 
    BTAG_lf_up       = btgSF_up[1]; BTAG_lf_dn       = btgSF_dn[1];
    BTAG_hf_up       = btgSF_up[2]; BTAG_hf_dn       = btgSF_dn[2]; 
    BTAG_hfstats1_up = btgSF_up[3]; BTAG_hfstats1_dn = btgSF_dn[3]; 
    BTAG_hfstats2_up = btgSF_up[4]; BTAG_hfstats2_dn = btgSF_dn[4]; 
    BTAG_lfstats1_up = btgSF_up[5]; BTAG_lfstats1_dn = btgSF_dn[5];
    BTAG_lfstats2_up = btgSF_up[6]; BTAG_lfstats2_dn = btgSF_dn[6];    
    BTAG_cferr1_up   = btgSF_up[7]; BTAG_cferr1_dn   = btgSF_dn[7];
    BTAG_cferr2_up   = btgSF_up[8]; BTAG_cferr2_dn   = btgSF_dn[8]; 
  }
  // Tau leptons //


 
 
  nTau = 0;

  for(const auto& tau1 : iEvent.get(tok_taus_) ) {

	if (tau1.pt()<mintauPt) continue;
    if(fabs(tau1.eta())>2.3) continue;

    Tau_pt[nTau] = tau1.pt();
    Tau_eta[nTau] = tau1.eta();
    Tau_phi[nTau] = tau1.phi();
    Tau_e[nTau] = tau1.energy();
    Tau_charge[nTau] = tau1.charge();
    Tau_isPF[nTau] = tau1.isPFTau();
    Tau_dxy[nTau] = tau1.dxy();
    //Tau_dz[nTau] = tau1.dz();

    if(!tau1.leadTrack().isNull()){
		Tau_leadtrkdxy[nTau] = tau1.leadTrack()->dxy(vertex.position());
        Tau_leadtrkdz[nTau] = tau1.leadTrack()->dz(vertex.position());
    }

    if(!tau1.leadChargedHadrCand().isNull()){

	//  Tau_dxy[nTau] = tau1.leadChargedHadrCand()->dxy(vertex.position());
	//  Tau_dz[nTau] = tau1.leadChargedHadrCand()->dz(vertex.position());
	
		Tau_leadtrkpt[nTau] = tau1.leadChargedHadrCand()->pt();
        Tau_leadtrketa[nTau] = tau1.leadChargedHadrCand()->eta();
        Tau_leadtrkphi[nTau] = tau1.leadChargedHadrCand()->phi();
    }
    
    // Id & iso variables //
    
    Tau_eiso2018_raw[nTau] = tau1.tauID("againstElectronMVA6Raw2018");
    Tau_eiso2018[nTau] = (0 + (int(tau1.tauID("againstElectronVLooseMVA62018"))) + (1<<(1*int(tau1.tauID("againstElectronLooseMVA62018")))) + (1<<(2*int(tau1.tauID("againstElectronMediumMVA62018")))) + (1<<(3*int(tau1.tauID("againstElectronTightMVA62018")))) + (1<<(4*int(tau1.tauID("againstElectronVTightMVA62018")))));
    
    Tau_jetiso_deeptau2017v2p1_raw[nTau] = tau1.tauID("byDeepTau2017v2p1VSjetraw");
    Tau_jetiso_deeptau2017v2p1[nTau] = (0 + (int(tau1.tauID("byVVVLooseDeepTau2017v2p1VSjet"))) + (1<<(1*int(tau1.tauID("byVVLooseDeepTau2017v2p1VSjet")))) + (1<<(2*int(tau1.tauID("byVLooseDeepTau2017v2p1VSjet")))) + (1<<(3*int(tau1.tauID("byLooseDeepTau2017v2p1VSjet")))) + (1<<(4*int(tau1.tauID("byMediumDeepTau2017v2p1VSjet")))) + (1<<(5*int(tau1.tauID("byTightDeepTau2017v2p1VSjet")))) + (1<<(6*int(tau1.tauID("byVTightDeepTau2017v2p1VSjet")))) + (1<<(7*int(tau1.tauID("byVVTightDeepTau2017v2p1VSjet")))) );

    Tau_eiso_deeptau2017v2p1_raw[nTau] = tau1.tauID("byDeepTau2017v2p1VSeraw");
    Tau_eiso_deeptau2017v2p1[nTau] = (0 + (int(tau1.tauID("byVVVLooseDeepTau2017v2p1VSe"))) + (1<<(1*int(tau1.tauID("byVVLooseDeepTau2017v2p1VSe")))) + (1<<(2*int(tau1.tauID("byVLooseDeepTau2017v2p1VSe")))) + (1<<(3*int(tau1.tauID("byLooseDeepTau2017v2p1VSe")))) + (1<<(4*int(tau1.tauID("byMediumDeepTau2017v2p1VSe")))) + (1<<(5*int(tau1.tauID("byTightDeepTau2017v2p1VSe")))) + (1<<(6*int(tau1.tauID("byVTightDeepTau2017v2p1VSe")))) + (1<<(7*int(tau1.tauID("byVVTightDeepTau2017v2p1VSe")))) );

    Tau_muiso_deeptau2017v2p1_raw[nTau] = tau1.tauID("byDeepTau2017v2p1VSmuraw");
    Tau_muiso_deeptau2017v2p1[nTau] = (0 + (int(tau1.tauID("byVLooseDeepTau2017v2p1VSmu"))) + (1<<(1*int(tau1.tauID("byLooseDeepTau2017v2p1VSmu")))) + (1<<(2*int(tau1.tauID("byMediumDeepTau2017v2p1VSmu")))) + (1<<(3*int(tau1.tauID("byTightDeepTau2017v2p1VSmu")))) );
    
    Tau_decayMode[nTau] = tau1.decayMode();
    Tau_decayModeinding[nTau] = tau1.tauID("decayModeFinding");
    Tau_decayModeindingNewDMs[nTau] = tau1.tauID("decayModeFindingNewDMs");

    Tau_rawiso[nTau] = tau1.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    Tau_rawisodR03[nTau] = (tau1.tauID("chargedIsoPtSumdR03")+TMath::Max(0.,tau1.tauID("neutralIsoPtSumdR03")-0.072*tau1.tauID("puCorrPtSum")));
    Tau_puCorr[nTau] = tau1.tauID("puCorrPtSum");

    if (++nTau>=njetmx) break;

  }
  
  // Photons //
  
  nPhoton = 0;
  edm::Handle<edm::View<pat::Photon>> photons;
  iEvent.getByToken(tok_photons_, photons);
  
  // for raw MVA score info //
  edm::Handle <edm::ValueMap <float> > mvaPhoID_FallV2_raw;
  iEvent.getByToken(tok_mvaPhoID_FallV2_raw, mvaPhoID_FallV2_raw);

  for(const auto& gamma1 : photons->ptrs() ) {

	if(gamma1->pt() < mingmPt) continue;

    Photon_pt[nPhoton] = gamma1->pt();
    Photon_e[nPhoton] = gamma1->energy();
    Photon_eta[nPhoton] = gamma1->eta();
    Photon_phi[nPhoton] = gamma1->phi();
    Photon_e1by9[nPhoton] = gamma1->maxEnergyXtal()/max(float(1),gamma1->e3x3());
    if (gamma1->hasConversionTracks()) { Photon_e1by9[nPhoton] *= -1; }
	Photon_e9by25[nPhoton] = gamma1->r9();
    Photon_hadbyem[nPhoton] = gamma1->hadronicOverEm();
    Photon_ietaieta[nPhoton] = gamma1->sigmaIetaIeta();

	// MVA id //

    Photon_mvaid_Fall17V2_WP90[nPhoton] =  gamma1->photonID(mPhoID_FallV2_WP90);
    Photon_mvaid_Fall17V2_WP80[nPhoton] = gamma1->photonID(mPhoID_FallV2_WP80);
//    Photon_mvaid_Fall17V2_raw[nPhoton] = (*mvaPhoID_FallV2_raw)[gamma1];
    Photon_mvaid_Fall17V2_raw[nPhoton] = gamma1->userFloat("PhotonMVAEstimatorRunIIFall17v2Values");
    Photon_mvaid_Spring16V1_WP90[nPhoton] = gamma1->photonID(mPhoID_SpringV1_WP90);
    Photon_mvaid_Spring16V1_WP80[nPhoton] = gamma1->photonID(mPhoID_SpringV1_WP80);

	// Isolation variables //

    Photon_trkiso[nPhoton] = gamma1->trkSumPtSolidConeDR04();
    Photon_emiso[nPhoton] = gamma1->ecalRecHitSumEtConeDR04();
    Photon_hadiso[nPhoton] = gamma1->hcalTowerSumEtConeDR04();
    Photon_phoiso[nPhoton] = gamma1->photonIso() ;
    Photon_chhadiso[nPhoton] = gamma1->chargedHadronIso();
    Photon_neuhadiso[nPhoton] = gamma1->neutralHadronIso();
    Photon_passEveto[nPhoton] = gamma1->passElectronVeto();
    Photon_PixelSeed[nPhoton] = gamma1->hasPixelSeed();

    Photon_energyScaleValue[nPhoton] = gamma1->userFloat("energyScaleValue");
    Photon_energySigmaValue[nPhoton] = gamma1->userFloat("energySigmaValue");
    Photon_energyScaleUp[nPhoton] = gamma1->userFloat("energyScaleUp");
    Photon_energyScaleDown[nPhoton] = gamma1->userFloat("energyScaleDown");
    Photon_energySigmaUp[nPhoton] = gamma1->userFloat("energySigmaUp");
    Photon_energySigmaDown[nPhoton] = gamma1->userFloat("energySigmaDown");
    Photon_energyScaleStatUp[nPhoton] = gamma1->userFloat("energyScaleStatUp");
    Photon_energyScaleStatDown[nPhoton] = gamma1->userFloat("energyScaleStatDown");
    Photon_energyScaleSystUp[nPhoton] = gamma1->userFloat("energyScaleSystUp");
    Photon_energyScaleSystDown[nPhoton] = gamma1->userFloat("energyScaleSystDown");
    Photon_energyScaleGainUp[nPhoton] = gamma1->userFloat("energyScaleGainUp");
    Photon_energyScaleGainDown[nPhoton] = gamma1->userFloat("energyScaleGainDown");
    Photon_energySigmaRhoUp[nPhoton] = gamma1->userFloat("energySigmaRhoUp");
    Photon_energySigmaRhoDown[nPhoton] = gamma1->userFloat("energySigmaRhoDown");
    Photon_energySigmaPhiUp[nPhoton] = gamma1->userFloat("energySigmaPhiUp");
    Photon_energySigmaPhiDown[nPhoton] = gamma1->userFloat("energySigmaPhiDown");
    Photon_energySmearNrSigma[nPhoton] = gamma1->userFloat("energySmearNrSigma");
    Photon_ecalEnergyPreCorr[nPhoton] = gamma1->userFloat("ecalEnergyPreCorr");
    Photon_ecalEnergyErrPreCorr[nPhoton] = gamma1->userFloat("ecalEnergyErrPreCorr");
    Photon_ecalEnergyPostCorr[nPhoton] = gamma1->userFloat("ecalEnergyPostCorr");
    Photon_ecalEnergyErrPostCorr[nPhoton] = gamma1->userFloat("ecalEnergyErrPostCorr");
    
    if (++nPhoton>=njetmx) break;

  }
       
  //cout<<"done!"<<endl;
  
  T2->Fill(); // filling the tree used to get sumofweights
  T1->Fill(); 

  // Skimming condition //
  
  //if(nPFJetAK8>=1 && (nMuon+nElectron)>=1){
  //	T1->Fill(); // filling the main tree
  //}
  
  // End of skimming 
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
Leptop::beginJob()
{
  
  Nevt = 0;
  
  ////JEC /////
  
  L1FastAK4       = new JetCorrectorParameters(mJECL1FastFileAK4.c_str());
  L2RelativeAK4   = new JetCorrectorParameters(mJECL2RelativeFileAK4.c_str());
  L3AbsoluteAK4   = new JetCorrectorParameters(mJECL3AbsoluteFileAK4.c_str());
  L2L3ResidualAK4 = new JetCorrectorParameters(mJECL2L3ResidualFileAK4.c_str());
  
  vecL1FastAK4.push_back(*L1FastAK4);
  vecL2RelativeAK4.push_back(*L2RelativeAK4);
  vecL3AbsoluteAK4.push_back(*L3AbsoluteAK4);
  vecL2L3ResidualAK4.push_back(*L2L3ResidualAK4);
  
  jecL1FastAK4       = new FactorizedJetCorrector(vecL1FastAK4);
  jecL2RelativeAK4   = new FactorizedJetCorrector(vecL2RelativeAK4);
  jecL3AbsoluteAK4   = new FactorizedJetCorrector(vecL3AbsoluteAK4);
  jecL2L3ResidualAK4 = new FactorizedJetCorrector(vecL2L3ResidualAK4);
  
  L1FastAK8       = new JetCorrectorParameters(mJECL1FastFileAK8.c_str());
  L2RelativeAK8   = new JetCorrectorParameters(mJECL2RelativeFileAK8.c_str());
  L3AbsoluteAK8   = new JetCorrectorParameters(mJECL3AbsoluteFileAK8.c_str());
  L2L3ResidualAK8 = new JetCorrectorParameters(mJECL2L3ResidualFileAK8.c_str());
  
  vecL1FastAK8.push_back(*L1FastAK8);
  vecL2RelativeAK8.push_back(*L2RelativeAK8);
  vecL3AbsoluteAK8.push_back(*L3AbsoluteAK8);
  vecL2L3ResidualAK8.push_back(*L2L3ResidualAK8);
  
  jecL1FastAK8       = new FactorizedJetCorrector(vecL1FastAK8);
  jecL2RelativeAK8   = new FactorizedJetCorrector(vecL2RelativeAK8);
  jecL3AbsoluteAK8   = new FactorizedJetCorrector(vecL3AbsoluteAK8);
  jecL2L3ResidualAK8 = new FactorizedJetCorrector(vecL2L3ResidualAK8);
  
  //BTag calibration
  BTagCalibration *calib = new BTagCalibration("DeepJet",mBtagSF_DeepFlav_itr.c_str()); 
  reader = BTagCalibrationReader(BTagEntry::OP_RESHAPING,"central",
                                                   {"up_jes","down_jes","up_lf","down_lf",
                                                        "up_hf","down_hf",
                                                        "up_hfstats1","down_hfstats1",
                                                        "up_hfstats2","down_hfstats2",
                                                        "up_lfstats1","down_lfstats1",
                                                        "up_lfstats2","down_lfstats2",
                                                        "up_cferr1","down_cferr1",
                                                        "up_cferr2","down_cferr2"});

  reader.load(*calib,BTagEntry::FLAV_B,"iterativefit");
  reader.load(*calib,BTagEntry::FLAV_C,"iterativefit");
  reader.load(*calib,BTagEntry::FLAV_UDSG,"iterativefit");
 
  for (int isrc = 0; isrc < nsrc; isrc++) {
    const char *name = jecsrcnames[isrc];
    JetCorrectorParameters *pAK4 = new JetCorrectorParameters(mJECUncFileAK4.c_str(), name) ;
    JetCorrectionUncertainty *uncAK4 = new JetCorrectionUncertainty(*pAK4);
    vsrc.push_back(uncAK4);
    JetCorrectorParameters *pAK8 = new JetCorrectorParameters(mJECUncFileAK8.c_str(), name) ;
    JetCorrectionUncertainty *uncAK8 = new JetCorrectionUncertainty(*pAK8);
    vsrcAK8.push_back(uncAK8);
  }
  
  if(read_btagSF){
	calib_deepcsv = BTagCalibration("DeepCSV", mBtagSF_DeepCSV.c_str());
	reader_deepcsv = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
	reader_deepcsv.load(calib_deepcsv, BTagEntry::FLAV_B, "comb");
	reader_deepcsv.load(calib_deepcsv, BTagEntry::FLAV_C, "comb");
	reader_deepcsv.load(calib_deepcsv, BTagEntry::FLAV_UDSG, "incl");
  
	calib_deepflav = BTagCalibration("DeepJet", mBtagSF_DeepFlav.c_str());
	reader_deepflav = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
	reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
	reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
	reader_deepflav.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");
  }
  
  if(isUltraLegacy)
  {
	if(year==2018){
		roch_cor.init((mRochcorFolder+"RoccoR2018UL.txt").c_str()); 
	}
	if(year==2017){
		roch_cor.init((mRochcorFolder+"RoccoR2017UL.txt").c_str()); 
	}
	if(year==2016){
		roch_cor.init((mRochcorFolder+"RoccoR2016aUL.txt").c_str()); 
	}
  }
  else{
		roch_cor.init((mRochcorFolder+"RoccoR2017.txt").c_str()); 
	  }
  
  //**Important**//
  //For precision top physics, change "comb" to "mujets" in BTagCalibrationReader above //
  //https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18#Additional_information
  
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Leptop::endJob() 
{
  theFile->cd();
  theFile->Write();
  theFile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
Leptop::beginRun(edm::Run const& iRun, edm::EventSetup const& pset)
{	
  bool changed(true);
  if(!isFastSIM){
	hltPrescaleProvider_.init(iRun,pset,theHLTTag,changed);
	HLTConfigProvider const&  hltConfig_ = hltPrescaleProvider_.hltConfigProvider();
  }
  //hltConfig_.dump("Triggers");
  //hltConfig_.dump("PrescaleTable");
}

// ------------ method called when ending the processing of a run  ------------
void 
Leptop::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Leptop::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Leptop::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Leptop::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Leptop);
