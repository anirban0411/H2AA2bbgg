// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>
#include <algorithm>

// ROOT includes
#include <TTree.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TPRegexp.h>

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

// CMSSW data formats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <DataFormats/Candidate/interface/ShallowCloneCandidate.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

//L1
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TH2D.h"

//BTag related
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondTools/BTau/interface/BTagCalibrationReader.h"

//JEC related
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
class TH1F;
class TH2F;
class TH2D;
class TStyle;
class TTree;

class TreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns> {

public:
  explicit TreeMaker(const edm::ParameterSet&);
  ~TreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:


  bool applyJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi);
  bool applyPileupJetID(const pat::Jet & jet, const std::string & level);

  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  //BTagCalibrationReader reader;
  TFile *f_vetomap;
  TH2D  *h_vetomap;
  TH2D  *h_vetomap_eep;

  void initializeBranches();

  // MC information
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken;
  const edm::EDGetTokenT<reco::GenParticleCollection > gensToken;
  const edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;
  const edm::EDGetTokenT<LHEEventProduct> lheInfoToken;
  // Trigger
  const edm::InputTag triggerResultsTag;  
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
  const edm::EDGetTokenT<edm::TriggerResults> metfilterspatLabel_;
  const edm::EDGetTokenT<edm::TriggerResults> metfiltersrecoLabel_;
  // Filter result
  const edm::InputTag filterResultsTag;  
  const edm::EDGetTokenT<edm::TriggerResults> filterResultsToken;
  const edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsToken;


  // Vertices from miniAOD
  const edm::EDGetTokenT<reco::VertexCollection > primaryVerticesToken;
  const edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> secondaryVerticesToken;
  edm::EDGetTokenT< double > rhoToken;
  edm::EDGetTokenT< double > rho_trk_Token;
  edm::EDGetTokenT< double > rho_calo_Token;
  // Pat objects from miniAOD
  const edm::EDGetTokenT<pat::MuonCollection> muonsToken;
  const edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
  const edm::EDGetTokenT<pat::PhotonCollection> photonsToken;
  const edm::EDGetTokenT<pat::TauCollection> tausToken;
  const edm::EDGetTokenT<pat::JetCollection> jetsToken;
  const edm::EDGetTokenT<pat::JetCollection> jetsAK8Token;
  std::string subJetCollectionName;
  const edm::EDGetTokenT<pat::METCollection> metToken;
  const edm::EDGetTokenT<reco::GenJetCollection > genJetsToken;
  const edm::EDGetTokenT<reco::GenJetCollection > genJetsWithNuToken;
  const edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
  //L1
  const edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1AlgosToken;
  const edm::ESGetToken<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd> l1GtMenuToken_;
  const edm::EDGetToken l1GtToken_;

  // trigger filters
  std::map<std::string,int> triggerPathsMap;

  // MET filters
  std::vector<std::string>   filterPathsVector;
  std::map<std::string, int> filterPathsMap;

  // Flag for the analyzer
  float muonPtMin, muonEtaMax;
  float electronPtMin, electronEtaMax;
  float photonPtMin, photonEtaMax;
  float tauPtMin, tauEtaMax;
  float jetPtMin;
  float jetAK8PtMin;
  float jetEtaMax;
  float dRJetGenMatch;
  bool  dumpOnlyJetMatchedToGen;
  bool  usePuppiJets;
  bool  saveLHEObjects;
  bool  isMC;
  std::vector<std::string> pnetDiscriminatorLabels;

  //jet veto map 
  bool event_veto_map = false;
  bool event_veto_map_eep = false; 

  Bool_t Flag_goodVertices_;
  Bool_t Flag_globalSuperTightHalo2016Filter_;
  Bool_t Flag_EcalDeadCellTriggerPrimitiveFilter_;
  Bool_t Flag_BadPFMuonFilter_;
  Bool_t Flag_BadPFMuonDzFilter_;
  Bool_t Flag_hfNoisyHitsFilter_;
  Bool_t Flag_eeBadScFilter_;
  Bool_t Flag_ecalBadCalibFilter_;

  // Event coordinates
  unsigned event, run, lumi;
  // Pileup information
  unsigned int putrue;
  // Flags for various event filters
  unsigned int flags;
  // Cross-section and event weight information for MC events
  std::string jetvetomapfile_;
  std::string mJECUncFile_;
  std::vector<JetCorrectionUncertainty*> vsrc;
  float xsec, wgt;
  float wgt_isr_up, wgt_isr_dn, wgt_fsr_up, wgt_fsr_dn;
  float wgt_scl_0, wgt_scl_1, wgt_scl_2, wgt_scl_3, wgt_scl_4, wgt_scl_5, wgt_scl_6, wgt_scl_7, wgt_scl_8;
  // Rho
  float rho;
  float rho_trk;
  float rho_calo;

  //Trigger
  std::vector<unsigned int> trigger_hlt_pass;
  std::vector<std::string>  trigger_hlt_path;

  // Generator-level information (Gen and LHE particles)
  std::vector<float>  lhe_particle_pt;
  std::vector<float>  lhe_particle_eta;
  std::vector<float>  lhe_particle_phi;
  std::vector<float>  lhe_particle_mass;
  std::vector<int>    lhe_particle_id;
  std::vector<unsigned int>  lhe_particle_status;

  std::vector<float>  gen_particle_pt;
  std::vector<float>  gen_particle_eta;
  std::vector<float>  gen_particle_phi;
  std::vector<float>  gen_particle_mass;
  std::vector<int>    gen_particle_id;
  std::vector<unsigned int>  gen_particle_status;
  std::vector<int>    gen_particle_daughters_id;
  std::vector<unsigned int> gen_particle_daughters_igen;
  std::vector<unsigned int> gen_particle_daughters_status;
  std::vector<float>  gen_particle_daughters_pt;
  std::vector<float>  gen_particle_daughters_eta;
  std::vector<float>  gen_particle_daughters_phi;
  std::vector<float>  gen_particle_daughters_mass;
  std::vector<int>    gen_particle_daughters_charge;
  // Vertex RECO
  unsigned int npv;
  unsigned int nsv;
  // Collection of muon 4-vectors, muon ID bytes, muon isolation values
  std::vector<float>  muon_pt;
  std::vector<float>  muon_eta;
  std::vector<float>  muon_phi;
  std::vector<float>  muon_mass;
  std::vector<unsigned int> muon_id;
  std::vector<int>    muon_charge;
  std::vector<unsigned int> muon_iso;
  std::vector<float>  muon_d0;
  std::vector<float>  muon_dz;
  // Collection of electron 4-vectors, electron ID bytes, electron isolation values
  std::vector<float>  electron_pt;
  std::vector<float>  electron_eta;
  std::vector<float>  electron_phi;
  std::vector<float>  electron_mass;
  std::vector<unsigned int> electron_id;
  std::vector<float>  electron_idscore;
  std::vector<int>    electron_charge;
  std::vector<float>  electron_d0;
  std::vector<float>  electron_dz;  
  // Collection of photon 4-vectors, photon ID bytes, photon isolation values                                                                                                                
  /*
  std::vector<float>  photon_pt;
  std::vector<float>  photon_pt_corr;
  std::vector<float>  photon_eta;
  std::vector<float>  photon_phi;
  std::vector<float>  photon_mass;
  std::vector<unsigned int> photon_id;
  std::vector<float>  photon_idscore;
  */
  // Collection of taus 4-vectors, tau ID bytes, tau isolation values
  /*
  std::vector<float>  tau_pt;
  std::vector<float>  tau_eta;
  std::vector<float>  tau_phi;
  std::vector<float>  tau_mass;
  std::vector<float>  tau_dxy;
  std::vector<float>  tau_dz;
  std::vector<unsigned int>  tau_decaymode;
  std::vector<unsigned int> tau_idjet_wp;
  std::vector<unsigned int> tau_idmu_wp;
  std::vector<unsigned int> tau_idele_wp;
  std::vector<float>  tau_idjet;
  std::vector<float>  tau_idele;
  std::vector<float>  tau_idmu;
  std::vector<int>    tau_charge;
  std::vector<float>  tau_genmatch_pt;
  std::vector<float>  tau_genmatch_eta;
  std::vector<float>  tau_genmatch_phi;
  std::vector<float>  tau_genmatch_mass;
  std::vector<int> tau_genmatch_decaymode;
  */
  // MET
  float met;
  float met_phi;
  float met_signif;
  float met_cov_00;
  float met_cov_01;
  float met_cov_10;
  float met_cov_11;
  // Collection of jet 4-vectors, jet ID bytes and b-tag discriminant values
  Int_t    njet;
  Int_t    njetak8;
  std::vector<float> jet_pt;
  std::vector<float> jet_eta;
  std::vector<float> jet_phi;
  std::vector<float> jet_mass;
  std::vector<float> jet_energy;
  std::vector<float> jet_pt_raw;
  std::vector<float> jet_mass_raw;
  std::vector<float> jet_energy_raw;
  std::vector<float> jet_chf;
  std::vector<float> jet_nhf;
  std::vector<float> jet_elf;
  std::vector<float> jet_phf;
  std::vector<float> jet_muf;
  std::vector<float> jet_deepjet_probb;
  std::vector<float> jet_deepjet_probbb;
  std::vector<float> jet_deepjet_problepb;
  std::vector<float> jet_deepjet_probc;
  std::vector<float> jet_deepjet_probg;
  std::vector<float> jet_deepjet_probuds;
  std::vector<float> jet_pnet_probb;
  std::vector<float> jet_pnet_probc;
  std::vector<float> jet_pnet_probuds;
  std::vector<float> jet_pnet_probg;
  std::vector<float> jet_pnet_probtauh;
  std::vector<float> jet_pnet_ptcorr;
  std::vector<float> jet_pnet_ptnu;
  std::vector<float> jet_pnet_ptres;
  std::vector<float> jet_pnet_probele;
  std::vector<float> jet_pnet_probmu;
  std::map<std::string,std::vector<float> > jet_pnetlast_score;
  std::vector<float> jet_breg_corr;
  std::vector<float> jet_breg_res;
  std::vector<float> jet_creg_corr;
  std::vector<float> jet_creg_res;
  std::vector<bool>  jet_veto;
  std::vector<bool>  jet_veto_eep;
  std::vector<unsigned int> jet_id;
  std::vector<unsigned int> jet_puid;
  std::vector<unsigned int> jet_ncand;
  std::vector<unsigned int> jet_nch;
  std::vector<unsigned int> jet_nnh;
  std::vector<unsigned int> jet_nel;
  std::vector<unsigned int> jet_nph;
  std::vector<unsigned int> jet_nmu;
  std::vector<unsigned int> jet_hflav;
  std::vector<int> jet_pflav;
  std::vector<unsigned int> jet_nbhad;
  std::vector<unsigned int> jet_nchad;
  std::vector<float> jet_genmatch_pt;
  std::vector<float> jet_genmatch_eta;
  std::vector<float> jet_genmatch_phi;
  std::vector<float> jet_genmatch_mass;
  std::vector<float> jet_genmatch_energy;
  std::vector<float> jet_genmatch_wnu_pt;
  std::vector<float> jet_genmatch_wnu_eta;
  std::vector<float> jet_genmatch_wnu_phi;
  std::vector<float> jet_genmatch_wnu_mass;
  std::vector<float> jet_genmatch_wnu_energy;

  //gen jet collection
  std::vector<float> genjet_pt;
  std::vector<float> genjet_eta;
  std::vector<float> genjet_phi;
  std::vector<float> genjet_energy;  
  std::vector<float> genjet_mass;
  std::vector<int> genjet_hflav;
  std::vector<int> genjet_pflav;
  std::vector<float> jet_jec_up;
  std::vector<float> jet_jec_dn;
  // AK8 jets
   
  std::vector<float> jetAK8_pt;
  std::vector<float> jetAK8_eta;
  std::vector<float> jetAK8_phi;
  std::vector<float> jetAK8_mass;
  std::vector<float> jetAK8_pt_raw;
  std::vector<float> jetAK8_mass_raw;
  std::vector<float> jetAK8_pnet_probHbb;
  std::vector<float> jetAK8_pnet_probHcc;
  std::vector<float> jetAK8_pnet_probHqq;
  std::vector<float> jetAK8_pnet_probHgg;
  std::vector<float> jetAK8_pnet_probHtt;
  std::vector<float> jetAK8_pnet_probHtm;
  std::vector<float> jetAK8_pnet_probHte;
  std::vector<float> jetAK8_pnet_probQCD2HF;
  std::vector<float> jetAK8_pnet_probQCD1HF;
  std::vector<float> jetAK8_pnet_probQCD0HF;
  std::vector<float> jetAK8_pnet_HbbVsQCD;
  std::vector<float> jetAK8_pnet_mass;
  std::vector<float> jetAK8_pnet_corr;
  std::vector<unsigned int> jetAK8_id;
  std::vector<unsigned int> jetAK8_ncand;
  std::vector<unsigned int> jetAK8_hflav;
  std::vector<int> jetAK8_pflav;
  std::vector<unsigned int> jetAK8_nbhad;
  std::vector<unsigned int> jetAK8_nchad;
  std::vector<float> jetAK8_softdrop_pt;
  std::vector<float> jetAK8_softdrop_eta;
  std::vector<float> jetAK8_softdrop_phi;
  std::vector<float> jetAK8_softdrop_mass;
  std::vector<float> jetAK8_softdrop_pt_raw;
  std::vector<float> jetAK8_softdrop_mass_raw;

  std::vector<float> jetAK8_softdrop_subjet_pt;
  std::vector<float> jetAK8_softdrop_subjet_pt_raw;
  std::vector<float> jetAK8_softdrop_subjet_eta;
  std::vector<float> jetAK8_softdrop_subjet_phi;
  std::vector<float> jetAK8_softdrop_subjet_mass;
  std::vector<float> jetAK8_softdrop_subjet_mass_raw;

  std::vector<float> jetAK8_legacy_pnet_mass;
  std::vector<float> jetAK8_legacy_pnet_probHbb;
  std::vector<float> jetAK8_legacy_pnet_QCD;

  //L1 object
  Int_t    nL1jet;
  std::vector<double> L1jet_pt;
  std::vector<double> L1jet_eta;
  std::vector<double> L1jet_phi;
  std::vector<double> L1jet_en;
  //Calo jet
  Int_t    nCalojet;
  std::vector<double> Calojet_pt;
  std::vector<double> Calojet_eta;
  std::vector<double> Calojet_phi;
  std::vector<double> Calojet_en;
 //Pixel jet
  Int_t   nPxljet;
  std::vector<double> Pxljet_pt;
  std::vector<double> Pxljet_eta;
  std::vector<double> Pxljet_phi;
  std::vector<double> Pxljet_en;
  //L3 jet
  Int_t    nL3jet;
  std::vector<double> L3jet_pt;
  std::vector<double> L3jet_eta;
  std::vector<double> L3jet_phi;
  std::vector<double> L3jet_en;
  //AK8 CALO-jet
  Int_t          nCaloAK8jet;
  vector<double> CaloAK8jet_pt;
  vector<double> CaloAK8jet_eta;
  vector<double> CaloAK8jet_phi;
  vector<double> CaloAK8jet_en;
  //AK8 L3-jet
  Int_t          nL3AK8jet1;
  vector<double> L3AK8jet1_pt;
  vector<double> L3AK8jet1_eta;
  vector<double> L3AK8jet1_phi;
  vector<double> L3AK8jet1_en;

  //L1 object
  Bool_t L1_QuadJet60er2p5= false;
  Bool_t L1_HTT280er= false;
  Bool_t L1_HTT320er= false;
  Bool_t L1_HTT360er= false;
  Bool_t L1_HTT400er= false;
  Bool_t L1_HTT450er= false;
  Bool_t L1_HTT280er_QuadJet_70_55_40_35_er2p5= false;
  Bool_t L1_HTT280er_QuadJet_70_55_40_40_er2p5= false;
  Bool_t L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3= false;
  Bool_t L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3= false;
  Bool_t L1_Mu6_HTT240er= false;
  Bool_t L1_SingleJet60= false;


  int idx_L1_HTT280er, idx_L1_QuadJet60er2p5, idx_L1_HTT320er, idx_L1_HTT360er, idx_L1_HTT400er, idx_L1_HTT450er;
  int idx_L1_HTT280er_QuadJet_70_55_40_35_er2p5, idx_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;
  int idx_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3, idx_L1_Mu6_HTT240er, idx_L1_SingleJet60;
  int idx_L1_HTT280er_QuadJet_70_55_40_40_er2p5;
 
  double Generator_weight;

  uint nPDFsets = 103;
  static const int nlhescalemax = 9;
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

  double LHE_weight;
  float Generator_qscale, Generator_x1, Generator_x2, Generator_xpdf1, Generator_xpdf2, Generator_scalePDF;
  int Generator_id1, Generator_id2;


  // TTree carrying the event weight information
  TTree* tree;
  TTree* tree1;

  // Sorters to order object collections in decreasing order of pT
  template<typename T> 
  class PatPtSorter {
  public:
    bool operator()(const T& i, const T& j) const {
      return (i.pt() > j.pt());
    }
  };
  PatPtSorter<pat::Muon>     muonSorter;
  PatPtSorter<pat::Electron> electronSorter;
  PatPtSorter<pat::Photon>   photonSorter;
  PatPtSorter<pat::Tau>      tauSorter;

  template<typename T> 
  class PatRefPtSorter {
  public:
    bool operator()(const T& i, const T& j) const {
      return (i->pt() > j->pt());
    }
  };
  PatRefPtSorter<pat::JetRef>      jetRefSorter;
  PatRefPtSorter<reco::GenJetRef>  genJetRefSorter;
  PatRefPtSorter<edm::Ptr<pat::Jet> > jetPtrSorter;

};

static const int nsrc = 1;
const char* srcnames[nsrc] = {"Total"};
const int njecmcmx = 2*nsrc + 1 ;

TreeMaker::TreeMaker(const edm::ParameterSet& iConfig): 
  pileupInfoToken          (mayConsume<std::vector<PileupSummaryInfo> >  (iConfig.getParameter<edm::InputTag>("pileUpInfo"))),
  gensToken                (mayConsume<reco::GenParticleCollection>   (iConfig.getParameter<edm::InputTag>("genParticles"))),
  genEvtInfoToken          (mayConsume<GenEventInfoProduct>  (iConfig.getParameter<edm::InputTag>("genEventInfo"))),
  lheInfoToken             (mayConsume<LHEEventProduct>               (iConfig.getParameter<edm::InputTag>("lheInfo"))),
  triggerResultsTag        (iConfig.getParameter<edm::InputTag>("triggerResults")),
  triggerResultsToken      (consumes<edm::TriggerResults>      (triggerResultsTag)),
  metfilterspatLabel_      (consumes<edm::TriggerResults>   (iConfig.getUntrackedParameter<edm::InputTag> ("metfilterspatLabel_"))),
  metfiltersrecoLabel_     (consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag> ("metfiltersrecoLabel_"))),
  filterResultsTag         (iConfig.getParameter<edm::InputTag>("filterResults")),
  filterResultsToken       (consumes<edm::TriggerResults>  (filterResultsTag)),
  triggerObjectsToken      (consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("objects"))),
  primaryVerticesToken     (consumes<reco::VertexCollection>  (iConfig.getParameter<edm::InputTag>("pVertices"))),
  secondaryVerticesToken   (consumes<reco::VertexCompositePtrCandidateCollection>  (iConfig.getParameter<edm::InputTag>("sVertices"))),
  rhoToken                 (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  rho_trk_Token            (consumes<double>(iConfig.getParameter<edm::InputTag>("rho_trk"))),
  rho_calo_Token           (consumes<double>(iConfig.getParameter<edm::InputTag>("rho_calo"))), 
  muonsToken               (consumes<pat::MuonCollection>  (iConfig.getParameter<edm::InputTag>("muons"))), 
  electronsToken           (consumes<pat::ElectronCollection>  (iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken             (consumes<pat::PhotonCollection>  (iConfig.getParameter<edm::InputTag>("photons"))), 
  tausToken                (consumes<pat::TauCollection>  (iConfig.getParameter<edm::InputTag>("taus"))), 
  jetsToken                (consumes<pat::JetCollection>  (iConfig.getParameter<edm::InputTag>("jets"))),
  jetsAK8Token             (consumes<pat::JetCollection>  (iConfig.getParameter<edm::InputTag>("jetsAK8"))),
  subJetCollectionName     (iConfig.getUntrackedParameter<string>("subjets")),
  metToken                 (consumes<pat::METCollection>  (iConfig.getParameter<edm::InputTag>("met"))),
  genJetsToken             (consumes<reco::GenJetCollection > (iConfig.getParameter<edm::InputTag>("genJets"))),
  genJetsWithNuToken       (consumes<reco::GenJetCollection > (iConfig.getParameter<edm::InputTag>("genJetsWithNu"))),
  jetFlavourInfosToken_    (consumes<reco::JetFlavourInfoMatchingCollection >    (iConfig.getParameter<edm::InputTag>("genJetsFlavour"))),
  l1GtMenuToken_           (esConsumes<L1TUtmTriggerMenu, L1TUtmTriggerMenuRcd>()),
  l1GtToken_               (consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1GtSrc"))),
  muonPtMin                (iConfig.existsAs<double>("muonPtMin")    ? iConfig.getParameter<double>("muonPtMin") : 15.),
  muonEtaMax               (iConfig.existsAs<double>("muonEtaMax")   ? iConfig.getParameter<double>("muonEtaMax") : 2.4),
  electronPtMin            (iConfig.existsAs<double>("electronPtMin")   ? iConfig.getParameter<double>("electronPtMin") : 15.),
  electronEtaMax           (iConfig.existsAs<double>("electronEtaMax")  ? iConfig.getParameter<double>("electronEtaMax") : 2.5),
  photonPtMin              (iConfig.existsAs<double>("photonPtMin")   ? iConfig.getParameter<double>("photonPtMin") : 15.),
  photonEtaMax             (iConfig.existsAs<double>("photonEtaMax")  ? iConfig.getParameter<double>("photonEtaMax") : 2.5),
  tauPtMin                 (iConfig.existsAs<double>("tauPtMin")   ? iConfig.getParameter<double>("tauPtMin") : 20.),
  tauEtaMax                (iConfig.existsAs<double>("tauEtaMax")  ? iConfig.getParameter<double>("tauEtaMax") : 2.5),
  jetPtMin                 (iConfig.existsAs<double>("jetPtMin")       ? iConfig.getParameter<double>("jetPtMin") : 25.),
  jetAK8PtMin              (iConfig.existsAs<double>("jetAK8PtMin")       ? iConfig.getParameter<double>("jetAK8PtMin") : 200.),
  jetEtaMax                (iConfig.existsAs<double>("jetEtaMax")      ? iConfig.getParameter<double>("jetEtaMax") : 2.5),
  dRJetGenMatch            (iConfig.existsAs<double>("dRJetGenMatch")        ? iConfig.getParameter<double>("dRJetGenMatch") : 0.4),
  dumpOnlyJetMatchedToGen  (iConfig.existsAs<bool>("dumpOnlyJetMatchedToGen")  ? iConfig.getParameter<bool>("dumpOnlyJetMatchedToGen") : false),
  usePuppiJets             (iConfig.existsAs<bool>("usePuppiJets") ? iConfig.getParameter<bool>("usePuppiJets") : false),
  saveLHEObjects           (iConfig.existsAs<bool>("saveLHEObjects")  ? iConfig.getParameter<bool>  ("saveLHEObjects") : true),
  isMC                     (iConfig.existsAs<bool>("isMC")   ? iConfig.getParameter<bool>  ("isMC") : true),
  pnetDiscriminatorLabels  (iConfig.existsAs<std::vector<std::string> > ("pnetDiscriminatorLabels") ? iConfig.getParameter<std::vector<std::string>>("pnetDiscriminatorLabels") : std::vector<std::string> ()),
  jetvetomapfile_          (iConfig.getParameter < std::string > ("jetvetomapfile")),
  mJECUncFile_             (iConfig.getParameter < std::string > ("JECUncFile")),
  xsec                     (iConfig.existsAs<double>("xsec") ? iConfig.getParameter<double>("xsec") * 1000.0 : 1.){
  
  usesResource("TFileService");

}

TreeMaker::~TreeMaker() {}

void TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // triggers
  edm::Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);

  edm::Handle<edm::TriggerResults>  filterResultsH;
  iEvent.getByToken(filterResultsToken, filterResultsH); 

  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
  iEvent.getByToken(triggerObjectsToken, triggerObjects);

  edm::Handle<reco::VertexCollection > primaryVerticesH;
  iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

  edm::Handle<reco::VertexCompositePtrCandidateCollection> secondaryVerticesH;
  iEvent.getByToken(secondaryVerticesToken, secondaryVerticesH);
  reco::VertexCompositePtrCandidateCollection svColl = *secondaryVerticesH;

  edm::Handle<double> rho_val;
  iEvent.getByToken(rhoToken,rho_val);  
  rho = *rho_val;
 
  edm::Handle<double> rho_trk_val;
  iEvent.getByToken(rho_trk_Token,rho_trk_val);
  rho_trk = *rho_trk_val;

  edm::Handle<double> rho_calo_val;
  iEvent.getByToken(rho_calo_Token,rho_calo_val);
  rho_calo = *rho_calo_val;

  edm::Handle<pat::MuonCollection> muonsH;
  iEvent.getByToken(muonsToken,muonsH);
  pat::MuonCollection muonsColl = *muonsH;

  edm::Handle<pat::ElectronCollection> electronsH;
  iEvent.getByToken(electronsToken,electronsH);
  pat::ElectronCollection electronsColl = *electronsH;

  edm::Handle<pat::PhotonCollection> photonsH;
  iEvent.getByToken(photonsToken,photonsH);
  pat::PhotonCollection photonsColl = *photonsH;

  edm::Handle<pat::TauCollection> tausH;
  iEvent.getByToken(tausToken,tausH);
  pat::TauCollection tausColl = *tausH;
    
  edm::Handle<pat::JetCollection> jetsH;
  iEvent.getByToken(jetsToken, jetsH);

  edm::Handle<pat::JetCollection> jetsAK8H;
  iEvent.getByToken(jetsAK8Token, jetsAK8H);

  edm::Handle<pat::METCollection> metH;
  iEvent.getByToken(metToken, metH);

  edm::ESHandle<L1TUtmTriggerMenu> menu;
  menu = iSetup.getHandle(l1GtMenuToken_);

  edm::Handle<BXVector<GlobalAlgBlk>> l1GtHandle;
  iEvent.getByToken(l1GtToken_, l1GtHandle);

  edm::Handle<std::vector<PileupSummaryInfo> > pileupInfoH;
  edm::Handle<GenEventInfoProduct> genEvtInfoH;
  edm::Handle<reco::GenParticleCollection> gensH;
  edm::Handle<LHEEventProduct> lheInfoH;
  edm::Handle<reco::GenJetCollection> genJetsH;
  edm::Handle<reco::GenJetCollection> genJetsWithNuH;
  edm::Handle<reco::JetFlavourInfoMatchingCollection> genJetsFlavourH;
  if(isMC){    
    iEvent.getByToken(pileupInfoToken, pileupInfoH);  
    iEvent.getByToken(genEvtInfoToken, genEvtInfoH);    
    iEvent.getByToken(lheInfoToken,lheInfoH);
    iEvent.getByToken(gensToken, gensH);
    iEvent.getByToken(genJetsToken, genJetsH);    
    iEvent.getByToken(genJetsWithNuToken, genJetsWithNuH);
    iEvent.getByToken(jetFlavourInfosToken_, genJetsFlavourH);
  }

  edm::Handle<edm::TriggerResults> METFilterResults;
  iEvent.getByToken(metfilterspatLabel_, METFilterResults);
  if(!(METFilterResults.isValid())) iEvent.getByToken(metfiltersrecoLabel_, METFilterResults);

  initializeBranches();
     
  event = iEvent.id().event();
  run   = iEvent.id().run();
  lumi  = iEvent.luminosityBlock();


  wgt = 1.0;
  if(genEvtInfoH.isValid())
    wgt = genEvtInfoH->weight(); 

  wgt_isr_up = 1.0, wgt_isr_dn = 1.0, wgt_fsr_up = 1.0, wgt_fsr_dn = 1.0;
  wgt_scl_0  = 1.0, wgt_scl_1 = 1.0, wgt_scl_2 = 1.0, wgt_scl_3 = 1.0;
  wgt_scl_4  = 1.0, wgt_scl_5 = 1.0, wgt_scl_6 = 1.0, wgt_scl_7 = 1.0, wgt_scl_8 = 1.0;
  if(isMC)
  {
    edm::Handle<LHEEventProduct>lheeventinfo ;
    iEvent.getByToken(lheInfoToken,lheeventinfo) ;
        nLHEScaleWeights = 0;
        nLHEPDFWeights = 0;
        nLHEAlpsWeights = 0;
        if(lheeventinfo.isValid()){
           LHE_weight = lheeventinfo->originalXWGTUP();
           for ( unsigned int index = 0; index < lheeventinfo->weights().size(); ++index ) {
                if(index<nlhescalemax && nLHEScaleWeights<nlhescalemax){
                        LHEScaleWeights[nLHEScaleWeights] = lheeventinfo->weights()[index].wgt/lheeventinfo->originalXWGTUP();
                        wgt_scl_0 = lheeventinfo->weights()[0].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_1 = lheeventinfo->weights()[1].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_2 = lheeventinfo->weights()[2].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_3 = lheeventinfo->weights()[3].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_4 = lheeventinfo->weights()[4].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_5 = lheeventinfo->weights()[5].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_6 = lheeventinfo->weights()[6].wgt/lheeventinfo->originalXWGTUP(); 
                        wgt_scl_7 = lheeventinfo->weights()[7].wgt/lheeventinfo->originalXWGTUP(); 
		        wgt_scl_8 = lheeventinfo->weights()[8].wgt/lheeventinfo->originalXWGTUP(); 
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
         }

      nLHEPSWeights = 8;
      if (genEvtInfoH.isValid()){
        Generator_weight = genEvtInfoH->weight();
        Generator_qscale = genEvtInfoH->qScale();
        Generator_x1 = (*genEvtInfoH->pdf()).x.first;
        Generator_x2 = (*genEvtInfoH->pdf()).x.second;
        Generator_id1 = (*genEvtInfoH->pdf()).id.first;
        Generator_id2 = (*genEvtInfoH->pdf()).id.second;
        Generator_xpdf1 = (*genEvtInfoH->pdf()).xPDF.first;
        Generator_xpdf2 = (*genEvtInfoH->pdf()).xPDF.second;
        Generator_scalePDF = (*genEvtInfoH->pdf()).scalePDF;
        if(genEvtInfoH->weights().size()>2){
                LHEPSWeights[0] = genEvtInfoH->weights()[2]/genEvtInfoH->weights()[1];
                LHEPSWeights[1] = genEvtInfoH->weights()[3]/genEvtInfoH->weights()[1];
                LHEPSWeights[2] = genEvtInfoH->weights()[4]/genEvtInfoH->weights()[1];
                LHEPSWeights[3] = genEvtInfoH->weights()[5]/genEvtInfoH->weights()[1];
                LHEPSWeights[4] = genEvtInfoH->weights()[24]/genEvtInfoH->weights()[1];
                LHEPSWeights[5] = genEvtInfoH->weights()[25]/genEvtInfoH->weights()[1];
                LHEPSWeights[6] = genEvtInfoH->weights()[26]/genEvtInfoH->weights()[1];
                LHEPSWeights[7] = genEvtInfoH->weights()[27]/genEvtInfoH->weights()[1];
                wgt_isr_up = genEvtInfoH->weights()[25]/genEvtInfoH->weights()[1]; wgt_isr_dn = genEvtInfoH->weights()[24]/genEvtInfoH->weights()[1];
                wgt_fsr_up = genEvtInfoH->weights()[3] /genEvtInfoH->weights()[1]; wgt_fsr_dn = genEvtInfoH->weights()[2] /genEvtInfoH->weights()[1];
                }
        }
       else
        {
                Generator_weight = Generator_qscale = Generator_x1 = Generator_x2 = Generator_id1 = Generator_id2 = Generator_xpdf1 = Generator_xpdf2 = Generator_scalePDF = -10000;
        }
     } //isMC  


  // Pileup information
  putrue = 0;
  if(pileupInfoH.isValid()){
    for (auto pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
      if (pileupInfo_iter->getBunchCrossing() == 0) 
	putrue = (unsigned int) pileupInfo_iter->getTrueNumInteractions();
    }
  }

  // Vertices RECO
  npv  = (unsigned int) primaryVerticesH->size();
  nsv  = (unsigned int) secondaryVerticesH->size();

  // MET filter info
  unsigned int flagvtx       = 1 * 1;
  unsigned int flaghalo      = 1 * 2;
  unsigned int flaghbhe      = 1 * 4;
  unsigned int flaghbheiso   = 1 * 8;
  unsigned int flagecaltp    = 1 * 16;
  unsigned int flagbadmuon   = 1 * 32; 
  unsigned int flagbadhad    = 1 * 64;
 
  // Which MET filters passed
  for (size_t i = 0; i < filterPathsVector.size(); i++) {
    //if (filterPathsMap[filterPathsVector[i]] == -1) continue;
    if (i == 0  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagvtx       = 0; // goodVertices
    if (i == 1  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghalo      = 0; // CSCTightHaloFilter
    if (i == 2  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbhe      = 0; // HBHENoiseFilter
    if (i == 3  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flaghbheiso   = 0; // HBHENoiseIsoFilter
    if (i == 4  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagecaltp    = 0; // EcalDeadCellTriggerPrimitiveFilter
    if (i == 5  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagbadmuon   = 0; // badmuon
    if (i == 6  and filterResultsH->accept(filterPathsMap[filterPathsVector[i]])) flagbadhad    = 0; // badhadrons
  }

  flags = flagvtx + flaghalo + flaghbhe + flaghbheiso + flagecaltp + flagbadmuon + flagbadhad;

  const edm::TriggerNames & metfilterName = iEvent.triggerNames(*METFilterResults);
  unsigned int goodVerticesIndex_ = metfilterName.triggerIndex("Flag_goodVertices");
  Flag_goodVertices_ = METFilterResults.product()->accept(goodVerticesIndex_);
  unsigned int globalSuperTightHalo2016FilterIndex_ = metfilterName.triggerIndex("Flag_globalSuperTightHalo2016Filter");
  Flag_globalSuperTightHalo2016Filter_ = METFilterResults.product()->accept(globalSuperTightHalo2016FilterIndex_);
  unsigned int EcalDeadCellTriggerPrimitiveFilterIndex_ = metfilterName.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter");
  Flag_EcalDeadCellTriggerPrimitiveFilter_ = METFilterResults.product()->accept(EcalDeadCellTriggerPrimitiveFilterIndex_);
  unsigned int BadPFMuonFilterIndex_ = metfilterName.triggerIndex("Flag_BadPFMuonFilter");
  Flag_BadPFMuonFilter_ = METFilterResults.product()->accept(BadPFMuonFilterIndex_);
  unsigned int BadPFMuonFilterDzIndex_ = metfilterName.triggerIndex("Flag_BadPFMuonDzFilter");
  Flag_BadPFMuonDzFilter_ = METFilterResults.product()->accept(BadPFMuonFilterDzIndex_);
  unsigned int hfNoisyHitsIndex_ = metfilterName.triggerIndex("Flag_hfNoisyHitsFilter");
  Flag_hfNoisyHitsFilter_ = METFilterResults.product()->accept(hfNoisyHitsIndex_);
  unsigned int eeBadScFilterIndex_ = metfilterName.triggerIndex("Flag_eeBadScFilter");
  Flag_eeBadScFilter_ = METFilterResults.product()->accept(eeBadScFilterIndex_);
  unsigned int ecalBadCalibFilterIndex_ = metfilterName.triggerIndex("Flag_ecalBadCalibFilter");
  Flag_ecalBadCalibFilter_ = METFilterResults.product()->accept(ecalBadCalibFilterIndex_); 


  // HLT triggers                                                                                                                                                                                     
  for(auto key : triggerPathsMap){
    trigger_hlt_path.push_back(key.first);
    trigger_hlt_pass.push_back(triggerResultsH->accept(key.second));
  }
        
  // LHE level information
  if(lheInfoH.isValid() and saveLHEObjects){
    const auto& hepeup = lheInfoH->hepeup();
    const auto& pup    = hepeup.PUP;
    for (unsigned int i = 0, n = pup.size(); i < n; ++i) {
      int status = hepeup.ISTUP[i];
      int id     = hepeup.IDUP[i];
      TLorentzVector p4(pup[i][0], pup[i][1], pup[i][2], pup[i][3]);
      if(abs(id) == 25 or abs(id) == 23 or abs(id) == 24 or abs(id)==6){ // Higgs, Z or W or top                                                                                             
        lhe_particle_pt.push_back(p4.Pt());
        lhe_particle_eta.push_back(p4.Eta());
        lhe_particle_phi.push_back(p4.Phi());
        lhe_particle_mass.push_back(p4.M());
        lhe_particle_status.push_back(status);
        lhe_particle_id.push_back(id);
      }
      else if(abs(id) >= 10 && abs(id) <= 17 and status == 1){ // final state leptons + neutrinos
        lhe_particle_pt.push_back(p4.Pt());
        lhe_particle_eta.push_back(p4.Eta());
        lhe_particle_phi.push_back(p4.Phi());
        lhe_particle_mass.push_back(p4.M());
        lhe_particle_status.push_back(status);
        lhe_particle_id.push_back(id);
      }
      else if(abs(id) >= 1 && abs(id) <= 5 and status == 1){ // final state quarks                                                                                                    
        lhe_particle_pt.push_back(p4.Pt());
        lhe_particle_eta.push_back(p4.Eta());
        lhe_particle_phi.push_back(p4.Phi());
        lhe_particle_mass.push_back(p4.M());
        lhe_particle_status.push_back(status);
	lhe_particle_id.push_back(id);
      }
      else if(abs(id) == 21 and status == 1){ // final state gluons                                                                                                               
        lhe_particle_pt.push_back(p4.Pt());
        lhe_particle_eta.push_back(p4.Eta());
        lhe_particle_phi.push_back(p4.Phi());
        lhe_particle_mass.push_back(p4.M());
        lhe_particle_status.push_back(status);
        lhe_particle_id.push_back(id);
      }
    }
  }

  // GEN particle info                                                                                                                                                                           
  if(gensH.isValid()){
    unsigned int igen = 0;
    for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
      // Higgs, Z and W beforethe decay and after the shower effects                                                                                                                            
      if((abs(gens_iter->pdgId()) == 25 or abs(gens_iter->pdgId()) == 24 or abs(gens_iter->pdgId()) == 23) and
         gens_iter->isLastCopy() and
         gens_iter->statusFlags().fromHardProcess()){
	
        gen_particle_pt.push_back(gens_iter->pt());
        gen_particle_eta.push_back(gens_iter->eta());
        gen_particle_phi.push_back(gens_iter->phi());
        gen_particle_mass.push_back(gens_iter->mass());
        gen_particle_id.push_back(gens_iter->pdgId());
        gen_particle_status.push_back(gens_iter->status());

        for(size_t idau = 0; idau < gens_iter->numberOfDaughters(); idau++){
          gen_particle_daughters_id.push_back(gens_iter->daughter(idau)->pdgId());
          gen_particle_daughters_igen.push_back(igen);
          gen_particle_daughters_pt.push_back(gens_iter->daughter(idau)->pt());
          gen_particle_daughters_eta.push_back(gens_iter->daughter(idau)->eta());
          gen_particle_daughters_phi.push_back(gens_iter->daughter(idau)->phi());
          gen_particle_daughters_mass.push_back(gens_iter->daughter(idau)->mass());
          gen_particle_daughters_status.push_back(gens_iter->daughter(idau)->status());
          gen_particle_daughters_charge.push_back(gens_iter->daughter(idau)->charge());
        }
        igen++;
      }
      // Final states Leptons (e,mu) and Neutrinos --> exclude taus. They need to be prompt or from Tau decay                                                                                        
      if (abs(gens_iter->pdgId()) > 10 and abs(gens_iter->pdgId()) < 17 and abs(gens_iter->pdgId()) != 15  and
          (gens_iter->isPromptFinalState() or
           gens_iter->isDirectPromptTauDecayProductFinalState())) {

        gen_particle_pt.push_back(gens_iter->pt());
        gen_particle_eta.push_back(gens_iter->eta());
        gen_particle_phi.push_back(gens_iter->phi());
        gen_particle_mass.push_back(gens_iter->mass());
        gen_particle_id.push_back(gens_iter->pdgId());
        gen_particle_status.push_back(gens_iter->status());
        igen++;
      }

      // Final state quarks or gluons from the hard process before the shower --> partons in which H/Z/W/top decay into                                                                              
      if (((abs(gens_iter->pdgId()) >= 1 and abs(gens_iter->pdgId()) <= 5) or abs(gens_iter->pdgId()) == 21) and
          gens_iter->statusFlags().fromHardProcess() and
          gens_iter->statusFlags().isFirstCopy()){
        gen_particle_pt.push_back(gens_iter->pt());
        gen_particle_eta.push_back(gens_iter->eta());
        gen_particle_phi.push_back(gens_iter->phi());
        gen_particle_mass.push_back(gens_iter->mass());
        gen_particle_id.push_back(gens_iter->pdgId());
        gen_particle_status.push_back(gens_iter->status());
        igen++;
      }

      // Special case of taus: last-copy, from hard process and, prompt and decayed                                                                                                             
      if(abs(gens_iter->pdgId()) == 15 and
         gens_iter->isLastCopy() and
         gens_iter->statusFlags().fromHardProcess() and
         gens_iter->isPromptDecayed()){ // hadronic taus                                                                                                                                              

        gen_particle_pt.push_back(gens_iter->pt());
        gen_particle_eta.push_back(gens_iter->eta());
        gen_particle_phi.push_back(gens_iter->phi());
        gen_particle_mass.push_back(gens_iter->mass());
        gen_particle_id.push_back(gens_iter->pdgId());
        gen_particle_status.push_back(gens_iter->status());

        // only store the final decay particles                                                                                                                                                     
        for(size_t idau = 0; idau < gens_iter->numberOfDaughters(); idau++){
          if(not dynamic_cast<const reco::GenParticle*>(gens_iter->daughter(idau))->statusFlags().isPromptTauDecayProduct()) continue;
          gen_particle_daughters_id.push_back(gens_iter->daughter(idau)->pdgId());
          gen_particle_daughters_igen.push_back(igen);
          gen_particle_daughters_pt.push_back(gens_iter->daughter(idau)->pt());
          gen_particle_daughters_eta.push_back(gens_iter->daughter(idau)->eta());
          gen_particle_daughters_phi.push_back(gens_iter->daughter(idau)->phi());
          gen_particle_daughters_mass.push_back(gens_iter->daughter(idau)->mass());
          gen_particle_daughters_status.push_back(gens_iter->daughter(idau)->status());
          gen_particle_daughters_charge.push_back(gens_iter->daughter(idau)->charge());
        }
        igen++;
      }
    }
  }

  // Sorting muons based on pT
  sort(muonsColl.begin(),muonsColl.end(),muonSorter);

  for (size_t i = 0; i < muonsColl.size(); i++) {

    if(muonsColl[i].pt() < muonPtMin) continue;
    if(fabs(muonsColl[i].eta()) > muonEtaMax) continue;

    muon_pt.push_back(muonsColl[i].pt());
    muon_eta.push_back(muonsColl[i].eta());
    muon_phi.push_back(muonsColl[i].phi());
    muon_mass.push_back(muonsColl[i].mass());

    // Muon isolation
    int isoval = 0;
    if(muonsColl[i].passed(reco::Muon::PFIsoLoose))
      isoval += 1;
    if(muonsColl[i].passed(reco::Muon::PFIsoMedium))
      isoval += 2;
    if(muonsColl[i].passed(reco::Muon::PFIsoTight))
      isoval += 4;
    if(muonsColl[i].passed(reco::Muon::PFIsoVeryTight))
      isoval += 8;
    if(muonsColl[i].passed(reco::Muon::MiniIsoLoose))
      isoval += 16;
    if(muonsColl[i].passed(reco::Muon::MiniIsoMedium))
      isoval += 32;
    if(muonsColl[i].passed(reco::Muon::MiniIsoTight))
      isoval += 64;
    
    muon_iso.push_back(isoval);
    
    // Muon id
    int midval = 0;
    if(muonsColl[i].passed(reco::Muon::CutBasedIdLoose))
       midval += 1;
    if(muonsColl[i].passed(reco::Muon::CutBasedIdMedium))
       midval += 2;
    if(muonsColl[i].passed(reco::Muon::CutBasedIdTight))
       midval += 4;
    if(muonsColl[i].passed(reco::Muon::MvaLoose))
       midval += 8;
    if(muonsColl[i].passed(reco::Muon::MvaMedium))
       midval += 16;
    if(muonsColl[i].passed(reco::Muon::MvaTight))
       midval += 32;
    
    muon_id.push_back(midval);
    muon_charge.push_back(muonsColl[i].charge());
    muon_d0.push_back(muonsColl[i].muonBestTrack()->dxy(primaryVerticesH->at(0).position()));
    muon_dz.push_back(muonsColl[i].muonBestTrack()->dz(primaryVerticesH->at(0).position()));
  }
  
  // Sorting electrons based on pT
  sort(electronsColl.begin(),electronsColl.end(),electronSorter);

  for (size_t i = 0; i < electronsColl.size(); i++) {

    if(electronsColl[i].pt() < electronPtMin) continue;
    if(fabs(electronsColl[i].eta()) > electronEtaMax) continue;

    electron_pt.push_back(electronsColl[i].pt());
    electron_eta.push_back(electronsColl[i].eta());
    electron_phi.push_back(electronsColl[i].phi());
    electron_mass.push_back(electronsColl[i].mass());

    int eidval = 0;   
    if(electronsColl[i].electronID("cutBasedElectronID-Fall17-94X-V2-loose"))
      eidval += 1;
    if(electronsColl[i].electronID("cutBasedElectronID-Fall17-94X-V2-medium"))
      eidval += 2;
    if(electronsColl[i].electronID("cutBasedElectronID-Fall17-94X-V2-tight"))
      eidval += 4;
    if(electronsColl[i].electronID("mvaEleID-Fall17-iso-V2-wpLoose"))
      eidval += 8;
    if(electronsColl[i].electronID("mvaEleID-Fall17-iso-V2-wp90"))
      eidval += 16;
    if(electronsColl[i].electronID("mvaEleID-Fall17-iso-V2-wp80"))
      eidval += 32;
    electron_id.push_back(eidval);
    electron_idscore.push_back(electronsColl[i].userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"));
    electron_charge.push_back(electronsColl[i].charge());
    electron_d0.push_back(electronsColl[i].gsfTrack()->dxy(primaryVerticesH->at(0).position()));
    electron_dz.push_back(electronsColl[i].gsfTrack()->dz(primaryVerticesH->at(0).position()));
  }

  /*
  // Sorting photons based on pT                                                                                                                                                                   
  sort(photonsColl.begin(),photonsColl.end(),photonSorter);

  for (size_t i = 0; i < photonsColl.size(); i++) {
    
    if(photonsColl[i].pt() < photonPtMin) continue;
    if(fabs(photonsColl[i].eta()) > photonEtaMax) continue;

    photon_pt.push_back(photonsColl[i].pt());
    photon_pt_corr.push_back(photonsColl[i].pt());
    //photon_pt_corr.push_back(photonsColl[i].pt()*photonsColl[i].userFloat("ecalEnergyPostCorr")/photonsColl[i].energy());
    photon_eta.push_back(photonsColl[i].eta());
    photon_phi.push_back(photonsColl[i].phi());
    photon_mass.push_back(photonsColl[i].mass());

    int pidval = 0;
    if(photonsColl[i].photonID("cutBasedPhotonID-Fall17-94X-V2-loose"))
      pidval += 1;
    if(photonsColl[i].photonID("cutBasedPhotonID-Fall17-94X-V2-medium"))
      pidval += 2;
    if(photonsColl[i].photonID("cutBasedPhotonID-Fall17-94X-V2-tight"))
      pidval += 4;
    if(photonsColl[i].photonID("mvaPhoID-RunIIFall17-v2-wp90"))
      pidval += 8;
    if(photonsColl[i].photonID("mvaPhoID-RunIIFall17-v2-wp80"))
      pidval += 16;

    photon_id.push_back(pidval);
    photon_idscore.push_back(photonsColl[i].userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));
  }

  // Sorting taus based on pT
  sort(tausColl.begin(),tausColl.end(),tauSorter);
  
  for (size_t i = 0; i < tausColl.size(); i++) {
    
    if(tausColl[i].pt() < tauPtMin) continue;
    if(fabs(tausColl[i].eta()) > tauEtaMax) continue;

    tau_pt.push_back(tausColl[i].pt());
    tau_eta.push_back(tausColl[i].eta());
    tau_phi.push_back(tausColl[i].phi());
    tau_mass.push_back(tausColl[i].mass());
    tau_dxy.push_back(tausColl[i].dxy());
    tau_dz.push_back((tausColl[i].leadChargedHadrCand().get() ? dynamic_cast<const pat::PackedCandidate*>(tausColl[i].leadChargedHadrCand().get())->dz() : 0.));
    tau_decaymode.push_back(tausColl[i].decayMode());
    tau_charge.push_back(tausColl[i].charge());    

    if(tausColl[i].isTauIDAvailable("byDeepTau2018v2p5VSjetraw"))
      tau_idjet.push_back(tausColl[i].tauID("byDeepTau2018v2p5VSjetraw"));
    else if(tausColl[i].isTauIDAvailable("byDeepTau2017v2p1VSjetraw"))
      tau_idjet.push_back(tausColl[i].tauID("byDeepTau2017v2p1VSjetraw"));

    if(tausColl[i].isTauIDAvailable("byDeepTau2018v2p5VSmuraw"))
      tau_idmu.push_back(tausColl[i].tauID("byDeepTau2018v2p5VSmuraw"));
    else if(tausColl[i].isTauIDAvailable("byDeepTau2017v2p1VSmuraw"))
      tau_idmu.push_back(tausColl[i].tauID("byDeepTau2017v2p1VSmuraw"));

    if(tausColl[i].isTauIDAvailable("byDeepTau2018v2p5VSeraw"))
      tau_idele.push_back(tausColl[i].tauID("byDeepTau2018v2p5VSeraw"));
    else if(tausColl[i].isTauIDAvailable("byDeepTau2017v2p1VSeraw"))
      tau_idele.push_back(tausColl[i].tauID("byDeepTau2017v2p1VSeraw"));

    int tauvsjetid = 0;
    if(tausColl[i].isTauIDAvailable("byVVVLooseDeepTau2017v2p1VSjet") and tausColl[i].tauID("byVVVLooseDeepTau2017v2p1VSjet"))
      tauvsjetid += 1;
    if(tausColl[i].isTauIDAvailable("byVVLooseDeepTau2017v2p1VSjet") and tausColl[i].tauID("byVVLooseDeepTau2017v2p1VSjet"))
      tauvsjetid += 2;
    if(tausColl[i].isTauIDAvailable("byVLooseDeepTau2017v2p1VSjet") and tausColl[i].tauID("byVLooseDeepTau2017v2p1VSjet"))
      tauvsjetid += 4;
    if(tausColl[i].isTauIDAvailable("byLooseDeepTau2017v2p1VSjet") and tausColl[i].tauID("byLooseDeepTau2017v2p1VSjet"))
      tauvsjetid += 8;
    if(tausColl[i].isTauIDAvailable("byMediumDeepTau2017v2p1VSjet") and tausColl[i].tauID("byMediumDeepTau2017v2p1VSjet"))
      tauvsjetid += 16;
    if(tausColl[i].isTauIDAvailable("byTightDeepTau2017v2p1VSjet") and tausColl[i].tauID("byTightDeepTau2017v2p1VSjet"))
      tauvsjetid += 32;

    tau_idjet_wp.push_back(tauvsjetid);

    int tauvsmuid = 0;
    if(tausColl[i].isTauIDAvailable("byVLooseDeepTau2017v2p1VSmu") and tausColl[i].tauID("byVLooseDeepTau2017v2p1VSmu"))
      tauvsmuid += 1;
    if(tausColl[i].isTauIDAvailable("byLooseDeepTau2017v2p1VSmu") and tausColl[i].tauID("byLooseDeepTau2017v2p1VSmu"))
      tauvsmuid += 2;
    if(tausColl[i].isTauIDAvailable("byMediumDeepTau2017v2p1VSmu") and tausColl[i].tauID("byMediumDeepTau2017v2p1VSmu"))
      tauvsmuid += 4;
    if(tausColl[i].isTauIDAvailable("byTightDeepTau2017v2p1VSmu") and tausColl[i].tauID("byTightDeepTau2017v2p1VSmu"))
      tauvsmuid += 8;

    tau_idmu_wp.push_back(tauvsmuid);

    int tauvselid = 0;
    if(tausColl[i].isTauIDAvailable("byVVVLooseDeepTau2017v2p1VSe") and tausColl[i].tauID("byVVVLooseDeepTau2017v2p1VSe"))
      tauvselid += 1;
    if(tausColl[i].isTauIDAvailable("byVVLooseDeepTau2017v2p1VSe") and tausColl[i].tauID("byVVLooseDeepTau2017v2p1VSe"))
      tauvselid += 2;
    if(tausColl[i].isTauIDAvailable("byVLooseDeepTau2017v2p1VSe") and tausColl[i].tauID("byVLooseDeepTau2017v2p1VSe"))
      tauvselid += 4;
    if(tausColl[i].isTauIDAvailable("byLooseDeepTau2017v2p1VSe") and tausColl[i].tauID("byLooseDeepTau2017v2p1VSe"))
      tauvselid += 8;
    if(tausColl[i].isTauIDAvailable("byMediumDeepTau2017v2p1VSe") and tausColl[i].tauID("byMediumDeepTau2017v2p1VSe"))
      tauvselid += 16;
    if(tausColl[i].isTauIDAvailable("byTightDeepTau2017v2p1VSe") and tausColl[i].tauID("byTightDeepTau2017v2p1VSe"))
      tauvselid += 32;

    tau_idele_wp.push_back(tauvselid);

    if(tausColl[i].genJet()){
      tau_genmatch_pt.push_back(tausColl[i].genJet()->pt());
      tau_genmatch_eta.push_back(tausColl[i].genJet()->eta());
      tau_genmatch_phi.push_back(tausColl[i].genJet()->phi());
      tau_genmatch_mass.push_back(tausColl[i].genJet()->mass());
      // reconstruct the decay mode of the gen-jet                                                                                                                                                     
      unsigned int tau_ch = 0;
      unsigned int tau_ph = 0;
      unsigned int tau_nh = 0;
      auto gen_constituents = tausColl[i].genJet()->getGenConstituents();
      for(size_t iconst = 0; iconst < gen_constituents.size(); iconst++){
        auto part = gen_constituents[iconst];
        if(part->status() != 1) continue;
        if(part->charge() == 0 and abs(part->pdgId()) == 22) tau_ph++;
        if(part->charge() == 0 and abs(part->pdgId()) != 22) tau_nh++;
        if(part->charge() != 0 and abs(part->pdgId()) != 11 and abs(part->pdgId()) != 13) tau_ch++;
      }
      tau_genmatch_decaymode.push_back(5*(tau_ch-1)+tau_ph/2+tau_nh);
    }
    else{
      tau_genmatch_pt.push_back(-1);
      tau_genmatch_eta.push_back(-1);
      tau_genmatch_phi.push_back(-1);
      tau_genmatch_mass.push_back(-1);
      tau_genmatch_decaymode.push_back(-1);
    }
  }
  */
  // Standard gen-jets excluding the neutrinos
   std::vector<reco::GenJetRef> jetv_gen;   
   if(genJetsH.isValid()){
     for (auto jets_iter = genJetsH->begin(); jets_iter != genJetsH->end(); ++jets_iter) {                                                                                                   
       reco::GenJetRef jref (genJetsH, jets_iter - genJetsH->begin());                                                                                                                      
       jetv_gen.push_back(jref);                                                                                                                                                              
     }
     sort(jetv_gen.begin(), jetv_gen.end(), genJetRefSorter);
   }

   std::vector<reco::GenJetRef> jetv_gen_wnu;   
   if(genJetsWithNuH.isValid()){
     std::vector<int> partonFlavour_AK4;
     std::vector<int> hadronFlavour_AK4;
     for (auto jets_iter = genJetsWithNuH->begin(); jets_iter != genJetsWithNuH->end(); ++jets_iter) {                                                                                                
       reco::GenJetRef jref (genJetsWithNuH, jets_iter - genJetsWithNuH->begin());                                                                                                                   
       jetv_gen_wnu.push_back(jref);                                                                                                                                                              
     }
      
     for (const reco::GenJet & jet : *genJetsWithNuH){
         bool matched = false;
	 for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *genJetsFlavourH) {
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
         //sort(jetv_gen_wnu.begin(), jetv_gen_wnu.end(), genJetRefSorter);

     for(size_t igen = 0; igen < jetv_gen_wnu.size(); igen++){
       genjet_pt.push_back(jetv_gen_wnu[igen]->pt());
       genjet_eta.push_back(jetv_gen_wnu[igen]->eta()); 
       genjet_phi.push_back(jetv_gen_wnu[igen]->phi());
       genjet_energy.push_back(jetv_gen_wnu[igen]->energy());
       genjet_mass.push_back(jetv_gen_wnu[igen]->mass());
       genjet_hflav.push_back(hadronFlavour_AK4[igen]); //   (*genJetsFlavourH)[edm::RefToBase<reco::Jet>(jetv_gen_wnu[igen])].getHadronFlavour());
       genjet_pflav.push_back(partonFlavour_AK4[igen]);   //     (*genJetsFlavourH)[edm::RefToBase<reco::Jet>(jetv_gen_wnu[igen])].getPartonFlavour());
     }
  }

  std::vector<pat::JetRef> jetv;   
  //JetVetoMap in the event
  Bool_t jetveto_map = false;
  Bool_t jetveto_eep_map = false;
  for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
    pat::JetRef jref(jetsH, jets_iter - jetsH->begin());
    if (jref->pt() > 15.0 && (jref->neutralEmEnergyFraction() + jref->chargedEmEnergyFraction()) < 0.9 && applyJetID(*jref,"tight",usePuppiJets)) 
      {
         TLorentzVector jet4v;
         jet4v.SetPtEtaPhiM(jref->pt(),jref->eta(),jref->phi(),jref->mass());
         for (size_t i = 0; i < muonsColl.size(); i++) {
             TLorentzVector muon4v; 
             muon4v.SetPtEtaPhiM(muonsColl[i].pt(), muonsColl[i].eta(), muonsColl[i].phi(), muonsColl[i].mass());
             if( muon4v.DeltaR(jet4v) > 0.2 )
             {
                  Int_t bin_eta = h_vetomap -> GetXaxis() -> FindBin(jref -> eta());
                  Int_t bin_phi = h_vetomap -> GetYaxis() -> FindBin(jref -> phi());
                  if(h_vetomap ->  GetBinContent(bin_eta, bin_phi) > 0 ) jetveto_map = true;
                  if(h_vetomap_eep -> GetBinContent(bin_eta, bin_phi) > 0 ) jetveto_eep_map = true;
             }
         }//ii   
      }    
  }
  event_veto_map = jetveto_map;
  event_veto_map_eep = jetveto_eep_map;
 
  //if( jetveto_map || jetveto_eep_map ){ std::cout << "The event is jet veto maped : " << std::endl; }
  //else {std::cout << "We are fine" << std::endl;}

  for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {                                                                                                                     
    pat::JetRef jref(jetsH, jets_iter - jetsH->begin());                                                                                                                                            
    if (jref->pt() < jetPtMin and jref->correctedJet("Uncorrected").pt() < jetPtMin) continue;
    if (fabs(jref->eta()) > jetEtaMax) continue;                 
    if (isMC and dumpOnlyJetMatchedToGen and not jref->genJet()) continue;
    jetv.push_back(jref);                                                                                                                                                                           
  }   

  sort(jetv.begin(), jetv.end(), jetRefSorter);
  /*
  double btgSF_temp = 1;
  for (size_t i = 0; i < jetv.size(); i++) {
   if( fabs(jetv[i]->eta()) < 2.5) {
       double btag_disc_val = jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probb")/(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probb") + jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probc") + jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probuds") + jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probg")) ;
       double temp = 1.0;
       std::cout << fabs(jetv[i]->eta()) << " " << jetv[i]->pt() << " " << btag_disc_val << std::endl;
       if(fabs(jetv[i]->hadronFlavour())==5)
       {
          //temp  = reader.eval(BTagEntry::FLAV_B, fabs(jetv[i]->eta()),jetv[i]->pt(),btag_disc_val );
          temp  = reader.eval_auto_bounds("central",BTagEntry::FLAV_B, fabs(jetv[i]->eta()),jetv[i]->pt(),btag_disc_val);
       }
       else if(fabs(jetv[i]->hadronFlavour())==4)
       {
          //temp  = reader.eval(BTagEntry::FLAV_C, fabs(jetv[i]->eta()),jetv[i]->pt(),btag_disc_val );
          temp  = reader.eval_auto_bounds("central",BTagEntry::FLAV_C, fabs(jetv[i]->eta()),jetv[i]->pt(),btag_disc_val);
       }
       else
       {
          temp  = reader.eval_auto_bounds("central",BTagEntry::FLAV_UDSG, fabs(jetv[i]->eta()),jetv[i]->pt(),btag_disc_val);
          //temp  = reader.eval(BTagEntry::FLAV_C, fabs(jetv[i]->eta()),jetv[i]->pt(),btag_disc_val );
       }
       if (temp!=0)  btgSF_temp = btgSF_temp * temp; 
    }
  } //i  
  std::cout << btgSF_temp << std::endl;
  */
  for (size_t i = 0; i < jetv.size(); i++) {
    jet_pt.push_back(jetv[i]->pt());
    jet_eta.push_back(jetv[i]->eta());
    jet_phi.push_back(jetv[i]->phi());
    jet_mass.push_back(jetv[i]->mass());
    jet_energy.push_back(jetv[i]->energy());
    jet_pt_raw.push_back(jetv[i]->correctedJet("Uncorrected").pt());
    jet_mass_raw.push_back(jetv[i]->correctedJet("Uncorrected").mass());
    jet_energy_raw.push_back(jetv[i]->correctedJet("Uncorrected").energy());

    jet_deepjet_probb.push_back(jetv[i]->bDiscriminator("pfDeepFlavourJetTags:probb"));
    jet_deepjet_probbb.push_back(jetv[i]->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    jet_deepjet_problepb.push_back(jetv[i]->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    jet_deepjet_probc.push_back(jetv[i]->bDiscriminator("pfDeepFlavourJetTags:probc"));
    jet_deepjet_probg.push_back(jetv[i]->bDiscriminator("pfDeepFlavourJetTags:probg"));
    jet_deepjet_probuds.push_back(jetv[i]->bDiscriminator("pfDeepFlavourJetTags:probuds"));

    if (usePuppiJets) {
        jet_pnet_probb.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probb"));
        jet_pnet_probc.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probc"));
        jet_pnet_probuds.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probuds"));
        jet_pnet_probg.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probg"));
        jet_pnet_probtauh.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h1p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h2p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup3h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup3h1p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h1p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h2p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum3h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum3h1p"));
        jet_pnet_ptcorr.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptcorr"));
        jet_pnet_ptnu.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptnu"));
        jet_pnet_ptres.push_back( 0.5*(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreshigh") - jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreslow")));
        jet_pnet_probele.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probele"));
        jet_pnet_probmu.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probmu"));
    }
    else {
        jet_pnet_probb.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probb"));
        jet_pnet_probc.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probc"));
        jet_pnet_probuds.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probuds"));
        jet_pnet_probg.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probg"));
        jet_pnet_probtauh.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probtaup1h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h1p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup1h2p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup3h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaup3h1p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h1p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum1h2p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum3h0p")+jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:probtaum3h1p"));
        jet_pnet_ptcorr.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:ptcorr"));
        jet_pnet_ptnu.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:ptnu"));
        jet_pnet_ptres.push_back( 0.5*(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:ptreshigh") - jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4PuppiCentralJetTags:ptreslow")));
        jet_pnet_probele.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probele"));
        jet_pnet_probmu.push_back(jetv[i]->bDiscriminator("pfParticleNetFromMiniAODAK4CHSCentralJetTags:probmu"));
        
    }

    //for(auto & ilabel : pnetDiscriminatorLabels) {
    //  //jet_pnetlast_score[ilabel].push_back(jetv[i]->bDiscriminator(("pfParticleNetAK4LastJetTags:"+ilabel).c_str()));
    //  std::cout<<"PNet label: "<<ilabel.c_str()<<std::endl;
    //  jet_pnetlast_score[ilabel].push_back(jetv[i]->bDiscriminator(ilabel.c_str()));
   // }
 
    if (not usePuppiJets){
      jet_breg_corr.push_back(jetv[i]->userFloat("bjetRegressionNN:corr"));
      jet_breg_res.push_back(jetv[i]->userFloat("bjetRegressionNN:res"));
      jet_creg_corr.push_back(jetv[i]->userFloat("cjetRegressionNN:corr"));
      jet_creg_res.push_back(jetv[i]->userFloat("cjetRegressionNN:res"));
    }

    int jetid = 0;
    if(applyJetID(*jetv[i],"tightLepVeto",usePuppiJets)) jetid += 1;
    jet_id.push_back(jetid);

    //storing jet veto maps for all the jets...................
    bool jetveto_map_ = false;
    bool jetveto_eep_map_ = false; 
    
    if ( (jetv[i]->neutralEmEnergyFraction() + jetv[i]->chargedEmEnergyFraction()) < 0.9 && applyJetID(*jetv[i],"tight",usePuppiJets))
      {
         TLorentzVector jet4v;
         jet4v.SetPtEtaPhiM(jetv[i]->pt(),jetv[i]->eta(),jetv[i]->phi(),jetv[i]->mass());
         for (size_t k = 0; k < muonsColl.size(); k++) {
             TLorentzVector muon4v;
             muon4v.SetPtEtaPhiM(muonsColl[k].pt(), muonsColl[k].eta(), muonsColl[k].phi(), muonsColl[k].mass());
             if( muon4v.DeltaR(jet4v) > 0.2 )
             {
                  Int_t bin_eta = h_vetomap -> GetXaxis() -> FindBin(jetv[i] -> eta());
                  Int_t bin_phi = h_vetomap -> GetYaxis() -> FindBin(jetv[i] -> phi());
                  if(h_vetomap ->  GetBinContent(bin_eta, bin_phi) > 0 ) jetveto_map_ = true;
                  if(h_vetomap_eep -> GetBinContent(bin_eta, bin_phi) > 0 ) jetveto_eep_map_ = true;
             }
         }//ii   
      }
    
    jet_veto.push_back(jetveto_map_);
    jet_veto_eep.push_back(jetveto_eep_map_);
    
    int jetpuid = 0;
    if(applyPileupJetID(*jetv[i],"loose"))  jetpuid += 1;
    if(applyPileupJetID(*jetv[i],"medium")) jetpuid += 2;
    if(applyPileupJetID(*jetv[i],"tight"))  jetpuid += 4;
    jet_puid.push_back(jetpuid);

    // Energy fractions
    jet_chf.push_back(jetv[i]->chargedHadronEnergyFraction());
    jet_nhf.push_back(jetv[i]->neutralHadronEnergyFraction());
    jet_elf.push_back(jetv[i]->electronEnergyFraction());
    jet_phf.push_back(jetv[i]->photonEnergyFraction());
    jet_muf.push_back(jetv[i]->muonEnergyFraction());

    // PF components
    jet_ncand.push_back(jetv[i]->chargedMultiplicity()+jetv[i]->neutralMultiplicity());
    jet_nch.push_back(jetv[i]->chargedHadronMultiplicity());
    jet_nnh.push_back(jetv[i]->neutralHadronMultiplicity());
    jet_nel.push_back(jetv[i]->electronMultiplicity());
    jet_nph.push_back(jetv[i]->photonMultiplicity());
    jet_nmu.push_back(jetv[i]->muonMultiplicity());
    jet_hflav.push_back(jetv[i]->hadronFlavour());
    jet_pflav.push_back(jetv[i]->partonFlavour());
    jet_nbhad.push_back(jetv[i]->jetFlavourInfo().getbHadrons().size());
    jet_nchad.push_back(jetv[i]->jetFlavourInfo().getcHadrons().size());
    // Matching with gen-jets
    int pos_matched = -1;
    float minDR = dRJetGenMatch;
    for(size_t igen = 0; igen < jetv_gen.size(); igen++){
      if(reco::deltaR(jetv_gen[igen]->p4(),jetv[i]->p4()) < minDR){
        pos_matched = igen;
        minDR = reco::deltaR(jetv_gen[igen]->p4(),jetv[i]->p4());
      }
    }

    if(pos_matched >= 0){
      jet_genmatch_pt.push_back(jetv_gen[pos_matched]->pt());
      jet_genmatch_eta.push_back(jetv_gen[pos_matched]->eta());
      jet_genmatch_phi.push_back(jetv_gen[pos_matched]->phi());
      jet_genmatch_mass.push_back(jetv_gen[pos_matched]->mass());
      jet_genmatch_energy.push_back(jetv_gen[pos_matched]->energy());
    }
    else{
      jet_genmatch_pt.push_back(0);
      jet_genmatch_eta.push_back(0);
      jet_genmatch_phi.push_back(0);
      jet_genmatch_mass.push_back(0);
      jet_genmatch_energy.push_back(0);
     }    
    //////////  
    pos_matched = -1;
    minDR = dRJetGenMatch;
    for(size_t igen = 0; igen < jetv_gen_wnu.size(); igen++){
      if(reco::deltaR(jetv_gen_wnu[igen]->p4(),jetv[i]->p4()) < minDR){
        pos_matched = igen;
        minDR = reco::deltaR(jetv_gen_wnu[igen]->p4(),jetv[i]->p4());
      }
    }

    if(pos_matched >= 0){
      jet_genmatch_wnu_pt.push_back(jetv_gen_wnu[pos_matched]->pt());
      jet_genmatch_wnu_eta.push_back(jetv_gen_wnu[pos_matched]->eta());
      jet_genmatch_wnu_phi.push_back(jetv_gen_wnu[pos_matched]->phi());
      jet_genmatch_wnu_mass.push_back(jetv_gen_wnu[pos_matched]->mass());
      jet_genmatch_wnu_energy.push_back(jetv_gen_wnu[pos_matched]->energy());
    }
    else{
      jet_genmatch_wnu_pt.push_back(0);
      jet_genmatch_wnu_eta.push_back(0);
      jet_genmatch_wnu_phi.push_back(0);
      jet_genmatch_wnu_mass.push_back(0);
      jet_genmatch_wnu_energy.push_back(0);
    }

    if(isMC)
    {
    double jecpt = jetv[i]->correctedJet("Uncorrected").pt();
    double JEC_up = 0.0; double JEC_dn = 0.0;
    for(int isrc =0 ; isrc<njecmcmx; isrc++)
       {
        double sup = 1.0 ;
        if((isrc>0)&&(isrc<=nsrc))
         {
            JetCorrectionUncertainty *jecUnc = vsrc[isrc-1];
            jecUnc->setJetEta(jetv[i]->eta());
            jecUnc->setJetPt(jecpt);
            sup += jecUnc->getUncertainty(true);
            if(isrc==1){ JEC_up = sup; }
         }
        else if(isrc>nsrc)
         {
             JetCorrectionUncertainty *jecUnc = vsrc[isrc-1-nsrc];
             jecUnc->setJetEta(jetv[i]->eta());
             jecUnc->setJetPt(jecpt);
             sup -= jecUnc->getUncertainty(false);
             if(isrc==2){ JEC_dn = sup; }
         }
        }
     jet_jec_up.push_back(JEC_up);
     jet_jec_dn.push_back(JEC_dn); 
     //std::cout << jetv[i]->pt() << " " << jetv[i]->correctedJet("Uncorrected").pt() << " " <<  jetv[i]->pt()* JEC_up << "  " << jetv[i]->pt()* JEC_dn << std::endl;
    }	    
    else
    {
     jet_jec_up.push_back(1.0);
     jet_jec_dn.push_back(1.0); 
    }
    njet++;
  }
  //std::cout << "end" << std::endl;  
  // AK8 jets
  std::vector<pat::JetRef> jetAK8v;   
  for (auto jetsAK8_iter = jetsAK8H->begin(); jetsAK8_iter != jetsAK8H->end(); ++jetsAK8_iter) {                                                                                                                     
    pat::JetRef jref(jetsAK8H, jetsAK8_iter - jetsAK8H->begin());                                                                                                                                            
    if (jref->pt() < jetAK8PtMin and jref->correctedJet("Uncorrected").pt() < jetAK8PtMin) continue;
    if (fabs(jref->eta()) > jetEtaMax) continue;                 
    if (isMC and dumpOnlyJetMatchedToGen and not jref->genJet()) continue;
    jetAK8v.push_back(jref);
                                             
  }    
  sort(jetAK8v.begin(), jetAK8v.end(), jetRefSorter);
  //std::cout << "start" << std::endl; 
  for (size_t i = 0; i < jetAK8v.size(); i++) {
    jetAK8_pt.push_back(jetAK8v[i]->pt());
    jetAK8_eta.push_back(jetAK8v[i]->eta());
    jetAK8_phi.push_back(jetAK8v[i]->phi());
    jetAK8_mass.push_back(jetAK8v[i]->mass());
    jetAK8_pt_raw.push_back(jetAK8v[i]->correctedJet("Uncorrected").pt());
    jetAK8_mass_raw.push_back(jetAK8v[i]->correctedJet("Uncorrected").mass());
        
    jetAK8_pnet_probHbb.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHbb"));
    jetAK8_pnet_probHcc.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHcc"));
    jetAK8_pnet_probHqq.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHqq"));
    jetAK8_pnet_probHgg.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHgg"));
    jetAK8_pnet_probHtt.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHtt"));
    jetAK8_pnet_probHtm.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHtm"));
    jetAK8_pnet_probHte.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHte"));
    jetAK8_pnet_probQCD2HF.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD2hf"));
    jetAK8_pnet_probQCD1HF.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD1hf"));
    jetAK8_pnet_probQCD0HF.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD0hf"));
    jetAK8_pnet_HbbVsQCD.push_back( jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHbb") / (jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHbb") + jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD2hf") + jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD1hf") + jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD0hf"))  );

    jetAK8_pnet_mass.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:masscorr")*jetAK8v[i]->mass());
    jetAK8_pnet_corr.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:masscorr")); 
    int jetid = 0;
    //std::cout << "starting jet id" << std::endl;
    if(jetAK8v[i]->isPFJet())
    {
      if(applyJetID(*jetAK8v[i],"tightLepVeto",true)) jetid += 1;
    }
    //std::cout << "ending jet id" << std::endl;
    jetAK8_id.push_back(jetid);
    if(jetAK8v[i]->isPFJet())
    {
    jetAK8_ncand.push_back(jetAK8v[i]->chargedHadronMultiplicity()+jetAK8v[i]->neutralHadronMultiplicity()+jetAK8v[i]->electronMultiplicity()+jetAK8v[i]->photonMultiplicity()+jetAK8v[i]->muonMultiplicity());   
    jetAK8_nbhad.push_back(jetAK8v[i]->jetFlavourInfo().getbHadrons().size());
    jetAK8_nchad.push_back(jetAK8v[i]->jetFlavourInfo().getcHadrons().size());
    jetAK8_hflav.push_back(jetAK8v[i]->hadronFlavour());
    jetAK8_pflav.push_back(jetAK8v[i]->partonFlavour());
    }
    else
    { 
       jetAK8_ncand.push_back(-10);
       jetAK8_nbhad.push_back(-10);
       jetAK8_nchad.push_back(-10);
       jetAK8_hflav.push_back(-10);
       jetAK8_pflav.push_back(-10);
    }
    TLorentzVector jetAK8_softdrop_4V;
    TLorentzVector jetAK8_softdrop_4V_raw;
    pat::JetPtrCollection subjets;
    if(jetAK8v[i]->nSubjetCollections() > 0){
    subjets = jetAK8v[i]->subjets(subJetCollectionName);
    for(size_t isub = 0; isub < subjets.size(); isub++){
             jetAK8_softdrop_subjet_pt.push_back(subjets.at(isub)->pt());
             jetAK8_softdrop_subjet_pt_raw.push_back(subjets.at(isub)->correctedJet("Uncorrected").pt());
             jetAK8_softdrop_subjet_eta.push_back(subjets.at(isub)->eta());
             jetAK8_softdrop_subjet_phi.push_back(subjets.at(isub)->phi());
             jetAK8_softdrop_subjet_mass.push_back(subjets.at(isub)->mass());
             jetAK8_softdrop_subjet_mass_raw.push_back(subjets.at(isub)->correctedJet("Uncorrected").mass());
             TLorentzVector subjet_4v, subjet_raw_4v;
             subjet_4v.SetPtEtaPhiM(subjets.at(isub)->pt(),subjets.at(isub)->eta(),subjets.at(isub)->phi(),subjets.at(isub)->mass());
             subjet_raw_4v.SetPtEtaPhiM(subjets.at(isub)->correctedJet("Uncorrected").pt(),subjets.at(isub)->correctedJet("Uncorrected").eta(),
                                   subjets.at(isub)->correctedJet("Uncorrected").phi(),subjets.at(isub)->correctedJet("Uncorrected").mass());
             jetAK8_softdrop_4V += subjet_4v;
             jetAK8_softdrop_4V_raw += subjet_raw_4v;
           }
    }
    jetAK8_softdrop_pt.push_back(jetAK8_softdrop_4V.Pt());
    jetAK8_softdrop_eta.push_back(jetAK8_softdrop_4V.Eta());
    jetAK8_softdrop_phi.push_back(jetAK8_softdrop_4V.Phi());
    jetAK8_softdrop_mass.push_back(jetAK8_softdrop_4V.M());
    jetAK8_softdrop_pt_raw.push_back(jetAK8_softdrop_4V_raw.Pt());
    jetAK8_softdrop_mass_raw.push_back(jetAK8_softdrop_4V_raw.M());
    //old score storage
    //
    jetAK8_legacy_pnet_mass.push_back(jetAK8v[i]->bDiscriminator("pfParticleNetMassRegressionJetTags:mass"));
    jetAK8_legacy_pnet_probHbb.push_back(jetAK8v[i]->bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probXbb"));
    jetAK8_legacy_pnet_QCD.push_back(  jetAK8v[i]->bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDbb") +  jetAK8v[i]->bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDcc") + jetAK8v[i]->bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDb") + jetAK8v[i]->bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDc") + jetAK8v[i]->bDiscriminator("pfMassDecorrelatedParticleNetJetTags:probQCDothers")  );

    //std::cout << jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8DiscriminatorsJetTags:HbbvsQCD") << " " << jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHbb")/ ( jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probHbb") + jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD2hf") + jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD1hf") + jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:probQCD0hf")  ) <<  std::endl;
    if( (jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8DiscriminatorsJetTags:HbbvsQCD") > 0.5) and (jetAK8v[i]->bDiscriminator("pfParticleNetFromMiniAODAK8JetTags:masscorr")*jetAK8v[i]->mass())>50  )
     {
       njetak8++;
     }
  } 
  //std::cout << "end" << std::endl;
  //Trigger level informations
  //std::cout << "start" << std::endl;
  const edm::TriggerNames &names_ = iEvent.triggerNames(*triggerResultsH);
  for(pat::TriggerObjectStandAlone obj : *triggerObjects)
  {
    obj.unpackPathNames(names_);
    //L1
    if(obj.collection()=="hltGtStage2Digis:Jet:HLT")
      {
        L1jet_pt.push_back(obj.pt());
        L1jet_eta.push_back(obj.eta());
        L1jet_phi.push_back(obj.phi());
        L1jet_en.push_back(obj.energy());
        nL1jet++;
      }
     //Calo
     if(obj.collection()=="hltAK4CaloJetsCorrectedIDPassed::HLT")
      {
        Calojet_pt.push_back(obj.pt());
        Calojet_eta.push_back(obj.eta());
        Calojet_phi.push_back(obj.phi());
        Calojet_en.push_back(obj.energy());
        nCalojet++;
      }
     //Pixel 
     if(obj.collection()=="hltAK4PixelOnlyPFJetsTightIDCorrected::HLT") 
      {
        Pxljet_pt.push_back(obj.pt());
        Pxljet_eta.push_back(obj.eta());
        Pxljet_phi.push_back(obj.phi());
        Pxljet_en.push_back(obj.energy());
        nPxljet++;
      }
      //L3-HLT
      if(obj.collection()=="hltAK4PFJetsTightIDCorrected::HLT")
        {
          L3jet_pt.push_back(obj.pt());
          L3jet_eta.push_back(obj.eta());
          L3jet_phi.push_back(obj.phi());
          L3jet_en.push_back(obj.energy());
          nL3jet++;
        }
       //AK8 calo object collection
       if(obj.collection()=="hltAK8CaloJetsCorrectedIDPassed::HLT")
          {
            CaloAK8jet_pt.push_back(obj.pt());
            CaloAK8jet_eta.push_back(obj.eta());
            CaloAK8jet_phi.push_back(obj.phi());
            CaloAK8jet_en.push_back(obj.energy());
            nCaloAK8jet++;
          } 
       if(obj.collection()=="hltAK8PFJets250SoftDropMass40::HLT")
          {
            L3AK8jet1_pt.push_back(obj.pt());
            L3AK8jet1_eta.push_back(obj.eta());
            L3AK8jet1_phi.push_back(obj.phi());
            L3AK8jet1_en.push_back(obj.energy());
            nL3AK8jet1++;
          } 
  }

  /// MET
  met      = metH->front().corPt();
  met_phi  = metH->front().corPhi();
  met_signif = metH->front().significance();
  met_cov_00 = metH->front().getSignificanceMatrix()(0,0);
  met_cov_01 = metH->front().getSignificanceMatrix()(0,1);
  met_cov_10 = metH->front().getSignificanceMatrix()(1,0);
  met_cov_11 = metH->front().getSignificanceMatrix()(1,1);

  for(auto const & keyval: menu->getAlgorithmMap())
       {
          if(keyval.second.getName() == "L1_HTT280er") idx_L1_HTT280er = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_QuadJet60er2p5") idx_L1_QuadJet60er2p5 = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT320er") idx_L1_HTT320er = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT360er") idx_L1_HTT360er = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT400er") idx_L1_HTT400er = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT450er") idx_L1_HTT450er = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT280er_QuadJet_70_55_40_35_er2p5") idx_L1_HTT280er_QuadJet_70_55_40_35_er2p5 = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3") idx_L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3") idx_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_Mu6_HTT240er") idx_L1_Mu6_HTT240er = keyval.second.getIndex();
          if(keyval.second.getName() == "L1_SingleJet60") idx_L1_SingleJet60 = keyval.second.getIndex();
       }

      Bool_t isL1_HTT280er = false;
      Bool_t isL1_QuadJet60er2p5 = false;
      Bool_t isL1_HTT320er = false;
      Bool_t isL1_HTT360er = false;
      Bool_t isL1_HTT400er = false;
      Bool_t isL1_HTT450er = false;
      Bool_t isL1_HTT280er_QuadJet_70_55_40_35_er2p5 = false;
      Bool_t isL1_HTT280er_QuadJet_70_55_40_40_er2p5 = false;
      Bool_t isL1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 = false;
      Bool_t isL1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = false;
      Bool_t isL1_Mu6_HTT240er = false;
      Bool_t isL1_SingleJet60 = false;

      if(l1GtHandle.isValid()){
         int ibx = 0;
         for(auto itr = l1GtHandle->begin(ibx); itr != l1GtHandle->end(ibx); ++itr){
          if(itr->getAlgoDecisionFinal(idx_L1_HTT280er)){ isL1_HTT280er = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_QuadJet60er2p5)){ isL1_QuadJet60er2p5 = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT320er)){isL1_HTT320er = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT360er)){isL1_HTT360er = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT400er)){isL1_HTT400er = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT450er)){isL1_HTT450er = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT280er_QuadJet_70_55_40_35_er2p5)){isL1_HTT280er_QuadJet_70_55_40_35_er2p5 = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT280er_QuadJet_70_55_40_40_er2p5)){isL1_HTT280er_QuadJet_70_55_40_40_er2p5 = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3)){ isL1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_Mu6_HTT240er)){isL1_Mu6_HTT240er = true;}
          if(itr->getAlgoDecisionFinal(idx_L1_SingleJet60)){isL1_SingleJet60 = true;}
         }
       }
      L1_HTT280er  =  isL1_HTT280er; L1_QuadJet60er2p5 = isL1_QuadJet60er2p5; L1_HTT320er = isL1_HTT320er;
      L1_HTT360er = isL1_HTT360er; L1_HTT400er = isL1_HTT400er; L1_HTT450er = isL1_HTT450er;
      L1_HTT280er_QuadJet_70_55_40_35_er2p5 = isL1_HTT280er_QuadJet_70_55_40_35_er2p5;
      L1_HTT280er_QuadJet_70_55_40_40_er2p5 = isL1_HTT280er_QuadJet_70_55_40_40_er2p5;
      L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 = isL1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3;
      L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 = isL1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3;
      L1_Mu6_HTT240er = isL1_Mu6_HTT240er; L1_SingleJet60 = isL1_SingleJet60;

  //std::cout << "end" << std::endl;
  if(isMC)
  {
    tree1->Fill();
  }      
   
  if(njet > 3 or njetak8>0)
  {
     tree->Fill();
  }
  
  //**************************************************************************//
  // This portion is only for the sorting of jets with b-tag cond.            //
  //**************************************************************************//  
  
  /*
  Float_t leadbtag_score = -999.9;
  Float_t subbtag_score = -999.9;
  Float_t temp_leadbtag_score = -999.9;
  Float_t temp_subbtag_score = -999.9;
  if(njet > 1 )
  {
     // Lead b-tag
     for (int i = 0; i < njet; i++)
     {
       temp_leadbtag_score = jet_pnet_probb.at(i) / ( jet_pnet_probb.at(i) + jet_pnet_probc.at(i) + jet_pnet_probuds.at(i) + jet_pnet_probg.at(i) );
       if( temp_leadbtag_score > leadbtag_score )
         {
            leadbtag_score = temp_leadbtag_score;
         }
      }
      //sublead b-tag
      for (int i = 0; i < njet; i++)
      {
        temp_subbtag_score = jet_pnet_probb.at(i) / ( jet_pnet_probb.at(i) + jet_pnet_probc.at(i) + jet_pnet_probuds.at(i) + jet_pnet_probg.at(i) );
        if( temp_subbtag_score > subbtag_score && temp_subbtag_score < leadbtag_score )
      {
            subbtag_score = temp_subbtag_score;
        }
      }
      if(leadbtag_score > 0.2605 && subbtag_score > 0.2605  )
      {
           tree->Fill();
      }
  }
  */
}

void TreeMaker::beginJob() {

  f_vetomap          = TFile::Open(jetvetomapfile_.c_str()); //afs/cern.ch/work/m/mukherje/private/RUN3/Analysis_UCSD/CMSSW_13_0_7/src/hh4banalysis/TreeMaker/test/Summer22EE_23Sep2023_RunEFG_v1.root");
  h_vetomap = (TH2D * ) f_vetomap -> Get("jetvetomap");
  h_vetomap_eep = (TH2D * ) f_vetomap -> Get("jetvetomap");

  for (int isrc = 0; isrc < nsrc; isrc++)
  {     
        const char *name = srcnames[isrc];
        JetCorrectorParameters *p = new JetCorrectorParameters(mJECUncFile_.c_str(), name) ;
        JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
        vsrc.push_back(unc);
  }

  // Access the TFileService
  edm::Service<TFileService> fs;  

  //Tree to store basic information without any cut
  tree1 = fs->make<TTree>("allinfo","allinfo");
  
  tree1->Branch("event", &event, "event/i");
  tree1->Branch("run", &run, "run/i");
  tree1->Branch("lumi", &lumi, "lumi/i");
  tree1->Branch("npv", &npv, "npv/i");
  //tree1->Branch("rho", &rho, "rho/F");
  //tree1->Branch("rho_trk", &rho_trk, "rho_trk/F");
  //tree1->Branch("rho_calo", &rho_calo, "rho_calo/F");

  if(isMC){
      tree1->Branch("xsec", &xsec, "xsec/F");
      tree1->Branch("wgt" , &wgt , "wgt/F");
      tree1->Branch("putrue", &putrue, "putrue/i");   
      //Total weight PS 
      tree1->Branch("wgt_isr_up", &wgt_isr_up, "wgt_isr_up/F");
      tree1->Branch("wgt_isr_dn", &wgt_isr_dn, "wgt_isr_dn/F");
      tree1->Branch("wgt_fsr_up", &wgt_fsr_up, "wgt_fsr_up/F");
      tree1->Branch("wgt_fsr_dn", &wgt_fsr_dn, "wgt_fsr_dn/F");
      //Total weight QCDscale
      tree1->Branch("wgt_scl_0", &wgt_scl_0, "wgt_scl_0/F");
      tree1->Branch("wgt_scl_1", &wgt_scl_1, "wgt_scl_1/F");
      tree1->Branch("wgt_scl_2", &wgt_scl_2, "wgt_scl_2/F");
      tree1->Branch("wgt_scl_3", &wgt_scl_3, "wgt_scl_3/F");
      tree1->Branch("wgt_scl_4", &wgt_scl_4, "wgt_scl_4/F");
      tree1->Branch("wgt_scl_5", &wgt_scl_5, "wgt_scl_5/F");
      tree1->Branch("wgt_scl_6", &wgt_scl_6, "wgt_scl_6/F");
      tree1->Branch("wgt_scl_7", &wgt_scl_7, "wgt_scl_7/F");
      tree1->Branch("wgt_scl_8", &wgt_scl_8, "wgt_scl_8/F");
  }
  
  tree = fs->make<TTree>("tree","tree");

  tree->Branch("event", &event, "event/i");
  tree->Branch("run", &run, "run/i");
  tree->Branch("lumi", &lumi, "lumi/i");
  tree->Branch("xsec", &xsec, "xsec/F");
  tree->Branch("wgt" , &wgt , "wgt/F");
  tree->Branch("flags", &flags, "flags/i");
  tree->Branch("putrue", &putrue, "putrue/i");
  tree->Branch("rho", &rho, "rho/F");
  tree->Branch("rho_trk", &rho_trk, "rho_trk/F");
  tree->Branch("rho_calo", &rho_calo, "rho_calo/F");
  tree->Branch("npv", &npv, "npv/i");
  tree->Branch("nsv", &nsv, "nsv/i");

  tree -> Branch("event_veto_map", &event_veto_map);
  tree -> Branch("event_veto_map_eep", &event_veto_map_eep);

  tree -> Branch("Flag_goodVertices", &Flag_goodVertices_);
  tree -> Branch("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter_);
  tree -> Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter_);
  tree -> Branch("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter_);
  tree -> Branch("Flag_BadPFMuonDzFilter", &Flag_BadPFMuonDzFilter_);
  tree -> Branch("Flag_hfNoisyHitsFilter", &Flag_hfNoisyHitsFilter_);
  tree -> Branch("Flag_eeBadScFilter", &Flag_eeBadScFilter_);
  tree -> Branch("Flag_ecalBadCalibFilter_", &Flag_ecalBadCalibFilter_);


  // GEN particles
  if(isMC){
    tree->Branch("gen_particle_pt","std::vector<float>", &gen_particle_pt);
    tree->Branch("gen_particle_eta","std::vector<float>", &gen_particle_eta);
    tree->Branch("gen_particle_phi","std::vector<float>", &gen_particle_phi);
    tree->Branch("gen_particle_mass","std::vector<float>", &gen_particle_mass);
    tree->Branch("gen_particle_id","std::vector<int>", &gen_particle_id);
    tree->Branch("gen_particle_status","std::vector<unsigned int>", &gen_particle_status);
    tree->Branch("gen_particle_daughters_id","std::vector<int>", &gen_particle_daughters_id);
    tree->Branch("gen_particle_daughters_igen","std::vector<unsigned int>", &gen_particle_daughters_igen);
    tree->Branch("gen_particle_daughters_pt","std::vector<float>", &gen_particle_daughters_pt);
    tree->Branch("gen_particle_daughters_eta","std::vector<float>", &gen_particle_daughters_eta);
    tree->Branch("gen_particle_daughters_phi","std::vector<float>", &gen_particle_daughters_phi);
    tree->Branch("gen_particle_daughters_mass","std::vector<float>", &gen_particle_daughters_mass);
    tree->Branch("gen_particle_daughters_status","std::vector<unsigned int>", &gen_particle_daughters_status);
    tree->Branch("gen_particle_daughters_charge","std::vector<int>", &gen_particle_daughters_charge);

    tree->Branch("lhe_particle_pt", "std::vector<float>", &lhe_particle_pt);
    tree->Branch("lhe_particle_eta", "std::vector<float>", &lhe_particle_eta);
    tree->Branch("lhe_particle_phi", "std::vector<float>", &lhe_particle_phi);
    tree->Branch("lhe_particle_mass", "std::vector<float>", &lhe_particle_mass);
    tree->Branch("lhe_particle_id", "std::vector<int>", &lhe_particle_id);
    tree->Branch("lhe_particle_status", "std::vector<unsigned int>", &lhe_particle_status);

  }

  tree->Branch("muon_pt", "std::vector<float>", &muon_pt);
  tree->Branch("muon_eta", "std::vector<float>", &muon_eta);
  tree->Branch("muon_phi", "std::vector<float>", &muon_phi);
  tree->Branch("muon_mass", "std::vector<float>", &muon_mass);
  tree->Branch("muon_id", "std::vector<unsigned int>" , &muon_id);
  tree->Branch("muon_iso", "std::vector<unsigned int>" , &muon_iso);
  tree->Branch("muon_charge", "std::vector<int>" , &muon_charge);
  tree->Branch("muon_d0", "std::vector<float>" , &muon_d0);
  tree->Branch("muon_dz", "std::vector<float>" , &muon_dz);

  tree->Branch("electron_pt", "std::vector<float>", &electron_pt);
  tree->Branch("electron_eta", "std::vector<float>", &electron_eta);
  tree->Branch("electron_phi", "std::vector<float>", &electron_phi);
  tree->Branch("electron_mass", "std::vector<float>", &electron_mass);
  tree->Branch("electron_id", "std::vector<unsigned int>" , &electron_id);
  tree->Branch("electron_idscore", "std::vector<float>" , &electron_idscore);
  tree->Branch("electron_charge", "std::vector<int>" , &electron_charge);
  tree->Branch("electron_d0", "std::vector<float>" , &electron_d0);
  tree->Branch("electron_dz", "std::vector<float>" , &electron_dz);
 
  /*
  tree->Branch("photon_pt", "std::vector<float>", &photon_pt);
  tree->Branch("photon_pt_corr", "std::vector<float>", &photon_pt_corr);
  tree->Branch("photon_eta", "std::vector<float>", &photon_eta);
  tree->Branch("photon_phi", "std::vector<float>", &photon_phi);
  tree->Branch("photon_mass", "std::vector<float>", &photon_mass);
  tree->Branch("photon_id", "std::vector<unsigned int>" , &photon_id);
  tree->Branch("photon_idscore", "std::vector<float>" , &photon_idscore);

  tree->Branch("tau_pt", "std::vector<float>" , &tau_pt);
  tree->Branch("tau_eta", "std::vector<float>" , &tau_eta);
  tree->Branch("tau_phi", "std::vector<float>" , &tau_phi);
  tree->Branch("tau_mass", "std::vector<float>" , &tau_mass);
  tree->Branch("tau_dxy", "std::vector<float>" , &tau_dxy);
  tree->Branch("tau_dz", "std::vector<float>" , &tau_dz);
  tree->Branch("tau_charge", "std::vector<int>" , &tau_charge);
  tree->Branch("tau_decaymode", "std::vector<unsigned int>" , &tau_decaymode);
  tree->Branch("tau_idjet_wp", "std::vector<unsigned int>" , &tau_idjet_wp);
  tree->Branch("tau_idmu_wp", "std::vector<unsigned int>" , &tau_idmu_wp);
  tree->Branch("tau_idele_wp", "std::vector<unsigned int>" , &tau_idele_wp);
  tree->Branch("tau_idjet", "std::vector<float>" , &tau_idjet);
  tree->Branch("tau_idele", "std::vector<float>" , &tau_idele);
  tree->Branch("tau_idmu", "std::vector<float>" , &tau_idmu);
  if(isMC){
    tree->Branch("tau_genmatch_pt", "std::vector<float>" , &tau_genmatch_pt);
    tree->Branch("tau_genmatch_eta", "std::vector<float>" , &tau_genmatch_eta);
    tree->Branch("tau_genmatch_phi", "std::vector<float>" , &tau_genmatch_phi);
    tree->Branch("tau_genmatch_mass", "std::vector<float>" , &tau_genmatch_mass);
  }
  */
  tree->Branch("met",&met,"met/F");
  tree->Branch("met_phi", &met_phi, "met_phi/F");
  tree->Branch("met_signif", &met_signif, "met_signif/F");
  tree->Branch("met_cov_00", &met_cov_00, "met_cov_00/F");
  tree->Branch("met_cov_01", &met_cov_01, "met_cov_01/F");
  tree->Branch("met_cov_10", &met_cov_10, "met_cov_10/F");
  tree->Branch("met_cov_11", &met_cov_11, "met_cov_11/F");

  tree->Branch("trigger_hlt_path","std::vector<std::string>",&trigger_hlt_path);
  tree->Branch("trigger_hlt_pass","std::vector<unsigned int>",&trigger_hlt_pass);


  tree->Branch("L1_QuadJet60er2p5" , &L1_QuadJet60er2p5 );
  tree->Branch("L1_HTT280er" , &L1_HTT280er );
  tree->Branch("L1_HTT320er" , &L1_HTT320er );
  tree->Branch("L1_HTT360er" , &L1_HTT360er );
  tree->Branch("L1_HTT400er" , &L1_HTT400er );
  tree->Branch("L1_HTT450er" , &L1_HTT450er );
  tree->Branch("L1_HTT280er_QuadJet_70_55_40_35_er2p5" ,       &L1_HTT280er_QuadJet_70_55_40_35_er2p5 );
  tree->Branch("L1_HTT280er_QuadJet_70_55_40_40_er2p5" ,       &L1_HTT280er_QuadJet_70_55_40_40_er2p5 );
  tree->Branch("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3" , &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3 );
  tree->Branch("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3" , &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3 );
  tree->Branch("L1_Mu6_HTT240er" , &L1_Mu6_HTT240er );
  tree->Branch("L1_SingleJet60",   &L1_SingleJet60);
  
  tree->Branch("njet", &njet, "njet/I");
  tree->Branch("njetak8", &njetak8, "njetak8/I");
  tree->Branch("jet_pt", "std::vector<float>" , &jet_pt);
  tree->Branch("jet_eta", "std::vector<float>" , &jet_eta);
  tree->Branch("jet_phi", "std::vector<float>" , &jet_phi);
  tree->Branch("jet_mass", "std::vector<float>" , &jet_mass);
  tree->Branch("jet_energy", "std::vector<float>" , &jet_energy);
  tree->Branch("jet_pt_raw", "std::vector<float>" , &jet_pt_raw);
  tree->Branch("jet_mass_raw", "std::vector<float>" , &jet_mass_raw);
  tree->Branch("jet_energy_raw", "std::vector<float>" , &jet_energy_raw);
  tree->Branch("jet_chf", "std::vector<float>" , &jet_chf);
  tree->Branch("jet_nhf", "std::vector<float>" , &jet_nhf);
  tree->Branch("jet_elf", "std::vector<float>" , &jet_elf);
  tree->Branch("jet_phf", "std::vector<float>" , &jet_phf);
  tree->Branch("jet_muf", "std::vector<float>" , &jet_muf);
  tree->Branch("jet_deepjet_probb", "std::vector<float>" , &jet_deepjet_probb);
  tree->Branch("jet_deepjet_probbb", "std::vector<float>" , &jet_deepjet_probbb);
  tree->Branch("jet_deepjet_probc", "std::vector<float>" , &jet_deepjet_probc);
  tree->Branch("jet_deepjet_problepb", "std::vector<float>" , &jet_deepjet_problepb);
  tree->Branch("jet_deepjet_probg", "std::vector<float>" , &jet_deepjet_probg);
  tree->Branch("jet_deepjet_probuds", "std::vector<float>" , &jet_deepjet_probuds);
  tree->Branch("jet_pnet_probb", "std::vector<float>" , &jet_pnet_probb);
  tree->Branch("jet_pnet_probc", "std::vector<float>" , &jet_pnet_probc);
  tree->Branch("jet_pnet_probuds", "std::vector<float>" , &jet_pnet_probuds);
  tree->Branch("jet_pnet_probg", "std::vector<float>" , &jet_pnet_probg);
  tree->Branch("jet_pnet_probtauh", "std::vector<float>" , &jet_pnet_probtauh);
  tree->Branch("jet_pnet_ptcorr", "std::vector<float>" , &jet_pnet_ptcorr);
  tree->Branch("jet_pnet_ptnu", "std::vector<float>" , &jet_pnet_ptnu);
  tree->Branch("jet_pnet_ptres", "std::vector<float>" , &jet_pnet_ptres);
  tree->Branch("jet_pnet_probele", "std::vector<float>" , &jet_pnet_probele);
  tree->Branch("jet_pnet_probmu", "std::vector<float>" , &jet_pnet_probmu);
  tree->Branch("jet_breg_corr","std::vector<float>",&jet_breg_corr);
  tree->Branch("jet_breg_res","std::vector<float>",&jet_breg_res);
  tree->Branch("jet_creg_corr","std::vector<float>",&jet_creg_corr);
  tree->Branch("jet_creg_res","std::vector<float>",&jet_creg_res);

  tree->Branch("jet_veto", "std::vector<bool>" , &jet_veto);
  tree->Branch("jet_veto_epp", "std::vector<bool>" , &jet_veto_eep);
  tree->Branch("jet_id", "std::vector<unsigned int>" , &jet_id);
  tree->Branch("jet_puid", "std::vector<unsigned int>" , &jet_puid);
  tree->Branch("jet_ncand", "std::vector<unsigned int>" , &jet_ncand);
  tree->Branch("jet_nch", "std::vector<unsigned int>" , &jet_nch);
  tree->Branch("jet_nnh", "std::vector<unsigned int>" , &jet_nnh);
  tree->Branch("jet_nel", "std::vector<unsigned int>" , &jet_nel);
  tree->Branch("jet_nph", "std::vector<unsigned int>" , &jet_nph);
  tree->Branch("jet_nmu", "std::vector<unsigned int>" , &jet_nmu);
  tree->Branch("jet_hflav", "std::vector<unsigned int>" , &jet_hflav);
  tree->Branch("jet_pflav", "std::vector<int>" , &jet_pflav);
  tree->Branch("jet_nbhad", "std::vector<unsigned int>" , &jet_nbhad);
  tree->Branch("jet_nchad", "std::vector<unsigned int>" , &jet_nchad);  
  if(isMC){
    tree->Branch("jet_jec_up","std::vector<float>", &jet_jec_up);
    tree->Branch("jet_jec_dn","std::vector<float>", &jet_jec_dn);    
    tree->Branch("jet_genmatch_pt","std::vector<float>" , &jet_genmatch_pt);
    tree->Branch("jet_genmatch_eta","std::vector<float>" , &jet_genmatch_eta);
    tree->Branch("jet_genmatch_phi","std::vector<float>" , &jet_genmatch_phi);
    tree->Branch("jet_genmatch_mass","std::vector<float>" , &jet_genmatch_mass);
    tree->Branch("jet_genmatch_energy","std::vector<float>" , &jet_genmatch_energy);
    tree->Branch("jet_genmatch_wnu_pt","std::vector<float>" , &jet_genmatch_wnu_pt);
    tree->Branch("jet_genmatch_wnu_eta","std::vector<float>" , &jet_genmatch_wnu_eta);
    tree->Branch("jet_genmatch_wnu_phi","std::vector<float>" , &jet_genmatch_wnu_phi);
    tree->Branch("jet_genmatch_wnu_mass","std::vector<float>" , &jet_genmatch_wnu_mass);
    tree->Branch("jet_genmatch_wnu_energy","std::vector<float>" , &jet_genmatch_wnu_energy);
    tree->Branch("genjet_pt","std::vector<float>" , &genjet_pt);
    tree->Branch("genjet_eta","std::vector<float>" , &genjet_eta);
    tree->Branch("genjet_phi","std::vector<float>" , &genjet_phi);
    tree->Branch("genjet_mass","std::vector<float>" , &genjet_mass);
    tree->Branch("genjet_energy","std::vector<float>" , &genjet_energy);
    tree->Branch("genjet_hflav","std::vector<int>" , &genjet_hflav);
    tree->Branch("genjet_pflav","std::vector<int>" , &genjet_pflav);
  }
  tree->Branch("jetAK8_pt", "std::vector<float>" , &jetAK8_pt);
  tree->Branch("jetAK8_eta", "std::vector<float>" , &jetAK8_eta);
  tree->Branch("jetAK8_phi", "std::vector<float>" , &jetAK8_phi);
  tree->Branch("jetAK8_mass", "std::vector<float>" , &jetAK8_mass);
  tree->Branch("jetAK8_pt_raw", "std::vector<float>" , &jetAK8_pt_raw);
  tree->Branch("jetAK8_mass_raw", "std::vector<float>" , &jetAK8_mass_raw);
  tree->Branch("jetAK8_pnet_probHbb", "std::vector<float>" , &jetAK8_pnet_probHbb);
  tree->Branch("jetAK8_pnet_probHcc", "std::vector<float>" , &jetAK8_pnet_probHcc);
  tree->Branch("jetAK8_pnet_probHqq", "std::vector<float>" , &jetAK8_pnet_probHqq);
  tree->Branch("jetAK8_pnet_probHgg", "std::vector<float>" , &jetAK8_pnet_probHgg);
  tree->Branch("jetAK8_pnet_probHtt", "std::vector<float>" , &jetAK8_pnet_probHtt);
  tree->Branch("jetAK8_pnet_probHtm", "std::vector<float>" , &jetAK8_pnet_probHtm);
  tree->Branch("jetAK8_pnet_probHte", "std::vector<float>" , &jetAK8_pnet_probHte);
  tree->Branch("jetAK8_pnet_probQCD2HF", "std::vector<float>" , &jetAK8_pnet_probQCD2HF);
  tree->Branch("jetAK8_pnet_probQCD1HF", "std::vector<float>" , &jetAK8_pnet_probQCD1HF);
  tree->Branch("jetAK8_pnet_probQCD0HF", "std::vector<float>" , &jetAK8_pnet_probQCD0HF);
  tree->Branch("jetAK8_pnet_HbbVsQCD", "std::vector<float>" , &jetAK8_pnet_HbbVsQCD);
  tree->Branch("jetAK8_pnet_mass", "std::vector<float>" , &jetAK8_pnet_mass);
  tree->Branch("jetAK8_pnet_corr", "std::vector<float>" , &jetAK8_pnet_corr);
  tree->Branch("jetAK8_id", "std::vector<unsigned int>" , &jetAK8_id);
  tree->Branch("jetAK8_ncand", "std::vector<unsigned int>" , &jetAK8_ncand);
  tree->Branch("jetAK8_hflav", "std::vector<unsigned int>" , &jetAK8_hflav);
  tree->Branch("jetAK8_pflav", "std::vector<int>" , &jetAK8_pflav);
  tree->Branch("jetAK8_nbhad", "std::vector<unsigned int>" , &jetAK8_nbhad);
  tree->Branch("jetAK8_nchad", "std::vector<unsigned int>" , &jetAK8_nchad);
  tree->Branch("jetAK8_softdrop_pt", "std::vector<float>" , &jetAK8_softdrop_pt);
  tree->Branch("jetAK8_softdrop_eta", "std::vector<float>" , &jetAK8_softdrop_eta);
  tree->Branch("jetAK8_softdrop_phi", "std::vector<float>" , &jetAK8_softdrop_phi);
  tree->Branch("jetAK8_softdrop_mass", "std::vector<float>" , &jetAK8_softdrop_mass);
  tree->Branch("jetAK8_softdrop_pt_raw", "std::vector<float>" , &jetAK8_softdrop_pt_raw);
  tree->Branch("jetAK8_softdrop_mass_raw", "std::vector<float>" , &jetAK8_softdrop_mass_raw);
  tree->Branch("jetAK8_legacy_pnet_mass", "std::vector<float>", &jetAK8_legacy_pnet_mass);
  tree->Branch("jetAK8_legacy_pnet_probHbb", "std::vector<float>", &jetAK8_legacy_pnet_probHbb);
  tree->Branch("jetAK8_legacy_pnet_QCD", "std::vector<float>", &jetAK8_legacy_pnet_QCD);

  tree-> Branch ("nL1jet", &nL1jet);
  tree-> Branch ("L1jet_pt", &L1jet_pt);
  tree-> Branch ("L1jet_eta", &L1jet_eta);
  tree-> Branch ("L1jet_phi", &L1jet_phi);
  tree-> Branch ("L1jet_en", &L1jet_en);

  tree-> Branch ("nCalojet", &nCalojet);
  tree-> Branch ("Calojet_pt", &Calojet_pt);
  tree-> Branch ("Calojet_eta", &Calojet_eta);
  tree-> Branch ("Calojet_phi", &Calojet_phi);
  tree-> Branch ("Calojet_en", &Calojet_en);

  tree-> Branch ("nCalojet", &nCalojet);
  tree-> Branch ("Calojet_pt", &Calojet_pt);
  tree-> Branch ("Calojet_eta", &Calojet_eta);
  tree-> Branch ("Calojet_phi", &Calojet_phi);
  tree-> Branch ("Calojet_en", &Calojet_en);

  tree-> Branch ("nPxljet",    &nPxljet);
  tree-> Branch ("Pxljet_pt",  &Pxljet_pt);
  tree-> Branch ("Pxljet_eta", &Pxljet_eta);
  tree-> Branch ("Pxljet_phi", &Pxljet_phi);
  tree-> Branch ("Pxljet_en",  &Pxljet_en);  

  tree-> Branch ("nL3jet", &nL3jet);
  tree-> Branch ("L3jet_pt", &L3jet_pt);
  tree-> Branch ("L3jet_eta", &L3jet_eta);
  tree-> Branch ("L3jet_phi", &L3jet_phi);
  tree-> Branch ("L3jet_en", &L3jet_en);

  tree-> Branch ("nCaloAK8jet",    &nCaloAK8jet);
  tree-> Branch ("CaloAK8jet_pt",  &CaloAK8jet_pt);
  tree-> Branch ("CaloAK8jet_eta", &CaloAK8jet_eta);
  tree-> Branch ("CaloAK8jet_phi", &CaloAK8jet_phi);
  tree-> Branch ("CaloAK8jet_en",  &CaloAK8jet_en);

  tree-> Branch ("nL3AK8jet1", &nL3AK8jet1);
  tree-> Branch ("L3AK8jet1_pt", &L3AK8jet1_pt);
  tree-> Branch ("L3AK8jet1_eta", &L3AK8jet1_eta);
  tree-> Branch ("L3AK8jet1_phi", &L3AK8jet1_phi);
  tree-> Branch ("L3AK8jet1_en", &L3AK8jet1_en);

  if(isMC)
  {

    //LHE
    tree->Branch("LHE_weight",&LHE_weight, "LHE_weight/D");
    tree->Branch("nLHEScaleWeights",&nLHEScaleWeights, "nLHEScaleWeights/I");
    tree->Branch("LHEScaleWeights",LHEScaleWeights,"LHEScaleWeights[nLHEScaleWeights]/F");
    tree->Branch("nLHEPDFWeights",&nLHEPDFWeights, "nLHEPDFWeights/I");
    tree->Branch("LHEPDFWeights",LHEPDFWeights,"LHEPDFWeights[nLHEPDFWeights]/F");
    tree->Branch("nLHEAlpsWeights",&nLHEAlpsWeights, "nLHEAlpsWeights/I");
    tree->Branch("LHEAlpsWeights",LHEAlpsWeights,"LHEAlpsWeights[nLHEAlpsWeights]/F");
    tree->Branch("nLHEPSWeights",&nLHEPSWeights, "nLHEPSWeights/I");
    tree->Branch("LHEPSWeights",LHEPSWeights,"LHEPSWeights[nLHEPSWeights]/F");

    //GEN 
    tree->Branch("Generator_weight", &Generator_weight, "Generator_weight/D") ;
    tree->Branch("Generator_qscale",&Generator_qscale,"Generator_qscale/F");
    tree->Branch("Generator_x1",&Generator_x1,"Generator_x1/F");
    tree->Branch("Generator_x2",&Generator_x2,"Generator_x2/F");
    tree->Branch("Generator_xpdf1",&Generator_xpdf1,"Generator_xpdf1/F");
    tree->Branch("Generator_xpdf2",&Generator_xpdf2,"Generator_xpdf2/F");
    tree->Branch("Generator_id1",&Generator_id1,"Generator_id1/I");
    tree->Branch("Generator_id2",&Generator_id2,"Generator_id2/I");
    tree->Branch("Generator_scalePDF",&Generator_scalePDF,"Generator_scalePDF/F");
  }
}  

void TreeMaker::endJob() {
}

void TreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {


  HLTConfigProvider fltrConfig;
  bool flag = false;
  fltrConfig.init(iRun, iSetup, filterResultsTag.process(), flag);
  
  // MET filter Paths
  filterPathsMap.clear();
  filterPathsVector.clear();
  filterPathsVector.push_back("Flag_goodVertices");
  filterPathsVector.push_back("Flag_globalSuperTightHalo2016Filter");
  filterPathsVector.push_back("Flag_HBHENoiseFilter");
  filterPathsVector.push_back("Flag_HBHENoiseIsoFilter");
  filterPathsVector.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  filterPathsVector.push_back("Flag_BadPFMuonFilter");
  filterPathsVector.push_back("Flag_BadChargedCandidateFilter");
  
  for (size_t i = 0; i < filterPathsVector.size(); i++) {
    filterPathsMap[filterPathsVector[i]] = -1;
  }
  
  for(size_t i = 0; i < filterPathsVector.size(); i++){
    TPRegexp pattern(filterPathsVector[i]);
    for(size_t j = 0; j < fltrConfig.triggerNames().size(); j++){
      std::string pathName = fltrConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	    filterPathsMap[filterPathsVector[i]] = j;
      }
    }
  }

  // trigger names
  triggerPathsMap.clear();
  HLTConfigProvider hltConfig;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), flag);
  for(size_t i = 0; i < hltConfig.triggerNames().size(); i++){
    std::string pathName = hltConfig.triggerNames()[i];
    if((not (TString(pathName).Contains("ParticleNet"))) && (not (TString(pathName).Contains("PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet"))) && (not (TString(pathName).Contains("AK8PFJet400_SoftDropMass40"))) && (not (TString(pathName).Contains("PNet")))  ) continue;
     triggerPathsMap[pathName] = i;
  }
}

void TreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void TreeMaker::initializeBranches(){

  event  = 0;
  run    = 0;
  lumi   = 0;
  putrue = 0;
  flags  = 0;
  wgt    = 0.;
  rho    = 0.;
 
  trigger_hlt_path.clear();
  trigger_hlt_pass.clear();

  lhe_particle_pt.clear();
  lhe_particle_eta.clear();
  lhe_particle_phi.clear();
  lhe_particle_mass.clear();
  lhe_particle_id.clear();
  lhe_particle_status.clear();
 
  gen_particle_pt.clear();
  gen_particle_eta.clear();
  gen_particle_phi.clear();
  gen_particle_mass.clear();
  gen_particle_id.clear();
  gen_particle_status.clear();
  gen_particle_daughters_id.clear();
  gen_particle_daughters_igen.clear();
  gen_particle_daughters_pt.clear();
  gen_particle_daughters_eta.clear();
  gen_particle_daughters_phi.clear();
  gen_particle_daughters_mass.clear();
  gen_particle_daughters_status.clear();
  gen_particle_daughters_charge.clear();

  muon_pt.clear();
  muon_eta.clear();
  muon_phi.clear();
  muon_mass.clear();
  muon_id.clear();
  muon_charge.clear();
  muon_iso.clear();
  muon_d0.clear();
  muon_dz.clear();

  electron_pt.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_mass.clear();
  electron_id.clear();
  electron_idscore.clear();
  electron_charge.clear();
  electron_d0.clear();
  electron_dz.clear();
  /*
  photon_pt.clear();
  photon_pt_corr.clear();
  photon_eta.clear();
  photon_phi.clear();
  photon_mass.clear();
  photon_id.clear();
  photon_idscore.clear();

  tau_pt.clear();
  tau_eta.clear();
  tau_phi.clear();
  tau_mass.clear();
  tau_dxy.clear();
  tau_dz.clear();
  tau_decaymode.clear();
  tau_idjet.clear();
  tau_idele.clear();
  tau_idmu.clear();
  tau_idjet_wp.clear();
  tau_idmu_wp.clear();
  tau_idele_wp.clear();
  tau_charge.clear();
  tau_genmatch_pt.clear();
  tau_genmatch_eta.clear();
  tau_genmatch_phi.clear();
  tau_genmatch_mass.clear();
  tau_genmatch_decaymode.clear();
  */
  met = 0.;
  met_phi = 0;
  met_signif = 0;
  met_cov_00 = 0;
  met_cov_10 = 0;
  met_cov_01 = 0;
  met_cov_11 = 0;

  npv = 0;
  nsv = 0;
  
  njet = 0;
  njetak8 = 0;
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_mass.clear();
  jet_energy.clear();
  jet_pt_raw.clear();
  jet_mass_raw.clear();
  jet_energy_raw.clear();
  jet_chf.clear();
  jet_nhf.clear();
  jet_elf.clear();
  jet_phf.clear();
  jet_muf.clear();
  jet_deepjet_probb.clear();
  jet_deepjet_probbb.clear();
  jet_deepjet_problepb.clear();
  jet_deepjet_probc.clear();
  jet_deepjet_probg.clear();
  jet_deepjet_probuds.clear();
  jet_pnet_probb.clear();
  jet_pnet_probc.clear();
  jet_pnet_probg.clear();
  jet_pnet_probuds.clear();
  jet_pnet_probtauh.clear();
  jet_pnet_probele.clear();
  jet_pnet_probmu.clear();
  jet_pnet_ptcorr.clear();
  jet_pnet_ptnu.clear();
  jet_pnet_ptres.clear();
  jet_jec_up.clear();
  jet_jec_dn.clear();
  for(auto & imap : jet_pnetlast_score)
    imap.second.clear();
  jet_breg_corr.clear();
  jet_breg_res.clear();
  jet_creg_corr.clear();
  jet_creg_res.clear();
  jet_veto.clear();
  jet_veto_eep.clear();
  jet_id.clear();
  jet_puid.clear();
  jet_ncand.clear();
  jet_nch.clear();
  jet_nnh.clear();
  jet_nel.clear();
  jet_nph.clear();
  jet_nmu.clear();
  jet_hflav.clear();
  jet_pflav.clear();
  jet_nbhad.clear();
  jet_nchad.clear();
  jet_genmatch_pt.clear();
  jet_genmatch_eta.clear();
  jet_genmatch_phi.clear();
  jet_genmatch_mass.clear();
  jet_genmatch_energy.clear();
  jet_genmatch_wnu_pt.clear();
  jet_genmatch_wnu_eta.clear();
  jet_genmatch_wnu_phi.clear();
  jet_genmatch_wnu_mass.clear();
  jet_genmatch_wnu_energy.clear();
  genjet_pt.clear();
  genjet_eta.clear();
  genjet_phi.clear();
  genjet_energy.clear();
  genjet_mass.clear();
  genjet_hflav.clear();
  genjet_pflav.clear();
  jetAK8_pt.clear();
  jetAK8_eta.clear();
  jetAK8_phi.clear();
  jetAK8_mass.clear();
  jetAK8_pt_raw.clear();
  jetAK8_mass_raw.clear();
  jetAK8_pnet_probHbb.clear();
  jetAK8_pnet_probHcc.clear();
  jetAK8_pnet_probHqq.clear();
  jetAK8_pnet_probHgg.clear();
  jetAK8_pnet_probHtt.clear();
  jetAK8_pnet_probHtm.clear();
  jetAK8_pnet_probHte.clear();
  jetAK8_pnet_probQCD2HF.clear();
  jetAK8_pnet_probQCD1HF.clear();
  jetAK8_pnet_probQCD0HF.clear();
  jetAK8_pnet_HbbVsQCD.clear();
  jetAK8_pnet_mass.clear();
  jetAK8_pnet_corr.clear();
  jetAK8_id.clear();
  jetAK8_ncand.clear();
  jetAK8_hflav.clear();
  jetAK8_pflav.clear();
  jetAK8_nbhad.clear();
  jetAK8_nchad.clear();
  jetAK8_softdrop_pt.clear();
  jetAK8_softdrop_eta.clear();
  jetAK8_softdrop_phi.clear();
  jetAK8_softdrop_mass.clear();
  jetAK8_softdrop_pt_raw.clear();
  jetAK8_softdrop_mass_raw.clear();
  jetAK8_softdrop_subjet_pt.clear();
  jetAK8_softdrop_subjet_pt_raw.clear();
  jetAK8_softdrop_subjet_eta.clear();
  jetAK8_softdrop_subjet_phi.clear();
  jetAK8_softdrop_subjet_mass.clear();
  jetAK8_softdrop_subjet_mass_raw.clear();
  jetAK8_legacy_pnet_mass.clear();
  jetAK8_legacy_pnet_probHbb.clear();
  jetAK8_legacy_pnet_QCD.clear();
  nL1jet=0; nCalojet=0; nL3jet = 0; nPxljet = 0;
  L1jet_pt.clear(); L1jet_eta.clear(); L1jet_phi.clear(); L1jet_en.clear();
  Calojet_pt.clear(); Calojet_eta.clear(); Calojet_phi.clear(); Calojet_en.clear();
  Pxljet_pt.clear(); Pxljet_eta.clear(); Pxljet_phi.clear(); Pxljet_en.clear();
  L3jet_pt.clear(); L3jet_eta.clear(); L3jet_phi.clear(); L3jet_en.clear();
  
  nCaloAK8jet = 0; CaloAK8jet_pt.clear(); CaloAK8jet_eta.clear(); CaloAK8jet_phi.clear(); CaloAK8jet_en.clear();
  nL3AK8jet1 = 0;  L3AK8jet1_pt.clear(); L3AK8jet1_eta.clear(); L3AK8jet1_phi.clear(); L3AK8jet1_en.clear();
  //
  L1_QuadJet60er2p5= false;
  L1_HTT280er= false;
  L1_HTT320er= false;
  L1_HTT360er= false;
  L1_HTT400er= false;
  L1_HTT450er= false;
  L1_HTT280er_QuadJet_70_55_40_35_er2p5= false;
  L1_HTT280er_QuadJet_70_55_40_40_er2p5= false;
  L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3= false;
  L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3= false;
  L1_Mu6_HTT240er= false;
  L1_SingleJet60= false; 
   
}

void TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// to apply jet ID: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL                                                                                                                             
bool TreeMaker::applyJetID(const pat::Jet & jet, const std::string & level, const bool & isPuppi){

  if(level != "tight" and level != "tightLepVeto")
    return true;


  double eta  = jet.eta();
  double nhf  = jet.neutralHadronEnergyFraction();
  double nemf = jet.neutralEmEnergyFraction();
  double chf  = jet.chargedHadronEnergyFraction();
  double muf  = jet.muonEnergyFraction();
  double cemf = jet.chargedEmEnergyFraction();
  int    np   = jet.chargedMultiplicity()+jet.neutralMultiplicity();
  int    chm  = jet.chargedMultiplicity();


  int jetid  = 0;
  
  if(isPuppi){
    if (fabs(eta) <= 2.6 and (nhf < 0.99 and nemf < 0.90 and np > 1 and chf > 0.01 and chm > 0)) jetid += 1;
    if (fabs(eta) <= 2.6 and (nhf < 0.99 and nemf < 0.90 and np > 1 and chf > 0.01 and chm > 0 and muf < 0.8 and cemf < 0.8)) jetid += 2;
  }
  else{
    if (fabs(eta) <= 2.6 and (nhf < 0.99 and nemf < 0.90 and np > 1 and chf > 0.01 and chm > 0)) jetid += 1;
    if (fabs(eta) <= 2.6 and (nhf < 0.99 and nemf < 0.90 and np > 1 and chf > 0.01 and chm > 0 and muf < 0.8 and cemf < 0.8)) jetid += 2;
  }

  if(level == "tight" and jetid > 0) return true;
  else if(level == "tightLepVeto" and jetid > 1) return true;
  else return false;

}

bool TreeMaker::applyPileupJetID(const pat::Jet & jet, const std::string & level){
  bool passpuid = false;
  if(jet.hasUserInt("pileupJetIdUpdated:fullId")){
    if(level == "loose"  and (bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 0)) or jet.pt() > 50)) passpuid = true;
    if(level == "medium" and (bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 1)) or jet.pt() > 50)) passpuid = true;
    if(level == "tight"  and (bool(jet.userInt("pileupJetIdUpdated:fullId") & (1 << 2)) or jet.pt() > 50)) passpuid = true;
  }
  else if (jet.hasUserInt("pileupJetId:fullId")){
    if(level == "loose"  and (bool(jet.userInt("pileupJetId:fullId") & (1 << 0)) or jet.pt() > 50)) passpuid = true;
    if(level == "medium" and (bool(jet.userInt("pileupJetId:fullId") & (1 << 1)) or jet.pt() > 50)) passpuid = true;
    if(level == "tight"  and (bool(jet.userInt("pileupJetId:fullId") & (1 << 2)) or jet.pt() > 50)) passpuid = true;
  }
  return passpuid;
}


DEFINE_FWK_MODULE(TreeMaker);

