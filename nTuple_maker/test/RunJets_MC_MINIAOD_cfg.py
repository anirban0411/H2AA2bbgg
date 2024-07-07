# -*- coding: utf-8 -*-

import FWCore.ParameterSet.Config as cms
import sys

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff import *
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.patSequences_cff import *
from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.slimming.metFilterPaths_cff import *


## Modified version of jetToolBox from https://github.com/cms-jet/jetToolbox
## Options for PUMethod: Puppi, CS, SK, CHS

# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms

process = cms.Process("Combined")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Conditions
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.load('RecoJets.JetProducers.TrackJetParameters_cfi')
process.load('RecoJets.JetProducers.PileupJetIDParams_cfi')

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


from RecoJets.Configuration.GenJetParticles_cff import *


process.GlobalTag.globaltag = "106X_upgrade2018_realistic_v15_L1v1"
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

from PhysicsTools.PatAlgos.tools.coreTools import *
process.load("PhysicsTools.PatAlgos.patSequences_cff")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Input
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inFiles = cms.untracked.vstring(
'root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18MiniAODv2/WHToAA_AToBB_AToGG_M-35_TuneCP5_13TeV_madgraph_pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/40000/134A4EF0-09F3-404F-8728-A63DF0E2B24D.root'
#'root://xrootd-cms.infn.it//store/mc/RunIISummer20UL18MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/00000/00435839-09ED-FA45-A52D-BEA0C9838058.root'
#'root://xrootd-cms.infn.it//store/user/abala/HAHM_ZHToAATo2B2G_MA-20GeV_TuneCP5_13TeV-madgraph_pythia8/HAHM_ZHToAATo2B2G_MA-20GeV_TuneCP5_13TeV-madgraph_pythia8_MINIAODSIM/220201_124019/0000/MINIAOD_7.root'
# 'file://afs/cern.ch/work/m/mukherje/public/134A4EF0-09F3-404F-8728-A63DF0E2B24D.root'
   )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2000))

#process.firstEvent = cms.untracked.PSet(input = cms.untracked.int32(1000))
process.source = cms.Source("PoolSource", fileNames = inFiles )

RealData = bool(False)
FastSIM = bool(False)

process.p = cms.Path()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#! Services
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.load('CommonTools.UtilAlgos.TFileService_cfi')

process.TFileService = cms.Service("TFileService",
fileName = cms.string('hist.root')             #largest data till April5,2016 
)

process.patJets.addTagInfos = True

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# Electron IDs for AOD/MINIAOD
switchOnVIDElectronIdProducer(process, DataFormat.MiniAOD)

# define which IDs to produce
el_id_modules = [
##    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff"
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff",
    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff"
##    "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff"
]
# Add them to the VID producer
for iModule in el_id_modules:

	setupAllVIDIdsInModule(process, iModule, setupVIDElectronSelection)

switchOnVIDPhotonIdProducer(process, DataFormat.MiniAOD)

pho_id_modules = [
	"RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff",
	"RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff"
]

for iModule in pho_id_modules:

	setupAllVIDIdsInModule(process, iModule, setupVIDPhotonSelection)

# E-Gamma scale and smearing corrections 
from Analysis.NTuplizer.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,era='2018-Prompt')  

#DNN-based AK8 jet tagger names
deep_discriminators = ["pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD",
                       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD",
                       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD",
		       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD",
		       "pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight",
		       "pfMassDecorrelatedParticleNetJetTags:probXbb",
		       "pfMassDecorrelatedParticleNetJetTags:probQCDbb",
		       "pfMassDecorrelatedParticleNetJetTags:probQCDcc",
		       "pfMassDecorrelatedParticleNetJetTags:probQCDb",
		       "pfMassDecorrelatedParticleNetJetTags:probQCDc",
		       "pfMassDecorrelatedParticleNetJetTags:probQCDothers",
		       "pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XbbvsQCD",
		       "pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XccvsQCD",
                       "pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XqqvsQCD",
		       "pfParticleNetDiscriminatorsJetTags:TvsQCD",
		       "pfParticleNetDiscriminatorsJetTags:WvsQCD",
		       "pfParticleNetDiscriminatorsJetTags:ZvsQCD",
		       "pfDeepCSVJetTags:probb", "pfDeepCSVJetTags:probbb"               
]
#from RecoBTag.MXNet.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
deep_discriminators += pfParticleNetJetTagsAll

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask, addToProcessAndTask
def _addProcessAndTask(proc, label, module):
	task = getPatAlgosToolsTask(proc)
    	addToProcessAndTask(label, module, proc, task)

#New PF-collection:
def producePF(process) :
        from CommonTools.PileupAlgos.Puppi_cff import puppi
        puppi.useExistingWeights = True
        puppi.candName = "packedPFCandidates"
        puppi.vertexName = "offlineSlimmedPrimaryVertices"
        from LeptonLessPFProducer_cff import leptonLessPFProducer
        _addProcessAndTask(process,"leptonLessPFProducer",leptonLessPFProducer.clone())
        _addProcessAndTask(process,"leptonLesspuppi",puppi.clone(candName = cms.InputTag('leptonLessPFProducer')))

## Define AK8 jet sequence:

def ak8JetSequences(process):
	
# jetToolbox

	JETCorrLevels = ['L1FastJet','L2Relative','L3Absolute']
	if RealData: JETCorrLevels.append('L2L3Residual') 

	subjetBTagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb']

	from jetToolbox_cff import jetToolbox

	jetToolbox( process, 'ak8', 'ak8JetSubs', 'noOutput', PUMethod='Puppi', dataTier="miniAOD", runOnMC=False, JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels, addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels, addNsub=True, bTagDiscriminators=['None'], subjetBTagDiscriminators=subjetBTagDiscriminators)

	jetToolbox( process, 'ak8', 'ak8JetSubs', 'noOutput', PUMethod='Puppi', dataTier="miniAOD", runOnMC=False, postFix='NoLep', newPFCollection=True, nameNewPFCollection='leptonLesspuppi', JETCorrPayload='AK8PFPuppi', JETCorrLevels=JETCorrLevels, addSoftDrop=True, addSoftDropSubjets=True, subJETCorrPayload='AK4PFPuppi', subJETCorrLevels=JETCorrLevels, addNsub=True, bTagDiscriminators=['None'], subjetBTagDiscriminators=subjetBTagDiscriminators)

#update collection

	from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

	updateJetCollection(
		process,
		jetSource=cms.InputTag('packedPatJetsAK8PFPuppiSoftDrop'),
		#jetSource=cms.InputTag('selectedPatJetsAK8PFPuppi'),
		pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
		svSource = cms.InputTag('slimmedSecondaryVertices'),
		rParam=0.8,
		jetCorrections = ('AK8PFPuppi', cms.vstring(JETCorrLevels), 'None'),
		btagDiscriminators = deep_discriminators, #pfParticleNetJetTagsAll, #deep_discriminators,
		postfix='AK8wLepWithPuppiDaughters', 
		printWarning = False
	)
	#new collection -> selectedUpdatedPatJetsAK8wLepWithPuppiDaughters

	updateJetCollection(
		process,
		jetSource=cms.InputTag('packedPatJetsAK8PFPuppiNoLepSoftDrop'),
		pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
		svSource = cms.InputTag('slimmedSecondaryVertices'),
		rParam=0.8,
		jetCorrections = ('AK8PFPuppi', cms.vstring(JETCorrLevels), 'None'),
		btagDiscriminators = deep_discriminators,
		postfix='AK8NoLepWithPuppiDaughters',
		printWarning = False
	)
	#new collection -> selectedUpdatedPatJetsAK8NoLepWithPuppiDaughters

def defaultJetSequences(process):
	producePF(process)
	ak8JetSequences(process)

defaultJetSequences(process)
#-> gives jet collection from jetToolbox+updateJetCollection

# Prefire correction #

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
	TheJets = cms.InputTag("slimmedJets"), #"updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs 
	DataEraECAL = cms.string("None"), 
	DataEraMuon = cms.string("20172018"), 
	UseJetEMPt = cms.bool(False),
	PrefiringRateSystematicUnctyECAL = cms.double(0.2),
	PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)


########## b Jet Energy Corrections ################


# Modules from NanoAOD
from PhysicsTools.NanoAOD.jets_cff import bJetVars,bjetNN 
#### update modules parameters
# create variables for regression
process.bJetVars = bJetVars.clone( 
    pvsrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    src = cms.InputTag("slimmedJets"),
    svsrc = cms.InputTag("slimmedSecondaryVertices"),
    gpsrc = cms.InputTag("prunedGenParticles"),
)


# add variables to the jet collection
# (using only bJetVars, which are variables relevant to regression)
process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag("slimmedJets"),
    userFloats = cms.PSet(
         leadTrackPt = cms.InputTag("bJetVars:leadTrackPt"),
         leptonPtRel = cms.InputTag("bJetVars:leptonPtRel"),
         leptonPtRatio = cms.InputTag("bJetVars:leptonPtRatio"),
         leptonPtRelInv = cms.InputTag("bJetVars:leptonPtRelInv"),
         leptonPtRelv0 = cms.InputTag("bJetVars:leptonPtRelv0"),
         leptonPtRatiov0 = cms.InputTag("bJetVars:leptonPtRatiov0"),
         leptonPtRelInvv0 = cms.InputTag("bJetVars:leptonPtRelInvv0"),
         leptonDeltaR = cms.InputTag("bJetVars:leptonDeltaR"),
         leptonPt = cms.InputTag("bJetVars:leptonPt"),
         vtxPt = cms.InputTag("bJetVars:vtxPt"),
         vtxMass = cms.InputTag("bJetVars:vtxMass"),
         vtx3dL = cms.InputTag("bJetVars:vtx3dL"),
         vtx3deL = cms.InputTag("bJetVars:vtx3deL"),
         ptD = cms.InputTag("bJetVars:ptD"),
         genPtwNu = cms.InputTag("bJetVars:genPtwNu"),
#         pileupJetIdDiscriminant = cms.InputTag("pileupJetIdUpdated:fullDiscriminant"),

         ),
    userInts = cms.PSet(
        vtxNtrk = cms.InputTag("bJetVars:vtxNtrk"),
        leptonPdgId = cms.InputTag("bJetVars:leptonPdgId"),
#        pileupJetId = cms.InputTag("pileupJetIdUpdated:fullId"),
    ),
)

# obtain the b jet regression corrections
process.bjetNN = bjetNN.clone(
    src = cms.InputTag("updatedJetsWithUserData"),
)

# add the b jet regression corrections to the jet collection
process.selectedUpdatedPatJetsBRegPileupJetID = cms.EDProducer("PATJetUserDataEmbedder",
     src = cms.InputTag("updatedJetsWithUserData"),
     userFloats = cms.PSet(
         bJetRegCorr = cms.InputTag("bjetNN:corr"),
         bJetRegRes = cms.InputTag("bjetNN:res"),
#         pileupJetIdDiscriminant = cms.InputTag("pileupJetIdUpdated:fullDiscriminant"),
         ),
#     userInts = cms.PSet(
#         pileupJetId = cms.InputTag("pileupJetIdUpdated:fullId"),
#         ),
)

process.BJetRegression=cms.Sequence( 
				process.bJetVars*
				process.updatedJetsWithUserData*
				process.bjetNN*
				process.selectedUpdatedPatJetsBRegPileupJetID 
			    )

# Analyzer #

process.mcjets =  cms.EDAnalyzer('Leptop',

	 Data =  cms.untracked.bool(RealData),
	 MonteCarlo =  cms.untracked.bool(not RealData),
	 FastSIM =  cms.untracked.bool(FastSIM),
         YEAR = cms.untracked.int32(2018),
         UltraLegacy =  cms.untracked.bool(True),                        
	 isReco = cms.untracked.bool(True),
 	 ReRECO = cms.untracked.bool(True),
	 SoftDrop_ON =  cms.untracked.bool(True),
	 add_prefireweights =  cms.untracked.bool(True),
	 store_electron_scalnsmear =  cms.untracked.bool(True),
	 store_electron_addvariab = cms.untracked.bool(False),
         store_fatjet_constituents = cms.untracked.bool(False),
	 Read_btagging_SF = cms.untracked.bool(False),
	 Subtract_Lepton_fromAK8 = cms.untracked.bool(False),
         Subtract_Lepton_fromAK4 = cms.untracked.bool(True),

 	 RootFileName = cms.untracked.string('rootuple.root'),  #largest data till April5,2016
	
#	 PFJetsAK8 = cms.InputTag("slimmedJetsAK8"),
#        PFJetsAK8 = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropPacked","SubJets","Combined"),
#        PFJetsAK8 = cms.InputTag("selectedPatJetsAK8PFPuppi","","Combined"),+ssDecorrelatedParticleNetJetTagsSlimmedJetsAK8
#	 PFJetsAK8 = cms.InputTag("updatedPatJetsTransientCorrectedSlimmedJetsAK8"),#"","PAT"),
	 PFJetsAK8 = cms.InputTag("selectedUpdatedPatJetsAK8NoLepWithPuppiDaughters"),
#	 PFJetsAK8 = cms.InputTag("selectedUpdatedPatJetsAK8wLepWithPuppiDaughters"),
#	 PFJetsAK8 = cms.InputTag("selectedPatJetsAK8PFPuppi"),
	 AK8PtCut = cms.untracked.double(0.),
         AK8GenPtCut = cms.untracked.double(0.),
	 softdropmass  = cms.untracked.string("ak8PFJetsSoftDropMass"),#ak8PFJetsPuppiSoftDropMass"),#('ak8PFJetsPuppiSoftDropMass'),
	 tau1  = cms.untracked.string("NjettinessAK8PuppiNoLep:tau1"),#'NjettinessAK8Puppi:tau1'),
	 tau2  = cms.untracked.string("NjettinessAK8PuppiNoLep:tau2"),#'NjettinessAK8Puppi:tau2'),
	 tau3  = cms.untracked.string("NjettinessAK8PuppiNoLep:tau3"),#'NjettinessAK8Puppi:tau3'),
	 subjets  = cms.untracked.string('SoftDrop'),#("SoftDrop"),#'SoftDropPuppi'),
#	 subjets = cms.untracked.string("slimmedJetsAK8PFPuppiSoftDropPacked"),
	 toptagger_DAK8 = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"),
	 Wtagger_DAK8 = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"),
	 Ztagger_DAK8 = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZvsQCD"),
	 Htagger_DAK8 = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"),
	 bbtagger_DAK8 = cms.untracked.string("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight"),	
	 toptagger_PNet = cms.untracked.string("pfParticleNetDiscriminatorsJetTags:TvsQCD"),
	 Wtagger_PNet = cms.untracked.string("pfParticleNetDiscriminatorsJetTags:WvsQCD"),
	 Ztagger_PNet = cms.untracked.string("pfParticleNetDiscriminatorsJetTags:ZvsQCD"),
	 Xbbtagger_PNet = cms.untracked.string("pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XbbvsQCD"),
	 Xcctagger_PNet = cms.untracked.string("pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XccvsQCD"),
	 Xqqtagger_PNet = cms.untracked.string("pfMassDecorrelatedParticleNetDiscriminatorsJetTags:XqqvsQCD"),
	 QCDtagger_PNet = cms.untracked.string("pfMassDecorrelatedParticleNetDiscriminatorsJetTags:QCD"),

	 PFJetsAK4 = cms.InputTag("slimmedJets"), 
         PFAK4JetsB = cms.InputTag("selectedUpdatedPatJetsBRegPileupJetID"),
	 minjPt = cms.untracked.double(10.),
	 maxEta = cms.untracked.double(3.),

	 GENJetAK8 = cms.InputTag("slimmedGenJetsAK8"),
	 GENJetAK4 = cms.InputTag("slimmedGenJets"),
         minGenPt = cms.untracked.double(10.),      
	 maxGenEta = cms.untracked.double(5.),         
         
	 Muons = cms.InputTag("slimmedMuons"),
	 minmuPt  = cms.untracked.double(10.),
	 EAFile_MuonMiniIso = cms.FileInPath("PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_94X.txt"),

	 Electrons = cms.InputTag("slimmedElectrons"),
	 minePt   = cms.untracked.double(10.),
	 electronID_isowp90        = cms.string('mvaEleID-Fall17-iso-V2-wp90'),
         electronID_noisowp90      = cms.string('mvaEleID-Fall17-noIso-V2-wp90'),
         electronID_isowp80        = cms.string('mvaEleID-Fall17-iso-V2-wp80'),
         electronID_noisowp80      = cms.string('mvaEleID-Fall17-noIso-V2-wp80'),
	 EAFile_EleMiniIso = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
	 EAFile_ElePFIso = cms.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),

	 Photons = cms.InputTag("slimmedPhotons"),
	 mingmPt  = cms.untracked.double(10.),
	 PhoID_FallV2_WP90 = cms.string("mvaPhoID-RunIIFall17-v2-wp90"),
         PhoID_FallV2_WP80 = cms.string("mvaPhoID-RunIIFall17-v2-wp80"),
         PhoID_SpringV1_WP90 = cms.string("mvaPhoID-Spring16-nonTrig-V1-wp90"),
         PhoID_SpringV1_WP80 = cms.string("mvaPhoID-Spring16-nonTrig-V1-wp80"),
	 label_mvaPhoID_FallV2_Value = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRunIIFall17v2Values"),
	 #label_mvaPhoID_FallV2_WP90 = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp90"),
         #label_mvaPhoID_FallV2_WP80 = cms.InputTag("egmPhotonIDs:mvaPhoID-RunIIFall17-v2-wp80"),

	 Taus = cms.InputTag("slimmedTausUpdated"),
	 mintauPt = cms.untracked.double(25.),

	 PFMet = cms.InputTag("slimmedMETs"),
         PuppiMet = cms.InputTag("slimmedMETsPuppi"),
         GENMet  = cms.InputTag("genMetTrue","","SIM"),

	 GenParticles = cms.InputTag("prunedGenParticles"),#("prunedGenParticles"),#("packedGenParticles"),
	 jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
	 #jetFlavourInfos = cms.InputTag("genJetAK8FlavourAssociation"),

	 pfCands = cms.InputTag("packedPFCandidates"),

         Beamspot = cms.InputTag("offlineBeamSpot"),
         PrimaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
         SecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
	 slimmedAddPileupInfo = cms.InputTag("slimmedAddPileupInfo"),
         PFRho = cms.InputTag("fixedGridRhoFastjetAll"),

	 Generator = cms.InputTag("generator"),
         LHEEventProductInputTag = cms.InputTag('externalLHEProducer'),
         GenEventProductInputTag = cms.InputTag('generator'),
	 nPDFsets = cms.untracked.uint32(103),

	 jecL1FastFileAK4          = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.txt'),
         jecL1FastFileAK8          = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L1FastJet_AK8PFPuppi.txt'),
         jecL2RelativeFileAK4      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK4PFchs.txt'),
         jecL2RelativeFileAK8      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2Relative_AK8PFPuppi.txt'),
         jecL3AbsoluteFileAK4      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.txt'),
         jecL3AbsoluteFileAK8      = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L3Absolute_AK8PFPuppi.txt'),
         jecL2L3ResidualFileAK4    = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2L3Residual_AK4PFchs.txt'),
         jecL2L3ResidualFileAK8    = cms.string('Summer19UL18_V5_MC/Summer19UL18_V5_MC_L2L3Residual_AK8PFPuppi.txt'),

	 PtResoFileAK4  = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt'),
         PtResoFileAK8  = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt'),
         PtSFFileAK4 = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt'),
         PtSFFileAK8 = cms.string('Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt'),

         JECUncFileAK4 = cms.string("Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.txt"),
	 JECUncFileAK8 = cms.string("Summer19UL18_V5_MC/Summer19UL18_V5_MC_UncertaintySources_AK8PFPuppi.txt"),

	 BtagSFFile_DeepCSV = cms.string("BtagRecommendation106XUL18/DeepCSV_106XUL18SF_V1p1.csv"),
	 BtagSFFile_DeepFlav = cms.string("BtagRecommendation106XUL18/DeepJet_106XUL18SF_V1p1.csv"),
         BtagSFFile_DeepFlav_Iter = cms.string("BtagRecommendation106XUL18/DeepJet_106XUL18SF.csv"),
	 RochcorFolder = cms.string("roccor.Run2.v5/"),

	 bits = cms.InputTag("TriggerResults","","HLT"),
         prescales = cms.InputTag("patTrigger","","RECO"),
         objects = cms.InputTag("slimmedPatTrigger")
)

#===== MET Filters ==

process.load('RecoMET.METFilters.primaryVertexFilter_cfi')
process.primaryVertexFilter.vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
if not FastSIM:
	process.load('RecoMET.METFilters.globalSuperTightHalo2016Filter_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')
process.HBHENoiseFilterResultProducerNoMinZ = process.HBHENoiseFilterResultProducer.clone(minZeros = cms.int32(99999))
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices") 
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadPFMuonDzFilter_cfi')
process.BadPFMuonDzFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonDzFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonDzFilter.vtx = cms.InputTag("offlineSlimmedPrimaryVertices")
process.BadPFMuonDzFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.EERecHitSource = cms.InputTag('reducedEgamma','reducedEERecHits')

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('eventoutput.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

if FastSIM:
	process.allMetFilterPaths=cms.Sequence(process.primaryVertexFilter*
	process.HBHENoiseFilter*
	process.HBHENoiseIsoFilter*
	process.EcalDeadCellTriggerPrimitiveFilter*
	process.BadPFMuonFilter*
	process.BadPFMuonDzFilter*
	process.ecalBadCalibFilter)
else:
	process.allMetFilterPaths=cms.Sequence(process.primaryVertexFilter*
	process.globalSuperTightHalo2016Filter*
	process.HBHENoiseFilter*
	process.HBHENoiseIsoFilter*
	process.EcalDeadCellTriggerPrimitiveFilter*
	process.BadPFMuonFilter*
	process.BadPFMuonDzFilter*
	process.ecalBadCalibFilter)


import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
updatedTauName = "slimmedTausUpdated"
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = True, updatedTauName = updatedTauName,
	toKeep = [ "2017v2",
        #       "dR0p32017v2", 
        #       "newDM2017v2", 
        "deepTau2017v2p1",
        "againstEle2018"\
        ]
)
tauIdEmbedder.runTauID()

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           isData=False,
                           )
from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True )
runMetCorAndUncFromMiniAOD(process,
                           isData=False,
                           metType="Puppi",
                           postfix="Puppi",
                           jetFlavor="AK4PFPuppi",
                           )
process.puppi.useExistingWeights = True

process.p = cms.Path(process.egmPhotonIDSequence 
 		     *process.HBHENoiseFilterResultProducer*process.HBHENoiseFilterResultProducerNoMinZ
		     *process.allMetFilterPaths
#		     *process.egmGsfElectronIDSequence*
		     *process.rerunMvaIsolationSequence*getattr(process,updatedTauName) #process.slimmedTausUpdated*  # this also works
		     *process.egammaPostRecoSeq
		     *process.prefiringweight
		     *process.fullPatMetSequence 
		     *process.puppiMETSequence
		     *process.fullPatMetSequencePuppi
		     *process.BJetRegression
		     *process.mcjets
		     )

process.schedule = cms.Schedule(process.p)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.options.numberOfThreads=cms.untracked.uint32(4)
#process.options.numberOfStreams=cms.untracked.uint32(0)
