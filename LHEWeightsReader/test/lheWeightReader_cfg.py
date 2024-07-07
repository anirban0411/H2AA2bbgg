import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')
options.parseArguments()

process = cms.Process("Analyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames  = cms.untracked.vstring(options.inputFiles)
    #fileNames = cms.untracked.vstring(
    #    '/store/mc/RunIIFall17MiniAODv2/PairVectorLQ_InclusiveDecay_M-1000_TuneCP5_13TeV-madgraph-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/110000/8E879E65-16AE-E811-A17C-44A84225CDA4.root'
    #)
)

#process.load("LHEWeightsReader.lheWeightReader_cfi")
process.analyzer = cms.EDAnalyzer('LHEWeightsReader',
	lheProduct = cms.InputTag('externalLHEProducer')
)

process.p = cms.Path(process.analyzer)
