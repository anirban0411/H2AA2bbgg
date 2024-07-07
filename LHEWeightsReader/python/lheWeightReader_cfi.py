import FWCore.ParameterSet.Config as cms

analyzer = cms.EDAnalyzer('LHEWeightsReader',
       lheProduct = cms.InputTag('externalLHEProducer')
)

