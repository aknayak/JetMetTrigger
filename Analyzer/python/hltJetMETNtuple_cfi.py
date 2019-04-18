import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq

isData = True
#isData = False
runOnData=isData #data/MC switch

def configureJetMetNtuple(process):
  #  Options parsing  #
  ###################
  
  from FWCore.ParameterSet.VarParsing import VarParsing
  
  options = VarParsing('analysis')
  options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')
  options.register('runJets',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Store Jet Variables')
  '''
  process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
  process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
  process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

  process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
  process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
  process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
  '''

  # Electron ID ==========================================================================================

  setupEgammaPostRecoSeq(process,
                         era='2017-Nov17ReReco')  


  ### END Electron ID ====================================================================================


  process.hltJetMetNtuple = cms.EDAnalyzer('HLTJetMETNtupleProducer',
                                           runJets = cms.bool(options.runJets),
                                           runMets = cms.untracked.bool(False),
                                           isData = cms.bool(isData),
                                           PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           MetCollectionTag = cms.InputTag('slimmedMETs'),
                                           applyMETFilters = cms.bool(options.applyMETFilters),
                                           BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
                                           BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
                                           MuonCollectionTag = cms.InputTag('slimmedMuons'),
                                           ElectronCollectionTag = cms.InputTag('slimmedElectrons'),
                                           PFJetCollectionTag = cms.InputTag('slimmedJets'),
                                           HLTPFJetCollectionTag = cms.InputTag('hltAK4PFJetsCorrected'),
                                           HLTCaloJetCollectionTag = cms.InputTag('hltAK4CaloJetsCorrected'),
                                           hltprocess = cms.InputTag('TriggerResults::HLT'),
                                           triggerObjects = cms.InputTag("slimmedPatTrigger"),
                                           triggerPaths = cms.untracked.vstring('HLT_PFMET200_NotCleaned_v',					 
                                                                                'HLT_PFMET200_HBHECleaned_v',
                                                                                'HLT_PFMET250_HBHECleaned_v',
                                                                                'HLT_PFMET300_HBHECleaned_v',
                                                                                'HLT_PFMET200_HBHE_BeamHaloCleaned_v',
                                                                                'HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v', 
                                                                                #'HLT_PFMET110_PFMHT110_IDTight_v', 
                                                                                'HLT_Ele25_WPTight_Gsf_v',
                                                                                'HLT_Ele25_eta2p1_WPTight_Gsf_v',
                                                                                'HLT_Ele27_WPTight_Gsf_v',
                                                                                'HLT_Ele27_eta2p1_WPTight_Gsf_v',
                                                                                'HLT_Ele30_WPTight_Gsf_v',
                                                                                'HLT_Ele30_eta2p1_WPTight_Gsf_v',
                                                                                'HLT_Ele32_WPTight_Gsf_v',
                                                                                'HLT_Ele32_eta2p1_WPTight_Gsf_v',
                                                                                'HLT_Ele35_WPTight_Gsf_v',
                                                                                'HLT_IsoMu24_v',
                                                                                'HLT_IsoMu24_eta2p1_v',
                                                                                'HLT_IsoMu27_v',
                                                                                'HLT_IsoMu30_v',
                                                                                'HLT_PFJet40_v',
                                                                                'HLT_PFJet60_v',
                                                                                'HLT_PFJet80_v',
                                                                                'HLT_PFJet140_v',
                                                                                'HLT_PFJet200_v',
                                                                                'HLT_PFJet260_v',
                                                                                'HLT_PFJet320_v',
                                                                                'HLT_PFJet400_v',
                                                                                'HLT_PFJet450_v',
                                                                                'HLT_PFJet500_v',
                                                                                'HLT_PFJet550_v',
                                                                                'HLT_PFJetFwd40_v', 
                                                                                'HLT_PFJetFwd60_v',
                                                                                'HLT_PFJetFwd80_v',
                                                                                'HLT_PFJetFwd140_v',
                                                                                'HLT_PFJetFwd200_v',
                                                                                'HLT_PFJetFwd260_v',
                                                                                'HLT_PFJetFwd320_v',
                                                                                'HLT_PFJetFwd400_v',
                                                                                'HLT_PFJetFwd450_v',
                                                                              'HLT_PFJetFwd500_v'
                                                                              )
                                             )
    
  process.JetMetNtupleSequence = cms.Sequence(process.egammaPostRecoSeq*
                                              #process.BadChargedCandidateFilter*
                                              #process.BadPFMuonFilter* 
                                              process.hltJetMetNtuple
  )
  
