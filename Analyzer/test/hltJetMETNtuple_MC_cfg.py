import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList


isData = False
#isData = False
runOnData=isData #data/MC switch

#####################
#  Options parsing  #
#####################

from FWCore.ParameterSet.VarParsing import VarParsing
import os, sys

options = VarParsing('analysis')
options.register('applyMETFilters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')

process = cms.Process("hltJetMET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
if runOnData:
  process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'
else:
  process.GlobalTag.globaltag = '100X_upgrade2018_realistic_v10'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_10_0_2/RelValADDMonoJet_d3MD3_13/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/10000/7E85E4B4-050F-E811-B4CE-0CC47A78A3EC.root',
    '/store/relval/CMSSW_10_0_2/RelValADDMonoJet_d3MD3_13/MINIAODSIM/100X_upgrade2018_realistic_v10-v1/10000/96E029B6-050F-E811-A73E-0CC47A7C3444.root'

    #'/store/relval/CMSSW_10_0_0/RelValADDMonoJet_d3MD3_13/MINIAODSIM/100X_upgrade2018_realistic_v6_mahiON-v1/10000/80FDD09F-8100-E811-98A3-0025905A6080.root',
    #'/store/relval/CMSSW_10_0_0/RelValADDMonoJet_d3MD3_13/MINIAODSIM/100X_upgrade2018_realistic_v6_mahiON-v1/10000/D807F941-8200-E811-B66D-0CC47A4D7690.root'

    #'/store/relval/CMSSW_10_0_0/RelValADDMonoJet_d3MD3_13/MINIAODSIM/100X_upgrade2018_realistic_v6_mahiOFF-v1/10000/101FFE1F-EDFE-E711-8001-0CC47A4D7602.root',
    #'/store/relval/CMSSW_10_0_0/RelValADDMonoJet_d3MD3_13/MINIAODSIM/100X_upgrade2018_realistic_v6_mahiOFF-v1/10000/98F20620-EDFE-E711-AD93-0025905B8592.root'
                                )
)

# Lumi filter ==========================================================================================

###process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_13TeV_2017_HCAL_DCS_GOOD_post_JEC_bugfix.txt').getVLuminosityBlockRange()

# Electron ID ==========================================================================================

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era='2017-Nov17ReReco') 

### END Electron ID ====================================================================================


process.hltJetMetNtuple = cms.EDAnalyzer('HLTJetMETNtupleProducer',
                                         runJets = cms.bool(False),
                                         runMets = cms.untracked.bool(False),
                                         isData = cms.bool(False),
                                         PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                         MetCollectionTag = cms.InputTag('slimmedMETs'),
                                         GenMetCollectionTag = cms.InputTag('genMetTrue'),
                                         applyMETFilters = cms.bool(options.applyMETFilters),
                                         BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
                                         BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
                                         MuonCollectionTag = cms.InputTag('slimmedMuons'),
                                         ElectronCollectionTag = cms.InputTag('slimmedElectrons'),
                                         PFJetCollectionTag = cms.InputTag('slimmedJets'),
                                         CaloJetCollectionTag = cms.InputTag('slimmedCaloJets'),
                                         GenJetCollectionTag = cms.InputTag('slimmedGenJets'),
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
									      'HLT_PFMET110_PFMHT110_IDTight_v', 
                                                                              'HLT_CaloMET80_NotCleaned_v',
                                                                              'HLT_CaloMET90_NotCleaned_v',
                                                                              'HLT_CaloMET100_NotCleaned_v',
                                                                              'HLT_CaloMET110_NotCleaned_v',
                                                                              'HLT_CaloMET250_NotCleaned_v',
                                                                              'HLT_CaloMET70_HBHECleaned_v',
                                                                              'HLT_CaloMET80_HBHECleaned_v',
                                                                              'HLT_CaloMET90_HBHECleaned_v',
                                                                              'HLT_CaloMET100_HBHECleaned_v',
                                                                              'HLT_CaloMET250_HBHECleaned_v',
                                                                              'HLT_CaloMET300_HBHECleaned_v',
                                                                              'HLT_CaloMET350_HBHECleaned_v',
                                                                              'HLT_Ele25_WPTight_Gsf_v',
                                                                              'HLT_Ele25_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_Ele27_WPTight_Gsf_v',
                                                                              'HLT_Ele27_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_Ele30_WPTight_Gsf_v',
                                                                              'HLT_Ele30_eta2p1_WPTight_Gsf_v',
                                                                              'HLT_Ele32_WPTight_Gsf_v',
                                                                              'HLT_Ele32_eta2p1_WPTight_Gsf_v',
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

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)
process.p = cms.Path(process.egammaPostRecoSeq*
		     process.BadChargedCandidateFilter*
		     process.BadPFMuonFilter* 
		     process.hltJetMetNtuple
		     )
