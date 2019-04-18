import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList


isData = True
runOnData=isData #data/MC switch
runOnElectronPD=False
runOnMuonPD=True

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
#Change the GT according to sample
if runOnData:
  process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'
else:
  process.GlobalTag.globaltag = '92X_upgrade2017_TSG_For90XSamples_V2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(

    #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/105643A1-328F-E711-943D-02163E014641.root',
    #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/1AE1B9B1-2F8F-E711-A3F6-02163E01441F.root',
    #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/247DE91B-318F-E711-BC5C-02163E012AFE.root'

    '/store/data/Run2017F/SingleElectron/MINIAOD/17Nov2017-v1/50000/F2BC1274-82E0-E711-9FEB-0CC47A78A340.root'
                                )
)

# Lumi filter ==========================================================================================

###process.source.lumisToProcess = LumiList.LumiList(filename = 'Cert_13TeV_2017_HCAL_DCS_GOOD_post_JEC_bugfix.txt').getVLuminosityBlockRange()

# Electron ID ==========================================================================================

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       era='2017-Nov17ReReco')

### END Electron ID ====================================================================================

### Add Event filters to reduce ntuple size ============================================================

process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlineSlimmedPrimaryVertices"),
   cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && abs(position.Rho) <= 2"), 
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.triggerSelection = cms.EDFilter("TriggerResultsFilter",
                                        triggerConditions = cms.vstring('HLT_IsoMu24_v*', 'HLT_IsoMu24_eta2p1_v*', 'HLT_IsoMu27_v*', 'HLT_IsoMu30_v*', 'HLT_Ele25_WPTight_Gsf_v*', 'HLT_Ele25_eta2p1_WPTight_Gsf_v*', 'HLT_Ele27_WPTight_Gsf_v*', 'HLT_Ele27_eta2p1_WPTight_Gsf_v*', 'HLT_Ele30_WPTight_Gsf_v*', 'HLT_Ele30_eta2p1_WPTight_Gsf_v*', 'HLT_Ele32_WPTight_Gsf_v*', 'HLT_Ele32_eta2p1_WPTight_Gsf_v*', 'HLT_Ele35_WPTight_Gsf_v*'),
                                        hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                                        l1tResults = cms.InputTag( "" ),
                                        throw = cms.bool(False)
                                        )

process.SelectedMuons = cms.EDFilter("PATMuonSelector",
                             src = cms.InputTag("slimmedMuons"),
                             cut = cms.string("pt > 25 && " + "abs(eta) < 2.4" +
                                              " && isGlobalMuon"+ 
                                              " && passed('CutBasedIdMediumPrompt')"+
                                              " && (pfIsolationR04.sumChargedHadronPt+max(0.,(pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt-0.5*pfIsolationR04.sumPUPt)))/pt < 0.15"
                                              )
                             )

process.minMuonFilter = cms.EDFilter("PATCandViewCountFilter",
                                     src = cms.InputTag('SelectedMuons'),
                                     minNumber = cms.uint32(1),
                                     maxNumber = cms.uint32(999999)
                                     )

process.selectedElectrons = cms.EDFilter("PATElectronSelector",
                                         src = cms.InputTag("slimmedElectrons"),
                                         cut = cms.string("pt > 25 && " + "abs(eta) < 2.4"+
                                                          " && (electronID('mvaEleID-Spring16-GeneralPurpose-V1-wp90')"+
                                                          " || electronID('mvaEleID-Fall17-iso-V2-wp90')"+
                                                          " || electronID('cutBasedElectronID-Summer16-80X-V1-medium')"+
                                                          " || electronID('cutBasedElectronID-Fall17-94X-V2-medium'))"+
                                                          " && (pfIsolationVariables.sumChargedHadronPt+max(0.,(pfIsolationVariables.sumNeutralHadronEt+pfIsolationVariables.sumPhotonEt-0.5*pfIsolationVariables.sumPUPt)))/pt < 0.15"
                                                          )
                                         )

process.minElectronFilter = cms.EDFilter("PATCandViewCountFilter",
                                         src = cms.InputTag('selectedElectrons'),
                                         minNumber = cms.uint32(1),
                                         maxNumber = cms.uint32(999999)
                                         )

if runOnMuonPD:
  process.leptonFilterSequence = cms.Sequence(process.SelectedMuons*process.minMuonFilter)
elif runOnElectronPD:
  process.leptonFilterSequence = cms.Sequence(process.selectedElectrons*process.minElectronFilter)
else:
  process.leptonFilterSequence = cms.Sequence()

######============================================================================================

process.hltJetMetNtuple = cms.EDAnalyzer('HLTJetMETNtupleProducer',
                                         runJets = cms.bool(False),
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

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)
process.p = cms.Path(process.primaryVertexFilter*
                     process.triggerSelection*
                     process.egammaPostRecoSeq*
                     process.leptonFilterSequence*
		     process.hltJetMetNtuple
		     )
