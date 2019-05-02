import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

isData = True
runOnData=isData #data/MC switch
runOnElectronPD=False
runOnMuonPD=True

def configureJetMetNtuple(process):
  #  Options parsing  #
  ###################
  
  from FWCore.ParameterSet.VarParsing import VarParsing
  
  #options = VarParsing('analysis')
  #options.register('applyMETFilters',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Apply MET filters')
  #options.register('runJets',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Store Jet Variables')

  process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
  process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
  process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

  process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
  process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
  process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

  # Electron ID ==========================================================================================

  #from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
  # turn on VID producer, indicate data format  to be
  # DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
  useAOD = False
  
  if useAOD == True :
    dataFormat = DataFormat.AOD
  else :
    dataFormat = DataFormat.MiniAOD
    
    switchOnVIDElectronIdProducer(process, dataFormat)
    
  # define which IDs we want to produce
  my_id_modules = [#'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
    #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff',
    #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff',
    'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff'
    ]


  #add them to the VID producer
  for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

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
                                                        " && globalTrack.normalizedChi2 < 3" +
                                                        " && combinedQuality.chi2LocalPosition < 12" +
                                                        " && combinedQuality.trkKink < 20" +
                                                        " && innerTrack.validFraction > 0.49" +
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
                                           runJets = cms.bool(True),
                                           runMets = cms.untracked.bool(True),
                                           isData = cms.bool(isData),
                                           PVCollectionTag = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           MetCollectionTag = cms.InputTag('slimmedMETs'),
                                           applyMETFilters = cms.bool(False),
                                           BadMuonFilter              = cms.InputTag("BadPFMuonFilter",""),
                                           BadChargedCandidateFilter = cms.InputTag("BadChargedCandidateFilter",""),
                                           MuonCollectionTag = cms.InputTag('slimmedMuons'),
                                           ElectronCollectionTag = cms.InputTag('slimmedElectrons'),
                                           PFJetCollectionTag = cms.InputTag('slimmedJets'),
                                           HLTPFJetCollectionTag = cms.InputTag('hltAK4PFJetsCorrected'),
                                           HLTCaloJetCollectionTag = cms.InputTag('hltAK4CaloJetsCorrected'),
                                           HLTPFMetTag = cms.InputTag('hltPFMETProducer'),
                                           HLTPFMetType1Tag = cms.InputTag('hltPFMETTypeOne'),
                                           HLTCaloMetTag = cms.InputTag('hltMet'),
                                           eleMvaSpring16WPMediumMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
                                           eleMvaSpring16WPTightMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
                                           mvaSpring16ValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                                           mvaSpring16CategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                                           eleSummer16VetoIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
                                           eleSummer16LooseIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                                           eleSummer16MediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
                                           eleSummer16TightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
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
    
  process.JetMetNtupleSequence = cms.Sequence(process.primaryVertexFilter*
                                              process.triggerSelection*
                                              process.leptonFilterSequence*
                                              process.egmGsfElectronIDSequence*
                                              process.BadChargedCandidateFilter*
                                              process.BadPFMuonFilter* 
                                              process.hltJetMetNtuple
                                              )
