import FWCore.ParameterSet.Config as cms

import FWCore.PythonUtilities.LumiList as LumiList

from JetMetTrigger.Analyzer.hltJetMETNtuple_cfi import *

process = cms.Process("hltJetMET")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
if runOnData:
  process.GlobalTag.globaltag = '92X_dataRun2_HLT_v7'
else:
  process.GlobalTag.globaltag = '92X_upgrade2017_TSG_For90XSamples_V2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
    #'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/C0CE5CF7-2671-E711-A65F-02163E014218.root',
    #'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/C81001BD-2C71-E711-BFAA-02163E01418D.root',
    #'/store/data/Run2017C/SingleElectron/MINIAOD/PromptReco-v1/000/299/594/00000/E03B5A6E-2F71-E711-B4D5-02163E01A4AD.root',
	
    '/store/data/Run2017F/SingleMuon/MINIAOD/PromptReco-v1/000/305/376/00000/2CFF0A26-D7BA-E711-BF4E-02163E0128ED.root'
 
    
                                )
)

configureJetMetNtuple(process)
process.ntuple = cms.Path(process.JetMetNtupleSequence)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)

