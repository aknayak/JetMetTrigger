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
#change GT according to sample
if runOnData:
  process.GlobalTag.globaltag = '94X_dataRun2_ReReco_EOY17_v2'
else:
  process.GlobalTag.globaltag = '92X_upgrade2017_TSG_For90XSamples_V2'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(

                              #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/105643A1-328F-E711-943D-02163E014641.root',
                              #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/1AE1B9B1-2F8F-E711-A3F6-02163E01441F.root',
                              #'/store/data/Run2017D/SingleElectron/MINIAOD/PromptReco-v1/000/302/031/00000/247DE91B-318F-E711-BC5C-02163E012AFE.root'
                              
                              '/store/data/Run2017F/SingleElectron/MINIAOD/17Nov2017-v1/50000/F2BC1274-82E0-E711-9FEB-0CC47A78A340.root'                              
	
                            )
                          )

configureJetMetNtuple(process)
process.ntuple = cms.Path(process.JetMetNtupleSequence)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("hltJetMetNtuple.root")
)

