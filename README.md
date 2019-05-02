# JetMetTrigger, Re-running of Jet Trigger sequence and production of trigger ntuple

For data in 2016 (07Aug17 Rereco), 2017 (prompt Reco), and 2018 (prompt reco)

Use the following release and instructions for setting up working area, for the re-running of trigger in data. 

https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGlobalHLT#CMSSW_10_1_X_2018_pp_data_taking

The following dataset is planned to be used. Remember to use miniaod along with RAW data (use useParent=1 or secondaryInputDataset option in crab)

For 2016:

/SingleMuon/Run2016H-07Aug17-v1/MINIAOD, /SingleMuon/Run2016H-v1/RAW


For 2017:

/SingleMuon/Run2017F-PromptReco-v1/MINIAOD, /SingleMuon/Run2017F-v1/RAW

For 2018:

/SingleMuon/Run2018D-PromptReco-v2/MINIAOD, /SingleMuon/Run2018D-v1/RAW


