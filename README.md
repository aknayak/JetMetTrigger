# JetMetTrigger

Instructions to get Fall17V2 electron ID, taken from https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2

2018 MiniAOD
-------------

cmsrel CMSSW_10_2_10
cd CMSSW_10_2_10/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8

From CMSSW_10_2_10, our ID PR has been merged so it is no longer necessary to do "git cms-merge-topic cms-egamma:EgammaID_1023". If you work in CMSSW < 10_2_10, then you have to merge that PR still

To get the IDs embed in the pat::Electron/pat::Photon, just do


from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #as energy corrections are not yet availible for 2018
                       era='2018-Prompt')  
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)


For 2017 MiniAOD V2
---------------------

for
data: dataset=/*/Run2017*31Mar2018*/MINIAOD.
MC: dataset=/*/*12Apr2018*/MINIAODSIM

First setup up an area, following the instructions below. We will add a egamma helper function to this release, its simply one python file to simplify the process.

cmsrel CMSSW_9_4_13
cd CMSSW_9_4_13/src
cmsenv
git cms-init
git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier
scram b -j 8

Then to run it, copy and paste this config snippet to your config file. From CMSSW_9_4_13, our ID PR has been merged so it is no longer necessary to do "git cms-merge-topic cms-egamma:EgammaID_949". If you work in CMSSW < 9_4_13, then you have to merge that PR still. The results are the same, its just one less thing to merge.

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2017-Nov17ReReco')  
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

2016 MiniAOD "V2" aka V3
------------------------
As noted above, in das this is the V3 but is conceptually the same as the V2 2017 miniAOD hence why its refered to as V2.

for
data: dataset=/*/Run2016*17Jul2018*/MINIAOD
MC: being produced

The 2016 "V2" miniAOD can be used out of the box and can be used in the same release area as used for the 2017 reminiAOD. If you dont want the Fall17V2 IDs, the "git cms-merge-topic cms-egamma:EgammaPostRecoTools_940" step strictly speaking is unnecessary as you dont need the script to apply any corrections, however it will do no harm if its in the release area so you might as well have it, particularly as its needed for reading the V1 miniAOD which is still needed for MC.

To apply the Fall17V2 IDs simply do

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False, #corrections by default are fine so no need to re-run
                       era='2016-Legacy')  
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)


Look at the instructions here to convert 2016 V1 miniAOD to V2 miniAOD on fly. 
https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Converting_V1_MiniAOD_to_V2_Mini

Accessing ID result
To determine if the electron or photon passes or fails and ID, simply do

 
pat::Electron::electronID("<name of id>");
The electron method returns a float which is either 0 (fail) or 1 (pass). For some out dated IDs, this can also be other values but for all modern recommended IDs, this is 0 or 1. 

The supported IDs are for electrons:
cutBasedElectronID-Fall17-94X-V1-loose cutBasedElectronID-Fall17-94X-V1-medium cutBasedElectronID-Fall17-94X-V1-tight cutBasedElectronID-Fall17-94X-V1-veto 
cutBasedElectronID-Summer16-80X-V1-loose cutBasedElectronID-Summer16-80X-V1-medium cutBasedElectronID-Summer16-80X-V1-tight cutBasedElectronID-Summer16-80X-V1-veto
heepElectronID-HEEPV70
mvaEleID-Fall17-iso-V1-wp80 mvaEleID-Fall17-iso-V1-wp90 mvaEleID-Fall17-iso-V1-wpLoose
mvaEleID-Fall17-noIso-V1-wp80,mvaEleID-Fall17-noIso-V1-wp90,mvaEleID-Fall17-noIso-V1-wpLoose
mvaEleID-Spring16-GeneralPurpose-V1-wp80 mvaEleID-Spring16-GeneralPurpose-V1-wp90
mvaEleID-Spring16-HZZ-V1-wpLoose



