### Instruction to run the L2TauNNTag test

## Part 1 : extract di-tau path from confDB
1. To get the CMSSW distribution:
\n\t export SCRAM_ARCH=slc7_amd64_gcc900
\n\t cmsrel CMSSW_12_1_0_pre3
\n\t cd CMSSW_12_1_0_pre3/src
\n\t cmsenv
\n\t git cms-init
\n\t git cms-addpkg HLTrigger/Configuration

1. You also need to have this repository in CMSSW_12_1_0_pre3/src. If you don't have it, the command to run is:
 cd CMSSW_12_1_0_pre3/src
\n\t cmsenv  
\n\t git cms-addpkg RecoTauTag/HLTProducers
\n\t git checkout -b L2TauTagNN
\n\t git push -u my-cmssw HEAD:L2TauTagNN
\n\t git checkout -b L2TauTagNN-debug
\n\t scram b -j8

1. You also need the RecoTauTag-Training files. If you don't have it, run the command:
\n\t git clone git@github.com:cms-data/RecoTauTag-TrainingFiles.git RecoTauTag/TrainingFiles

1. To get a user configuration file with the paths you need run the command: HLTriggerFirstPath,HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4,HLTGlobalPFTauHPSSequence,HLTriggerFinalPath --unprescale --cff > HLTrigger/Configuration/python/HLT_User_cff.py
\n \t please, check that all the paths have been correctly saved in the file CMSSW_12_1_0_pre3/src/HLTrigger/Configuration/HLT_User_cff.py \n\t otherwise you have to manually add all the missing paths (and the dependent sequences)

1. You need to authenticate to run the test command since it exploits files that are saved on store. So you need to run the command:
voms-proxy-init --rfc --voms cms

1. So you are now ready to run the test for L2TauTagNN
\n\tcd CMSSW_12_1_0_pre3/src
\n\tcmsenv
\n\tscram b -j 10
\n\tcmsRun RecoTauTag/HLTProducers/test/testL2TauTagNN.py


\n\tif you want to run on a limited number of events :
\n\tcmsRun RecoTauTag/HLTProducers/test/testL2TauTagNN.py maxEvents=20
