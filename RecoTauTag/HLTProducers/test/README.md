## Tests for tau HLT producers


### Example on how to run di-tau trigger on H->TauTau vbf sample - Get the path

To get the CMSSW distribution:
```
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_12_1_0_pre3
cd CMSSW_12_1_0_pre3/src
cmsenv
git cms-init
git cms-addpkg HLTrigger/Configuration
```
To get a user configuration file with the paths you need run the command:
```
hltGetConfiguration HLTriggerFirstPath,HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4,HLTGlobalPFTauHPSSequence,HLTriggerFinalPath --unprescale --cff > HLTrigger/Configuration/python/HLT_User_cff.py
```

please, check that all the paths have been correctly saved in the file CMSSW_12_1_0_pre3/src/HLTrigger/Configuration/HLT_User_cff.py *  otherwise you have to manually add all the missing paths (and the dependent sequences). I had to manually inser the HLTGlobalPFTauHPSSequence.

### Installing the RecoTauTag repository
You also need to have this repository in CMSSW_12_1_0_pre3/src. If you don't have it, the command to run is:
```
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_12_1_0_pre3
cd CMSSW_12_1_0_pre3/src
cmsenv  
git cms-merge-topic valeriadamante:L2TauTagNN
git cms-addpkg RecoTauTag/HLTProducers
scram b -j8
```

You also need the RecoTauTag-Training files. If you don't have it, run the command:
```
git clone git@github.com:cms-data/RecoTauTag-TrainingFiles.git RecoTauTag/TrainingFiles
```

### L2TauNNTag test
You need to authenticate to run the test command since it exploits files that are saved on store. So you need to run the command:
```
voms-proxy-init --rfc --voms cms
```

So you are now ready to run the test for L2TauTagNN
```
cd CMSSW_12_1_0_pre3/src
cmsenv
scram b -j 10
cmsRun RecoTauTag/HLTProducers/test/testL2TauTagNN.py
```

if you want to run on a limited number of events :
```
cmsRun RecoTauTag/HLTProducers/test/testL2TauTagNN.py maxEvents=20
```
